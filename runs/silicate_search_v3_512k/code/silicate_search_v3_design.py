"""Declared profile/column search with ice fixed by measured v2 column targets."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from silicate_pilot_design import TEMPLATE_SHA256, digest, require
from silicate_structure_design import FIXED, column_normalization
from silicate_search_v2_design import (DUST as V2_DUST, NUMERICS as V2_NUMERICS,
                                       SPECIES as V2_SPECIES, SUCCESS as V2_SUCCESS,
                                       sample)


NAME = "silicate_search_v3"
RUN_NAME = "silicate_search_v3_512k"
REFERENCE_MASS = 2.25e-4
MATERIALS = ("draine", "pyroxene_mg50")
CARBON_FRACTIONS = (.15, .20, .25)
CONFIGURATIONS = (
    ("anchor", -1.5, 1.),
    ("shallow", -1.25, .86),
    ("flat_high", -1., .83),
    ("flat_low", -1., .66),
)
SPECIES = {key: dict(V2_SPECIES[key]) for key in (*MATERIALS, "carbon", "ice")}
NUMERICS = dict(V2_NUMERICS)
DUST = {**V2_DUST, "species": SPECIES,
        "ice_fraction_policy": "Derived from the exact per-silicate v2 median bracketed ice-column crossing"}
SUCCESS = {**V2_SUCCESS, "ice_success_quantity": "raw_flux",
           "seed_margin_policy": "Success within the measured seed difference of a tolerance is marginal",
           "seed_difference_warning_dex": .03}
ICE_FRACTION_FORMULA = (
    "f_ice=T*A/(M_ref*c*k_ext_ref-T*(k_ice-k_sil)); "
    "A=(1-f_C)*k_sil+f_C*k_C; all k at 2.2 um"
)


def _finite_number(value, name, *, positive=False):
    require(isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value), f"{name} must be a finite number")
    require(not positive or value > 0, f"{name} must be positive")
    return float(value)


def _freeze(value):
    """Match prior campaign's stable cross-platform 13-significant-digit inputs."""
    return float(format(value, ".13g"))


def validate_targets(targets):
    require(isinstance(targets, dict) and set(targets) == set(MATERIALS),
            "Require exact Draine and pyroxene ice-column targets")
    for material, value in targets.items():
        _finite_number(value, f"{material} ice-column target", positive=True)
    return targets


def declared_cases():
    cases = []
    for configuration, exponent, factor in CONFIGURATIONS:
        for material in MATERIALS:
            for carbon in CARBON_FRACTIONS:
                cases.append(dict(index=len(cases), configuration_id=configuration,
                    density_exponent=exponent, column_factor=factor,
                    material_id=material, carbon_mass_fraction=carbon,
                    inclination_deg=70., random_seed=43001, replica_of=None,
                    role="primary"))
    for original in (7, 10):
        cases.append({**cases[original], "index": len(cases), "random_seed": 43002,
                      "replica_of": original, "role": "seed_replica"})
    require(len(cases) == 26, "Wrong v3 run count")
    return cases


def model_numerics(record):
    seed = record["random_seed"]
    require(type(seed) is int and seed in (43001, 43002), "Undeclared v3 photon seed")
    return {**NUMERICS, "random_seed": seed}


def validate_opacities(opacities):
    require(set(opacities) == set(SPECIES), "Require all four measured DHS species")
    first = None
    for row in opacities.values():
        wave = np.asarray(row["wavelength_um"], dtype=float)
        require(wave.shape == (99,) and np.all(np.isfinite(wave))
                and np.all(wave > 0) and np.all(np.diff(wave) > 0),
                "DHS opacities require 99 increasing positive wavelengths")
        if first is None:
            first = wave
        require(np.array_equal(wave, first), "Pure species use different wavelength grids")
        arrays = {key: np.asarray(row[key], dtype=float) for key in
                  ("kappa_ext_cm2_g", "kappa_abs_cm2_g", "kappa_sca_cm2_g", "albedo")}
        require(all(array.shape == wave.shape and np.all(np.isfinite(array)) and np.all(array >= 0)
                    for array in arrays.values()), "Nonfinite or negative DHS opacity")
        require(np.all(arrays["albedo"] <= 1), "Invalid DHS albedo")
        require(np.allclose(arrays["kappa_abs_cm2_g"] + arrays["kappa_sca_cm2_g"],
                            arrays["kappa_ext_cm2_g"], rtol=1e-8, atol=1e-12),
                "DHS absorption/scattering do not sum to extinction")
        require(np.allclose(arrays["kappa_sca_cm2_g"], arrays["kappa_ext_cm2_g"]*arrays["albedo"],
                            rtol=1e-8, atol=1e-12), "DHS albedo disagrees with scattering/extinction")
        for wavelength in (2.2, 9.7, 18.):
            require(np.count_nonzero(np.isclose(wave, wavelength, rtol=2e-6, atol=0)) == 1,
                    f"Missing or ambiguous DHS normalization wavelength {wavelength}")
    return opacities


def matched_composition(opacities, material, carbon, exponent, column_factor, target):
    require(material in MATERIALS, "V3 supports Draine and pyroxene only")
    carbon = _finite_number(carbon, "Carbon mass fraction")
    require(0 <= carbon < 1, "Carbon mass fraction must lie in [0, 1)")
    factor = _finite_number(column_factor, "Column factor", positive=True)
    target = _finite_number(target, "Ice-column target", positive=True)
    exponent = _finite_number(exponent, "Density exponent")
    require(-3 < exponent < 0, "Density exponent must lie strictly between -3 and 0")
    ext = {name: sample(opacities, name, "kappa_ext_cm2_g", 2.2) for name in SPECIES}
    reference_ext = .96*ext["draine"] + .04*ext["ice"]
    reference_abs = (.96*sample(opacities, "draine", "kappa_abs_cm2_g", 9.7)
                     + .04*sample(opacities, "ice", "kappa_abs_cm2_g", 9.7))
    require(reference_ext > 0 and reference_abs > 0, "Invalid reference DHS matching opacity")
    base = (1-carbon)*ext[material] + carbon*ext["carbon"]
    denominator = REFERENCE_MASS*factor*reference_ext - target*(ext["ice"]-ext[material])
    require(math.isfinite(denominator) and denominator > 0, "Ice-column target has no positive denominator")
    ice = _freeze(target*base/denominator)
    silicate = _freeze(1-carbon-ice)
    require(ice > 0 and silicate > 0 and carbon+ice < 1,
            "Ice-column target cannot be reached with positive silicate and physical ice fractions")
    fractions = {material: silicate, "carbon": carbon, "ice": ice}

    def mixture(field, wavelength):
        return sum(fraction*sample(opacities, species, field, wavelength)
                   for species, fraction in fractions.items())

    extinction = mixture("kappa_ext_cm2_g", 2.2)
    absorption = mixture("kappa_abs_cm2_g", 9.7)
    long_absorption = mixture("kappa_abs_cm2_g", 18.)
    require(all(math.isfinite(v) and v > 0 for v in (extinction, absorption)),
            "Invalid mixture DHS matching opacity")
    profile = column_normalization(exponent, REFERENCE_MASS)["mass_ratio"]
    mass = _freeze(REFERENCE_MASS*factor*reference_ext/extinction*profile)
    radial_ice = ice*mass/profile
    require(math.isclose(radial_ice, target, rel_tol=2e-12, abs_tol=0),
            "Frozen ice fraction and dust mass do not preserve the declared ice column")
    return dict(reference_mass_msun=REFERENCE_MASS, normalization_wavelength_um=2.2,
        reference_extinction_cm2_g=reference_ext, mixture_extinction_cm2_g=extinction,
        mixture_absorption_9p7_cm2_g=absorption, mixture_absorption_18_cm2_g=long_absorption,
        S=absorption/extinction, S_over_reference=(absorption/extinction)/(reference_abs/reference_ext),
        absorption_18_over_9p7=long_absorption/absorption, composition_mass_factor=reference_ext/extinction,
        profile_mass_factor=profile, column_factor=factor, actual_mass_msun=mass,
        ice_mass_fraction=ice, silicate_mass_fraction=silicate, absolute_ice_mass_msun=mass*ice,
        ice_column_target_msun=target, radial_ice_column_msun=radial_ice,
        ice_fraction_formula=ICE_FRACTION_FORMULA, ice_fraction_significant_digits=13,
        mass_significant_digits=13, discrete_optical_depth_match_certified=False,
        method="mass-weighted DHS extinction at 2.2 um times declared column and analytic profile factors; ice from fixed radial-column target")


def catalogue(opacities, targets):
    validate_opacities(opacities)
    validate_targets(targets)
    records = declared_cases()
    for record in records:
        match = matched_composition(opacities, record["material_id"], record["carbon_mass_fraction"],
            record["density_exponent"], record["column_factor"], targets[record["material_id"]])
        record.update(opacity_matching=match, reference_mass_msun=REFERENCE_MASS,
            ice_mass_fraction=match["ice_mass_fraction"], silicate_mass_fraction=match["silicate_mass_fraction"],
            numerics=model_numerics(record), parameters={**FIXED,
                "envelope_dust_mass_msun": match["actual_mass_msun"],
                "envelope_silicate_file": SPECIES[record["material_id"]]["filename"],
                "envelope_carbon_mass_fraction": record["carbon_mass_fraction"],
                "envelope_ice_mass_fraction": match["ice_mass_fraction"],
                "envelope_density_exponent": record["density_exponent"],
                "inclination_deg": 70., "envelope_amax_um": .4})
    return records


def validate_frozen_campaign(bundle, experiment, manifest, index=None):
    import sys
    bundle = Path(bundle)
    sys.path.insert(0, str(bundle/"code/src"))
    from mcfost_grid.physics import render_parameter
    require(experiment["experiment_id"] == NAME and len(manifest["models"]) == 26,
            "Wrong silicate v3 experiment")
    opacities = json.loads((bundle/"inputs/stage0_opacities.json").read_text())
    wave_plan = json.loads((bundle/"inputs/stage0_wavelengths.json").read_text())
    required_wave = sorted(set([anchor["wavelength_um"] for anchor in manifest["anchors"]]+[2.2, 9.7, 18.]))
    require(len(manifest["anchors"]) == 98 and len(required_wave) == 99
            and wave_plan["wavelengths_um"] == required_wave
            and all(row["wavelength_um"] == required_wave for row in opacities.values()),
            "Frozen DHS opacities must cover all 99 production/screen wavelengths")
    require(experiment["stage0_opacity"] == opacities, "Reported Stage 0 opacities changed")
    require(digest(bundle/"inputs/stage0_receipt.json") == experiment["stage0_receipt_sha256"],
            "Copied Stage 0 receipt changed")
    receipt = json.loads((bundle/"inputs/stage0_receipt.json").read_text())
    require(receipt["status"] == "passed" and receipt["wavelengths_um"] == required_wave
            and [row["id"] for row in receipt["species"]] == list(SPECIES),
            "Frozen Stage 0 receipt is incomplete")
    reference = json.loads((bundle/"inputs/v2_reference.json").read_text())
    targets = experiment["ice_column_targets_msun"]
    require(targets == reference["ice_column_targets_msun"], "Exact v2 ice-column targets changed")
    expected = catalogue(opacities, targets)
    require(experiment["catalogue"] == expected, "Measured-opacity/ice-column matched catalogue changed")
    require(experiment["opacity_matching"] == [dict(index=row["index"], **row["opacity_matching"]) for row in expected],
            "Reported opacity matching changed")
    require(experiment["dust_prescription"] == DUST and experiment["success_criteria"] == SUCCESS,
            "V3 dust prescription or criteria changed")
    require(experiment["numerics"] == manifest["configuration"]["numerics"] == NUMERICS,
            "V3 base numerics changed")
    require(manifest["temperature_policy"] == "fresh equilibrium for every physical model",
            "Fresh temperatures required for every run, including seed replicas")
    require(digest(bundle/"inputs/template.para") == TEMPLATE_SHA256, "V02 template changed")
    for species in SPECIES.values():
        require(digest(bundle/"inputs/utils/Dust"/species["filename"]) == species["sha256"],
                "Optical constants changed")
    models = manifest["models"]
    require(len({model["id"] for model in models}) == 26, "V3 seed replicas need distinct model IDs")
    require(len({json.dumps(model["parameters"], sort_keys=True) for model in models}) == 24,
            "V3 requires 24 physical models and two exact seed replicas")
    for model, row in zip(models, expected):
        require(model["index"] == row["index"] and model["parameters"] == row["parameters"],
                "V3 physical model changed")
        require(model["numerics"] == model_numerics(row), "V3 per-model numerical seed changed")
    for replica, original in ((24, 7), (25, 10)):
        require(models[replica]["parameters"] == models[original]["parameters"]
                and models[replica]["numerics"]["random_seed"] == 43002
                and models[original]["numerics"]["random_seed"] == 43001,
                "V3 seed replica identity changed")
    if index is not None:
        require(type(index) is int and 0 <= index < 26, "V3 model index out of range")
    template = (bundle/"inputs/template.para").read_text()
    for model in (models if index is None else [models[index]]):
        folder = bundle/"models"/model["id"]
        numerics = model["numerics"]
        require((folder/"temperature.para").read_text() == render_parameter(template, model["parameters"], numerics, "temperature"),
                "V3 temperature input changed")
        for anchor in manifest["anchors"]:
            for stage in ("image", "coeval"):
                require((folder/"anchors"/anchor["id"]/(stage+".para")).read_text()
                        == render_parameter(template, model["parameters"], numerics, stage, anchor["wavelength_um"]),
                        "V3 image input changed")
