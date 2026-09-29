"""Predeclared composition search; physical masses require measured DHS opacities."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from silicate_structure_design import MATERIALS, NUMERICS, FIXED, column_normalization
from silicate_pilot_design import ICE_NAME, LAB_SHA256, TEMPLATE_SHA256, digest, require

NAME = "silicate_search_v2"
RUN_NAME = "silicate_search_v2_512k"
REFERENCE_MASS = 2.25e-4
SPECIES = {key: dict(MATERIALS[key]) for key in ("draine", "pyroxene_mg50", "olivine_mg50")}
SPECIES.update(carbon=dict(filename="ac_opct.dat", density_g_cm3=1.8,
    sha256="ef8f09849580a0859f29f61216486e7e249103759f851bcfa7abe11c09a44de2"),
    ice=dict(filename=ICE_NAME, density_g_cm3=.94, sha256=LAB_SHA256))
DUST = dict(family="separate_DHS_silicate_carbon_H2O30K", species=SPECIES,
    amin_um=.03, amax_um=.4, size_exponent=2.75, grain_bins=50,
    vmax=.1, porosity=0., ice_laboratory_temperature_k=30.,
    temperature_dependent_ice_survival=False, optical_tables_modified=False,
    native_out_of_range_treatment_accepted=True, tail_physical_accuracy_certified=False,
    carbon_size_distribution_assumption="Same as envelope silicate; disk remains 70/30 silicate/carbon Mie")
SUCCESS = dict(silicate_log_shape_abs_dex_max=.10, ice_raw_log_residual_abs_dex_max=.15,
    long_miri_raw_log_residual_abs_dex_max=.15, nir_continuum_log_residual_rms_dex_max=.20,
    ice_bands=["h02", "h03"], long_miri_bands=["s05", "s06"],
    ice_normalized_shapes_also_reported=True, h04_in_success_criterion=False,
    statistical_acceptance_region=False)


def declared_cases():
    cases, seen = [], {}
    def add(stage, material, carbon, ice, inclination, exponent):
        key = (material, carbon, ice, inclination, exponent)
        if key in seen:
            cases[seen[key]]["stages"].append(stage)
            return
        seen[key] = len(cases)
        cases.append(dict(index=len(cases), stages=[stage],
            role="composition" if stage == 1 else "material_structure",
            material_id=material, carbon_mass_fraction=carbon, ice_mass_fraction=ice,
            silicate_mass_fraction=float(format(1-carbon-ice, ".13g")),
            inclination_deg=inclination, density_exponent=exponent))
    for carbon in (0., .1, .2, .3):
        for ice in (.04, .08, .12):
            for inclination in (60., 70.):
                add(1, "draine", carbon, ice, inclination, -1.5)
    for material in ("draine", "pyroxene_mg50", "olivine_mg50"):
        for exponent in (-1.5, -1.75):
            for ice in (.08, .12):
                add(2, material, .2, ice, 60., exponent)
    require(len(cases) == 34 and sum(2 in r["stages"] for r in cases) == 12,
            "Wrong predeclared v2 catalogue")
    return cases


def validate_opacities(opacities):
    require(set(opacities) == set(SPECIES), "Require all five measured DHS species")
    first = None
    for species, row in opacities.items():
        wave = np.asarray(row["wavelength_um"], dtype=float)
        require(wave.ndim == 1 and len(wave) >= 3 and np.all(np.isfinite(wave))
                and np.all(wave > 0) and np.all(np.diff(wave) > 0), "Invalid DHS wavelength grid")
        if first is None:
            first = wave
        require(np.array_equal(wave, first), "Pure species use different wavelength grids")
        arrays = {k: np.asarray(row[k], dtype=float) for k in
            ("kappa_ext_cm2_g", "kappa_abs_cm2_g", "kappa_sca_cm2_g", "albedo")}
        require(all(v.shape == wave.shape and np.all(np.isfinite(v)) and np.all(v >= 0)
                    for v in arrays.values()), "Nonfinite or negative DHS opacity")
        require(np.all(arrays["albedo"] <= 1), "Invalid DHS albedo")
        require(np.allclose(arrays["kappa_abs_cm2_g"]+arrays["kappa_sca_cm2_g"], arrays["kappa_ext_cm2_g"], rtol=1e-8, atol=1e-12),
                "DHS absorption/scattering do not sum to extinction")
        require(np.allclose(arrays["kappa_sca_cm2_g"], arrays["kappa_ext_cm2_g"]*arrays["albedo"], rtol=1e-8, atol=1e-12),
                "DHS albedo disagrees with scattering/extinction")
        for wavelength in (2.2, 9.7, 18.):
            require(np.count_nonzero(np.isclose(wave, wavelength, rtol=2e-6, atol=0)) == 1,
                    f"Missing or ambiguous DHS normalization wavelength {wavelength}")
    return opacities


def sample(opacities, species, field, wavelength):
    row = opacities[species]
    indices = np.flatnonzero(np.isclose(row["wavelength_um"], wavelength, rtol=2e-6, atol=0))
    require(len(indices) == 1, "Missing or ambiguous measured opacity wavelength")
    return float(row[field][indices[0]])


def matched_composition(opacities, material, carbon, ice, exponent=-1.5):
    fractions = {material: float(format(1-carbon-ice, ".13g")), "carbon": carbon, "ice": ice}
    require(all(math.isfinite(v) and v >= 0 for v in fractions.values())
            and fractions[material] > 0 and abs(sum(fractions.values())-1) < 1e-12,
            "Invalid mixture mass fractions")
    def mixture(field, wavelength):
        return sum(fraction*sample(opacities, species, field, wavelength) for species, fraction in fractions.items())
    reference_ext = .96*sample(opacities, "draine", "kappa_ext_cm2_g", 2.2)+.04*sample(opacities, "ice", "kappa_ext_cm2_g", 2.2)
    reference_abs = .96*sample(opacities, "draine", "kappa_abs_cm2_g", 9.7)+.04*sample(opacities, "ice", "kappa_abs_cm2_g", 9.7)
    ext, absorption, long_absorption = mixture("kappa_ext_cm2_g", 2.2), mixture("kappa_abs_cm2_g", 9.7), mixture("kappa_abs_cm2_g", 18.)
    require(min(reference_ext, reference_abs, ext, absorption) > 0, "Invalid DHS matching opacity")
    profile_factor = column_normalization(exponent, REFERENCE_MASS)["mass_ratio"]
    mass = float(format(REFERENCE_MASS*reference_ext/ext*profile_factor, ".13g"))
    return dict(reference_mass_msun=REFERENCE_MASS, normalization_wavelength_um=2.2,
        reference_extinction_cm2_g=reference_ext, mixture_extinction_cm2_g=ext,
        mixture_absorption_9p7_cm2_g=absorption, mixture_absorption_18_cm2_g=long_absorption,
        S=absorption/ext, S_over_reference=(absorption/ext)/(reference_abs/reference_ext),
        absorption_18_over_9p7=long_absorption/absorption,
        composition_mass_factor=reference_ext/ext, profile_mass_factor=profile_factor,
        actual_mass_msun=mass, absolute_ice_mass_msun=mass*ice,
        mass_significant_digits=13, discrete_optical_depth_match_certified=False,
        method="mass-weighted pure-species DHS extinction at 2.2 um, then analytic radial-column profile factor")


def catalogue(opacities):
    validate_opacities(opacities)
    records = declared_cases()
    for record in records:
        match = matched_composition(opacities, record["material_id"], record["carbon_mass_fraction"],
                                    record["ice_mass_fraction"], record["density_exponent"])
        record.update(opacity_matching=match, reference_mass_msun=REFERENCE_MASS,
            parameters={**FIXED, "envelope_dust_mass_msun": match["actual_mass_msun"],
                "envelope_silicate_file": SPECIES[record["material_id"]]["filename"],
                "envelope_carbon_mass_fraction": record["carbon_mass_fraction"],
                "envelope_ice_mass_fraction": record["ice_mass_fraction"],
                "envelope_density_exponent": record["density_exponent"],
                "inclination_deg": record["inclination_deg"], "envelope_amax_um": .4})
    return records


def validate_frozen_campaign(bundle, experiment, manifest, index=None):
    import sys
    bundle = Path(bundle)
    sys.path.insert(0, str(bundle/"code/src"))
    from mcfost_grid.physics import render_parameter
    require(experiment["experiment_id"] == NAME and len(manifest["models"]) == 34,
            "Wrong silicate v2 experiment")
    opacities = json.loads((bundle/"inputs/stage0_opacities.json").read_text())
    wave_plan = json.loads((bundle/"inputs/stage0_wavelengths.json").read_text())
    required_wave = sorted(set([a["wavelength_um"] for a in manifest["anchors"]]+[2.2, 9.7, 18.]))
    require(len(required_wave) == 99 and wave_plan["wavelengths_um"] == required_wave
            and all(row["wavelength_um"] == required_wave for row in opacities.values()),
            "Frozen DHS opacities must cover all 99 production/screen wavelengths")
    require(experiment["stage0_opacity"] == opacities, "Reported Stage 0 opacities changed")
    require(digest(bundle/"inputs/stage0_receipt.json") == experiment["stage0_receipt_sha256"],
            "Copied Stage 0 receipt changed")
    receipt = json.loads((bundle/"inputs/stage0_receipt.json").read_text())
    require(receipt["status"] == "passed" and receipt["wavelengths_um"] == required_wave
            and [row["id"] for row in receipt["species"]] == list(SPECIES),
            "Frozen Stage 0 receipt is incomplete")
    expected = catalogue(opacities)
    require(experiment["catalogue"] == expected, "Measured-opacity matched catalogue changed")
    require(experiment["opacity_matching"] == [dict(index=row["index"], **row["opacity_matching"]) for row in expected],
            "Reported opacity matching changed")
    require(experiment["dust_prescription"] == DUST and experiment["success_criteria"] == SUCCESS,
            "V2 dust prescription or criteria changed")
    require(experiment["numerics"] == manifest["configuration"]["numerics"] == NUMERICS, "V2 numerics changed")
    require(manifest["temperature_policy"] == "fresh equilibrium for every physical model", "Fresh temperatures required")
    require(digest(bundle/"inputs/template.para") == TEMPLATE_SHA256, "V02 template changed")
    for spec in SPECIES.values():
        require(digest(bundle/"inputs/utils/Dust"/spec["filename"]) == spec["sha256"], "Optical constants changed")
    for model, row in zip(manifest["models"], expected):
        require(model["index"] == row["index"] and model["parameters"] == row["parameters"], "V2 physical model changed")
    template = (bundle/"inputs/template.para").read_text()
    selected = manifest["models"] if index is None else [manifest["models"][index]]
    for model in selected:
        folder = bundle/"models"/model["id"]
        require((folder/"temperature.para").read_text() == render_parameter(template, model["parameters"], NUMERICS, "temperature"), "V2 temperature input changed")
        for anchor in manifest["anchors"]:
            for stage in ("image", "coeval"):
                require((folder/"anchors"/anchor["id"]/(stage+".para")).read_text() == render_parameter(template, model["parameters"], NUMERICS, stage, anchor["wavelength_um"]), "V2 image input changed")
