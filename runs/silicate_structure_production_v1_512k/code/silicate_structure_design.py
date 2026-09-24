"""Predeclared mass, silicate and radial-envelope contrasts; no fit selection."""
from __future__ import annotations

import json
import math
from pathlib import Path

from silicate_pilot_design import (ICE_NAME, LAB_SHA256, DRAINE_SHA256, TEMPLATE_SHA256,
                                  FIXED as PILOT_FIXED, NUMERICS as PILOT_NUMERICS,
                                  DUST as PILOT_DUST, digest, optical_table, require,
                                  validate_ice)

NAME = "silicate_structure_production_v1"
RUN_NAME = "silicate_structure_production_v1_512k"
NUMERICS = dict(PILOT_NUMERICS)
MATERIAL_MANIFEST_SHA256 = "5d18ab8f7013f9eb30dbbc179c6aae286ef74ce47c9269901d1d04b2a4fd49f6"
MATERIALS = {
    "draine": dict(filename="Draine_Si_sUV.dat", sha256=DRAINE_SHA256,
                   density_g_cm3=3.5, wavelength_min_um=.001, wavelength_max_um=1000.,
                   rows=1201, label="Draine astronomical silicate"),
    "olivine_mg50": dict(filename="Olivine_MgFeSiO4_Dorschner1995_mcfost.dat",
                         sha256="41e34acb22986f83aeb1113d31e73acee19707add259198036ca51baa25f3965",
                         density_g_cm3=3.71, wavelength_min_um=.2, wavelength_max_um=500.,
                         rows=109, label="Amorphous MgFeSiO4 glass"),
    "pyroxene_mg50": dict(filename="Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat",
                          sha256="05c1ebf12963ff3c9de9dd32e168135e62986e527e60343da5e432ff018b44a2",
                          density_g_cm3=3.2, wavelength_min_um=.2, wavelength_max_um=500.,
                          rows=109, label="Amorphous Mg0.5Fe0.5SiO3 glass"),
}
FIXED = {k: v for k, v in PILOT_FIXED.items() if k != "envelope_dust_mass_msun"}
BASELINE_EXPONENT = -1.5
RIN_AU, ROUT_AU = 1., 3000.
INCLINATIONS = (50., 60., 70.)
MASSES_MSUN = (.6e-4, 1.e-4, 1.5e-4, 2.25e-4)
MATCHED_REFERENCE_MASSES_MSUN = MASSES_MSUN[1:]
DUST = {**PILOT_DUST,
        "silicate": "Selected per envelope model; disk always retains Draine silicate/carbon",
        "silicate_materials": MATERIALS,
        "silicate_native_extrapolation_accepted": True,
        "silicate_coverage_policy": "Exact supplied tables; native simulator extrapolation outside tabulated coverage; no generated optical tails",
        "laboratory_silicate_temperature_k": None,
        "laboratory_silicate_density_policy": "Retain each supplied material's bulk density; equal dust mass is not equal optical depth",
        "tail_physical_accuracy_certified": False}


def power_integral(exponent, lower=RIN_AU, upper=ROUT_AU):
    """Integral of r**exponent dr over positive finite radial bounds."""
    values = (exponent, lower, upper)
    require(all(isinstance(v, (int, float)) and not isinstance(v, bool)
                and math.isfinite(v) for v in values), "Power-integral inputs must be finite numbers")
    require(0 < lower < upper, "Power-integral radii must obey 0 < Rin < Rout")
    q = exponent + 1.
    if abs(q) < 1e-12:
        return math.log(upper / lower)
    return lower**q * math.expm1(q * math.log(upper / lower)) / q


def column_mass_ratio(exponent, reference_exponent=BASELINE_EXPONENT):
    """M(p)/M(p0) at equal one-sided central radial column, same cavity.

    For rho=A*r**p, M=Omega*A*J(p+2) and Sigma=A*J(p).
    The common angular factor Omega cancels. This is an analytic continuous
    envelope normalization, not a certification of MCFOST's discrete grid.
    It is also not an equal-aperture-mass or equal emergent-flux condition.
    """
    return (power_integral(exponent + 2) / power_integral(exponent)
            * power_integral(reference_exponent) / power_integral(reference_exponent + 2))


def column_normalization(exponent, reference_mass):
    # Freeze 13 significant digits so platform libm rounding cannot alter
    # catalogue equality between preparation on macOS and checks on Linux.
    ratio = float(format(column_mass_ratio(exponent), ".13g"))
    return dict(method="analytic_continuous_radial_column", density_law="rho=A*r**p",
                exponent=exponent, reference_exponent=BASELINE_EXPONENT,
                rin_au=RIN_AU, rout_au=ROUT_AU,
                reference_mass_msun=reference_mass, actual_mass_msun=reference_mass * ratio,
                mass_ratio=ratio, cavity_half_opening_deg=FIXED["cavity_half_opening_deg"],
                mass_ratio_significant_digits=13,
                formula="M(p)=M(p0)*[J(p+2)/J(p)]/[J(p0+2)/J(p0)]; J(q)=integral_Rin^Rout r^q dr",
                angular_factor="Same conical cavity; common angular factor cancels",
                scope="Envelope-only one-sided central radial mass column outside the cavity; disk unchanged",
                discrete_mcfost_column_certified=False)


def catalogue():
    records = []

    def append(role, material, reference_mass, inclination, amax=.4, exponent=BASELINE_EXPONENT):
        normalization = column_normalization(exponent, reference_mass)
        parent = next((r["index"] for r in records
                       if r["role"] == "mass_size" and r["material_id"] == "draine"
                       and r["reference_mass_msun"] == reference_mass
                       and r["parameters"]["inclination_deg"] == inclination
                       and r["parameters"]["envelope_amax_um"] == .4), None)
        records.append(dict(index=len(records), role=role, material_id=material,
                            reference_mass_msun=reference_mass,
                            column_mass_ratio=normalization["mass_ratio"],
                            column_normalization=normalization,
                            parent_index=parent,
                            contrast_group_id=f"mass{reference_mass:.8g}_i{inclination:g}",
                            ice_mass_fraction=.04, dust_family="disjoint_DHS_H2O30K",
                            parameters={**FIXED, "envelope_dust_mass_msun": normalization["actual_mass_msun"],
                                        "inclination_deg": inclination, "envelope_amax_um": amax,
                                        "envelope_density_exponent": exponent,
                                        "envelope_silicate_file": MATERIALS[material]["filename"]}))

    for mass in MASSES_MSUN:
        for inclination in INCLINATIONS:
            for amax in (.4, 1.):
                append("mass_size", "draine", mass, inclination, amax)
    for material in ("olivine_mg50", "pyroxene_mg50"):
        for mass in MATCHED_REFERENCE_MASSES_MSUN:
            for inclination in INCLINATIONS:
                append("material", material, mass, inclination)
    for exponent in (-1.25, -1.75):
        for mass in MATCHED_REFERENCE_MASSES_MSUN:
            for inclination in INCLINATIONS:
                append("density_profile", "draine", mass, inclination, exponent=exponent)
    return records


def plan():
    return dict(schema_version=1, experiment_id=NAME, default_run_name=RUN_NAME,
                model_count=60, catalogue=catalogue(), dust_prescription=DUST, numerics=NUMERICS,
                arm_counts=dict(mass_size=24, material=18, density_profile=18),
                temperature_solves=60, probes_per_model=98, observable_count=19, image_requests=5880,
                cpus_per_task=64, max_concurrent_tasks=16, peak_cpus=1024,
                launch_ready=False, input_policy_ready=True,
                selection_rule="Predeclared matched contrasts; all 60 configurations analysed, no fit-selected subgrid",
                ice_coverage_policy="historical_native_extrapolation",
                pending="Prepare frozen bundle; configure cluster; validate exact DHS material initialization before submission",
                temperature_policy="Fresh equilibrium for each physical model; no previous temperature reuse",
                structure_normalization="Analytic central radial envelope column at fixed Rin/Rout and cavity; discrete column not certified",
                fixed_ice_fraction_is_ice_abundance_constraint=False,
                formal_likelihood=False, posterior_inference=False)


def validate_materials(dust_directory, provenance_path):
    """Require unchanged supplied material tables and their source metadata."""
    dust_directory, provenance_path = Path(dust_directory), Path(provenance_path)
    require(provenance_path.is_file() and digest(provenance_path) == MATERIAL_MANIFEST_SHA256,
            "Laboratory-silicate provenance manifest changed")
    provenance = json.loads(provenance_path.read_text())
    result = {}
    for material, expected in MATERIALS.items():
        path = dust_directory / expected["filename"]
        require(path.is_file() and not path.is_symlink() and digest(path) == expected["sha256"],
                f"Supplied silicate constants changed or absent: {expected['filename']}")
        if material != "draine":
            header, rows = optical_table(path)
            require(header == [expected["density_g_cm3"], 1500.]
                    and len(rows) == expected["rows"]
                    and rows[0][0] == expected["wavelength_min_um"]
                    and rows[-1][0] == expected["wavelength_max_um"],
                    f"Laboratory silicate metadata disagrees with its table: {material}")
            metadata = provenance["materials"][material]
            require(all(metadata[key] == expected[target] for key, target in (
                ("filename", "filename"), ("output_sha256", "sha256"),
                ("density_g_cm3", "density_g_cm3"), ("rows", "rows"),
                ("wavelength_min_um", "wavelength_min_um"), ("wavelength_max_um", "wavelength_max_um"))),
                f"Laboratory material provenance differs: {material}")
            require(metadata["n_k_modified"] is False and metadata["measurement_temperature_k"] is None,
                    "Laboratory silicate samples or measurement-temperature assumptions changed")
        result[material] = {**expected, "native_extrapolation_required": True,
                            "native_extrapolation_accepted": True, "n_k_modified": False,
                            "tail_physical_accuracy_certified": False}
    return dict(materials=result, provenance_sha256=digest(provenance_path),
                coverage_policy="exact_tables_native_extrapolation", policy_accepted=True)


def validate_frozen_campaign(bundle, experiment, manifest, index=None):
    """Gate this exact physical design and re-render selected frozen inputs."""
    import sys
    sys.path.insert(0, str(Path(bundle) / "code/src"))
    from mcfost_grid.physics import render_parameter
    bundle = Path(bundle)
    expected_catalogue = catalogue()
    require(experiment.get("catalogue") == expected_catalogue, "60-model silicate/structure design changed")
    require(experiment.get("dust_prescription") == DUST, "Fixed silicate/H2O material prescription changed")
    require(manifest["configuration"]["numerics"] == NUMERICS, "Campaign numerics changed")
    require(manifest.get("temperature_policy") == "fresh equilibrium for every physical model", "Fresh temperatures required")
    require(len(manifest["models"]) == 60, "Campaign requires exactly 60 models")
    for model, record in zip(manifest["models"], expected_catalogue):
        require(model["index"] == record["index"] and model["parameters"] == record["parameters"],
                "Campaign model parameters disagree with the declared design")
    receipt = validate_ice(bundle / "inputs/utils/Dust" / ICE_NAME,
                           bundle / "inputs/ice_input_provenance.json",
                           bundle / "inputs/H2O30K_laboratory_reference.dat")
    ice_provenance = json.loads((bundle / "inputs/ice_input_provenance.json").read_text())
    require(NAME in ice_provenance.get("accepted_experiment_ids", []),
            "Historical H2O extrapolation policy must explicitly cover this production campaign")
    require(receipt == experiment.get("ice_input_check"), "Supplied H2O input receipt changed")
    receipt = validate_materials(bundle / "inputs/utils/Dust", bundle / "inputs/laboratory_silicates_manifest.json")
    require(receipt == experiment.get("material_input_check"), "Supplied silicate input receipt changed")
    require(digest(bundle / "inputs/template.para") == TEMPLATE_SHA256,
            "Archived v02 thermal-grid/geometry template changed")
    template = (bundle / "inputs/template.para").read_text()
    if index is not None:
        require(type(index) is int and 0 <= index < 60, "Campaign model index out of range")
    selected = manifest["models"] if index is None else [manifest["models"][index]]
    for model in selected:
        directory = bundle / "models" / model["id"]
        expected = render_parameter(template, model["parameters"], NUMERICS, "temperature")
        require((directory / "temperature.para").read_text() == expected, "Campaign temperature parameters changed")
        for anchor in manifest["anchors"]:
            for stage in ("image", "coeval"):
                expected = render_parameter(template, model["parameters"], NUMERICS, stage, anchor["wavelength_um"])
                require((directory / "anchors" / anchor["id"] / f"{stage}.para").read_text() == expected,
                        f"Campaign {stage} parameters changed: {model['id']}/{anchor['id']}")
