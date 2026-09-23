"""Fixed four-model design and supplied-ice input checks (no simulations)."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

NAME = "silicate_size_pilot_v1"
RUN_NAME = "silicate_size_pilot_v1_512k"
ICE_NAME = "H2O_30K_Leiden_mcfost.dat"
LAB_SHA256 = "29b572f8907993870d0b15823e1182a9c4813b8cdb8f4fbdd40ef39cf28fb00c"
DRAINE_SHA256 = "3a6af751429f6f8dc948e60675b13be01ffcbed0662eec62934ec6d9c61dea80"
TEMPLATE_SHA256 = "28e8ce08b312e24af72cf45cb679b895e70d6ea7d2dd1e03b1cab2169f03b9ca"
FIXED = dict(envelope_dust_mass_msun=.000225, envelope_size_exponent=2.75,
             envelope_ice_mass_fraction=.04, cavity_half_opening_deg=17.5,
             distance_pc=140., stellar_temperature_k=4000., stellar_radius_rsun=2.5,
             stellar_mass_msun=.5, accretion_rate_msun_per_year=7.3011756855e-8)
NUMERICS = dict(photons_temperature=512000, photons_image=512000, photons_sed=512000,
                random_seed=43001, image_npix=6001, image_size_au=6000.,
                grains=50, grid_nr=100, grid_ntheta=70, grid_n_inner=20)
DUST = dict(family="disjoint_DHS", silicate="Draine_Si_sUV.dat", ice=ICE_NAME,
            silicate_mass_fraction=.96, ice_mass_fraction=.04, vmax=.1, porosity=0.,
            amin_um=.03, size_exponent=2.75, grain_bins=50, ice_amax_um=.4,
            silicate_amax_um=[.4, 1.], ice_density_g_cm3=.94,
            ice_laboratory_temperature_k=30., H2O30K_used=True,
            mass_convention="Total envelope dust mass includes both disjoint species",
            size_change_scope="Silicate only; H2O grain sizes stay fixed",
            temperature_dependent_ice_survival=False,
            draine_long_wavelength_policy="Retain the existing Draine table (ends at 1000 um) and simulator treatment outside it through the 3000-um thermal grid; no new silicate tail or tail-convergence claim.")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def catalogue():
    return [dict(index=i, role="silicate_size_pair", parent_index=None,
                 contrast_group_id=f"i{incl:g}", ice_mass_fraction=.04,
                 dust_family="disjoint_DHS_Draine_H2O30K",
                 parameters={**FIXED, "inclination_deg": incl, "envelope_amax_um": amax})
            for i, (incl, amax) in enumerate(( (50., .4), (50., 1.), (70., .4), (70., 1.) ))]


def plan():
    return dict(schema_version=1, experiment_id=NAME, default_run_name=RUN_NAME, model_count=4,
                catalogue=catalogue(), dust_prescription=DUST, numerics=NUMERICS,
                temperature_solves=4, probes_per_model=98, observable_count=19,
                image_requests=392, cpus_per_task=64, max_concurrent_tasks=4,
                launch_ready=False,
                input_policy_ready=True,
                ice_coverage_policy="historical_native_extrapolation",
                pending="Prepare frozen bundle and configure the shared cluster environment",
                temperature_policy="Fresh equilibrium for each physical model, no old temperature reuse")


def optical_table(path):
    """Read the MCFOST density/T header followed by wavelength, n, k rows."""
    rows = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        try:
            values = [float(v.replace("D", "e").replace("d", "e")) for v in line.split()]
        except ValueError as exc:
            raise ValueError(f"Invalid optical table numbers at line {number}: {path}") from exc
        require(all(math.isfinite(v) for v in values), f"Nonfinite optical constants: {path}:{number}")
        rows.append(values)
    require(len(rows) >= 4 and len(rows[0]) == 2 and all(len(r) == 3 for r in rows[1:]),
            "Expected density/sublimation-temperature header, then three columns: wavelength_um n k")
    header, samples = rows[0], rows[1:]
    require(header[0] > 0 and header[1] > 0, "Nonpositive optical table density/temperature")
    require(all(w > 0 and n > 0 and k >= 0 for w, n, k in samples), "Invalid wavelength, n or k")
    require(all(a[0] < b[0] for a, b in zip(samples, samples[1:])),
            "Supplied H2O wavelengths must be strictly increasing without duplicates")
    return header, samples


def validate_ice(path, provenance_path, laboratory_path):
    """Check the explicitly chosen coverage policy and measured material.

    The historical policy accepts only the exact archived table. A separately
    supplied full-range extension still has to meet its stricter input contract.
    Neither policy validates the physical accuracy of thermal tails.
    No extrapolation, clipping, splicing or optical-constant generation occurs.
    """
    path, provenance_path, laboratory_path = map(Path, (path, provenance_path, laboratory_path))
    require(path.is_file(), f"Missing H2O optical-constants table: {path}")
    header, samples = optical_table(path)
    require(provenance_path.is_file(), f"Missing H2O input provenance: {provenance_path}")
    provenance = json.loads(provenance_path.read_text())
    policy = provenance.get("coverage_policy", "supplied_full_range")
    require(policy in {"supplied_full_range", "historical_native_extrapolation"}, "Unknown H2O coverage policy")
    coverage_passed = samples[0][0] <= .1 and samples[-1][0] >= 3000.
    if policy == "supplied_full_range":
        require(coverage_passed,
                f"H2O coverage is {samples[0][0]:.12g}-{samples[-1][0]:.12g} um; the supplied_full_range policy requires 0.1-3000 um. The original table is supported only with explicit historical provenance.")
    else:
        require(digest(path) == LAB_SHA256, "Historical H2O policy requires the exact archived table; modified constants are not accepted")
        require(provenance.get("native_extrapolation_accepted_for_pilot") is True,
                "Historical H2O policy requires explicit acceptance of native extrapolation for this pilot")
        require(provenance.get("zero_k_policy") == "preserve_original_samples_without_floor",
                "Historical H2O policy must preserve original zero-k samples without a floor")
    require(header == [.94, 150.], "Retain the v02 H2O density/sublimation assumptions: 0.94 g/cm3, 150 K")
    require(provenance.get("schema_version") == 1 and provenance.get("table_sha256") == digest(path),
            "H2O provenance must have schema_version=1 and the exact table_sha256")
    require(digest(laboratory_path) == LAB_SHA256 and provenance.get("laboratory_sha256") == LAB_SHA256,
            "H2O laboratory reference hash differs")
    require(provenance.get("laboratory_temperature_k") == 30., "H2O provenance must identify the 30 K laboratory region")
    require(provenance.get("measured_region_treatment") == "preserved_original_samples",
            "This pilot requires preserved original H2O30K samples; a resampled/reprocessed table needs a reviewed input contract")
    for field in ("sources", "extension_method", "join_description", "tail_temperature_assumptions"):
        value = provenance.get(field)
        items = value if isinstance(value, list) else [value]
        require(items and all(isinstance(v, str) and v.strip() and "REQUIRED" not in v and "TODO" not in v for v in items),
                f"Document H2O {field}; placeholder metadata is not launchable")
    # Keep the archived measured spectrum fixed. An extension can add samples
    # outside it, but cannot silently replace or smooth its spectral features.
    _, lab = optical_table(laboratory_path)
    inside = [r for r in samples if lab[0][0] <= r[0] <= lab[-1][0]]
    require(len(inside) == len(lab) and all(
        math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-13)
        for row, old in zip(inside, lab) for a, b in zip(row, old)),
        "Extended H2O table changes the archived measured region; supply the preserving extension or review the changed material first")
    return dict(table_sha256=digest(path), provenance_sha256=digest(provenance_path),
                laboratory_sha256=LAB_SHA256, wavelength_min_um=samples[0][0],
                wavelength_max_um=samples[-1][0], rows=len(samples), density_g_cm3=header[0],
                measured_region_preserved=True, coverage_passed=coverage_passed,
                coverage_policy=policy, policy_accepted=True,
                native_extrapolation_required=not coverage_passed,
                zero_k_rows=sum(row[2] == 0 for row in samples),
                zero_k_floor_applied=False,
                tail_physical_accuracy_certified=False)


def validate_frozen_pilot(bundle, experiment, manifest, index=None):
    """Gate shared supplied inputs and re-render the exact selected model files."""
    import sys
    sys.path.insert(0, str(Path(bundle) / "code/src"))
    from mcfost_grid.physics import render_parameter
    bundle = Path(bundle)
    require(experiment.get("catalogue") == catalogue(), "Four-model silicate design changed")
    require(experiment.get("dust_prescription") == DUST, "Fixed v02 dust prescription changed")
    require(manifest["configuration"]["numerics"] == NUMERICS, "Silicate pilot numerics changed")
    require(manifest.get("temperature_policy") == "fresh equilibrium for every physical model", "Fresh temperatures required")
    receipt = validate_ice(bundle / "inputs/utils/Dust" / ICE_NAME,
                           bundle / "inputs/ice_input_provenance.json",
                           bundle / "inputs/H2O30K_laboratory_reference.dat")
    require(receipt == experiment.get("ice_input_check"), "Supplied H2O input receipt changed")
    require(digest(bundle / "inputs/utils/Dust/Draine_Si_sUV.dat") == DRAINE_SHA256,
            "Existing Draine constants changed")
    require(digest(bundle / "inputs/template.para") == TEMPLATE_SHA256,
            "Archived v02 structure or 0.1-3000 um thermal-grid template changed")
    template = (bundle / "inputs/template.para").read_text()
    selected = manifest["models"] if index is None else [manifest["models"][index]]
    for model in selected:
        directory = bundle / "models" / model["id"]
        expected = render_parameter(template, model["parameters"], NUMERICS, "temperature")
        require((directory / "temperature.para").read_text() == expected, "Pilot temperature parameters changed")
        for anchor in manifest["anchors"]:
            for stage in ("image", "coeval"):
                expected = render_parameter(template, model["parameters"], NUMERICS, stage, anchor["wavelength_um"])
                require((directory / "anchors" / anchor["id"] / f"{stage}.para").read_text() == expected,
                        f"Pilot {stage} parameters changed: {model['id']}/{anchor['id']}")
