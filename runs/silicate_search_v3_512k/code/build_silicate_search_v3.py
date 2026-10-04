#!/usr/bin/env python3
"""Freeze a portable Stage-0 package, then materialize measured-opacity models."""
from __future__ import annotations

import argparse
import fcntl
import json
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent/"src" if (Path(__file__).parent/"src").is_dir() else ROOT/"src"))
from mcfost_grid.configuration import atomic_json, canonical_hash, now, sha256
from mcfost_grid.physics import render_parameter
from silicate_search_v3_design import (NAME, RUN_NAME, SPECIES, DUST, SUCCESS, NUMERICS,
    TEMPLATE_SHA256, catalogue, declared_cases, model_numerics, require)

COPIES = {
    "build_silicate_search_v3.py": "build_silicate_search_v3.py",
    "silicate_search_v3_design.py": "silicate_search_v3_design.py",
    "silicate_search_v2_design.py": "silicate_search_v2_design.py",
    "silicate_search_v3_stage0.py": "silicate_search_v3_stage0.py",
    "silicate_search_v3_task.py": "silicate_search_v3_task.py",
    "configure_silicate_search_v3.py": "configure_silicate_search_v3.py",
    "analyze_silicate_search_v3.py": "analyze_silicate_search_v3.py",
    "extinction_ice_production_task.py": "production_task.py",
    "analyze_extinction_ice_production.py": "analyze_extinction_ice_production.py",
    "configure_extinction_ice_numerics_v2.py": "configure_extinction_ice_numerics_v2.py",
    "analyze_silicate_structure_production.py": "analyze_silicate_structure_production.py",
    "silicate_structure_design.py": "silicate_structure_design.py",
    "silicate_pilot_design.py": "silicate_pilot_design.py",
}


def prepare_reference(inputs):
    """Freeze the measured evidence; planning masses never become run inputs."""
    import silicate_search_v3_plan as planning
    previous = ROOT/"runs/silicate_search_v2_512k"
    report = json.loads((previous/"results/summary.json").read_text())
    campaign = json.loads((previous/"results/campaign_summary.json").read_text())
    require(sha256(previous/"results/summary.json") == campaign["production_summary_sha256"], "V2 summary changed")
    require(sha256(inputs/"production_observation_contract.json") == report["observation_contract_sha256"], "V2 observation contract changed")
    plan = planning.plan()
    require(plan["configurations"] == [(-1.5, 1.), (-1.25, .86), (-1., .83), (-1., .66)], "V3 column levels no longer reproduce the design")
    require(plan["ice_column_target_crossings"] == {"draine": 7, "pyroxene_mg50": 2}, "V3 bracketed ice crossings changed")
    old_opacities = report["campaign_context"]["stage0_opacity"]
    # Use the exact recorded Stage 0 values as the reproducibility reference.
    # Least-squares recovery gives -1.8e-13 for the truly zero ice 9.7 opacity;
    # that roundoff is retained only as a recovery diagnostic, never as opacity.
    screens = {}
    for species in SPECIES:
        saved = old_opacities[species]
        screens[species] = {}
        for key, field, wavelength in (("ext_2p2", "kappa_ext_cm2_g", 2.2),
                ("abs_9p7", "kappa_abs_cm2_g", 9.7), ("abs_18", "kappa_abs_cm2_g", 18.)):
            indices = np.flatnonzero(np.isclose(saved["wavelength_um"], wavelength, rtol=2e-6, atol=0))
            require(len(indices) == 1, "Ambiguous v2 reference wavelength")
            value = float(saved[field][indices[0]])
            recovered = plan["species_opacities_cm2_g"][species][key]
            require(value >= 0 and np.isclose(value, recovered, rtol=1e-10, atol=1e-8), "Recovered opacities disagree with v2 Stage 0")
            screens[species][key] = value
    source_paths = {
        "v1_band_predictions.csv": ROOT/"runs/silicate_structure_production_v1_512k/results/band_predictions.csv",
        "v1_ranking.csv": ROOT/"runs/silicate_structure_production_v1_512k/results/ranking.csv",
        "v2_opacity_matched_catalogue.csv": previous/"results/opacity_matched_catalogue.csv",
        "ice_column_crossings.csv": ROOT/"validation/silicate_search_v2_512k/ice_column_crossings.csv",
        "v2_review_source.py": ROOT/"scripts/review_silicate_search_v2.py",
        "v3_planning_source.py": ROOT/"scripts/silicate_search_v3_plan.py",
    }
    provenance_dir = inputs/"provenance"
    provenance_dir.mkdir()
    for name, source in source_paths.items():
        shutil.copy2(source, provenance_dir/name)
    provenance = dict(source_hashes={str(p.relative_to(ROOT)): sha256(p) for p in source_paths.values()},
        frozen_source_hashes={"inputs/provenance/"+name: sha256(path) for name, path in source_paths.items()},
        v2_summary_sha256=sha256(previous/"results/summary.json"),
        v2_campaign_summary_sha256=sha256(previous/"results/campaign_summary.json"),
        v2_reported_production_manifest_sha256=report["manifest_sha256"],
        v2_raw_stage0_revalidated=False,
        column_factor_derivation=plan["column_factor_derivation"], configurations=plan["configurations"],
        column_factor_method="v1 log-linear band interpolation on a 2001-point log-column grid; round optimum/reference to .86; bracket .86 squared by +/-0.05 dex then round to .83/.66",
        ice_column_target_crossings=plan["ice_column_target_crossings"],
        ice_column_crossings_sha256=sha256(source_paths["ice_column_crossings.csv"]),
        ice_target_method="Per-material median column-equivalent ice mass of bracketed v2 raw-ice zero crossings; no extrapolated crossings",
        ice_fraction_formula="f_ice = T*((1-f_C)*k_sil+f_C*k_C)/(M_ref*c*k_ref-T*(k_ice-k_sil))",
        opacity_recovery=plan["species_opacities_cm2_g"], opacity_recovery_max_relative_fit_residual=plan["max_relative_fit_residual"],
        miri_reduction_status=dict(corrected_reduction_available=False, collaborators_level_change_known=False,
            user_update="No corrected reduction or update yet", recorded_date="2026-10-04",
            spectrum_policy="Retain the existing frozen observations; a corrected 12-20 um level shift greater than 0.1 dex requires new residuals/projections and a new frozen package before launch"),
        interpretation="User-provided v3 design is tested conditionally; introductory projections are not simulated outcomes")
    reference = dict(schema_version=1, species_opacities_cm2_g=screens,
        ice_column_targets_msun=plan["ice_column_targets_msun"],
        observation_contract_sha256=report["observation_contract_sha256"],
        bridge_models=[m for m in report["models"] if m["model_index"] in (15, 17)],
        design_provenance=provenance)
    require([m["model_index"] for m in reference["bridge_models"]] == [15, 17], "Missing v2 bridge models")
    atomic_json(inputs/"v2_reference.json", reference)
    return reference


def validate_bootstrap(bundle):
    bundle = Path(bundle).resolve()
    path = bundle/"bootstrap.json"
    require(sha256(path) == (bundle/"bootstrap.sha256").read_text().split()[0], "Bootstrap manifest changed")
    plan = json.loads(path.read_text())
    require(plan["experiment_id"] == NAME and plan["cases"] == declared_cases(), "V3 declaration changed")
    require(plan["model_count"] == 26 and plan["comparison_count"] == 54, "Wrong v3 size")
    require(plan["dust_prescription"] == DUST and plan["numerics"] == NUMERICS
            and plan["success_criteria"] == SUCCESS, "Bootstrap physical or numerical policy changed")
    for relative, checksum in plan["input_hashes"].items():
        p = bundle/relative
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts
                and p.resolve().is_relative_to(bundle) and not p.is_symlink()
                and p.is_file() and sha256(p) == checksum, f"Bootstrap input changed: {relative}")
    for species in SPECIES.values():
        require(sha256(bundle/"inputs/utils/Dust"/species["filename"]) == species["sha256"], "Material table changed")
    reference = json.loads((bundle/"inputs/v2_reference.json").read_text())
    require(plan["ice_column_targets_msun"] == reference["ice_column_targets_msun"]
            and plan["design_provenance"] == reference["design_provenance"], "V3 source provenance changed")
    return plan


def prepare_bootstrap(destination):
    destination = Path(destination).resolve()
    require(not destination.exists(), f"Refusing to replace existing package: {destination}")
    previous = ROOT/"runs/silicate_search_v2_512k"
    require(sha256(previous/"inputs/template.para") == TEMPLATE_SHA256, "Reference template changed")
    for species in SPECIES.values():
        require(sha256(previous/"inputs/utils/Dust"/species["filename"]) == species["sha256"], "Reference dust changed")
    for name in COPIES:
        require((ROOT/"scripts"/name).is_file(), f"Missing package code {name}")
    for name in ("SILICATE_SEARCH_V3.md", "SILICATE_SEARCH_V3_RUN.md"):
        require((ROOT/"docs"/name).is_file(), f"Missing package documentation {name}")
    destination.mkdir(parents=True)
    inputs, code = destination/"inputs", destination/"code"
    inputs.mkdir(); code.mkdir()
    for name in ("template.para", "observed_spectrum.ecsv", "production_observation_contract.json",
                 "production_probes.ecsv", "H2O30K_laboratory_reference.dat", "laboratory_silicates_manifest.json"):
        shutil.copy2(previous/"inputs"/name, inputs/name)
    shutil.copytree(previous/"inputs/utils", inputs/"utils")
    provenance = json.loads((previous/"inputs/ice_input_provenance.json").read_text())
    provenance["accepted_experiment_ids"] = sorted(set(provenance.get("accepted_experiment_ids", [])) | {NAME})
    provenance["policy_scope"] = "silicate_search_v3; unchanged supplied H2O30K and native out-of-range treatment explicitly retained"
    provenance["user_direction"] = "Prepare docs/SILICATE_SEARCH_V3.md with shallower profiles, column factors and derived ice fractions"
    atomic_json(inputs/"ice_input_provenance.json", provenance)
    anchors = json.loads((previous/"inputs/production_anchors.json").read_text())
    require(len(anchors) == 98, "Expected inherited 98 image probes")
    atomic_json(inputs/"production_anchors.json", anchors)
    waves = sorted(set([a["wavelength_um"] for a in anchors]+[2.2, 9.7, 18.]))
    require(len(waves) == 99, "Expected 99 distinct Stage-0 wavelengths")
    atomic_json(inputs/"stage0_wavelengths.json", dict(wavelengths_um=waves,
        production_wavelengths_um=[a["wavelength_um"] for a in anchors],
        extra_normalization_wavelengths_um=[2.2, 9.7, 18.], count_header=False))
    reference = prepare_reference(inputs)
    for source, target in COPIES.items():
        shutil.copy2(ROOT/"scripts"/source, code/target)
    shutil.copytree(ROOT/"src/mcfost_grid", code/"src/mcfost_grid", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(ROOT/"workflow.py", code/"workflow.py")
    shutil.copy2(ROOT/"requirements-cluster-recovery.txt", destination/"requirements-cluster.txt")
    shutil.copy2(ROOT/"docs/SILICATE_SEARCH_V3.md", destination/"DESIGN.md")
    shutil.copy2(ROOT/"docs/SILICATE_SEARCH_V3_RUN.md", destination/"README.md")
    (destination/".gitignore").write_text("__pycache__/\n*.pyc\n.DS_Store\nlogs/\nresults/\nstage0/\nproduction/\nbuild_attempts/\nsubmissions/\n.mplconfig/\nmachine*.cluster.json\ncluster_environment.json\njob_*.sh\nsubmit.sh\n*.lock\n")
    plan = dict(schema_version=1, experiment_id=NAME, created_utc=now(),
        model_count=26, comparison_count=54, cases=declared_cases(),
        physical_models=24, seed_replicas=2, derived_ice_fractions_available=False,
        matched_masses_available=False, stage0_required=True, production_parameters_rendered=False,
        dust_prescription=DUST, numerics=NUMERICS, success_criteria=SUCCESS,
        ice_column_targets_msun=reference["ice_column_targets_msun"], design_provenance=reference["design_provenance"],
        selection_rule="All configurations fixed before RT outputs; Stage 0 sets mass and ice fraction from frozen column targets",
        parent_bootstrap_sha256=sha256(previous/"bootstrap.json"),
        input_hashes={str(p.relative_to(destination)):sha256(p) for p in sorted(destination.rglob("*")) if p.is_file()})
    atomic_json(destination/"bootstrap.json", plan)
    (destination/"bootstrap.sha256").write_text(sha256(destination/"bootstrap.json")+"  bootstrap.json\n")
    atomic_json(destination/"machine.template.json", dict(schema_version=1, mcfost_executable="mcfost",
        mcfost_utils="$MCFOST_UTILS", backend="image_method2", threads=64, timeout_seconds=14400,
        max_memory_gb=112, slurm=dict(time="24:00:00", memory="160G", max_parallel=16,
        python="/SET/BY/CONFIGURE", setup_lines=[], analysis_cpus=2, analysis_memory="8G", analysis_time="02:00:00")))
    validate_bootstrap(destination)
    return dict(bundle=str(destination), model_count=26, comparison_count=54,
        stage0_required=True, mass_estimates_used=False, jobs_submitted=False)


def materialize(bundle):
    from silicate_search_v3_stage0 import validate_stage0_receipt
    from analyze_silicate_search_v3 import predeclared_contrasts
    if (Path(__file__).parent/"production_task.py").is_file():
        from production_task import validate_package
    else:
        from extinction_ice_production_task import validate_package
    bundle = Path(bundle).resolve()
    validate_bootstrap(bundle)
    with (bundle/".materialize.lock").open("a+") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        receipt = validate_stage0_receipt(bundle)
        receipt_hash = sha256(bundle/"stage0/receipt.json")
        production = bundle/"production"
        if production.exists():
            existing = validate_package(production)
            require(existing["stage0_receipt_sha256"] == receipt_hash,
                    "Existing production was frozen against different Stage-0 outputs")
            return dict(status="cached", production=str(production), models=26)
        opacities = {row["id"]: json.loads((bundle/row["opacity_file"]).read_text()) for row in receipt["species"]}
        reference = json.loads((bundle/"inputs/v2_reference.json").read_text())
        records = catalogue(opacities, reference["ice_column_targets_msun"])
        require(len(predeclared_contrasts(records)) == 54, "Wrong comparison count")
        models = [dict(index=row["index"],
                       id=f"m{row['index']:04d}_{canonical_hash(dict(parameters=row['parameters'], random_seed=row['random_seed']))[:8]}",
                       parameters=row["parameters"], numerics=model_numerics(row)) for row in records]
        require(len({m["id"] for m in models}) == 26, "Duplicate v3 model")
        attempts = bundle/"build_attempts"
        attempts.mkdir(exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix="build-", dir=attempts))/"production"
        staging.mkdir()
        shutil.copytree(bundle/"inputs", staging/"inputs")
        shutil.copytree(bundle/"code", staging/"code", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        shutil.copy2(bundle/"README.md", staging/"README.md")
        shutil.copy2(bundle/"stage0/receipt.json", staging/"inputs/stage0_receipt.json")
        atomic_json(staging/"inputs/stage0_opacities.json", opacities)
        anchors = json.loads((bundle/"inputs/production_anchors.json").read_text())
        template = (bundle/"inputs/template.para").read_text()
        for model in models:
            directory = staging/"models"/model["id"]
            directory.mkdir(parents=True)
            (directory/"temperature.para").write_text(render_parameter(template, model["parameters"], model["numerics"], "temperature"))
            for anchor in anchors:
                folder = directory/"anchors"/anchor["id"]
                folder.mkdir(parents=True)
                for stage in ("image", "coeval"):
                    (folder/(stage+".para")).write_text(render_parameter(template, model["parameters"], model["numerics"], stage, anchor["wavelength_um"]))
                (folder/"wavelength.lambda").write_text(f"{anchor['wavelength_um']:.17g}\n")
        config = dict(schema_version=1, run_name=RUN_NAME, output_dir="..", numerical_only=True,
            smoke_test=False, max_models=26, template="template.para", models=[m["parameters"] for m in models],
            numerics=NUMERICS, dust_files={p.name:"utils/Dust/"+p.name for p in (staging/"inputs/utils/Dust").iterdir()},
            observations=dict(anchors_file="production_probes.ecsv", spectrum_file="observed_spectrum.ecsv",
                aperture_radius_arcsec=1., target_distance_pc=147., aperture_subpixels=64, quality_policy="aperture_v2"))
        atomic_json(staging/"inputs/configuration.json", config)
        experiment = dict(schema_version=1, experiment_id=NAME, catalogue=records,
            tasks=[dict(index=m["index"], model_id=m["id"]) for m in models],
            model_count=26, temperature_solves=26, image_requests=2548, probes_per_model=98, observable_count=19,
            physical_models=24, seed_replicas=2, fresh_temperatures_for_replicas=True,
            dust_prescription=DUST, numerics=NUMERICS, success_criteria=SUCCESS,
            ice_column_targets_msun=reference["ice_column_targets_msun"],
            design_provenance=reference["design_provenance"],
            observation_contract_sha256=reference["observation_contract_sha256"],
            stage0_receipt_sha256=receipt_hash, stage0_opacity=opacities,
            opacity_matching=[dict(index=r["index"], **r["opacity_matching"]) for r in records],
            predeclared_contrasts=predeclared_contrasts(records), bootstrap_sha256=sha256(bundle/"bootstrap.json"),
            posterior_inference=False, physical_cause_uniquely_identified=False,
            selection_rule="All 26 declared cases; no selection by fit or by Stage-0 S ratios",
            resources=dict(cpus_per_task=64, max_concurrent_tasks=16, memory="160G"))
        for name in ("experiment.json", "production_experiment.json"):
            atomic_json(staging/name, experiment)
        (staging/"experiment.sha256").write_text(sha256(staging/"experiment.json")+"  experiment.json\n")
        manifest = dict(schema_version=1, run_id=RUN_NAME, created_utc=now(), configuration=config,
            models=models, anchors=anchors, measurement=dict(aperture_radius_arcsec=1., target_distance_pc=147., aperture_subpixels=64, quality_policy="aperture_v2"),
            observation_spectrum="inputs/observed_spectrum.ecsv", model_directory="models/{id}",
            estimator="signed_direct_method2_image_total_I", formal_likelihood=False,
            temperature_policy="fresh equilibrium for every physical model",
            input_hashes={str(p.relative_to(staging)):sha256(p) for p in sorted(staging.rglob("*")) if p.is_file()})
        atomic_json(staging/"manifest.json", manifest)
        (staging/"manifest.sha256").write_text(sha256(staging/"manifest.json")+"  manifest.json\n")
        validate_package(staging)
        require(sha256(bundle/"stage0/receipt.json") == receipt_hash, "Stage 0 changed during build")
        staging.rename(production)
    return dict(status="prepared", production=str(production), models=26, comparisons=54, fresh_temperatures=26)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--materialize", type=Path)
    mode.add_argument("--validate", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT/"runs"/RUN_NAME)
    args = parser.parse_args()
    if args.materialize:
        result = materialize(args.materialize)
    elif args.validate:
        plan = validate_bootstrap(args.validate)
        result = dict(bootstrap_valid=True, models=plan["model_count"], stage0_required=True)
    else:
        result = prepare_bootstrap(args.output)
    print(json.dumps(result, indent=2, allow_nan=False))
