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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent/"src" if (Path(__file__).parent/"src").is_dir() else ROOT/"src"))
from mcfost_grid.configuration import atomic_json, canonical_hash, now, sha256
from mcfost_grid.physics import render_parameter
from silicate_search_v2_design import (NAME, RUN_NAME, SPECIES, DUST, SUCCESS, NUMERICS,
    TEMPLATE_SHA256, catalogue, declared_cases, require)

COPIES = {
    "build_silicate_search_v2.py": "build_silicate_search_v2.py",
    "silicate_search_v2_design.py": "silicate_search_v2_design.py",
    "silicate_search_v2_stage0.py": "silicate_search_v2_stage0.py",
    "silicate_search_v2_task.py": "silicate_search_v2_task.py",
    "configure_silicate_search_v2.py": "configure_silicate_search_v2.py",
    "analyze_silicate_search_v2.py": "analyze_silicate_search_v2.py",
    "extinction_ice_production_task.py": "production_task.py",
    "analyze_extinction_ice_production.py": "analyze_extinction_ice_production.py",
    "configure_extinction_ice_numerics_v2.py": "configure_extinction_ice_numerics_v2.py",
    "analyze_silicate_structure_production.py": "analyze_silicate_structure_production.py",
    "silicate_structure_design.py": "silicate_structure_design.py",
    "silicate_pilot_design.py": "silicate_pilot_design.py",
}


def validate_bootstrap(bundle):
    bundle = Path(bundle).resolve()
    path = bundle/"bootstrap.json"
    require(sha256(path) == (bundle/"bootstrap.sha256").read_text().split()[0], "Bootstrap manifest changed")
    plan = json.loads(path.read_text())
    require(plan["experiment_id"] == NAME and plan["cases"] == declared_cases(), "V2 declaration changed")
    require(plan["model_count"] == 34 and plan["comparison_count"] == 48, "Wrong v2 size")
    require(plan["dust_prescription"] == DUST and plan["numerics"] == NUMERICS
            and plan["success_criteria"] == SUCCESS, "Bootstrap physical or numerical policy changed")
    for relative, checksum in plan["input_hashes"].items():
        p = bundle/relative
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts
                and p.resolve().is_relative_to(bundle) and not p.is_symlink()
                and p.is_file() and sha256(p) == checksum, f"Bootstrap input changed: {relative}")
    for species in SPECIES.values():
        require(sha256(bundle/"inputs/utils/Dust"/species["filename"]) == species["sha256"], "Material table changed")
    return plan


def prepare_bootstrap(destination):
    destination = Path(destination).resolve()
    require(not destination.exists(), f"Refusing to replace existing package: {destination}")
    previous = ROOT/"runs/silicate_structure_production_v1_512k"
    require(sha256(previous/"inputs/template.para") == TEMPLATE_SHA256, "Reference template changed")
    for species in SPECIES.values():
        require(sha256(previous/"inputs/utils/Dust"/species["filename"]) == species["sha256"], "Reference dust changed")
    for name in COPIES:
        require((ROOT/"scripts"/name).is_file(), f"Missing package code {name}")
    for name in ("SILICATE_SEARCH_V2.md", "SILICATE_SEARCH_V2_RUN.md"):
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
    provenance["policy_scope"] = "silicate_search_v2; unchanged supplied H2O30K and native out-of-range treatment explicitly retained by the v2 design"
    provenance["user_direction"] = "Prepare docs/SILICATE_SEARCH_V2.md with unchanged optical constants and matched-extinction DHS composition grid"
    atomic_json(inputs/"ice_input_provenance.json", provenance)
    manifest = json.loads((previous/"manifest.json").read_text())
    anchors = manifest["anchors"]
    require(len(anchors) == 98, "Expected inherited 98 image probes")
    atomic_json(inputs/"production_anchors.json", anchors)
    waves = sorted(set([a["wavelength_um"] for a in anchors]+[2.2, 9.7, 18.]))
    require(len(waves) == 99, "Expected 99 distinct Stage-0 wavelengths")
    atomic_json(inputs/"stage0_wavelengths.json", dict(wavelengths_um=waves,
        production_wavelengths_um=[a["wavelength_um"] for a in anchors],
        extra_normalization_wavelengths_um=[2.2, 9.7, 18.], count_header=False))
    for source, target in COPIES.items():
        shutil.copy2(ROOT/"scripts"/source, code/target)
    shutil.copytree(ROOT/"src/mcfost_grid", code/"src/mcfost_grid", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(ROOT/"workflow.py", code/"workflow.py")
    shutil.copy2(ROOT/"requirements-cluster-recovery.txt", destination/"requirements-cluster.txt")
    shutil.copy2(ROOT/"docs/SILICATE_SEARCH_V2.md", destination/"DESIGN.md")
    shutil.copy2(ROOT/"docs/SILICATE_SEARCH_V2_RUN.md", destination/"README.md")
    (destination/".gitignore").write_text("__pycache__/\n*.pyc\n.DS_Store\nlogs/\nresults/\nstage0/\nproduction/\nbuild_attempts/\nsubmissions/\n.mplconfig/\nmachine*.cluster.json\ncluster_environment.json\njob_*.sh\nsubmit.sh\n*.lock\n")
    plan = dict(schema_version=1, experiment_id=NAME, created_utc=now(),
        model_count=34, comparison_count=48, cases=declared_cases(),
        matched_masses_available=False, stage0_required=True, production_parameters_rendered=False,
        dust_prescription=DUST, numerics=NUMERICS, success_criteria=SUCCESS,
        selection_rule="All compositions and comparisons fixed before opacity/RT outputs; Stage 0 sets only mass normalization",
        parent_manifest_sha256=sha256(previous/"manifest.json"),
        input_hashes={str(p.relative_to(destination)):sha256(p) for p in sorted(destination.rglob("*")) if p.is_file()})
    atomic_json(destination/"bootstrap.json", plan)
    (destination/"bootstrap.sha256").write_text(sha256(destination/"bootstrap.json")+"  bootstrap.json\n")
    atomic_json(destination/"machine.template.json", dict(schema_version=1, mcfost_executable="mcfost",
        mcfost_utils="$MCFOST_UTILS", backend="image_method2", threads=64, timeout_seconds=14400,
        max_memory_gb=112, slurm=dict(time="24:00:00", memory="160G", max_parallel=16,
        python="/SET/BY/CONFIGURE", setup_lines=[], analysis_cpus=2, analysis_memory="8G", analysis_time="02:00:00")))
    validate_bootstrap(destination)
    return dict(bundle=str(destination), model_count=34, comparison_count=48,
        stage0_required=True, mass_estimates_used=False, jobs_submitted=False)


def materialize(bundle):
    from silicate_search_v2_stage0 import validate_stage0_receipt
    from analyze_silicate_search_v2 import predeclared_contrasts
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
            return dict(status="cached", production=str(production), models=34)
        opacities = {row["id"]: json.loads((bundle/row["opacity_file"]).read_text()) for row in receipt["species"]}
        records = catalogue(opacities)
        require(len(predeclared_contrasts(records)) == 48, "Wrong comparison count")
        models = [dict(index=row["index"], id=f"m{row['index']:04d}_{canonical_hash(row['parameters'])[:8]}",
                       parameters=row["parameters"]) for row in records]
        require(len({m["id"] for m in models}) == 34, "Duplicate v2 model")
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
            (directory/"temperature.para").write_text(render_parameter(template, model["parameters"], NUMERICS, "temperature"))
            for anchor in anchors:
                folder = directory/"anchors"/anchor["id"]
                folder.mkdir(parents=True)
                for stage in ("image", "coeval"):
                    (folder/(stage+".para")).write_text(render_parameter(template, model["parameters"], NUMERICS, stage, anchor["wavelength_um"]))
                (folder/"wavelength.lambda").write_text(f"{anchor['wavelength_um']:.17g}\n")
        config = dict(schema_version=1, run_name=RUN_NAME, output_dir="..", numerical_only=True,
            smoke_test=False, max_models=34, template="template.para", models=[m["parameters"] for m in models],
            numerics=NUMERICS, dust_files={p.name:"utils/Dust/"+p.name for p in (staging/"inputs/utils/Dust").iterdir()},
            observations=dict(anchors_file="production_probes.ecsv", spectrum_file="observed_spectrum.ecsv",
                aperture_radius_arcsec=1., target_distance_pc=147., aperture_subpixels=64, quality_policy="aperture_v2"))
        atomic_json(staging/"inputs/configuration.json", config)
        experiment = dict(schema_version=1, experiment_id=NAME, catalogue=records,
            tasks=[dict(index=m["index"], model_id=m["id"]) for m in models],
            model_count=34, temperature_solves=34, image_requests=3332, probes_per_model=98, observable_count=19,
            dust_prescription=DUST, numerics=NUMERICS, success_criteria=SUCCESS,
            stage0_receipt_sha256=receipt_hash, stage0_opacity=opacities,
            opacity_matching=[dict(index=r["index"], **r["opacity_matching"]) for r in records],
            predeclared_contrasts=predeclared_contrasts(records), bootstrap_sha256=sha256(bundle/"bootstrap.json"),
            posterior_inference=False, physical_cause_uniquely_identified=False,
            selection_rule="All 34 declared cases; no selection by fit or by Stage-0 S ratios",
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
    return dict(status="prepared", production=str(production), models=34, comparisons=48, fresh_temperatures=34)


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
