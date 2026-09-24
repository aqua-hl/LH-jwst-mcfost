#!/usr/bin/env python3
"""Freeze the approved 60-model silicate/structure search without submitting jobs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.configuration import atomic_json, prepare_run, sha256
from production_observables import build_contract
from silicate_structure_design import (NAME, RUN_NAME, ICE_NAME, DUST, NUMERICS,
                                      TEMPLATE_SHA256, MATERIALS, catalogue, plan,
                                      validate_ice, validate_materials)


def build_package(destination):
    from astropy.table import Table
    from analyze_silicate_structure_production import predeclared_contrasts
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"Refusing to replace an existing run: {destination}")
    template = ROOT / "reference/parameters/ice_v02_dust.para"
    if sha256(template) != TEMPLATE_SHA256:
        raise ValueError("Archived v02 template changed; review the design before preparation")
    ice = ROOT / "reference/dust" / ICE_NAME
    provenance = ROOT / "silicate/constants/H2O_30K_production.provenance.json"
    material_manifest = ROOT / "laboratory_silicates_v1/manifest.json"
    material_dir = ROOT / "laboratory_silicates_v1/Dust"
    ice_check = validate_ice(ice, provenance, ice)
    material_check = validate_materials(material_dir, material_manifest)
    records, contract = catalogue(), build_contract(ROOT)
    copies = {
        "extinction_ice_production_task.py": "production_task.py",
        "analyze_extinction_ice_production.py": "analyze_extinction_ice_production.py",
        "analyze_silicate_structure_production.py": "analyze_silicate_structure_production.py",
        "silicate_structure_design.py": "silicate_structure_design.py",
        "silicate_pilot_design.py": "silicate_pilot_design.py",
        "configure_extinction_ice_numerics_v2.py": "configure_cluster.py",
        "check_silicate_structure_materials.py": "check_materials.py",
        "production_observables.py": "production_observables.py",
    }
    doc = ROOT / "docs/SILICATE_STRUCTURE_PRODUCTION_V1.md"
    sources = [template, ice, provenance, material_manifest, doc, Path(__file__).resolve(),
               ROOT / "requirements-cluster-recovery.txt", ROOT / "workflow.py",
               ROOT / "reference/observations/continuum_sed_R100.ecsv",
               ROOT / "reference/observations/tmc1a_sed_unstitched.ecsv",
               ROOT / "reference/observations/continuum_fit_mask_v1_manifest.json"]
    sources += [material_dir / m["filename"] for m in MATERIALS.values()]
    sources += [ROOT / "scripts" / p for p in copies]
    sources += sorted((ROOT / "src/mcfost_grid").glob("*.py"))
    before = {str(p.relative_to(ROOT)): sha256(p) for p in sources}
    dust_files = {m["filename"]: str(material_dir / m["filename"]) for m in MATERIALS.values()}
    dust_files[ICE_NAME] = str(ice)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="silicate-structure-inputs-", dir=destination.parent) as tmp:
        tmp = Path(tmp)
        probes = tmp / "probes.ecsv"
        Table(rows=contract["probes"]).write(probes, format="ascii.ecsv")
        config = dict(schema_version=1, run_name=destination.name, output_dir=str(destination.parent),
                      template=str(template), dust_files=dust_files, numerical_only=True,
                      smoke_test=False, max_models=60, models=[r["parameters"] for r in records],
                      numerics=NUMERICS,
                      observations=dict(anchors_file=str(probes),
                                        spectrum_file=str(ROOT / "reference/observations/continuum_sed_R100.ecsv"),
                                        aperture_radius_arcsec=1., target_distance_pc=147.,
                                        aperture_subpixels=64, quality_policy="aperture_v2"),
                      notes="Predeclared 24 mass-size, 18 material and 18 density-profile models; 4% supplied H2O30K separate DHS. Profile masses match analytic radial envelope columns, not verified discrete line-of-sight optical depths. Fresh temperature per model. 512k packets, 1 AU pixels, mandatory signed-flux checks. No cached pilot products imported.")
        atomic_json(tmp / "grid.json", config)
        run = prepare_run(tmp / "grid.json")
        shutil.copy2(probes, run / "inputs/production_probes.ecsv")
    for source, target in copies.items():
        shutil.copy2(ROOT / "scripts" / source, run / "code" / target)
    shutil.copy2(doc, run / "README.md")
    shutil.copy2(ROOT / "requirements-cluster-recovery.txt", run / "requirements-cluster.txt")
    shutil.copy2(material_manifest, run / "inputs/laboratory_silicates_manifest.json")
    shutil.copy2(ice, run / "inputs/H2O30K_laboratory_reference.dat")
    shutil.copy2(provenance, run / "inputs/ice_input_provenance.json")
    atomic_json(run / "inputs/production_observation_contract.json", contract)
    atomic_json(run / "inputs/search_plan.json", plan())
    manifest = json.loads((run / "manifest.json").read_text())
    frozen = manifest["configuration"]
    frozen.update(output_dir="../..", template="template.para",
                  dust_files={name: "utils/Dust/" + name for name in dust_files})
    frozen["observations"].update(anchors_file="production_probes.ecsv", spectrum_file="observed_spectrum.ecsv")
    atomic_json(run / "inputs/configuration.json", frozen)
    experiment = dict(schema_version=1, experiment_id=NAME, catalogue=records,
                      tasks=[dict(index=m["index"], model_id=m["id"]) for m in manifest["models"]],
                      model_count=60, temperature_solves=60, image_requests=60*len(contract["probes"]),
                      probes_per_model=len(contract["probes"]), observable_count=len(contract["bands"]),
                      selection_rule="All three arms and all matched comparisons fixed before outcomes; no selection on fit",
                      resources=dict(cpus_per_task=64, max_concurrent_tasks=16, peak_cpus=1024,
                                     slurm_memory_per_task="160G", mcfost_memory_gb=112),
                      numerics=NUMERICS, dust_prescription=DUST,
                      predeclared_contrasts=predeclared_contrasts(records),
                      ice_input_check=ice_check, material_input_check=material_check,
                      block_weights=contract["block_weights"], quadrature_reference_fraction=.01,
                      production_convergence_certified=False, formal_likelihood=False,
                      posterior_inference=False, source_sha256=before,
                      prior_products_imported=False,
                      column_normalization="Analytic radial envelope column at Rin=1/Rout=3000 AU and fixed conical cavity. Discrete numerical columns and total source optical depths are not certified equal.",
                      preflight="Bounded envelope-only dust initialization on the cluster; production/thermal/broad-coverage wavelengths must pass for every prescription; extra zero-k stress samples remain recorded diagnostics")
    for filename in ("experiment.json", "production_experiment.json"):
        atomic_json(run / filename, experiment)
    (run / "experiment.sha256").write_text(sha256(run / "experiment.json") + "  experiment.json\n")
    # Only simulator inputs and source are versioned. Mutable cluster products
    # are excluded to keep subsequent git pulls independent of runtime output.
    (run / ".gitignore").write_text("""__pycache__/
*.pyc
.DS_Store
.mplconfig/
logs/
results/
checks/
material_checks/
material_preflight/
submissions/
machine*.cluster.json
job_*.sh
submit.sh
runtime_binding.json
production_runtime_binding.json
material_preflight.json
*.lock
models/**/status.json
models/**/measurements.json
models/**/temperature_complete.json
models/**/runtime_binding.json
models/**/temperature/
models/**/attempts/
""")
    for p in sorted(run.rglob("*")):
        if p.is_file() and p.name not in {"manifest.json", "manifest.sha256", "machine.template.json"}:
            manifest["input_hashes"][str(p.relative_to(run))] = sha256(p)
    atomic_json(run / "manifest.json", manifest)
    (run / "manifest.sha256").write_text(sha256(run / "manifest.json") + "  manifest.json\n")
    atomic_json(run / "machine.template.json", dict(
        schema_version=1, mcfost_executable="mcfost", mcfost_utils="$MCFOST_UTILS",
        backend="image_method2", threads=64, timeout_seconds=14400, max_memory_gb=112,
        slurm=dict(time="24:00:00", memory="160G", max_parallel=16, python="/SET/BY/CONFIGURE_CLUSTER",
                   setup_lines=[], analysis_cpus=2, analysis_memory="8G", analysis_time="02:00:00"),
        notes="60 models at 512k, 64 CPUs each, up to 16 concurrent. Shared venv. Bounded material preflight then production array and automatic analysis."))
    if before != {p: sha256(ROOT / p) for p in before}:
        raise RuntimeError("Sources changed during preparation; do not launch this package")
    from extinction_ice_production_task import validate_package
    validate_package(run)
    return dict(bundle=str(run), models=60, temperature_solves=60,
                images=experiment["image_requests"], bands=len(contract["bands"]),
                probes_per_model=len(contract["probes"]), input_validation_passed=True,
                mcfost_invoked=False, jobs_submitted=False,
                cluster_material_preflight_required=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / RUN_NAME)
    args = parser.parse_args(argv)
    result = plan() if args.plan else build_package(args.output)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
