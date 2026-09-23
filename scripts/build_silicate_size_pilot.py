#!/usr/bin/env python3
"""Prepare four fresh-temperature silicate-size cases after supplied-ice checks."""
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
from silicate_pilot_design import (NAME, RUN_NAME, ICE_NAME, DUST, NUMERICS, TEMPLATE_SHA256,
                                  catalogue, plan, validate_ice)


def build_package(destination, ice_constants=None, ice_provenance=None):
    from astropy.table import Table
    ice_constants = ice_constants if ice_constants is not None else ROOT / "reference/dust" / ICE_NAME
    ice_provenance = ice_provenance if ice_provenance is not None else ROOT / "silicate/constants/H2O_30K_historical.provenance.json"
    destination, ice_constants, ice_provenance = map(lambda p: Path(p).resolve(),
                                                    (destination, ice_constants, ice_provenance))
    if destination.exists():
        raise FileExistsError(f"Refusing to replace an existing run: {destination}")
    lab = ROOT / "reference/dust" / ICE_NAME
    receipt = validate_ice(ice_constants, ice_provenance, lab)
    contract, records = build_contract(ROOT), catalogue()
    template = ROOT / "reference/parameters/ice_v02_dust.para"
    if sha256(template) != TEMPLATE_SHA256:
        raise ValueError("Archived v02 structure/thermal-grid template changed; review the pilot design first")
    copies = {"extinction_ice_production_task.py": "production_task.py",
              "analyze_extinction_ice_production.py": "analyze_extinction_ice_production.py",
              "analyze_silicate_size_pilot.py": "analyze_silicate_size_pilot.py",
              "silicate_pilot_design.py": "silicate_pilot_design.py",
              "configure_extinction_ice_numerics_v2.py": "configure_cluster.py",
              "production_observables.py": "production_observables.py"}
    sources = [template, lab, ice_constants, ice_provenance, Path(__file__).resolve(),
               ROOT / "silicate/screen_draine.json", ROOT / "silicate/FOUR_MODEL_PILOT.md",
               ROOT / "reference/observations/continuum_sed_R100.ecsv",
               ROOT / "reference/observations/tmc1a_sed_unstitched.ecsv",
               ROOT / "reference/observations/continuum_fit_mask_v1_manifest.json",
               ROOT / "requirements-cluster-recovery.txt", ROOT / "workflow.py"]
    sources += [ROOT / "scripts" / p for p in copies] + list((ROOT / "src/mcfost_grid").glob("*.py"))
    before = {str(p): sha256(p) for p in sources}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="silicate-pilot-inputs-", dir=destination.parent) as tmp:
        tmp = Path(tmp)
        probes = tmp / "probes.ecsv"
        Table(rows=contract["probes"]).write(probes, format="ascii.ecsv")
        config = dict(schema_version=1, run_name=destination.name, output_dir=str(destination.parent),
                      template=str(template), dust_files={ICE_NAME: str(ice_constants)},
                      numerical_only=True, smoke_test=False, max_models=4,
                      models=[r["parameters"] for r in records], numerics=NUMERICS,
                      observations=dict(anchors_file=str(probes),
                                        spectrum_file=str(ROOT / "reference/observations/continuum_sed_R100.ecsv"),
                                        aperture_radius_arcsec=1., target_distance_pc=147.,
                                        aperture_subpixels=64, quality_policy="aperture_v2"),
                      notes="Four fresh-temperature models. Only Draine silicate amax and inclination vary; H2O grains stay fixed at archived v02 settings. Signed image/aperture checks remain mandatory. Paired broad-band responses are primary; inherited weighted scores are descriptive.")
        atomic_json(tmp / "grid.json", config)
        run = prepare_run(tmp / "grid.json")
        shutil.copy2(probes, run / "inputs/production_probes.ecsv")
    for source, target in copies.items():
        shutil.copy2(ROOT / "scripts" / source, run / "code" / target)
    shutil.copy2(ROOT / "silicate/FOUR_MODEL_PILOT.md", run / "README.md")
    shutil.copy2(ROOT / "requirements-cluster-recovery.txt", run / "requirements-cluster.txt")
    shutil.copy2(ROOT / "silicate/screen_draine.json", run / "inputs/silicate_screen_draine.json")
    shutil.copy2(lab, run / "inputs/H2O30K_laboratory_reference.dat")
    shutil.copy2(ice_provenance, run / "inputs/ice_input_provenance.json")
    atomic_json(run / "inputs/production_observation_contract.json", contract)
    manifest = json.loads((run / "manifest.json").read_text())
    config = manifest["configuration"]
    config.update(output_dir="../..", template="template.para", dust_files={ICE_NAME: "utils/Dust/" + ICE_NAME})
    config["observations"].update(anchors_file="production_probes.ecsv", spectrum_file="observed_spectrum.ecsv")
    atomic_json(run / "inputs/configuration.json", config)
    experiment = dict(schema_version=1, experiment_id=NAME, catalogue=records,
                      tasks=[dict(index=m["index"], model_id=m["id"]) for m in manifest["models"]],
                      model_count=4, temperature_solves=4, image_requests=4*len(contract["probes"]),
                      probes_per_model=len(contract["probes"]), observable_count=len(contract["bands"]),
                      selection_rule="Predeclared full 2x2 contrast, not selection on fit",
                      resources=dict(cpus_per_task=64, max_concurrent_tasks=4, peak_cpus=256,
                                     slurm_memory_per_task="160G", mcfost_memory_gb=112),
                      numerics=NUMERICS, dust_prescription=DUST, ice_input_check=receipt,
                      pilot=dict(primary="Unscaled broad-band flux/residual changes within each inclination",
                                 axes=["silicate amax 0.4/1.0 um", "inclination 50/70 deg"],
                                 ice_grain_sizes_fixed=True, fresh_temperature_per_model=True,
                                 bare_Mie_screen_is_quantitative_prediction=False,
                                 normalization="Fixed total envelope dust mass; no adjustment to fix NIR optical depth",
                                 decisive_opacity_vs_structure_verdict=False),
                      block_weights=contract["block_weights"], quadrature_reference_fraction=.01,
                      production_convergence_certified=False, formal_likelihood=False,
                      posterior_inference=False, source_sha256=before)
    for filename in ("experiment.json", "production_experiment.json"):
        atomic_json(run / filename, experiment)
    (run / "experiment.sha256").write_text(sha256(run / "experiment.json") + "  experiment.json\n")
    for p in run.rglob("*"):
        if p.is_file() and p.name not in {"manifest.json", "manifest.sha256", "machine.template.json"}:
            manifest["input_hashes"][str(p.relative_to(run))] = sha256(p)
    atomic_json(run / "manifest.json", manifest)
    (run / "manifest.sha256").write_text(sha256(run / "manifest.json") + "  manifest.json\n")
    atomic_json(run / "machine.template.json", dict(
        schema_version=1, mcfost_executable="mcfost", mcfost_utils="$MCFOST_UTILS",
        backend="image_method2", threads=64, timeout_seconds=14400, max_memory_gb=112,
        slurm=dict(time="24:00:00", memory="160G", max_parallel=4, python="/SET/BY/CONFIGURE_CLUSTER",
                   setup_lines=[], analysis_cpus=2, analysis_memory="8G", analysis_time="02:00:00"),
        notes="Four models, up to 256 CPUs total; 1 AU/pixel; fresh temperatures. Shared venv required. Automatic paired analysis afterany; no automatic material/structure verdict."))
    if before != {p: sha256(p) for p in before}:
        raise RuntimeError("Sources changed during preparation; do not launch this package")
    from extinction_ice_production_task import validate_package
    validate_package(run)
    return dict(bundle=str(run), models=4, temperature_solves=4, images=experiment["image_requests"],
                bands=len(contract["bands"]), probes_per_model=len(contract["probes"]),
                input_validation_passed=True, mcfost_invoked=False, jobs_submitted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true", help="Print the four-model design without requiring or generating optical constants")
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / RUN_NAME)
    parser.add_argument("--ice-constants", type=Path, default=ROOT / "reference/dust" / ICE_NAME,
                        help="Defaults to the unchanged archived H2O30K table")
    parser.add_argument("--ice-provenance", type=Path, default=ROOT / "silicate/constants/H2O_30K_historical.provenance.json",
                        help="Defaults to the explicit historical/native-extrapolation assumption for this pilot")
    args = parser.parse_args()
    result = plan() if args.plan else build_package(args.output, args.ice_constants, args.ice_provenance)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
