#!/usr/bin/env python3
"""Freeze the 20-replicate V2 numerical experiment, with an optional archive.

No MCFOST process, scheduler command or network connection is invoked.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.configuration import prepare_run, sha256
from mcfost_grid.photometry import gaussian_psf_fwhm_arcsec

NAME = "extinction_ice_numerics_v2"


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def build_package(destination, archive=None):
    from astropy.table import Table
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"Refusing to replace an existing numerical package: {destination}")
    if archive is not None and Path(archive).exists():
        raise FileExistsError(f"Archive already exists: {archive}")
    design_path = ROOT / "config/extinction_ice_identifiability_v2.design.json"
    design = json.loads(design_path.read_text())
    sources = [design_path, ROOT / "reference/parameters/continuum_original_dust.para",
               ROOT / "reference/observations/continuum_sed_R100.ecsv",
               ROOT / "scripts/extinction_ice_numerical_task.py",
               ROOT / "scripts/analyze_extinction_ice_numerics_v2.py",
               ROOT / "scripts/configure_extinction_ice_numerics_v2.py",
               ROOT / "docs/EXTINCTION_ICE_NUMERICS_V2.md", Path(__file__).resolve()]
    for path in sources:
        if not path.is_file():
            raise FileNotFoundError(path)
    before = {str(p.relative_to(ROOT)): sha256(p) for p in sources}
    for folder in ("inputs", "config", "code", "subruns", "logs"):
        (destination / folder).mkdir(parents=True, exist_ok=True)
    shutil.copy2(sources[1], destination / "inputs/template.para")
    shutil.copy2(sources[2], destination / "inputs/observed_spectrum_reference.ecsv")
    shutil.copy2(design_path, destination / "inputs/design_reference.json")
    for source, target in (("extinction_ice_numerical_task.py", "numerical_task.py"),
                           ("analyze_extinction_ice_numerics_v2.py", "analyze_extinction_ice_numerics_v2.py"),
                           ("configure_extinction_ice_numerics_v2.py", "configure_cluster.py")):
        shutil.copy2(ROOT / "scripts" / source, destination / "code" / target)
    shutil.copy2(ROOT / "docs/EXTINCTION_ICE_NUMERICS_V2.md", destination / "README.md")
    shutil.copy2(ROOT / "requirements-cluster-recovery.txt", destination / "requirements-cluster.txt")
    probes = []
    waves = design["numerical_gate"]["wavelengths_um"]
    for idx, wave in enumerate(waves):
        instrument = "NIRSpec" if wave < 4.9 else "MIRI"
        # Preserve established continuum map scales; these are probes of
        # sampling noise, not a spatial-resolution or source-centre validation.
        pixels, size = ((2401, 6000.) if wave < 2.7 else
                        (2401, 3003.7468767198525) if wave < 4.9 else
                        (1201, 1502.4989583259237))
        probes.append(dict(id=f"n{idx+1:03d}", wavelength_um=wave, instrument=instrument,
                           region="numerical_probe", score=False, image_npix=pixels,
                           image_size_au=size, psf_fwhm_arcsec=gaussian_psf_fwhm_arcsec(instrument, wave)))
    table = Table(rows=probes)
    table.meta = {"purpose": "Unscored numerical probes; no observed fluxes or uncertainty placeholders",
                  "aperture_radius_arcsec": 1., "aperture_centre": "model source centre",
                  "not_an_h2o_optical_depth_or_observational_fit": True}
    table.write(destination / "inputs/numerical_probes.ecsv")
    physical = [design["models"][i]["parameters"] for i in design["numerical_gate"]["case_design_indices"]]
    assert [p["inclination_deg"] for p in physical] == [50., 70.]
    tasks = []
    for photons, seed in itertools.product(design["numerical_gate"]["photon_budgets"], design["numerical_gate"]["seed_values"]):
        subname = f"p{photons:07d}_s{seed}"
        config = {"schema_version": 1, "run_name": subname, "output_dir": "../subruns",
                  "template": "../inputs/template.para", "numerical_only": True, "smoke_test": False,
                  "max_models": 2, "models": physical,
                  "numerics": {"photons_temperature": photons, "photons_image": photons,
                               "photons_sed": photons, "random_seed": seed,
                               "grid_nr": 100, "grid_ntheta": 70, "grid_n_inner": 20, "grains": 50},
                  "observations": {"anchors_file": "../inputs/numerical_probes.ecsv",
                                   "spectrum_file": "../inputs/observed_spectrum_reference.ecsv",
                                   "aperture_radius_arcsec": 1., "target_distance_pc": 147.,
                                   "aperture_subpixels": 64, "quality_policy": "aperture_v2"},
                  "notes": "Numerical-only replicate, fresh temperature per model/seed/budget. Source-centred 1-arcsec probe. No observational score or H2O optical depth. Boundary warnings remain diagnostic under aperture_v2; finite support and flux checks remain mandatory."}
        config_path = destination / "config" / f"{subname}.json"
        write_json(config_path, config)
        run = prepare_run(config_path)
        manifest = json.loads((run / "manifest.json").read_text())
        for model in manifest["models"]:
            tasks.append(dict(index=len(tasks), run_path=str(run.relative_to(destination)),
                              model_index=model["index"], model_id=model["id"],
                              inclination_deg=model["parameters"]["inclination_deg"],
                              photon_packets=photons, seed=seed,
                              manifest_sha256=sha256(run / "manifest.json")))
    machine = {"schema_version": 1, "mcfost_executable": "mcfost", "mcfost_utils": "$MCFOST_UTILS",
               "backend": "image_method2", "threads": 64, "timeout_seconds": 14400,
               "max_memory_gb": 12,
               "slurm": {"time": "24:00:00", "memory": "16G", "max_parallel": 16,
                         "python": "/SET/BY/CONFIGURE_CLUSTER", "setup_lines": [],
                         "analysis_cpus": 2, "analysis_memory": "4G", "analysis_time": "01:00:00"},
               "notes": "Edit partition/account/reservation/setup if required. Configure with the shared mcfost-v11 Python before submitting. 64 cores/task and16 concurrent;24h/task,4h/MCFOST call are initial bounds, not timing estimates. Numerical tests only; no dust/ice likelihood."}
    write_json(destination / "machine.template.json", machine)
    experiment = {"schema_version": 1, "experiment_id": NAME,
                  "purpose": "Sampling-noise and photon-budget test; no observational fit, material inference or resolution certification",
                  "tasks": tasks, "probes": probes,
                  "resources": {"cpus_per_task": 64, "max_concurrent_tasks": 16, "peak_cpus": 1024},
                  "temperature_solves": 20, "image_requests": 120,
                  "photon_budgets": [512000, 2048000], "seeds": list(range(42001, 42006)),
                  "acceptance": {"scatter_fraction": .01, "budget_change_fraction": .01, "confidence": .95},
                  "aperture_contract": {"radius_arcsec": 1., "centre": "model source centre",
                                        "distance_pc": 140., "target_distance_pc": 147.,
                                        "quality_policy": "aperture_v2"},
                  "scope_limits": ["No observed-data scoring; existing continuum-centre offset is not validated",
                                   "No 0.35-arcsec H2O/LSF/continuum operator validation",
                                   "No broad-band quadrature, dense SED or total-flux closure validation",
                                   "No material changes, H2O30K optical extension or dust-limit test",
                                   "Six probes and two backgrounds cannot certify untested regimes",
                                   "One-third-of-physical-contrast criterion remains unassessed",
                                   "Zero measured scatter is degenerate, not proof of zero true variance"],
                  "temperature_policy": "Independent fresh equilibrium for every seed, budget and inclination; cached resume only for identical replicate",
                  "baseline_dust": "Original generic coated Mie; no H2O30K inputs required",
                  "physical_inputs_reference": "inputs/design_reference.json",
                  "source_sha256": before,
                  "input_hashes": {str(p.relative_to(destination)): sha256(p) for p in sorted(destination.rglob("*"))
                                   if p.is_file() and p.name != "machine.template.json"}}
    write_json(destination / "experiment.json", experiment)
    (destination / "experiment.sha256").write_text(sha256(destination / "experiment.json") + "  experiment.json\n")
    if before != {str(p.relative_to(ROOT)): sha256(p) for p in sources}:
        raise RuntimeError("A source changed while preparing the package; do not launch this bundle")
    if archive is not None:
        archive = Path(archive).resolve()
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "w:gz") as stream:
            stream.add(destination, arcname=NAME)
        archive.with_name(archive.name + ".sha256").write_text(sha256(archive) + "  " + archive.name + "\n")
    return {"bundle": str(destination), "archive": str(archive) if archive else None,
            "tasks": len(tasks), "temperature_solves": 20, "images": 120,
            "mcfost_invoked": False, "jobs_submitted": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / NAME)
    delivery = parser.add_mutually_exclusive_group()
    delivery.add_argument("--archive", type=Path, default=ROOT / "transfer" / f"{NAME}.tar.gz")
    delivery.add_argument("--no-archive", action="store_true",
                          help="Prepare directly in this checkout after git pull, without creating a transfer archive")
    args = parser.parse_args()
    print(json.dumps(build_package(args.output, None if args.no_archive else args.archive), indent=2))


if __name__ == "__main__":
    main()
