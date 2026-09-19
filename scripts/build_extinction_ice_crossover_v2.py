#!/usr/bin/env python3
"""Prepare ten pixel-scale images with one saved temperature, without Slurm."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from extinction_ice_crossover_task import (
    DEFAULT_RUN_NAME, DIAGNOSTIC_ID, GEOMETRIES, IDENTITY_KEYS, NAME, SEEDS, TEMPERATURE_SEEDS, WAVELENGTH,
    contained, digest, read_json, require,
)
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.configuration import atomic_json, now
from mcfost_grid.physics import render_parameter
from mcfost_grid.runner import _valid_temperature


def source_models(source):
    """Validate the selected temperature source only; permit a partial download."""
    source = Path(source).resolve()
    exp = read_json(source / "experiment.json")
    require(exp.get("experiment_id") == "extinction_ice_numerics_v2" and exp.get("schema_version") == 1,
            "Source must be the original V2 numerical experiment")
    checksum = digest(source / "experiment.json")
    require(checksum == (source / "experiment.sha256").read_text().split()[0], "Source experiment checksum differs")
    global_binding = read_json(source / "runtime_binding.json")
    require(global_binding["experiment_sha256"] == checksum, "Source runtime refers to a different experiment")
    records = []
    for seed in TEMPERATURE_SEEDS:
        selected = [t for t in exp["tasks"] if t["photon_packets"] == 2048000
                    and t["seed"] == seed and t["inclination_deg"] == 50]
        require(len(selected) == 1, f"Missing or duplicate source case for temperature seed {seed}")
        task = selected[0]
        run = contained(source, task["run_path"])
        manifest_path = run / "manifest.json"
        manifest = read_json(manifest_path)
        require(digest(manifest_path) == task["manifest_sha256"]
                == exp["input_hashes"][f"{task['run_path']}/manifest.json"], "Source manifest checksum differs")
        model = manifest["models"][task["model_index"]]
        require(model["id"] == task["model_id"] and model["parameters"]["inclination_deg"] == 50,
                "Source model identity differs")
        numerics = manifest["configuration"]["numerics"]
        require(numerics["random_seed"] == seed and all(numerics[k] == 2048000 for k in
                ("photons_temperature", "photons_image", "photons_sed")), "Wrong source photon setting or seed")
        logical = Path("models") / model["id"]
        directory = run / logical
        if not (directory / "temperature_complete.json").is_file():
            directory = run / model["id"]
        require(directory.resolve().is_relative_to(source), "Source model escapes the source folder")

        def source_file(relative):
            path = Path(relative)
            return contained(directory, path.relative_to(logical)) if path.is_relative_to(logical) else contained(run, relative)

        verified = {}
        for relative, expected in manifest["input_hashes"].items():
            if relative.startswith((str(logical) + "/", "inputs/utils/", "code/", "inputs/template.para")):
                path = source_file(relative)
                require(path.is_file() and not path.is_symlink() and digest(path) == expected,
                        f"Source input changed or missing: {task['run_path']}/{relative}")
                verified[relative] = expected
        require(f"{logical}/anchors/n002/image.para" in verified and f"{logical}/temperature.para" in verified,
                "Source parameter files are not pinned")
        receipt = read_json(directory / "temperature_complete.json")
        temperature = source_file(receipt["path"])
        require(temperature.resolve().is_relative_to(directory.resolve()) and temperature.is_file()
                and digest(temperature) == receipt["sha256"], "Source temperature checksum/path differs")
        _valid_temperature(temperature)
        binding = read_json(directory / "runtime_binding.json")
        measured = read_json(directory / "measurements.json")
        require(receipt["fingerprint"] == binding["fingerprint"] == measured["fingerprint"]
                and receipt["fingerprint"]["manifest_sha256"] == task["manifest_sha256"], "Source model fingerprints differ")
        require(measured["complete"] is True and measured["model_id"] == model["id"]
                and measured["parameters"] == model["parameters"], "Source model is incomplete or has different parameters")
        identity = {k: binding["fingerprint"][k] for k in IDENTITY_KEYS if k != "threads"}
        identity["threads"] = binding["runtime"]["threads"]
        require(identity["threads"] == 64 and identity["backend"] == "image_method2"
                and all(global_binding["fingerprint"][k] == identity[k] for k in IDENTITY_KEYS), "Source runtimes differ")
        attempt = next((p for p in temperature.parents if p.is_relative_to(directory)
                        and (p / "command.json").is_file()), None)
        require(attempt is not None, "Source temperature lacks execution receipts")
        command = read_json(attempt / "command.json")
        execution = read_json(attempt / "execution.json")
        args = command["argv"]
        require(args.count("-seed") == 1 and args[args.index("-seed") + 1] == str(seed)
                and "-img" not in args and "-no_T" not in args and int(command["threads"]) == 64
                and execution["returncode"] == 0 and "Processing complete" in (attempt / "mcfost.log").read_text(),
                "Source temperature was not a successful independent solve with the required seed")
        anchor = next(a for a in manifest["anchors"] if a["id"] == "n002")
        require(anchor["wavelength_um"] == WAVELENGTH and anchor["image_npix"] == 2401
                and anchor["image_size_au"] == 6000., "Unexpected source n002 wavelength/geometry")
        provenance = [manifest_path, directory / "temperature_complete.json", directory / "runtime_binding.json",
                      directory / "measurements.json", attempt / "command.json", attempt / "execution.json", attempt / "mcfost.log"]
        records.append(dict(seed=seed, run=run, directory=directory, task=task, model=model, manifest=manifest,
                            temperature=temperature, temperature_sha256=receipt["sha256"],
                            image_parameter=directory / "anchors/n002/image.para", verified_inputs=verified,
                            runtime=identity, machine=binding["runtime"], anchor=anchor, provenance=provenance))
    return records


def build_package(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    require(not destination.is_relative_to(source) and not source.is_relative_to(destination),
            "Destination must be separate from the source experiment")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite an existing diagnostic: {destination}")
    records = source_models(source)
    files = {"scripts/extinction_ice_crossover_task.py": "code/crossover_task.py",
             "scripts/analyze_extinction_ice_crossover_v2.py": "code/analyze_extinction_ice_crossover_v2.py",
             "scripts/configure_extinction_ice_numerics_v2.py": "code/configure_cluster.py",
             "docs/EXTINCTION_ICE_CROSSOVER_V2.md": "README.md",
             "requirements-cluster-recovery.txt": "requirements-cluster.txt"}
    for path in sorted((ROOT / "src/mcfost_grid").glob("*.py")):
        files[str(path.relative_to(ROOT))] = f"code/src/mcfost_grid/{path.name}"
    source_hashes = {relative: digest(ROOT / relative) for relative in [*files, "scripts/build_extinction_ice_crossover_v2.py"]}
    destination.mkdir(parents=True)
    for relative, target in files.items():
        path = destination / target
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, path)
    for folder in ("inputs/temperatures", "inputs/provenance", "logs", "tasks"):
        (destination / folder).mkdir(parents=True, exist_ok=True)
    for filename in ("experiment.json", "experiment.sha256", "runtime_binding.json"):
        shutil.copy2(source / filename, destination / "inputs/provenance" / filename)
    first = records[0]
    shutil.copy2(first["image_parameter"], destination / "inputs/template_image.para")
    shutil.copytree(first["run"] / "inputs/utils", destination / "inputs/utils")
    require(digest(destination / "inputs/template_image.para") == first["verified_inputs"][
        f"models/{first['model']['id']}/anchors/n002/image.para"], "Source image template changed during copying")
    for relative, expected in first["verified_inputs"].items():
        if relative.startswith("inputs/utils/"):
            require(digest(destination / relative) == expected, f"Source utility changed during copying: {relative}")
    provenance = []
    tasks = []
    for record in records:
        t_seed = record["seed"]
        temp_path = f"inputs/temperatures/T{t_seed}.fits.gz"
        shutil.copy2(record["temperature"], destination / temp_path)
        provenance_dir = destination / "inputs/provenance" / f"T{t_seed}"
        provenance_dir.mkdir()
        for path in record["provenance"]:
            shutil.copy2(path, provenance_dir / path.name)
        provenance.append({"temperature_seed": t_seed, "source_run_path": record["task"]["run_path"],
                           "source_model_id": record["model"]["id"], "source_model_directory": str(record["directory"].relative_to(source)),
                           "temperature_sha256": record["temperature_sha256"],
                           "source_temperature_path": str(record["temperature"].relative_to(source)),
                           "verified_source_inputs": record["verified_inputs"]})
        for seed in SEEDS:
            index = len(tasks)
            task_path = f"tasks/t{index:03d}_T{t_seed}_s{seed}"
            task = dict(index=index, temperature_seed=t_seed, image_seed=seed, task_path=task_path,
                        temperature_path=temp_path, temperature_sha256=record["temperature_sha256"])
            tasks.append(task)
            for geometry in GEOMETRIES:
                work = destination / task_path / geometry["id"]
                work.mkdir(parents=True)
                text = render_parameter(record["image_parameter"].read_text(), {}, {
                    "image_npix": geometry["image_npix"], "image_size_au": geometry["image_size_au"]}, "image", WAVELENGTH)
                (work / "image.para").write_text(text)
    machine = {key: value for key, value in first["machine"].items() if key in
               {"schema_version", "mcfost_executable", "mcfost_utils", "backend", "threads", "timeout_seconds", "max_memory_gb", "slurm"}}
    machine.update(schema_version=1, backend="image_method2", threads=64)
    machine.setdefault("timeout_seconds", 14400)
    # A 4801², eight-plane, 64-thread MCFOST image array alone is about
    # 47 GB in the audited default-real build. Increase both limits for the
    # entire pair so resolution is the only changed simulation setting.
    machine["max_memory_gb"] = 64
    machine["slurm"] = {**machine.get("slurm", {}), "max_parallel": 16, "python": "/SET/BY/CONFIGURE_CLUSTER"}
    machine["slurm"].setdefault("time", "24:00:00")
    machine["slurm"]["memory"] = "96G"
    machine["slurm"].setdefault("analysis_cpus", 2)
    machine["slurm"]["analysis_memory"] = "8G"
    machine["slurm"].setdefault("analysis_time", "01:00:00")
    atomic_json(destination / "machine.template.json", machine)
    experiment = dict(schema_version=2, experiment_id=NAME, diagnostic_id=DIAGNOSTIC_ID, created_utc=now(),
        purpose="Pixel-scale diagnostic at one fixed temperature and field; image seeds are matched replicates, no observational fit",
        tasks=tasks, geometries=GEOMETRIES, temperature_seeds=list(TEMPERATURE_SEEDS), image_seeds=list(SEEDS),
        photon_packets=2048000, temperature_solves=0, image_requests=10, wavelength_um=WAVELENGTH,
        parameters=first["model"]["parameters"], psf_fwhm_arcsec=first["anchor"]["psf_fwhm_arcsec"],
        measurement={"aperture_radius_arcsec": 1., "target_distance_pc": 147., "aperture_subpixels": 64, "quality_policy": "aperture_v2"},
        resources={"cpus_per_task": 64, "max_concurrent_tasks": 16, "scheduled_max_concurrent_tasks": 5, "peak_cpus": 320},
        comparison={"confidence": .95, "reference_fraction": .01, "production_convergence_certified": False},
        source_experiment_sha256=digest(source / "experiment.json"), source_runtime=first["runtime"],
        source_models=provenance, source_sha256=source_hashes,
        limitations=["Only pixel count changes: both geometries retain a 6000-AU field, fixed physics, temperature and image-photon setting.",
                     "Two pixel scales do not certify spatial convergence.",
                     "The temperature arm is excluded by design; every image uses the saved 42004 temperature. No photon-budget contrast is made.",
                     "Same image seeds label paired comparisons but do not guarantee identical OpenMP random streams.",
                     "A resumed coarse/fine pair can span workers; analysis must report that fact.",
                     "No observational likelihood, broad-band, material or H2O inference is performed."],
        input_hashes={str(p.relative_to(destination)): digest(p) for p in sorted(destination.rglob("*"))
                      if p.is_file() and p.name != "machine.template.json"})
    atomic_json(destination / "experiment.json", experiment)
    (destination / "experiment.sha256").write_text(digest(destination / "experiment.json") + "  experiment.json\n")
    require(source_hashes == {relative: digest(ROOT / relative) for relative in source_hashes}, "Code changed during preparation")
    for record in records:
        require(digest(record["temperature"]) == record["temperature_sha256"], "Source temperature changed during preparation")
        require(digest(destination / f"inputs/temperatures/T{record['seed']}.fits.gz") == record["temperature_sha256"], "Temperature copy mismatch")
    return dict(bundle=str(destination), tasks=5, images=10, temperature_solves=0, peak_cpus=320,
                mcfost_invoked=False, jobs_submitted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Existing V2 numerical run with the 2.048M, 50-degree seed-42004 source model")
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / DEFAULT_RUN_NAME)
    args = parser.parse_args()
    print(json.dumps(build_package(args.source, args.output), indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
