#!/usr/bin/env python3
"""Validate and dispatch one frozen numerical replicate, or report its status."""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def contained(bundle, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Unsafe bundle path: {relative}")
    result = bundle / path
    if not result.resolve().is_relative_to(bundle.resolve()):
        raise ValueError(f"Bundle path escapes root: {relative}")
    return result


def validate_package(bundle):
    bundle = Path(bundle).resolve()
    expected = (bundle / "experiment.sha256").read_text().split()[0]
    if digest(bundle / "experiment.json") != expected:
        raise ValueError("Experiment manifest checksum differs")
    experiment = json.loads((bundle / "experiment.json").read_text())
    if experiment.get("schema_version") != 1 or experiment.get("experiment_id") != "extinction_ice_numerics_v2":
        raise ValueError("Unknown numerical experiment schema/identity")
    for relative, checksum in experiment["input_hashes"].items():
        path = contained(bundle, relative)
        if path.is_symlink() or not path.is_file() or digest(path) != checksum:
            raise ValueError(f"Frozen package file changed or missing: {relative}")
    tasks = experiment["tasks"]
    if [t["index"] for t in tasks] != list(range(20)):
        raise ValueError("Expected 20 contiguous numerical tasks")
    expected_cases = {(i, p, s) for i in (50., 70.) for p in (512000, 2048000) for s in range(42001, 42006)}
    if {(t["inclination_deg"], t["photon_packets"], t["seed"]) for t in tasks} != expected_cases:
        raise ValueError("Numerical cases do not cover the frozen factorial design")
    for task in tasks:
        run = contained(bundle, task["run_path"])
        manifest = json.loads((run / "manifest.json").read_text())
        if digest(run / "manifest.json") != task["manifest_sha256"]:
            raise ValueError("Task manifest checksum differs")
        model = manifest["models"][task["model_index"]]
        config = manifest["configuration"]
        if (model["id"] != task["model_id"] or model["parameters"]["inclination_deg"] != task["inclination_deg"]
                or config.get("numerical_only") is not True or any(a["score"] for a in manifest["anchors"])
                or config["numerics"]["random_seed"] != task["seed"]
                or any(config["numerics"][key] != task["photon_packets"] for key in
                       ("photons_temperature", "photons_image", "photons_sed"))):
            raise ValueError(f"Task metadata differs from frozen run: {task['index']}")
        actual = [(a["id"], a["wavelength_um"]) for a in manifest["anchors"]]
        if actual != [(a["id"], a["wavelength_um"]) for a in experiment["probes"]]:
            raise ValueError("Task probes differ from experiment")
    return experiment


def load_runner(bundle, experiment):
    code = contained(bundle, experiment["tasks"][0]["run_path"]) / "code/src"
    sys.path.insert(0, str(code))
    runner = importlib.import_module("mcfost_grid.runner")
    if not Path(runner.__file__).resolve().is_relative_to(code.resolve()):
        raise RuntimeError("A different workflow is already imported; launch this dispatcher in a fresh Python process")
    return runner


def run_task(bundle, index, machine):
    bundle = Path(bundle).resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(bundle / ".mplconfig"))
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    experiment = validate_package(bundle)
    if index < 0 or index >= len(experiment["tasks"]):
        raise ValueError("Task index must be 0..19")
    task = experiment["tasks"][index]
    runner = load_runner(bundle, experiment)
    runtime = runner.runtime_config(machine)
    if runtime["threads"] != experiment["resources"]["cpus_per_task"]:
        raise ValueError("Numerical replicas require the frozen 64-thread setting")
    if runtime["backend"] != "image_method2":
        raise ValueError("This numerical package requires image_method2")
    if runtime["threads"] > runner.cpu_capacity():
        raise ValueError("Allocated CPU count is smaller than the required thread count")
    packages = {name: importlib.import_module(name).__version__ for name in ("numpy", "scipy", "astropy", "matplotlib")}
    fingerprint = {"executable_sha256": digest(runtime["mcfost_executable"]),
                   "utilities_content_sha256": runner._utilities_hash(runtime),
                   "backend": runtime["backend"], "threads": runtime["threads"],
                   "python_version": sys.version.split()[0], "packages": packages}
    binding = bundle / "runtime_binding.json"
    with (bundle / ".runtime.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        identity = {"fingerprint": fingerprint, "experiment_sha256": digest(bundle / "experiment.json")}
        if binding.exists():
            previous = json.loads(binding.read_text())
            if any(previous.get(key) != value for key, value in identity.items()):
                raise RuntimeError("Numerical experiment runtime changed across replicas; use a new package")
        else:
            runner.atomic_json(binding, {**identity, "created_utc": runner.now()})
    print(json.dumps({"numerical_task": index, "inclination_deg": task["inclination_deg"],
                      "photons": task["photon_packets"], "seed": task["seed"], "runtime": fingerprint}), flush=True)
    return runner.run_model(contained(bundle, task["run_path"]), task["model_index"], Path(machine).resolve())


def status(bundle):
    bundle = Path(bundle).resolve()
    experiment = validate_package(bundle)
    rows = []
    for task in experiment["tasks"]:
        path = contained(bundle, task["run_path"]) / "models" / task["model_id"] / "status.json"
        value = json.loads(path.read_text()) if path.exists() else {"state": "not_started"}
        rows.append({**task, **{k: value[k] for k in ("state", "stage", "completed_anchors", "error") if k in value}})
    return {"experiment_id": experiment["experiment_id"], "counts": dict(Counter(t["state"] for t in rows)),
            "tasks": rows, "note": "File-based status; use squeue/sacct for the scheduler. Cached successes resume without new calculations."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--index", type=int)
    modes.add_argument("--status", action="store_true")
    modes.add_argument("--validate", action="store_true")
    parser.add_argument("--machine", type=Path)
    args = parser.parse_args()
    if args.status:
        result = status(args.bundle)
    elif args.validate:
        experiment = validate_package(args.bundle)
        result = {"package_valid": True, "tasks": len(experiment["tasks"]), "mcfost_invoked": False}
    else:
        if args.machine is None:
            parser.error("--index requires --machine")
        result = run_task(args.bundle, args.index, args.machine)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
