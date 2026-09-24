#!/usr/bin/env python3
"""Validate and run one frozen 1-AU production model; resume identical products."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys

NAME = "extinction_ice_production_1au_v1"
SILICATE_PILOT = "silicate_size_pilot_v1"
SILICATE_STRUCTURE = "silicate_structure_production_v1"


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_package(bundle, index=None):
    bundle = Path(bundle).resolve()
    for filename in ("experiment", "manifest"):
        require(digest(bundle / f"{filename}.json") == (bundle / f"{filename}.sha256").read_text().split()[0],
                f"{filename} checksum differs")
    experiment = json.loads((bundle / "experiment.json").read_text())
    require((bundle / "production_experiment.json").read_bytes() == (bundle / "experiment.json").read_bytes(),
            "Production experiment copies differ")
    require(experiment.get("schema_version") == 1 and experiment.get("experiment_id") in {NAME, SILICATE_PILOT, SILICATE_STRUCTURE}, "Wrong experiment")
    manifest = json.loads((bundle / "manifest.json").read_text())
    if index is not None:
        require(type(index) is int and 0 <= index < len(manifest["models"]), "Model index out of range")
    model_id = manifest["models"][index]["id"] if index is not None else None
    for relative, checksum in manifest["input_hashes"].items():
        parts = Path(relative).parts
        require(parts and not Path(relative).is_absolute() and ".." not in parts, f"Unsafe frozen path: {relative}")
        if model_id is not None and parts[0] == "models" and len(parts) > 1 and parts[1] != model_id:
            continue  # Full preparation/analysis checks all; workers check their own model plus every shared input.
        path = bundle / relative
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts
                and path.resolve().is_relative_to(bundle) and not path.is_symlink()
                and path.is_file() and digest(path) == checksum, f"Frozen input changed: {relative}")
    require(len(experiment["tasks"]) == len(manifest["models"]) == len(experiment["catalogue"]), "Wrong model count")
    require([r["index"] for r in experiment["tasks"]] == list(range(len(manifest["models"]))), "Task indices differ")
    for task, model, metadata in zip(experiment["tasks"], manifest["models"], experiment["catalogue"]):
        require(task["index"] == model["index"] == metadata["index"] and task["model_id"] == model["id"]
                and metadata["parameters"] == model["parameters"], "Task/physics identity differs")
    config = manifest["configuration"]
    require(config.get("numerical_only") is True and all(a["score"] is False for a in manifest["anchors"]),
            "Production probes must use their separate integrated-band analysis")
    photon_count = 512000 if experiment["experiment_id"] in {SILICATE_PILOT, SILICATE_STRUCTURE} else 2048000
    require(all(config["numerics"][k] == photon_count for k in ("photons_temperature", "photons_sed", "photons_image")),
            "Photon count differs")
    require(config["numerics"]["image_npix"] == 6001 and config["numerics"]["image_size_au"] == 6000.,
            "Expected approximately 1 AU/pixel at fixed 6000-AU field")
    require(manifest["measurement"] == dict(aperture_radius_arcsec=1., target_distance_pc=147.,
            aperture_subpixels=64, quality_policy="aperture_v2"), "Measurement policy differs")
    contract = json.loads((bundle / "inputs/production_observation_contract.json").read_text())
    require([(a["id"], a["wavelength_um"]) for a in manifest["anchors"]]
            == [(a["id"], a["wavelength_um"]) for a in contract["probes"]], "Probe contract differs")
    if experiment["experiment_id"] == SILICATE_PILOT:
        from silicate_pilot_design import validate_frozen_pilot
        validate_frozen_pilot(bundle, experiment, manifest, index)
    elif experiment["experiment_id"] == SILICATE_STRUCTURE:
        require(len(manifest["models"]) == 60, "Silicate structure campaign requires 60 models")
        from silicate_structure_design import validate_frozen_campaign
        validate_frozen_campaign(bundle, experiment, manifest, index)
    return experiment


def load_runner(bundle):
    code = Path(bundle).resolve() / "code/src"
    sys.path.insert(0, str(code))
    runner = importlib.import_module("mcfost_grid.runner")
    require(Path(runner.__file__).resolve().is_relative_to(code), "A different workflow is already imported; use a fresh process")
    return runner


def run_task(bundle, index, machine):
    bundle = Path(bundle).resolve()
    experiment = validate_package(bundle, index)
    require(0 <= index < len(experiment["tasks"]), "Model index out of range")
    runner = load_runner(bundle)
    runtime = runner.runtime_config(machine)
    require(runtime["threads"] == 64 and runtime["backend"] == "image_method2"
            and runtime["max_memory_gb"] == 112, "Production requires 64 threads, image_method2 and max_memory_gb=112")
    if experiment["experiment_id"] == SILICATE_STRUCTURE:
        sys.path.insert(0, str(bundle / "code"))
        from check_materials import validate_receipt
        validate_receipt(bundle, runtime)
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(bundle / ".mplconfig"))
    packages = {name: {"version": importlib.import_module(name).__version__,
                       "file": importlib.import_module(name).__file__}
                for name in ("numpy", "scipy", "astropy", "matplotlib")}
    identity = dict(python=os.path.abspath(sys.executable), python_version=sys.version.split()[0],
                    packages=packages, threads=runtime["threads"], max_memory_gb=runtime["max_memory_gb"])
    with (bundle / ".production_runtime.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        binding = bundle / "production_runtime_binding.json"
        if binding.exists():
            require(json.loads(binding.read_text())["identity"] == identity,
                    "Production Python/packages/thread/memory environment changed across workers")
        else:
            runner.atomic_json(binding, dict(identity=identity, created_utc=runner.now()))
    return runner.run_model(bundle, index, Path(machine).resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--index", type=int)
    modes.add_argument("--validate", action="store_true")
    modes.add_argument("--status", action="store_true")
    parser.add_argument("--machine", type=Path)
    args = parser.parse_args()
    if args.validate:
        experiment = validate_package(args.bundle)
        result = dict(package_valid=True, models=len(experiment["tasks"]), mcfost_invoked=False)
    elif args.status:
        validate_package(args.bundle)
        result = load_runner(args.bundle).status(args.bundle)
    else:
        if args.machine is None:
            parser.error("--index requires --machine")
        result = run_task(args.bundle, args.index, args.machine)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
