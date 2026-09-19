#!/usr/bin/env python3
"""Run paired image geometries with frozen temperatures; never solve temperature."""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import sys

NAME = "extinction_ice_crossover_v2"  # Retained tool family; schema 2 rejects the old temperature arm.
DIAGNOSTIC_ID = "fixed_temperature_pixel_scale_v2"
DEFAULT_RUN_NAME = "extinction_ice_pixel_scale_v2"
SEEDS = tuple(range(42001, 42006))
TEMPERATURE_SEEDS = (42004,)
WAVELENGTH = 2.546352514425122
GEOMETRIES = [{"id": "coarse", "image_npix": 2401, "image_size_au": 6000.0},
              {"id": "fine", "image_npix": 4801, "image_size_au": 6000.0}]
IDENTITY_KEYS = ("executable_sha256", "utilities_content_sha256", "backend", "threads")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    def reject(value):
        raise ValueError(f"Nonfinite JSON value: {value}")
    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=reject)


def contained(bundle, relative):
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ValueError(f"Unsafe bundle path: {relative}")
    path = Path(bundle) / relative
    if not path.resolve().is_relative_to(Path(bundle).resolve()):
        raise ValueError(f"Path escapes bundle: {relative}")
    return path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_modules(bundle):
    """Import only the workflow frozen inside this bundle, under a private name."""
    code = Path(bundle).resolve() / "code/src/mcfost_grid"
    name = "_crossover_" + hashlib.sha256(str(code).encode()).hexdigest()[:16]
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, code / "__init__.py",
                                                      submodule_search_locations=[str(code)])
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return importlib.import_module(name + ".runner"), importlib.import_module(name + ".photometry")


def validate_package(bundle):
    bundle = Path(bundle).resolve()
    require(digest(bundle / "experiment.json") == (bundle / "experiment.sha256").read_text().split()[0],
            "Experiment checksum differs")
    exp = read_json(bundle / "experiment.json")
    require(exp.get("schema_version") == 2 and exp.get("experiment_id") == NAME
            and exp.get("diagnostic_id") == DIAGNOSTIC_ID, "Expected schema-2 fixed-temperature pixel-scale experiment")
    hashes = exp.get("input_hashes", {})
    require(bool(hashes), "Missing frozen input hashes")
    for relative, expected in hashes.items():
        path = contained(bundle, relative)
        require(not path.is_symlink() and path.is_file() and digest(path) == expected,
                f"Frozen input changed or missing: {relative}")
    tasks = exp["tasks"]
    require(len(tasks) == 5 and [t["index"] for t in tasks] == list(range(5)), "Expected five paired tasks")
    require({(t["temperature_seed"], t["image_seed"]) for t in tasks} ==
            {(t, s) for t in TEMPERATURE_SEEDS for s in SEEDS}, "Expected fixed temperature 42004 and five image seeds")
    require(len({t["temperature_sha256"] for t in tasks}) == 1, "All images must reuse the same temperature bytes")
    require(exp["temperature_seeds"] == [42004] and exp["image_seeds"] == list(SEEDS), "Temperature variation is excluded from this design")
    require(len({t["task_path"] for t in tasks}) == 5, "Duplicate task destinations")
    require(exp["geometries"] == GEOMETRIES and exp["wavelength_um"] == WAVELENGTH,
            "Wrong geometries or wavelength")
    require(exp["photon_packets"] == 2048000 and exp["temperature_solves"] == 0
            and exp["image_requests"] == 10 and exp["resources"]["cpus_per_task"] == 64,
            "Wrong photon budget, task resources or calculation counts")
    require(exp["parameters"]["inclination_deg"] == 50 and exp["parameters"]["distance_pc"] == 140,
            "Wrong physical case/distance")
    require(exp["measurement"] == {"aperture_radius_arcsec": 1., "target_distance_pc": 147.,
                                  "aperture_subpixels": 64, "quality_policy": "aperture_v2"},
            "Changed measurement contract")
    required = {"code/crossover_task.py", "code/analyze_extinction_ice_crossover_v2.py",
                "code/configure_cluster.py", "inputs/template_image.para"}
    required.update(f"code/src/mcfost_grid/{name}.py" for name in
                    ("__init__", "configuration", "physics", "runner", "photometry", "quality"))
    for task in tasks:
        contained(bundle, task["task_path"])
        require(task["temperature_path"] == f"inputs/temperatures/T{task['temperature_seed']}.fits.gz",
                "Wrong fixed temperature path")
        require(hashes.get(task["temperature_path"]) == task["temperature_sha256"], "Unpinned temperature")
        required.add(task["temperature_path"])
        for geometry in GEOMETRIES:
            required.add(f"{task['task_path']}/{geometry['id']}/image.para")
    require(required <= set(hashes), "Required executable/temperature/parameter input is not pinned")
    require(any(p.startswith("inputs/utils/Dust/") for p in hashes)
            and any(p.startswith("inputs/utils/Stellar_Spectra/") for p in hashes), "Missing frozen utilities")
    runner, _ = load_modules(bundle)
    physics = importlib.import_module(runner.__package__ + ".physics")
    template = (bundle / "inputs/template_image.para").read_text()
    for task in tasks:
        runner._valid_temperature(contained(bundle, task["temperature_path"]))
        for geometry in GEOMETRIES:
            expected = physics.render_parameter(template, {}, {
                "image_npix": geometry["image_npix"], "image_size_au": geometry["image_size_au"]}, "image", WAVELENGTH)
            require((contained(bundle, task["task_path"]) / geometry["id"] / "image.para").read_text() == expected,
                    "Image parameter file differs from the frozen geometry")
    return exp


def image_measurement(image, exp, geometry, photometry):
    """Canonical total I and diagnostic components, with actual geometry checks."""
    import numpy as np
    from astropy.io import fits
    measured = photometry.measure_image(image, exp["wavelength_um"], "NIRSpec",
                                         model_distance_pc=exp["parameters"]["distance_pc"],
                                         psf_fwhm_arcsec=exp["psf_fwhm_arcsec"], **exp["measurement"])
    d = measured["diagnostics"]
    expected_pixel = geometry["image_size_au"] / geometry["image_npix"] / exp["parameters"]["distance_pc"]
    require(d["nx"] == d["ny"] == geometry["image_npix"]
            and math.isclose(d["model_pixel_arcsec"], expected_pixel, rel_tol=5e-6), "Image geometry differs from request")
    weights, _, _ = photometry._adjoint_weights((d["ny"], d["nx"]), d["centre_x_zero_based"],
        d["centre_y_zero_based"], d["aperture_radius_pixels"], d["psf_sigma_pixels"], d["aperture_subpixels"])
    factor = d["distance_flux_factor"] * exp["wavelength_um"] * 1e-6 / photometry.SPEED_OF_LIGHT_M_S * 1e26
    # Only five intensity planes are needed. Convert one at a time to avoid
    # an additional 1.5-GB float64 copy of the full fine-grid Stokes cube.
    flux = np.zeros(8, dtype=np.float64)
    with fits.open(image, memmap=False) as hdus:
        for i in (0, 4, 5, 6, 7):
            plane = np.asarray(hdus[0].data[i, 0, 0], dtype=np.float64)
            flux[i] = np.einsum("ij,ij->", plane, weights, optimize=False) * factor
    require(math.isclose(float(flux[0]), measured["flux_jy"], rel_tol=1e-10), "Component plane-zero mismatch")
    measured["components_aperture_jy"] = {name: float(flux[i]) for i, name in
        ((0, "total_i"), (4, "direct_star"), (5, "scattered_star"), (6, "direct_thermal"), (7, "scattered_thermal"))}
    measured["component_sum_minus_total_aperture_fraction"] = float((flux[4:8].sum() - flux[0]) / flux[0])
    return measured


def _identity(bundle, exp, task, runtime):
    runner, _ = load_modules(bundle)
    identity = {"executable_sha256": digest(runtime["mcfost_executable"]),
                "utilities_content_sha256": runner._utilities_hash(runtime),
                "backend": runtime["backend"], "threads": runtime["threads"]}
    require(all(identity[k] == exp["source_runtime"][k] for k in IDENTITY_KEYS),
            "Simulator/utilities/backend/threads differ from the source experiment")
    identity.update(max_memory_gb=runtime["max_memory_gb"],
                    python_version=sys.version.split()[0], packages={name: importlib.import_module(name).__version__
                    for name in ("numpy", "scipy", "astropy", "matplotlib")})
    fingerprint = {"experiment_sha256": digest(bundle / "experiment.json"), "runtime": identity}
    with (bundle / ".runtime.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = bundle / "runtime_binding.json"
        if path.exists():
            require(read_json(path)["fingerprint"] == fingerprint, "Crossover runtime changed; prepare a new bundle")
        else:
            runner.atomic_json(path, {"fingerprint": fingerprint, "created_utc": runner.now()})
    return fingerprint


def _checked_rows(bundle, exp, task, result, remeasure=False):
    runner, photometry = load_modules(bundle)
    fingerprint = result["fingerprint"]
    require(fingerprint["experiment_sha256"] == digest(bundle / "experiment.json"), "Cached experiment differs")
    binding = read_json(bundle / "runtime_binding.json")
    require(binding["fingerprint"] == fingerprint, "Cached runtime differs")
    require(all(fingerprint["runtime"][k] == exp["source_runtime"][k] for k in IDENTITY_KEYS), "Wrong source runtime")
    require(result["task_index"] == task["index"], "Wrong task result")
    seen, rows = set(), result["measurements"]
    for row in rows:
        geometry = next((g for g in GEOMETRIES if g["id"] == row["geometry"]), None)
        require(geometry is not None and geometry["id"] not in seen, "Duplicate or unknown geometry")
        seen.add(geometry["id"])
        require(all(row[k] == task[k] for k in ("temperature_seed", "image_seed", "temperature_sha256"))
                and row["task_index"] == task["index"] and row["threads"] == 64
                and row["image_size_au"] == geometry["image_size_au"]
                and row["image_npix"] == geometry["image_npix"] and row["wavelength_um"] == WAVELENGTH,
                "Cached measurement has wrong case metadata")
        error = importlib.import_module(runner.__package__ + ".quality").cached_policy_error(row, "aperture_v2")
        require(not error and row.get("quality_pass") is True and math.isfinite(row["flux_jy"])
                and row["flux_jy"] > 0, f"Invalid cached measurement: {error}")
        image = contained(bundle, row["image_path"])
        require(image.is_file() and digest(image) == row["image_sha256"], "Cached image changed or missing")
        attempt = contained(bundle, row["attempt_path"])
        work = contained(bundle, task["task_path"]) / geometry["id"]
        require(attempt.parent == work / "attempts" and image.is_relative_to(attempt / "output"),
                "Cached product belongs to another task/geometry")
        require(read_json(attempt / "measurement.json") == row, "Cached measurement differs from its attempt receipt")
        for name, expected in row["receipt_sha256"].items():
            require(name in ("command.json", "execution.json", "mcfost.log")
                    and digest(attempt / name) == expected, "Execution receipt changed")
        require(set(row["receipt_sha256"]) == {"command.json", "execution.json", "mcfost.log"}, "Missing execution receipts")
        command = read_json(attempt / "command.json")
        args = command["argv"]
        require(command["threads"] == "64" and read_json(attempt / "execution.json")["returncode"] == 0,
                "Failed or wrong-thread execution receipt")
        require(args == row["command_argv"], "Command identity differs")
        for flag, value in (("-seed", str(task["image_seed"])), ("-img", str(WAVELENGTH)),
                            ("-Tfile", "Temperature.fits.gz")):
            require(args.count(flag) == 1 and args[args.index(flag) + 1] == value, f"Wrong {flag} in command")
        require("-no_T" in args and "-rt2" in args and "Processing complete" in (attempt / "mcfost.log").read_text(),
                "Missing fixed-temperature Method-2 completion")
        require(digest(attempt / "output/Temperature.fits.gz") == task["temperature_sha256"], "Attempt temperature differs")
        if remeasure:
            fresh = image_measurement(image, exp, geometry, photometry)
            require(math.isclose(fresh["flux_jy"], row["flux_jy"], rel_tol=1e-10), "Stored flux differs from raw image")
            for name, value in fresh["components_aperture_jy"].items():
                require(math.isclose(value, row["components_aperture_jy"][name], rel_tol=1e-10, abs_tol=1e-25),
                        "Stored component differs from raw image")
    require(not result["complete"] or seen == {"coarse", "fine"}, "Complete task lacks a geometry")
    return rows


def inspect_task(bundle, experiment, task, remeasure=True):
    bundle = Path(bundle).resolve()
    path = contained(bundle, task["task_path"])
    status = {**task, "state": "not_started"}
    try:
        if (path / "status.json").is_file():
            marker = read_json(path / "status.json")
            status.update(state=marker.get("state", "not_started"), recorded_state=marker.get("state"), error=marker.get("error"))
        result_path = path / "measurements.json"
        if not result_path.is_file():
            require(status["state"] != "complete", "Complete marker without measurements")
            return status, []
        result = read_json(result_path)
        rows = _checked_rows(bundle, experiment, task, result, remeasure)
        status.update(state="complete" if result["complete"] else "partial", completed_images=len(rows))
        return status, rows
    except (ValueError, RuntimeError, OSError, KeyError, TypeError, IndexError) as exc:
        status.update(state="invalid", error=str(exc))
        return status, []


def run_task(bundle, index, machine):
    bundle = Path(bundle).resolve()
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(bundle / ".mplconfig"))
    exp = validate_package(bundle)
    require(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < 5, "Task index must be 0..4")
    task = exp["tasks"][index]
    runner, photometry = load_modules(bundle)
    runtime = runner.runtime_config(machine)
    require(runtime["threads"] == 64 and runtime["backend"] == "image_method2"
            and runner.cpu_capacity() >= 64, "Require 64 allocated CPUs and image_method2")
    directory = contained(bundle, task["task_path"])
    directory.mkdir(exist_ok=True)
    with (directory / ".task.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("This paired task is already running") from exc
        status_path = directory / "status.json"
        try:
            fingerprint = _identity(bundle, exp, task, runtime)
            result_path = directory / "measurements.json"
            result = read_json(result_path) if result_path.exists() else {
                "task_index": index, "fingerprint": fingerprint, "complete": False, "measurements": []}
            require(result["fingerprint"] == fingerprint, "Cached runtime differs")
            rows = _checked_rows(bundle, exp, task, result) if result_path.exists() else []
            if result["complete"]:
                runner.atomic_json(status_path, {"state": "complete", "completed_images": 2, "updated_utc": runner.now()})
                return {"index": index, "status": "cached", "images": 2}
            runner.check_backend(bundle, runtime)
            env = runner.environment(bundle, runtime)
            physics = importlib.import_module(runner.__package__ + ".physics")
            arguments = [*physics.physical_args(exp["parameters"]), *physics.numerical_args({"random_seed": task["image_seed"]})]
            temperature = contained(bundle, task["temperature_path"])
            for geometry in GEOMETRIES:
                if any(r["geometry"] == geometry["id"] for r in rows):
                    continue
                runner.atomic_json(status_path, {"state": "running", "geometry": geometry["id"],
                    "completed_images": len(rows), "host": socket.gethostname(), "pid": os.getpid(),
                    "slurm_job": os.environ.get("SLURM_JOB_ID"), "updated_utc": runner.now()})
                work = directory / geometry["id"]
                attempt = runner._new_attempt(work / "attempts")
                output = attempt / "output"
                output.mkdir()
                (output / "Temperature.fits.gz").symlink_to(os.path.relpath(temperature, output))
                command = [runtime["mcfost_executable"], "image.para", "-img", str(WAVELENGTH),
                           *arguments, "-rt2", "-no_T", "-Tfile", "Temperature.fits.gz",
                           "-root_dir", os.path.relpath(output, work), "-max_mem", str(runtime["max_memory_gb"]), "-no_backup"]
                elapsed, log = runner._invoke(command, work, env, runtime["timeout_seconds"], attempt)
                require(any(mark in log for mark in ("Using ray-tracing method 2", "RT source-function method = 2",
                            "Using scattering method 2")), "Method-2 completion marker missing")
                require(digest(temperature) == task["temperature_sha256"], "Frozen temperature was modified during image execution")
                image = runner._one_product(output, "RT.fits.gz")
                measured = image_measurement(image, exp, geometry, photometry)
                measured.update(task_index=index, temperature_seed=task["temperature_seed"], image_seed=task["image_seed"],
                    temperature_sha256=task["temperature_sha256"], temperature_path=task["temperature_path"],
                    geometry=geometry["id"], image_npix=geometry["image_npix"], image_size_au=geometry["image_size_au"],
                    photon_packets=2048000, inclination_deg=50., threads=64, host=socket.gethostname(),
                    slurm_job=os.environ.get("SLURM_JOB_ID"), simulation_seconds=elapsed,
                    image_path=str(image.relative_to(bundle)), image_sha256=digest(image),
                    attempt_path=str(attempt.relative_to(bundle)), command_argv=command,
                    receipt_sha256={name: digest(attempt / name) for name in ("command.json", "execution.json", "mcfost.log")})
                runner.atomic_json(attempt / "measurement.json", measured)
                rows.append(measured)
                result.update(measurements=rows, updated_utc=runner.now())
                runner.atomic_json(result_path, result)
            result["complete"] = True
            runner.atomic_json(result_path, result)
            runner.atomic_json(status_path, {"state": "complete", "completed_images": 2, "updated_utc": runner.now()})
            return {"index": index, "status": "complete", "images": 2}
        except Exception as exc:
            runner.atomic_json(status_path, {"state": "failed", "error": f"{type(exc).__name__}: {exc}", "updated_utc": runner.now()})
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--validate", action="store_true")
    modes.add_argument("--status", action="store_true")
    modes.add_argument("--index", type=int)
    parser.add_argument("--machine", type=Path)
    args = parser.parse_args()
    if args.validate or args.status:
        exp = validate_package(args.bundle)
        result = {"experiment_id": NAME, "diagnostic_id": DIAGNOSTIC_ID,
                  "package_valid": True, "tasks": 5, "images": 10, "temperature_solves": 0}
        if args.status:
            rows = [inspect_task(args.bundle, exp, task, remeasure=False)[0] for task in exp["tasks"]]
            result.update(tasks=rows, counts=dict(Counter(r["state"] for r in rows)),
                          note="File status; use squeue/sacct for live scheduler state.")
    else:
        if args.machine is None:
            parser.error("--index requires --machine")
        result = run_task(args.bundle, args.index, args.machine)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
