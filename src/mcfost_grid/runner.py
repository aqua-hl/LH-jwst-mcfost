"""One independent physical model per task; resumable, bounded local execution."""
from __future__ import annotations

import fcntl
import json
import math
import os
import re
import shlex
import shutil
import signal
import socket
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .configuration import atomic_json, canonical_hash, load_json, now, resolve, sha256, _strict_keys
from .physics import numerical_args, physical_args
from .quality import STRICT_POLICY, cached_policy_error


def runtime_config(path):
    path = Path(path).resolve()
    value = load_json(path)
    _strict_keys(value, {"schema_version", "mcfost_executable", "mcfost_utils", "backend", "threads",
                         "timeout_seconds", "max_memory_gb", "slurm", "notes"}, "machine")
    if value.get("schema_version") != 1:
        raise ValueError("machine schema_version must be 1")
    executable = os.path.expandvars(value.get("mcfost_executable", os.environ.get("MCFOST_EXE", "mcfost")))
    found = shutil.which(executable) if "/" not in executable else str(resolve(path.parent, executable))
    if not found or not Path(found).is_file() or not os.access(found, os.X_OK):
        raise ValueError(f"MCFOST executable not available: {executable}")
    utils = value.get("mcfost_utils", os.environ.get("MCFOST_UTILS"))
    if not utils:
        raise ValueError("Set mcfost_utils in the machine JSON or MCFOST_UTILS")
    utils = resolve(path.parent, os.path.expandvars(utils))
    if not all((utils / name).is_dir() for name in ("Dust", "Lambda", "Stellar_Spectra")):
        raise ValueError(f"Incomplete MCFOST utilities: {utils}")
    backend = value.get("backend", "image_method2")
    if backend not in ("image_method2", "coeval_method2"):
        raise ValueError("backend must be image_method2 or coeval_method2")
    threads = value.get("threads", 2)
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    for name, default in (("timeout_seconds", 3600), ("max_memory_gb", 4)):
        number = value.get(name, default)
        if (isinstance(number, bool) or not isinstance(number, (int, float))
                or not math.isfinite(number) or number <= 0):
            raise ValueError(f"{name} must be positive")
    return {**value, "mcfost_executable": str(Path(found).resolve()), "mcfost_utils": str(utils),
            "backend": backend, "threads": threads, "timeout_seconds": value.get("timeout_seconds", 3600),
            "max_memory_gb": value.get("max_memory_gb", 4)}


def cpu_capacity():
    if os.environ.get("SLURM_CPUS_PER_TASK"):
        return int(os.environ["SLURM_CPUS_PER_TASK"])
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def environment(run, runtime):
    env = os.environ.copy()
    env.update(MCFOST_AUTO_UPDATE="0", MCFOST_UTILS=runtime["mcfost_utils"],
               MY_MCFOST_UTILS=str(run / "inputs/utils"), OMP_NUM_THREADS=str(runtime["threads"]),
               OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    return env


def check_backend(run, runtime):
    result = subprocess.run([runtime["mcfost_executable"], "-help"], env=environment(run, runtime),
                            capture_output=True, text=True, timeout=30)
    output = result.stdout + result.stderr
    if "usage : mcfost" not in output.lower() or "ERROR:" in output:
        raise RuntimeError("MCFOST startup failed: " + output[-2000:])
    if runtime["backend"] == "coeval_method2" and "-rt-sed-method" not in output:
        raise RuntimeError("This binary lacks -rt-sed-method. Select image_method2 or install the patched build; no fallback is automatic.")
    return output.splitlines()[:3]


def validate_inputs(run, manifest, model_id=None):
    """Check frozen shared inputs and, optionally, one model's parameter files.

    Preparation/configuration and independent full audits keep the default
    all-model check. Array workers can avoid repeatedly opening every other
    model's parameter files while still checking all shared executable inputs.
    Unsafe manifest paths are rejected even when their model is not selected.
    """
    run = Path(run).resolve()
    known_models = None
    if model_id is not None:
        if not isinstance(model_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", model_id):
            raise ValueError("Selected model ID must be a safe filename token")
        known_models = {model["id"] for model in manifest.get("models", [])}
        if model_id not in known_models:
            raise ValueError(f"Unknown selected model ID: {model_id}")
    for relative, digest in manifest["input_hashes"].items():
        relative_path = Path(relative)
        parts = relative_path.parts
        if relative_path.is_absolute() or ".." in parts or not parts:
            raise RuntimeError(f"Unsafe frozen run input path: {relative}")
        if model_id is not None and len(parts) >= 3 and parts[0] == "models":
            if parts[1] not in known_models:
                raise RuntimeError(f"Frozen input refers to an unknown model: {relative}")
            if parts[1] != model_id:
                continue
        path = run / relative
        if (not path.resolve().is_relative_to(run) or path.is_symlink()
                or not path.is_file() or sha256(path) != digest):
            raise RuntimeError(f"Frozen run input changed or missing: {relative}. Prepare a new run instead of mixing configurations.")


def _utilities_hash(runtime):
    """Bind cached products to utility bytes, not just a mutable directory name.

    The full utility tree is deliberately covered: automatic stellar-spectrum
    selection can use files other than the nominal 4000-K reference spectrum.
    This once-per-model read is small compared with a radiative-transfer task.
    """
    root = Path(runtime["mcfost_utils"])
    records = {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*")) if p.is_file()}
    if not records:
        raise RuntimeError("MCFOST utility directory contains no files")
    return canonical_hash(records)


def _cached_measurements(previous, manifest, run):
    """Validate resumable anchor identities and their preserved raw products."""
    anchors = {row["id"]: row for row in manifest["anchors"]}
    done = set()
    for item in previous.get("measurements", []):
        name = item.get("anchor_id")
        policy_error = cached_policy_error(item, manifest.get("measurement", {}).get("quality_policy", STRICT_POLICY))
        if policy_error:
            raise RuntimeError(f"Cached {name}: {policy_error}")
        if name not in anchors or name in done or item.get("quality_pass") is not True:
            raise RuntimeError("Cached measurements have duplicate, unknown or failed anchors")
        try:
            valid_wave = math.isclose(float(item["wavelength_um"]), anchors[name]["wavelength_um"], rel_tol=1e-10)
            valid_flux = math.isfinite(item["flux_jy"]) and item["flux_jy"] > 0
        except (KeyError, TypeError, ValueError):
            valid_wave = valid_flux = False
        if not valid_wave or not valid_flux:
            raise RuntimeError(f"Cached measurement has invalid wavelength/flux: {name}")
        for kind in ("image", "coeval_sed"):
            if kind == "coeval_sed" and previous["fingerprint"]["backend"] != "coeval_method2":
                continue
            path = item.get(kind + "_path")
            digest = item.get(kind + "_sha256")
            if (not path or not digest or Path(path).is_absolute() or ".." in Path(path).parts
                    or not (run / path).is_file() or sha256(run / path) != digest):
                raise RuntimeError(f"Cached {kind} changed or missing for {name}; cannot resume")
        done.add(name)
    if previous.get("complete") and done != set(anchors):
        raise RuntimeError("Complete cache does not contain every requested anchor")
    return done


def _check_run_runtime(run, fingerprint, commit=False):
    """One simulator/backend/utility identity across every model in the grid."""
    with (run / ".runtime.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        binding = run / "runtime_binding.json"
        if binding.exists():
            if load_json(binding).get("fingerprint") != fingerprint:
                raise RuntimeError("Grid runtime or manifest changed; all models must use the same binary, backend and utilities")
        elif commit:
            atomic_json(binding, {"fingerprint": fingerprint, "created_utc": now()})


def _invoke(command, cwd, env, timeout, attempt):
    log = attempt / "mcfost.log"
    started = time.monotonic()
    atomic_json(attempt / "command.json", {"argv": command, "cwd": str(cwd), "started_utc": now(),
                                          "threads": env["OMP_NUM_THREADS"]})
    with log.open("wb") as stream:
        proc = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            returncode = proc.wait(timeout=timeout)
        except BaseException:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
            raise
    elapsed = time.monotonic() - started
    atomic_json(attempt / "execution.json", {"returncode": returncode, "elapsed_seconds": elapsed, "finished_utc": now()})
    text = log.read_text(errors="replace")
    if (returncode != 0 or "Processing complete" not in text
            or "fitsio error status" in text.casefold()
            or any(line.lstrip().startswith("ERROR:") for line in text.splitlines())):
        raise RuntimeError(f"MCFOST failed or incomplete (exit {returncode}); see {log}")
    return elapsed, text


def _new_attempt(parent):
    parent.mkdir(parents=True, exist_ok=True)
    index = 1
    while True:
        directory = parent / f"attempt_{index:03d}"
        try:
            directory.mkdir()
            return directory
        except FileExistsError:
            index += 1


def _one_product(root, name):
    found = list(root.rglob(name))
    if len(found) != 1:
        raise RuntimeError(f"Expected exactly one {name} in {root}; found {len(found)}")
    return found[0]


def _valid_temperature(path):
    import numpy as np
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdus:
        data = np.asarray(hdus[0].data)
        if data.size == 0 or not np.all(np.isfinite(data)) or not np.any(data > 0) or np.any(data < 0):
            raise RuntimeError("Temperature FITS is empty, negative, or nonfinite")


def run_model(run_dir, index, machine_path):
    from .photometry import measure_image
    run = Path(run_dir).resolve()
    manifest = load_json(run / "manifest.json")
    if index < 0 or index >= len(manifest["models"]):
        raise ValueError(f"Task index must be 0..{len(manifest['models']) - 1}")
    runtime = runtime_config(machine_path)
    if runtime["threads"] > cpu_capacity():
        raise ValueError(f"Requested {runtime['threads']} threads exceeds allocated/local CPU count {cpu_capacity()}")
    model = manifest["models"][index]
    validate_inputs(run, manifest, model_id=model["id"])
    directory = run / "models" / model["id"]
    with (directory / ".task.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"Model is already running: {model['id']}")
        manifest_hash = sha256(run / "manifest.json")
        fingerprint = {"manifest_sha256": manifest_hash, "executable_sha256": sha256(runtime["mcfost_executable"]),
                       "utilities_content_sha256": _utilities_hash(runtime), "backend": runtime["backend"]}
        _check_run_runtime(run, fingerprint)
        binding = directory / "runtime_binding.json"
        if binding.exists() and load_json(binding)["fingerprint"] != fingerprint:
            raise RuntimeError("Runtime or manifest changed; use a new run_name to avoid mixing builds/backends")
        status_path = directory / "status.json"
        result_path = directory / "measurements.json"
        temp_receipt = directory / "temperature_complete.json"
        if temp_receipt.exists():
            receipt = load_json(temp_receipt)
            temperature = run / receipt["path"]
            if (receipt["fingerprint"] != fingerprint or not temperature.is_file()
                    or sha256(temperature) != receipt["sha256"]):
                raise RuntimeError("Cached temperature changed; cannot resume")
        if result_path.exists():
            previous = load_json(result_path)
            if previous.get("fingerprint") != fingerprint:
                raise RuntimeError("Cached measurements have different inputs/runtime")
            _cached_measurements(previous, manifest, run)
            if not temp_receipt.exists():
                raise RuntimeError("Cached measurements lack their temperature receipt")
            if previous.get("complete"):
                atomic_json(status_path, {"state": "complete", "updated_utc": now(), "model_id": model["id"],
                                          "completed_anchors": len(previous["measurements"]),
                                          "total_anchors": len(manifest["anchors"]), "validated_cached_result": True})
                return {"model_id": model["id"], "status": "cached", "index": index}
        else:
            previous = {"model_id": model["id"], "model_index": index, "parameters": model["parameters"],
                        "fingerprint": fingerprint, "complete": False, "measurements": [],
                        "smoke_test": manifest["configuration"].get("smoke_test", False)}
        banner = check_backend(run, runtime)
        _check_run_runtime(run, fingerprint, commit=True)
        if not binding.exists():
            atomic_json(binding, {"fingerprint": fingerprint, "runtime": runtime, "banner": banner,
                                  "host": socket.gethostname(), "created_utc": now()})
        env = environment(run, runtime)
        args = [*physical_args(model["parameters"]),
                *numerical_args(manifest["configuration"].get("numerics", {}))]
        common = ["-max_mem", str(runtime["max_memory_gb"]), "-no_backup"]
        atomic_json(status_path, {"state": "running", "stage": "temperature", "started_utc": now(),
                                  "model_id": model["id"], "pid": os.getpid(), "host": socket.gethostname(),
                                  "slurm_job": os.environ.get("SLURM_JOB_ID"), "threads": runtime["threads"]})
        try:
            if temp_receipt.exists():
                receipt = load_json(temp_receipt)
                temperature = run / receipt["path"]
                if receipt["fingerprint"] != fingerprint or sha256(temperature) != receipt["sha256"]:
                    raise RuntimeError("Cached temperature changed; cannot resume")
            else:
                attempt = _new_attempt(directory / "temperature")
                output = attempt / "output"
                # MCFOST 4.1 prefixes root_dir with './' internally. Absolute
                # output paths can therefore silently produce FITSIO failures.
                command = [runtime["mcfost_executable"], "temperature.para", *args,
                           "-root_dir", os.path.relpath(output, directory), *common]
                elapsed, _ = _invoke(command, directory, env, runtime["timeout_seconds"], attempt)
                temperature = _one_product(output, "Temperature.fits.gz")
                _valid_temperature(temperature)
                atomic_json(temp_receipt, {"fingerprint": fingerprint, "path": str(temperature.relative_to(run)),
                                          "sha256": sha256(temperature), "elapsed_seconds": elapsed})
            done = {item["anchor_id"] for item in previous["measurements"] if item.get("quality_pass")}
            for anchor in manifest["anchors"]:
                if anchor["id"] in done:
                    continue
                atomic_json(status_path, {"state": "running", "stage": anchor["id"], "updated_utc": now(),
                                          "completed_anchors": len(done), "total_anchors": len(manifest["anchors"]),
                                          "model_id": model["id"], "host": socket.gethostname(), "pid": os.getpid()})
                work = directory / "anchors" / anchor["id"]
                attempt = _new_attempt(work / "attempts")
                output = attempt / "output"
                # MCFOST resolves -Tfile underneath root_dir, not underneath
                # cwd. A relative link keeps this attempt portable and avoids
                # duplicating the model's equilibrium-temperature product.
                output.mkdir()
                (output / "Temperature.fits.gz").symlink_to(os.path.relpath(temperature, output))
                if runtime["backend"] == "image_method2":
                    command = [runtime["mcfost_executable"], "image.para", "-img", str(anchor["wavelength_um"])]
                else:
                    command = [runtime["mcfost_executable"], "coeval.para", "-rt-sed-method", "2"]
                command += [*args, "-rt2", "-no_T", "-Tfile", "Temperature.fits.gz",
                            "-root_dir", os.path.relpath(output, work), *common]
                elapsed, log = _invoke(command, work, env, runtime["timeout_seconds"], attempt)
                if not any(marker in log for marker in ("Using ray-tracing method 2", "RT source-function method = 2",
                                                       "Using scattering method 2")):
                    raise RuntimeError(f"Method-2 completion marker missing: {attempt / 'mcfost.log'}")
                image = _one_product(output, "RT.fits.gz")
                sed = _one_product(output, "sed_rt.fits.gz") if runtime["backend"] == "coeval_method2" else None
                start = time.monotonic()
                measured = measure_image(image, wavelength_um=anchor["wavelength_um"], instrument=anchor["instrument"],
                                         model_distance_pc=model["parameters"]["distance_pc"],
                                         coeval_sed_path=sed, psf_fwhm_arcsec=anchor.get("psf_fwhm_arcsec"),
                                         **manifest["measurement"])
                measured.update(anchor_id=anchor["id"], wavelength_um=anchor["wavelength_um"],
                                simulation_seconds=elapsed, measurement_seconds=time.monotonic() - start,
                                image_path=str(image.relative_to(run)), image_sha256=sha256(image),
                                backend=runtime["backend"])
                if sed is not None:
                    measured.update(coeval_sed_path=str(sed.relative_to(run)), coeval_sed_sha256=sha256(sed))
                if not measured.get("quality_pass"):
                    raise RuntimeError(f"Image photometry quality failed: {anchor['id']}")
                atomic_json(attempt / "measurement.json", measured)
                previous["measurements"] = [m for m in previous["measurements"] if m["anchor_id"] != anchor["id"]] + [measured]
                done.add(anchor["id"])
                previous["updated_utc"] = now()
                atomic_json(result_path, previous)
            previous["complete"] = True
            atomic_json(result_path, previous)
            atomic_json(status_path, {"state": "complete", "updated_utc": now(), "model_id": model["id"],
                                      "completed_anchors": len(done), "total_anchors": len(manifest["anchors"])})
            return {"model_id": model["id"], "index": index, "status": "complete"}
        except BaseException as exc:
            atomic_json(status_path, {"state": "failed", "updated_utc": now(), "model_id": model["id"],
                                      "error": f"{type(exc).__name__}: {exc}",
                                      "failed_checks": list(getattr(exc, "failed_checks", [])),
                                      "diagnostics": getattr(exc, "diagnostics", {})})
            raise


def run_local(run, machine_path, workers=1):
    runtime = runtime_config(machine_path)
    if workers < 1 or workers * runtime["threads"] > cpu_capacity():
        raise ValueError("workers × threads exceeds available CPUs (or workers is not positive)")
    manifest = load_json(Path(run) / "manifest.json")
    start = time.monotonic()
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_model, str(run), item["index"], str(machine_path)): item for item in manifest["models"]}
        for future in as_completed(futures):
            item = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"model_id": item["id"], "index": item["index"], "status": "failed", "error": str(exc)}
            results.append(result)
            print(json.dumps(result), flush=True)
    return {"models": results, "elapsed_seconds": time.monotonic() - start,
            "failed": sum(item["status"] == "failed" for item in results)}


def status(run):
    run = Path(run)
    manifest = load_json(run / "manifest.json")
    counts, models = {}, []
    for item in manifest["models"]:
        path = run / "models" / item["id"] / "status.json"
        entry = load_json(path) if path.exists() else {"state": "not_started"}
        entry["model_id"] = item["id"]
        entry["index"] = item["index"]
        counts[entry["state"]] = counts.get(entry["state"], 0) + 1
        models.append(entry)
    return {"run": manifest["run_id"], "models_total": len(models), "counts": counts, "models": models,
            "note": "File-based status, not a live scheduler query; interrupted jobs can leave a running marker. Safe rerun acquires a per-model lock."}
