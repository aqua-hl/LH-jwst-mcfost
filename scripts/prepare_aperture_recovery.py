#!/usr/bin/env python3
"""Fork a finished strict run into aperture_v2, preserving all original science files.

Existing FITS stay in the parent run and are referenced by relative symlinks.
Successful stricter measurements retain their exact fluxes with migration
provenance. Audited boundary failures are remeasured; only uncomputed anchors
need MCFOST after preparation. This command never submits or runs MCFOST.
"""
from __future__ import annotations

from contextlib import ExitStack
from copy import deepcopy
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from mcfost_grid.configuration import atomic_json, load_json, now, sha256
from mcfost_grid.photometry import KERNEL_ID, measure_image
from mcfost_grid.quality import APERTURE_POLICY, BOUNDARY_CHECKS, STRICT_POLICY, migrate_strict_measurement
from mcfost_grid.runner import _cached_measurements, validate_inputs


def check(ok, message):
    if not ok:
        raise ValueError(message)


def parent_product(parent, relative, model_id):
    path = Path(relative)
    prefix = Path("models") / model_id
    check(not path.is_absolute() and ".." not in path.parts and path.is_relative_to(prefix),
          f"Invalid source product path: {relative}")
    full = parent / path
    check(full.is_file() and full.resolve().is_relative_to(parent), f"Missing/out-of-run source product: {relative}")
    return full


def inherited_path(relative, model_id):
    prefix = Path("models") / model_id
    return str(prefix / "parent_products" / Path(relative).relative_to(prefix))


def checked_audit(parent, manifest, audit_path):
    audit = load_json(audit_path)
    check(audit.get("audit") == "read_only_aperture_crop_invariance_v1", "Unsupported audit format")
    check(audit.get("run_id") == manifest["run_id"] and audit.get("manifest_sha256") == sha256(parent/"manifest.json"),
          "Audit does not identify this exact parent manifest")
    check(audit.get("audit_script_sha256") == sha256(ROOT/"scripts/audit_aperture_boundaries.py"),
          "Audit implementation differs from the reviewed checker")
    check(audit.get("original_acceptance_and_results_unchanged") is True, "Audit did not preserve the parent acceptance")
    rows = {}
    for row in audit["models"]:
        index = row["index"]
        check(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(manifest["models"]),
              "Invalid audited index")
        check(index not in rows and row["model_id"] == manifest["models"][index]["id"], "Duplicate/mismatched audited model")
        crop = row["crop_comparison"]
        check(row.get("original_diagnostics_reproduced") is True
              and crop.get("aperture_crop_invariance_verified") is True
              and math.isfinite(crop["relative_difference"]) and 0 <= crop["relative_difference"] <= 1e-10
              and crop["maximum_absolute_aperture_weight_outside_crop"] == 0,
              f"Model {index}: aperture audit did not pass")
        full, direct = crop["full_image_adjoint_aperture_w_m2"], crop["central_crop_direct_aperture_w_m2"]
        check(math.isfinite(full) and full > 0 and math.isfinite(direct) and direct > 0
              and abs(full-direct)/full <= 1e-10, f"Model {index}: invalid audited aperture fluxes")
        failures = set(row["original_failed_checks"])
        check(failures and failures <= BOUNDARY_CHECKS, f"Model {index}: audit includes another failure")
        rows[index] = row
    check(rows, "Audit contains no models")
    return audit, rows


def remeasure_audited(parent, manifest, model, row):
    """Require the audited image, successful simulator receipts and unchanged flux."""
    anchor = next(a for a in manifest["anchors"] if a["id"] == row["anchor_id"])
    check(row["wavelength_um"] == anchor["wavelength_um"], "Audited wavelength differs from manifest")
    image = parent_product(parent, row["image"], model["id"])
    check(image.is_relative_to(parent/"models"/model["id"]/"anchors"/anchor["id"]), "Image belongs to another anchor")
    check(sha256(image) == row["image_sha256"], f"Audited FITS changed: {image}")
    parts = Path(row["image"]).parts
    attempt = parent / Path(*parts[:parts.index("output")])
    execution = load_json(attempt/"execution.json")
    command = load_json(attempt/"command.json")
    log = (attempt/"mcfost.log").read_text(errors="replace")
    check(execution["returncode"] == 0 and "Processing complete" in log
          and "fitsio error status" not in log.casefold()
          and not any(line.lstrip().startswith("ERROR:") for line in log.splitlines())
          and any(s in log for s in ("Using ray-tracing method 2", "RT source-function method = 2", "Using scattering method 2")),
          "Audited image lacks successful Method-2 completion records")
    argv = command["argv"]
    check("-img" in argv and "-rt2" in argv and "-no_T" in argv and "-Tfile" in argv
          and math.isclose(float(argv[argv.index("-img")+1]), anchor["wavelength_um"], rel_tol=1e-12),
          "Audited image command does not match the requested image calculation")
    start = time.monotonic()
    result = measure_image(image, wavelength_um=anchor["wavelength_um"], instrument=anchor["instrument"],
                           psf_fwhm_arcsec=anchor.get("psf_fwhm_arcsec"),
                           model_distance_pc=model["parameters"]["distance_pc"],
                           **{**manifest["measurement"], "quality_policy": APERTURE_POLICY})
    for key in ("aperture_flux_at_model_distance_w_m2", "plane0_edge_positive_flux_fraction",
                "plane0_convolution_flux_loss_fraction"):
        check(math.isclose(result["diagnostics"][key], row["original_diagnostics"][key], rel_tol=1e-9, abs_tol=0),
              f"Audited diagnostic no longer reproduces: {key}")
    check(set(result["quality_warnings"]) == set(row["original_failed_checks"]), "Audited boundary flags changed")
    check(sha256(image) == row["image_sha256"], "Audited image changed during remeasurement")
    result.update(anchor_id=anchor["id"], image_path=row["image"], image_sha256=row["image_sha256"],
                  simulation_seconds=execution["elapsed_seconds"], measurement_seconds=time.monotonic()-start,
                  backend="image_method2", source_acceptance={"method": "remeasured_existing_audited_image",
                  "source_manifest_sha256": sha256(parent/"manifest.json"),
                  "execution_sha256": sha256(attempt/"execution.json"), "command_sha256": sha256(attempt/"command.json"),
                  "log_sha256": sha256(attempt/"mcfost.log")})
    return result


def prepare_recovery(parent, audit_path, name):
    parent, audit_path = Path(parent).resolve(), Path(audit_path).resolve()
    check(isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9_-]+(?:\.[A-Za-z0-9_-]+)*", name), "Unsafe recovery name")
    target = parent.parent / name
    if target.exists():
        raise FileExistsError(f"Recovery already exists: {target}; resume it instead of preparing again")
    manifest = load_json(parent/"manifest.json")
    check(manifest.get("measurement", {}).get("quality_policy", STRICT_POLICY) == STRICT_POLICY,
          "Recovery requires the original strict policy")
    validate_inputs(parent, manifest)
    check(sha256(parent/"code/src/mcfost_grid/physics.py") == sha256(ROOT/"src/mcfost_grid/physics.py"),
          "Simulation-parameter implementation changed; checkpoint migration needs a separate review")
    original_hash = sha256(parent/"manifest.json")
    audit, audited = checked_audit(parent, manifest, audit_path)
    fingerprint = load_json(parent/"runtime_binding.json")["fingerprint"]
    check(fingerprint.get("manifest_sha256") == original_hash and fingerprint.get("backend") == "image_method2",
          "Parent runtime does not identify this image_method2 run")
    snapshots, records, temperatures, states = {}, {}, {}, {}
    with ExitStack() as locks:
        # No model can write its checkpoints while the source snapshot is taken.
        for model in manifest["models"]:
            lock = locks.enter_context((parent/"models"/model["id"]/".task.lock").open("a+"))
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError(f"Parent model {model['index']} is still running; wait for the parent array") from None
        for model in manifest["models"]:
            mid, index = model["id"], model["index"]
            directory = parent/"models"/mid
            status = load_json(directory/"status.json")
            check(status["state"] in {"complete", "failed"}, f"Model {index}: parent work is not finished")
            check((status["state"] == "failed") == (index in audited),
                  f"Model {index}: every failed model needs a matching audit, and audited models must still be failed")
            if index in audited:
                check(set(status.get("failed_checks", [])) == set(audited[index]["original_failed_checks"]),
                      f"Model {index}: failure changed after audit")
            binding = load_json(directory/"runtime_binding.json")
            temperature = load_json(directory/"temperature_complete.json")
            check(binding["fingerprint"] == temperature["fingerprint"] == fingerprint, f"Model {index}: mixed source runtime")
            path = parent_product(parent, temperature["path"], mid)
            check(sha256(path) == temperature["sha256"], f"Model {index}: source temperature changed")
            payload_path = directory/"measurements.json"
            payload = load_json(payload_path) if payload_path.exists() else {
                "model_id": mid, "model_index": index, "parameters": model["parameters"], "measurements": [],
                "complete": False, "fingerprint": fingerprint, "smoke_test": manifest["configuration"].get("smoke_test", False)}
            check(payload["fingerprint"] == fingerprint and payload["model_id"] == mid
                  and payload["model_index"] == index and payload["parameters"] == model["parameters"],
                  f"Model {index}: source measurement identity differs")
            check(payload.get("complete") is (status["state"] == "complete"), f"Model {index}: source completion flags differ")
            _cached_measurements(payload, manifest, parent)
            migrated = []
            for item in payload["measurements"]:
                parent_product(parent, item["image_path"], mid)
                check(item.get("kernel_id") == KERNEL_ID, "Source aperture estimator differs from the recovery estimator")
                changed = migrate_strict_measurement(item)
                changed["source_acceptance"] = {"method": "previously_passed_stricter_policy_flux_unchanged",
                                                 "source_manifest_sha256": original_hash,
                                                 "source_measurements_sha256": sha256(payload_path)}
                migrated.append(changed)
            if index in audited:
                check(audited[index]["anchor_id"] not in {m["anchor_id"] for m in migrated}, "Audited anchor is already cached")
                migrated.append(remeasure_audited(parent, manifest, model, audited[index]))
            new_payload = deepcopy(payload)
            new_payload["measurements"] = migrated
            new_payload["complete"] = len(migrated) == len(manifest["anchors"])
            records[mid], temperatures[mid], states[mid] = new_payload, temperature, status
            snapshots[mid] = {p.name: p.read_bytes() for p in
                              [directory/n for n in ("status.json", "runtime_binding.json", "temperature_complete.json", "measurements.json")]
                              if p.is_file()}
        check(sha256(parent/"manifest.json") == original_hash, "Parent manifest changed during preparation")
        stage = Path(tempfile.mkdtemp(prefix=f".{name}.prepare-", dir=parent.parent))
        try:
            for relative in manifest["input_hashes"]:
                dest = stage/relative
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(parent/relative, dest)
            shutil.rmtree(stage/"code")
            shutil.copytree(ROOT/"src/mcfost_grid", stage/"code/src/mcfost_grid",
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            shutil.copy2(ROOT/"workflow.py", stage/"code/workflow.py")
            config = deepcopy(manifest["configuration"])
            config["run_name"] = name
            config.setdefault("observations", {})["quality_policy"] = APERTURE_POLICY
            atomic_json(stage/"inputs/configuration.json", config)
            shutil.copy2(parent/"manifest.json", stage/"inputs/parent_manifest.json")
            shutil.copy2(audit_path, stage/"inputs/aperture_boundary_audit.json")
            shutil.copy2(Path(__file__), stage/"inputs/recovery_preparer.py")
            for mid, files in snapshots.items():
                folder = stage/"inputs/source_metadata"/mid
                folder.mkdir(parents=True)
                for filename, data in files.items():
                    (folder/filename).write_bytes(data)
            updated = deepcopy(manifest)
            updated.update(run_id=name, created_utc=now(), configuration=config,
                           measurement={**manifest["measurement"], "quality_policy": APERTURE_POLICY},
                           recovery={"parent_run_relative": "../"+parent.name, "parent_manifest_sha256": original_hash,
                                     "policy": APERTURE_POLICY, "audited_indices": sorted(audited),
                                     "simulation_inputs_unchanged": True,
                                     "source_products": "Relative links to parent; keep both run directories together"})
            updated["input_hashes"] = {str(p.relative_to(stage)): sha256(p) for p in sorted(stage.rglob("*")) if p.is_file()}
            atomic_json(stage/"manifest.json", updated)
            new_fingerprint = {**fingerprint, "manifest_sha256": sha256(stage/"manifest.json")}
            atomic_json(stage/"runtime_binding.json", {"fingerprint": new_fingerprint, "created_utc": now(),
                                                       "inherited_simulation_identity": fingerprint})
            missing = 0
            for model in manifest["models"]:
                mid = model["id"]
                directory = stage/"models"/mid
                (directory/"parent_products").symlink_to(os.path.relpath(parent/"models"/mid, directory), target_is_directory=True)
                temperature = deepcopy(temperatures[mid])
                temperature.update(fingerprint=new_fingerprint, path=inherited_path(temperature["path"], mid),
                                   inherited_from_manifest_sha256=original_hash)
                atomic_json(directory/"temperature_complete.json", temperature)
                binding = json.loads(snapshots[mid]["runtime_binding.json"])
                binding.update(fingerprint=new_fingerprint, inherited_from_manifest_sha256=original_hash)
                atomic_json(directory/"runtime_binding.json", binding)
                payload = records[mid]
                payload.update(fingerprint=new_fingerprint, updated_utc=now(), quality_policy=APERTURE_POLICY,
                               inherited_from_manifest_sha256=original_hash)
                for item in payload["measurements"]:
                    item["image_path"] = inherited_path(item["image_path"], mid)
                atomic_json(directory/"measurements.json", payload)
                count = len(payload["measurements"])
                missing += len(manifest["anchors"])-count
                atomic_json(directory/"status.json", {"state": "complete" if payload["complete"] else "not_started",
                    "model_id": mid, "completed_anchors": count, "total_anchors": len(manifest["anchors"]),
                    "updated_utc": now(), "recovery_policy": APERTURE_POLICY})
            receipt = {"run": str(target), "parent_run": str(parent), "parent_manifest_sha256": original_hash,
                       "manifest_sha256": new_fingerprint["manifest_sha256"], "models": len(manifest["models"]),
                       "audited_images_remeasured": len(audited), "remaining_image_calculations": missing,
                       "new_temperature_calculations": 0, "submitted": False,
                       "note": "Original science inputs/products untouched; keep parent files for relative product links"}
            atomic_json(stage/"recovery_receipt.json", receipt)
            validate_inputs(stage, updated)
            check(not target.exists(), f"Recovery destination appeared during preparation: {target}")
            stage.rename(target)
            return receipt
        except BaseException:
            if stage.exists():
                shutil.rmtree(stage)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--run-name", default="continuum_production_v1.1_aperture_v2")
    args = parser.parse_args()
    print(json.dumps(prepare_recovery(args.parent, args.audit, args.run_name), indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError, StopIteration) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
