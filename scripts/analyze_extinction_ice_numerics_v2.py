#!/usr/bin/env python3
"""Audit the 20-task V2 numerical experiment; never fit observations or rank models.

This intentionally requires the complete downloaded bundle, including raw FITS
products. A status marker alone is insufficient evidence of a valid calculation.
Five independent seed draws are used for Gaussian sample-SD intervals; matching
seed labels pair the two photon budgets without assuming variance cancellation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import importlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import chi2, t


SEEDS = tuple(range(42001, 42006))
BUDGETS = (512000, 2048000)
INCLINATIONS = (50., 70.)
IDENTITY_KEYS = ("executable_sha256", "utilities_content_sha256", "backend")


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {value}")))


def confined(root, relative):
    """Permit only ordinary files under this portable bundle, never external paths."""
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path resolves outside bundle: {relative}")
    return path


def check_hashes(root, hashes):
    for relative, expected in hashes.items():
        path = confined(root, relative)
        if path.is_symlink() or not path.is_file() or digest(path) != expected:
            raise ValueError(f"Frozen input changed or missing: {relative}")


def _frozen_runner(run):
    package = run / "code/src/mcfost_grid"
    if not (package / "__init__.py").is_file():
        raise ValueError("Missing frozen workflow code/src/mcfost_grid")
    # Import under a private package name: an unrelated installed package or a
    # repository checkout must not silently replace the bundle's frozen code.
    name = "_v2_numerical_frozen_" + hashlib.sha256(str(package).encode()).hexdigest()[:16]
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, package / "__init__.py",
                                                      submodule_search_locations=[str(package)])
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return importlib.import_module(name + ".runner")


def validate_experiment(experiment):
    if experiment.get("schema_version") != 1 or experiment.get("experiment_id") != "extinction_ice_numerics_v2":
        raise ValueError("Expected extinction_ice_numerics_v2 experiment schema_version=1")
    if not isinstance(experiment.get("input_hashes"), dict) or not experiment["input_hashes"]:
        raise ValueError("Experiment must contain nonempty frozen input_hashes")
    tasks, probes = experiment["tasks"], experiment["probes"]
    expected = {(inclination, budget, seed) for inclination in INCLINATIONS
                for budget in BUDGETS for seed in SEEDS}
    actual = {(task["inclination_deg"], task["photon_packets"], task["seed"]) for task in tasks}
    if len(tasks) != 20 or actual != expected or {task["index"] for task in tasks} != set(range(20)):
        raise ValueError("Task catalogue must contain each of the two inclinations, two budgets and five seeds once")
    if len({(task["run_path"], task["model_id"]) for task in tasks}) != 20:
        raise ValueError("Duplicate task run/model destination")
    if len(probes) != 6 or {probe["id"] for probe in probes} != {f"n{i:03d}" for i in range(1, 7)}:
        raise ValueError("Expected six unique probes n001..n006")
    for probe in probes:
        if not math.isfinite(probe["wavelength_um"]) or probe["wavelength_um"] <= 0:
            raise ValueError("Probe wavelengths must be finite and positive")
    acceptance = experiment["acceptance"]
    for key in ("scatter_fraction", "budget_change_fraction", "confidence"):
        if not 0 < acceptance[key] < 1:
            raise ValueError(f"Invalid acceptance.{key}")
    threads = experiment.get("resources", {}).get("cpus_per_task", 64)
    if threads != 64:
        raise ValueError("This experiment requires 64 threads per task")


def _command_receipt(run, product, task, threads, source_hashes):
    attempt = next((parent for parent in product.parents
                    if parent.is_relative_to(run) and (parent / "command.json").is_file()), None)
    if attempt is None:
        raise ValueError(f"Missing command receipt for {product.relative_to(run)}")
    command_path, execution_path = attempt / "command.json", attempt / "execution.json"
    command, execution = read_json(command_path), read_json(execution_path)
    for path in (command_path, execution_path):
        source_hashes[str(path.relative_to(run))] = digest(path)
    if int(command.get("threads", 0)) != threads:
        raise ValueError(f"Wrong thread count in {command_path.relative_to(run)}")
    args = command.get("argv", [])
    if args.count("-seed") != 1 or args.index("-seed") + 1 >= len(args):
        raise ValueError("Command lacks one explicit seed")
    if str(args[args.index("-seed") + 1]) != str(task["seed"]):
        raise ValueError("Command seed differs from experiment")
    if execution.get("returncode") != 0:
        raise ValueError("MCFOST execution receipt is not successful")
    return args


def _inspect_task(bundle, experiment, task, run_cache):
    """Return task status and flux rows; failures remain explicit and unscored."""
    status = {key: task[key] for key in ("index", "run_path", "model_index", "model_id",
                                       "inclination_deg", "photon_packets", "seed")}
    status.update(state="invalid", reason=None, source_hashes={})
    rows = []
    try:
        run = confined(bundle, task["run_path"])
        if task["run_path"] not in run_cache:
            try:
                manifest_path = run / "manifest.json"
                manifest = read_json(manifest_path)
                manifest_hash = digest(manifest_path)
                if task.get("manifest_sha256") and manifest_hash != task["manifest_sha256"]:
                    raise ValueError("Subrun manifest differs from experiment's pinned hash")
                check_hashes(run, manifest["input_hashes"])
                runner = _frozen_runner(run)
                runner.validate_inputs(run, manifest)
                code_identity = {key: value for key, value in manifest["input_hashes"].items()
                                 if key.startswith("code/")}
                if not code_identity:
                    raise ValueError("Subrun does not pin its workflow code")
                run_cache[task["run_path"]] = (manifest, manifest_hash, runner, code_identity)
            except Exception as exc:
                run_cache[task["run_path"]] = exc
        entry = run_cache[task["run_path"]]
        if isinstance(entry, Exception):
            raise ValueError(str(entry))
        manifest, manifest_hash, runner, code_identity = entry
        if task.get("manifest_sha256") and manifest_hash != task["manifest_sha256"]:
            raise ValueError("Subrun manifest differs from experiment's pinned hash")
        status["manifest_sha256"] = manifest_hash
        status["frozen_code_identity"] = code_identity
        if manifest["configuration"].get("numerical_only") is not True:
            raise ValueError("Subrun must explicitly declare numerical_only=true")
        numerics = manifest["configuration"]["numerics"]
        if numerics.get("random_seed") != task["seed"]:
            raise ValueError("Frozen numerical seed differs from experiment")
        if any(numerics[key] != task["photon_packets"] for key in
               ("photons_temperature", "photons_image", "photons_sed")):
            raise ValueError("Frozen photon budgets differ from experiment")
        model = manifest["models"][task["model_index"]]
        if (model["id"] != task["model_id"] or model["index"] != task["model_index"]
                or model["parameters"]["inclination_deg"] != task["inclination_deg"]):
            raise ValueError("Frozen model identity differs from task catalogue")
        expected_probes = {p["id"]: p["wavelength_um"] for p in experiment["probes"]}
        anchors = manifest["anchors"]
        if len(anchors) != 6 or {a["id"]: a["wavelength_um"] for a in anchors} != expected_probes:
            raise ValueError("Frozen anchors do not exactly match the six probes")
        if any(anchor.get("score", True) for anchor in anchors):
            raise ValueError("Numerical probes must have score=false")
        directory = confined(run, f"models/{model['id']}")
        status_path, measurement_path = directory / "status.json", directory / "measurements.json"
        marker = read_json(status_path) if status_path.is_file() else {}
        if status_path.is_file():
            status["source_hashes"][str(status_path.relative_to(run))] = digest(status_path)
        status["recorded_state"] = marker.get("state", "not_started")
        if not measurement_path.is_file():
            status.update(state=marker.get("state", "not_started"), reason="measurements.json is absent")
            if status["state"] == "complete":
                status["state"] = "invalid"
            return status, rows
        measured = read_json(measurement_path)
        status["source_hashes"][str(measurement_path.relative_to(run))] = digest(measurement_path)
        if measured.get("complete") is not True:
            status.update(state="incomplete", reason="Measurements are incomplete")
            return status, rows
        if (measured.get("model_id") != model["id"] or measured.get("model_index") != model["index"]
                or measured.get("parameters") != model["parameters"]):
            raise ValueError("Measurements do not match frozen model identity/parameters")
        fingerprint = measured["fingerprint"]
        if fingerprint.get("manifest_sha256") != manifest_hash:
            raise ValueError("Measurement fingerprint has a different manifest")
        for key in IDENTITY_KEYS[:2]:
            if not isinstance(fingerprint.get(key), str) or len(fingerprint[key]) != 64:
                raise ValueError(f"Invalid runtime fingerprint {key}")
        if fingerprint.get("backend") != "image_method2":
            raise ValueError("Numerical experiment requires the image_method2 backend")
        threads = experiment.get("resources", {}).get("cpus_per_task", 64)
        for binding_path in (run / "runtime_binding.json", directory / "runtime_binding.json"):
            binding = read_json(binding_path)
            status["source_hashes"][str(binding_path.relative_to(run))] = digest(binding_path)
            if binding.get("fingerprint") != fingerprint:
                raise ValueError("Runtime binding does not match measurements")
            if binding_path.parent == directory and binding.get("runtime", {}).get("threads") != threads:
                raise ValueError("Model runtime binding has the wrong thread count")
        receipt_path = directory / "temperature_complete.json"
        receipt = read_json(receipt_path)
        status["source_hashes"][str(receipt_path.relative_to(run))] = digest(receipt_path)
        temperature = confined(run, receipt["path"])
        if (receipt.get("fingerprint") != fingerprint or not temperature.is_file()
                or temperature.is_symlink() or digest(temperature) != receipt["sha256"]):
            raise ValueError("Temperature receipt/fingerprint/raw hash mismatch")
        runner._valid_temperature(temperature)
        _command_receipt(run, temperature, task, threads, status["source_hashes"])
        runner._cached_measurements(measured, manifest, run)
        for item in measured["measurements"]:
            if item["wavelength_um"] != expected_probes[item["anchor_id"]]:
                raise ValueError("Measurement wavelength is not exactly the frozen probe")
            if item.get("backend") != fingerprint["backend"]:
                raise ValueError("Measurement backend differs from fingerprint")
            image = confined(run, item["image_path"])
            args = _command_receipt(run, image, task, threads, status["source_hashes"])
            if fingerprint["backend"] == "image_method2":
                if "-img" not in args or float(args[args.index("-img") + 1]) != item["wavelength_um"]:
                    raise ValueError("Image command wavelength differs from frozen probe")
            if "-rt2" not in args or "-no_T" not in args:
                raise ValueError("Image command lacks method-2/reused-temperature flags")
            rows.append({key: task[key] for key in ("index", "run_path", "model_id", "inclination_deg",
                                                   "photon_packets", "seed")}
                        | {"probe_id": item["anchor_id"], "wavelength_um": item["wavelength_um"],
                           "flux_jy": item["flux_jy"], "image_path": item["image_path"],
                           "image_sha256": item["image_sha256"], "temperature_sha256": receipt["sha256"]})
        status.update(state="complete", reason=None,
                      runtime_identity={key: fingerprint[key] for key in IDENTITY_KEYS} | {"threads": threads},
                      temperature_sha256=receipt["sha256"], completed_probes=len(rows))
        return status, rows
    except (OSError, ValueError, TypeError, KeyError, IndexError, RuntimeError) as exc:
        status.update(state="invalid", reason=f"{type(exc).__name__}: {exc}")
        return status, []


def scatter_statistics(values, confidence=.95):
    """Two-sided Gaussian SD confidence interval, normalised by the sample mean.

    The mean denominator is a plug-in estimate; these are not exact confidence
    intervals for the coefficient of variation. Five draws do not test normality.
    """
    values = np.asarray(values, dtype=float)
    count = len(values)
    if count < 2 or np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Scatter statistics require at least two finite positive values")
    identical = bool(np.ptp(values) == 0)
    mean = float(np.mean(values))
    sample_sd = 0. if identical else float(np.std(values, ddof=1))
    alpha, degrees = 1 - confidence, count - 1
    lower = sample_sd * math.sqrt(degrees / chi2.ppf(1 - alpha / 2, degrees))
    upper = sample_sd * math.sqrt(degrees / chi2.ppf(alpha / 2, degrees))
    return {"n_seeds": count, "mean_flux_jy": mean, "sample_sd_jy": sample_sd,
            "sample_sd_fraction": sample_sd / mean, "sd_lower_fraction": lower / mean,
            "sd_upper_fraction": upper / mean, "scatter_interval_degenerate": identical or sample_sd == 0}


def paired_statistics(low, high, confidence=.95):
    """Pair on explicit seed IDs; evaluate d_seed=(F_high-F_low)/F_low."""
    if set(low) != set(high) or len(low) < 2:
        raise ValueError("Paired budgets require the same seed IDs and at least two pairs")
    seeds = sorted(low)
    low_values, high_values = np.array([low[s] for s in seeds]), np.array([high[s] for s in seeds])
    if (np.any(~np.isfinite(low_values)) or np.any(~np.isfinite(high_values))
            or np.any(low_values <= 0) or np.any(high_values <= 0)):
        raise ValueError("Paired fluxes must be finite and positive")
    differences = (high_values - low_values) / low_values
    mean = float(np.mean(differences))
    sd = float(np.std(differences, ddof=1))
    half_width = float(t.ppf((1 + confidence) / 2, len(seeds) - 1) * sd / math.sqrt(len(seeds)))
    return {"n_pairs": len(seeds), "mean_fractional_change": mean,
            "sample_sd_fractional_change": sd, "ci_half_width": half_width,
            "ci_lower": mean - half_width, "ci_upper": mean + half_width,
            "absolute_change_upper_bound": abs(mean) + half_width}


def numerical_statistics(experiment, rows, eligible):
    scatter, changes = [], []
    confidence = experiment["acceptance"]["confidence"]
    for inclination in INCLINATIONS:
        for probe in experiment["probes"]:
            groups = {}
            base = {"inclination_deg": inclination, "probe_id": probe["id"],
                    "wavelength_um": probe["wavelength_um"]}
            for budget in BUDGETS:
                selected = {r["seed"]: r["flux_jy"] for r in rows if r["inclination_deg"] == inclination
                            and r["photon_packets"] == budget and r["probe_id"] == probe["id"]}
                groups[budget] = selected
                result = base | {"photon_packets": budget, "n_seeds": len(selected),
                                 "required_for_gate": budget == max(BUDGETS),
                                 "eligible_for_gate": eligible and set(selected) == set(SEEDS),
                                 "passes_scatter": None}
                if len(selected) >= 2:
                    result.update(scatter_statistics(list(selected.values()), confidence))
                if result["eligible_for_gate"] and not result["scatter_interval_degenerate"]:
                    result["passes_scatter"] = result["sd_upper_fraction"] <= experiment["acceptance"]["scatter_fraction"]
                scatter.append(result)
            paired_seeds = set(groups[BUDGETS[0]]) & set(groups[BUDGETS[1]])
            result = base | {"low_photon_packets": BUDGETS[0], "high_photon_packets": BUDGETS[1],
                             "n_pairs": len(paired_seeds), "eligible_for_gate": eligible and paired_seeds == set(SEEDS),
                             "passes_budget_change": None}
            if len(paired_seeds) >= 2:
                result.update(paired_statistics({s: groups[BUDGETS[0]][s] for s in paired_seeds},
                                                {s: groups[BUDGETS[1]][s] for s in paired_seeds}, confidence))
            if result["eligible_for_gate"]:
                result["passes_budget_change"] = result["absolute_change_upper_bound"] <= experiment["acceptance"]["budget_change_fraction"]
            changes.append(result)
    return scatter, changes


def _write_csv(path, records, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def plot_convergence(output, summary, scatter, changes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5), constrained_layout=True)
    for inclination, color in zip(INCLINATIONS, ("#1767a6", "#ba4c22")):
        for budget, style in zip(BUDGETS, ("--", "-")):
            subset = [r for r in scatter if r["inclination_deg"] == inclination
                      and r["photon_packets"] == budget and "sd_upper_fraction" in r]
            axes[0].plot([r["wavelength_um"] for r in subset], [100 * r["sd_upper_fraction"] for r in subset],
                         marker="o", linestyle=style, color=color,
                         label=f"{inclination:g}° / {budget / 1e6:g}M packets")
        subset = [r for r in changes if r["inclination_deg"] == inclination and "mean_fractional_change" in r]
        axes[1].errorbar([r["wavelength_um"] for r in subset], [100 * r["mean_fractional_change"] for r in subset],
                         yerr=[100 * r["ci_half_width"] for r in subset], fmt="o-", capsize=3,
                         color=color, label=f"{inclination:g}°")
    axes[0].axhline(100 * summary["acceptance"]["scatter_fraction"], color="black", ls=":", label="1% target")
    limit = 100 * summary["acceptance"]["budget_change_fraction"]
    axes[1].axhspan(-limit, limit, color="grey", alpha=.15)
    axes[1].axhline(0, color="grey", linewidth=.7)
    axes[0].set_ylabel("Upper 95% SD interval / seed mean (%)")
    axes[1].set_ylabel("Paired photon-budget change with 95% CI (%)")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Probe wavelength (µm)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.18)
    fig.suptitle(f"V2 numerical experiment: {summary['status']} — {summary['complete_tasks']}/20 tasks validated")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"numerical_convergence.{suffix}", dpi=170)
    plt.close(fig)


def analyze_bundle(bundle, output=None, make_plot=True):
    bundle = Path(bundle).resolve()
    experiment_path = bundle / "experiment.json"
    experiment = read_json(experiment_path)
    validate_experiment(experiment)
    output = Path(output).resolve() if output else bundle / "results"
    global_errors = []
    try:
        checksum = (bundle / "experiment.sha256").read_text().split()
        if len(checksum) != 2 or checksum[1] != "experiment.json" or checksum[0] != digest(experiment_path):
            raise ValueError("Experiment hash differs from experiment.sha256")
    except (OSError, ValueError) as exc:
        global_errors.append(str(exc))
    try:
        check_hashes(bundle, experiment.get("input_hashes", {}))
    except (OSError, ValueError) as exc:
        global_errors.append(str(exc))
    cache, statuses, rows = {}, [], []
    for task in sorted(experiment["tasks"], key=lambda item: item["index"]):
        status, task_rows = _inspect_task(bundle, experiment, task, cache)
        statuses.append(status)
        rows.extend(task_rows)
    identities = {json.dumps(item["runtime_identity"], sort_keys=True) for item in statuses if item["state"] == "complete"}
    code_identities = {json.dumps(item["frozen_code_identity"], sort_keys=True) for item in statuses
                       if "frozen_code_identity" in item}
    if len(identities) > 1:
        global_errors.append("Mixed executable, utilities, backend or thread identity across subruns; inference refused")
    if len(code_identities) > 1:
        global_errors.append("Mixed frozen workflow code across subruns; inference refused")
    # The array launcher must establish one bundle identity before any task.
    # An entirely unstarted bundle can still generate an incomplete report.
    bundle_binding = bundle / "runtime_binding.json"
    if bundle_binding.is_file():
        binding = read_json(bundle_binding)
        if binding.get("experiment_sha256") != digest(experiment_path):
            global_errors.append("Bundle runtime binding has a different experiment hash")
        bound_identity = binding.get("runtime_identity", binding.get("identity", binding.get("fingerprint")))
        if bound_identity is None:
            global_errors.append("Bundle runtime binding does not contain an identity")
        else:
            projected = {key: bound_identity.get(key) for key in (*IDENTITY_KEYS, "threads")}
            if identities and identities != {json.dumps(projected, sort_keys=True)}:
                global_errors.append("Bundle runtime binding differs from completed tasks")
    elif any(item["state"] == "complete" or item.get("recorded_state") == "complete" for item in statuses):
        global_errors.append("Bundle runtime_binding.json is required when completed tasks exist")
    count_complete = sum(item["state"] == "complete" for item in statuses)
    eligible = count_complete == 20 and not global_errors
    scatter, changes = numerical_statistics(experiment, rows, eligible)
    passed = bool(eligible and all(row["passes_scatter"] for row in scatter if row["required_for_gate"])
                  and all(row["passes_budget_change"] for row in changes))
    degenerate = any(row.get("scatter_interval_degenerate") for row in scatter if row["required_for_gate"])
    invalid = any(item["state"] == "invalid" for item in statuses)
    summary = {
        "schema_version": 1, "experiment_id": experiment["experiment_id"],
        "analysis": "numerical_convergence_only", "experiment_sha256": digest(experiment_path),
        "analyzer_sha256": digest(Path(__file__)), "total_tasks": 20, "complete_tasks": count_complete,
        "task_state_counts": dict(Counter(item["state"] for item in statuses)),
        "status": "integrity_error" if global_errors or invalid else ("incomplete" if not eligible else
                  ("passed" if passed else ("inconclusive" if degenerate else "not_passed"))),
        "numerical_gate_passed": passed, "inference_eligible": eligible,
        "acceptance": experiment["acceptance"], "errors": global_errors, "tasks": statuses,
        "runtime_identities": [json.loads(value) for value in sorted(identities)],
        "scatter": scatter, "budget_changes": changes,
        "assumptions": {
            "scatter_interval": "Two-sided Gaussian chi-square confidence interval for sample SD with n-1 degrees of freedom, divided by sample mean as a plug-in normalisation; not an exact coefficient-of-variation interval.",
            "budget_change": "Per seed d=(F_2048000-F_512000)/F_512000; report mean d and two-sided Student-t confidence interval with n-1 degrees of freedom. Match seed labels; do not assume variance cancellation.",
            "independent_draws": "Five separately simulated temperature and image solutions per budget/case. Five seeds do not establish Gaussianity or rule out rare outliers.",
            "zero_scatter": "Five identical fluxes produce a degenerate interval and are inconclusive, not a pass; check that the simulator uses independent seed streams before interpreting zero observed variance.",
            "multiple_tests": "Intervals are marginal per probe/case, not simultaneous familywise confidence intervals.",
            "scope": "Only these two physical cases and six monochromatic apertures; no extrapolation to all materials, broad-band quadrature, full spectra or production geometry.",
        },
        "physical_contrast_gate": "not_assessed", "observation_fit_performed": False,
        "model_ranking_performed": False, "production_authorized": False,
        "file_status_note": "File-based product validation, not a live scheduler query. Incomplete, mixed-runtime or corrupted tasks cannot pass.",
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "numerical_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    _write_csv(output / "seed_fluxes.csv", rows,
               ["index", "run_path", "model_id", "inclination_deg", "photon_packets", "seed", "probe_id", "wavelength_um",
                "flux_jy", "image_path", "image_sha256", "temperature_sha256"])
    _write_csv(output / "numerical_scatter.csv", scatter,
               ["inclination_deg", "probe_id", "wavelength_um", "photon_packets", "n_seeds", "mean_flux_jy", "sample_sd_jy",
                "sample_sd_fraction", "sd_lower_fraction", "sd_upper_fraction", "scatter_interval_degenerate",
                "required_for_gate", "eligible_for_gate", "passes_scatter"])
    _write_csv(output / "budget_changes.csv", changes,
               ["inclination_deg", "probe_id", "wavelength_um", "low_photon_packets", "high_photon_packets", "n_pairs",
                "mean_fractional_change", "sample_sd_fractional_change", "ci_half_width", "ci_lower", "ci_upper",
                "absolute_change_upper_bound", "eligible_for_gate", "passes_budget_change"])
    if make_plot:
        plot_convergence(output, summary, scatter, changes)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path, help="Complete extinction_ice_numerics_v2 run bundle")
    parser.add_argument("--output", type=Path, help="Report directory (default: BUNDLE/results)")
    args = parser.parse_args()
    try:
        result = analyze_bundle(args.bundle, args.output)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps({key: result[key] for key in ("status", "complete_tasks", "total_tasks", "numerical_gate_passed")}, indent=2))
    return 0 if result["numerical_gate_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
