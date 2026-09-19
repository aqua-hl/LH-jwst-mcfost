#!/usr/bin/env python3
"""Remeasure the fixed-temperature V2 pixel-scale test; never rank or certify.

Ten images use one saved temperature, five image seeds and two pixel scales
at the same field of view. Each coarse/fine pair shares one allocation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy.stats import chi2, t


TEMPERATURE_SEED = 42004
DIAGNOSTIC_ID = "fixed_temperature_pixel_scale_v2"
IMAGE_SEEDS = tuple(range(42001, 42006))
GEOMETRIES = ("coarse", "fine")
COMPONENTS = ("total_i", "direct_star", "scattered_star", "direct_thermal", "scattered_thermal")
EXPECTED_CELLS = {(TEMPERATURE_SEED, seed, geometry)
                  for seed in IMAGE_SEEDS for geometry in GEOMETRIES}
CONFIDENCE = .95


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

    def invalid(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")

    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)


def _frozen_dispatcher(bundle):
    """Check the manifest and its inputs before executing any frozen Python."""
    experiment = read_json(bundle / "experiment.json")
    sidecar = (bundle / "experiment.sha256").read_text().split()
    if sidecar != [digest(bundle / "experiment.json"), "experiment.json"]:
        raise ValueError("Experiment hash differs from experiment.sha256")
    if (experiment.get("schema_version") != 2
            or experiment.get("diagnostic_id") != DIAGNOSTIC_ID):
        raise ValueError("Expected schema 2 fixed-temperature pixel-scale diagnostic; old temperature-arm bundles are incompatible")
    hashes = experiment.get("input_hashes", {})
    if not isinstance(hashes, dict) or "code/crossover_task.py" not in hashes:
        raise ValueError("Frozen code/crossover_task.py must be pinned in input_hashes")
    for relative, expected in hashes.items():
        relative = Path(relative)
        path = bundle / relative
        if (relative.is_absolute() or ".." in relative.parts
                or not path.resolve().is_relative_to(bundle) or path.is_symlink()
                or not path.is_file() or digest(path) != expected):
            raise ValueError(f"Frozen input changed, missing or unsafe: {relative}")
    name = "_v2_crossover_" + hashlib.sha256(str(bundle).encode()).hexdigest()[:16]
    spec = importlib.util.spec_from_file_location(name, bundle / "code/crossover_task.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        validated = module.validate_package(bundle)
    except Exception:
        sys.modules.pop(name, None)
        raise
    if validated != experiment:
        raise ValueError("Frozen validator returned a different experiment")
    return module, experiment


def scatter_statistics(values, confidence=CONFIDENCE):
    """Gaussian SD CI / plug-in mean; this is not an exact CV interval."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError("Scatter statistics require at least two finite values")
    if not 0 < confidence < 1:
        raise ValueError("Confidence must be between zero and one")
    mean = float(np.mean(values))
    identical = bool(np.ptp(values) == 0)
    sd = 0. if identical else float(np.std(values, ddof=1))
    alpha, degrees = 1 - confidence, len(values) - 1
    lower = sd * math.sqrt(degrees / chi2.ppf(1 - alpha / 2, degrees))
    upper = sd * math.sqrt(degrees / chi2.ppf(alpha / 2, degrees))
    return {"n_seeds": len(values), "mean_flux_jy": mean, "sample_sd_jy": sd,
            "sample_sd_fraction": sd / mean if mean > 0 else None,
            "sd_lower_fraction": lower / mean if mean > 0 else None,
            "sd_upper_fraction": upper / mean if mean > 0 else None,
            "scatter_interval_degenerate": identical or sd == 0,
            "positive_mean_for_fraction": mean > 0}


def paired_statistics(reference, comparison, confidence=CONFIDENCE):
    """Keep seed labels; fractions use the reference flux as denominator."""
    if set(reference) != set(comparison) or len(reference) < 2:
        raise ValueError("Paired comparisons require the same seed IDs and at least two pairs")
    if not 0 < confidence < 1:
        raise ValueError("Confidence must be between zero and one")
    seeds = sorted(reference)
    first = np.array([reference[seed] for seed in seeds], dtype=float)
    second = np.array([comparison[seed] for seed in seeds], dtype=float)
    if np.any(~np.isfinite(first)) or np.any(~np.isfinite(second)):
        raise ValueError("Paired values must be finite")
    absolute = second - first
    mean_delta = float(np.mean(absolute))
    delta_sd = 0. if np.ptp(absolute) == 0 else float(np.std(absolute, ddof=1))
    critical = float(t.ppf((1 + confidence) / 2, len(seeds) - 1))
    delta_half = critical * delta_sd / math.sqrt(len(seeds))
    result = {"n_pairs": len(seeds), "paired_seeds": seeds, "mean_delta_jy": mean_delta,
              "delta_ci_lower_jy": mean_delta - delta_half,
              "delta_ci_upper_jy": mean_delta + delta_half,
              "fractional_comparison_defined": bool(np.all(first > 0)),
              "mean_fractional_change": None, "sample_sd_fractional_change": None,
              "ci_half_width": None, "ci_lower": None, "ci_upper": None,
              "paired_interval_degenerate": None}
    if np.all(first > 0):
        differences = absolute / first
        mean = float(np.mean(differences))
        sd = 0. if np.ptp(differences) == 0 else float(np.std(differences, ddof=1))
        half = critical * sd / math.sqrt(len(seeds))
        result.update(mean_fractional_change=mean, sample_sd_fractional_change=sd,
                      ci_half_width=half, ci_lower=mean - half, ci_upper=mean + half,
                      paired_interval_degenerate=sd == 0)
    return result


def _selection(rows, geometry):
    return {row["image_seed"]: row for row in rows
            if row["temperature_seed"] == TEMPERATURE_SEED and row["geometry"] == geometry}


def _component(row, component):
    return row["flux_jy"] if component == "total_i" else row["components_aperture_jy"][component]


def _comparison(first, second, metadata, confidence):
    seeds = sorted(set(first) & set(second))
    result = dict(metadata, n_pairs=len(seeds), paired_seeds=seeds,
                  complete_five_seed_comparison=seeds == list(IMAGE_SEEDS),
                  missing_reference_seeds=sorted(set(IMAGE_SEEDS) - set(first)),
                  missing_comparison_seeds=sorted(set(IMAGE_SEEDS) - set(second)),
                  same_host_pairs=sum(first[s].get("host") == second[s].get("host")
                                      and bool(first[s].get("host")) for s in seeds),
                  components={})
    if len(seeds) >= 2:
        for component in COMPONENTS:
            # Component values remain signed. Zero denominators give an absolute
            # difference interval only; they are never clipped or renormalized.
            if all(component == "total_i" or component in first[s].get("components_aperture_jy", {})
                   and component in second[s].get("components_aperture_jy", {}) for s in seeds):
                statistics = paired_statistics({s: _component(first[s], component) for s in seeds},
                                               {s: _component(second[s], component) for s in seeds}, confidence)
                if component == "total_i":
                    result.update(statistics)
                else:
                    result["components"][component] = statistics
    return result


def summarize(rows, confidence=CONFIDENCE):
    """Describe every pixel-scale cell; invalid/duplicate rows cannot look complete.

    Partial groups are explicitly labelled and still get descriptive statistics
    when two values exist. Every original row is retained in the returned record.
    No flux-based rejection or production-convergence decision is made.
    """
    rows = list(rows)
    cells, errors, invalid = {}, [], []
    for position, row in enumerate(rows):
        try:
            key = (row["temperature_seed"], row["image_seed"], row["geometry"])
            if key not in EXPECTED_CELLS:
                raise ValueError(f"Unexpected pixel-scale cell: {key}")
            flux = row["flux_jy"]
            if isinstance(flux, bool) or not math.isfinite(flux) or flux <= 0:
                raise ValueError("Total aperture flux must be finite and positive")
            if row.get("quality_pass", True) is not True:
                raise ValueError("Image did not pass aperture quality checks")
            for value in row.get("components_aperture_jy", {}).values():
                if isinstance(value, bool) or not math.isfinite(value):
                    raise ValueError("Component flux must be finite")
            cells.setdefault(key, []).append((position, row))
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(f"Row {position}: {exc}")
            invalid.append({"row_index": position, "reason": str(exc), "row": row})
    valid = []
    for key, entries in sorted(cells.items()):
        if len(entries) != 1:
            reason = f"Duplicate pixel-scale cell {key}: row indices {[p for p, _ in entries]}"
            errors.append(reason)
            invalid.extend({"row_index": position, "reason": reason, "row": row}
                           for position, row in entries)
        else:
            valid.append(entries[0][1])
    valid_cells = {(r["temperature_seed"], r["image_seed"], r["geometry"]) for r in valid}
    missing = [{"temperature_seed": temp, "image_seed": seed, "geometry": geometry}
               for temp, seed, geometry in sorted(EXPECTED_CELLS - valid_cells)]
    scatter = []
    for geometry in GEOMETRIES:
        selected = _selection(valid, geometry)
        result = {"temperature_seed": TEMPERATURE_SEED, "geometry": geometry,
                  "n_seeds": len(selected), "image_seeds": sorted(selected),
                  "complete_five_seed_group": set(selected) == set(IMAGE_SEEDS),
                  "hosts": sorted({r["host"] for r in selected.values() if r.get("host")}),
                  "components": {}}
        if len(selected) >= 2:
            result.update(scatter_statistics([r["flux_jy"] for r in selected.values()], confidence))
            for component in COMPONENTS[1:]:
                if all(component in r.get("components_aperture_jy", {}) for r in selected.values()):
                    result["components"][component] = scatter_statistics(
                        [_component(r, component) for r in selected.values()], confidence)
        scatter.append(result)
    spatial = [_comparison(_selection(valid, "coarse"), _selection(valid, "fine"),
                           {"temperature_seed": TEMPERATURE_SEED, "reference_geometry": "coarse",
                            "comparison_geometry": "fine"}, confidence)]
    complete = not errors and not missing and len(valid) == 10
    return {"status": "integrity_error" if errors else ("complete_diagnostic" if complete else "incomplete"),
            "design_complete": complete, "expected_images": 10, "validated_images": len(valid),
            "missing_cells": missing, "errors": errors, "invalid_rows": invalid,
            "rows": rows, "valid_rows": valid, "scatter": scatter, "spatial_changes": spatial,
            "confidence": confidence,
            "production_convergence_certified": False, "outliers_removed": 0,
            "observation_fit_performed": False, "model_ranking_performed": False}


ASSUMPTIONS = {
    "scatter": "Two-sided Gaussian chi-square SD intervals with n-1 degrees of freedom, divided by the sample mean as a plug-in normalization; not exact coefficient-of-variation intervals.",
    "paired": "Two-sided Student-t intervals for the mean of per-image-seed differences. Pixel-scale fractions are (fine-coarse)/coarse. Seed labels are matched without assuming variance cancellation.",
    "independence": "Five image seeds are separate draws conditional on one fixed temperature. Five values cannot establish Gaussianity, stream independence, or absence of rare outliers. Matching image seeds does not ensure identical random trajectories or bitwise OpenMP reproducibility.",
    "fixed_temperature": "Every image reuses the identical saved temperature file from seed 42004, the original anomalous image's temperature. Temperature is excluded as an experimental variable; this experiment does not estimate temperature uncertainty.",
    "hosts": "Coarse/fine images are requested within one task allocation; resuming a partial task may place its second geometry on another worker. Host matches are reported explicitly. Different image seeds may also run on different workers.",
    "degenerate": "Identical values produce a degenerate interval and do not establish numerical convergence. No production pass is defined, even when all measured differences are small.",
    "multiple_comparisons": "Confidence intervals are marginal, not simultaneous familywise intervals. Component comparisons are diagnostics and their sum need not exactly equal total I.",
    "scope": "One 50-degree model, one 2.546352514425122-micron wavelength, 2.048M image photons, one frozen temperature, and two pixel scales (2401 and 4801 pixels across the same 6000-AU field) with a source-centred 1-arcsecond aperture. No broad-band, ice-aperture, photon-budget, whole-field SED or production-grid certification.",
    "model_data_tension": "The reported 5.1% direct-thermal fraction at 2.546 microns belongs to the audited nominal 50-degree model. It is not a measured fraction for every v1.1 model and does not resolve the separately reported 27.46% observational residual at w005.",
}


def _csv_value(value):
    if isinstance(value, (list, dict)):
        return json.dumps(value, sort_keys=True, allow_nan=False)
    return value


def _write_csv(path, rows, preferred):
    fields = list(preferred)
    fields.extend(sorted({key for row in rows for key in row} - set(fields)))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _percent(value):
    return "—" if value is None else f"{100 * value:.4f}%"


def _review(summary):
    lines = ["# Fixed-temperature pixel-scale diagnostic — V2", "",
             f"**{summary['status']}**: {summary['complete_tasks']}/5 tasks and "
             f"{summary['validated_images']}/10 image cells validated and remeasured.", "",
             "All ten images use one frozen temperature (seed 42004), 2.048M image photons, "
             "and a 6000-AU field. Pixel count changes from 2401 to 4801; five paired image seeds "
             "measure repeatability at each scale. There is no temperature arm or new temperature solve.", "",
             "This is a diagnostic comparison. Production convergence remains uncertified; "
             "no observational fit, model ranking or flux-based outlier rejection is performed.", ""]
    if summary["errors"]:
        lines.extend(["## Integrity issues", "", *[f"- {error}" for error in summary["errors"]], ""])
    if summary.get("warnings"):
        lines.extend(["## Execution notes", "", *[f"- {warning}" for warning in summary["warnings"]], ""])
    incomplete = [r for r in summary["tasks"] if r["state"] != "complete"]
    if incomplete:
        lines.extend(["## Tasks needing attention", "", "| Index | State | Detail |", "|---|---|---|"])
        for item in incomplete:
            detail = str(item.get("error") or item.get("reason") or "").replace("|", "\\|").replace("\n", " ")
            lines.append(f"| {item['index']} | {item['state']} | {detail} |")
        lines.append("")
    lines.extend(["## Scatter at fixed temperature and geometry", "",
                  "| Temperature seed | Geometry | Seeds | Sample SD / mean | 95% SD interval / mean | Degenerate |",
                  "|---|---|---|---|---|---|"])
    for row in summary["scatter"]:
        bounds = f"{_percent(row.get('sd_lower_fraction'))} to {_percent(row.get('sd_upper_fraction'))}"
        lines.append(f"| {row['temperature_seed']} | {row['geometry']} | {row['n_seeds']}/5 | "
                     f"{_percent(row.get('sample_sd_fraction'))} | {bounds} | "
                     f"{row.get('scatter_interval_degenerate', 'unavailable')} |")
    lines.extend(["", "## Pixel-scale changes: (fine − coarse) / coarse", "",
                  "| Temperature seed | Seed pairs | Mean change | 95% interval | Same-host pairs |",
                  "|---|---|---|---|---|"])
    for row in summary["spatial_changes"]:
        lines.append(f"| {row['temperature_seed']} | {row['n_pairs']}/5 | {_percent(row.get('mean_fractional_change'))} | "
                     f"{_percent(row.get('ci_lower'))} to {_percent(row.get('ci_upper'))} | {row['same_host_pairs']} |")
    lines.extend(["", "Partial groups are labelled by their available seed count and cannot complete the design. "
                  "The field of view is identical at both pixel scales. The primary comparison is the same supported "
                  "central aperture, not a production-grid or whole-field SED certification. "
                  "All component fluxes remain signed and are recorded without renormalization.", "",
                  "Component scatter and paired differences are included in the JSON and CSV reports. "
                  "The total-I plane is the canonical measurement; component sums are diagnostic.", "",
                  "## Interpretation limits", "", *[f"- {text}" for text in ASSUMPTIONS.values()], "",
                  "## Files", "", "- `summary.json`: provenance, every task and image record, statistics and assumptions.",
                  "- `seed_fluxes.csv`: all validated image measurements, component fluxes, hashes and worker metadata.",
                  "- `scatter.csv`, `spatial_changes.csv`: two pixel-scale groups and one paired comparison.",
                  "- `crossover_diagnostics.png` and `.pdf`: seed realizations and matched comparisons (when plotting is enabled).", "",
                  "Status reflects local files, not a live Slurm query. Original image and temperature files are read only.", ""])
    return "\n".join(lines)


def plot_diagnostics(output, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    colors = {"coarse": "#1767a6", "fine": "#ba4c22"}
    rows = summary["valid_rows"]
    for geometry, style in (("coarse", "o--"), ("fine", "s-")):
        selected = _selection(rows, geometry)
        seeds = sorted(selected)
        axes[0].plot([s - 42000 for s in seeds], [selected[s]["flux_jy"] for s in seeds],
                     style, color=colors[geometry], label=geometry)
        if all("direct_thermal" in selected[s].get("components_aperture_jy", {}) for s in seeds):
            axes[1].plot([s - 42000 for s in seeds],
                         [selected[s]["components_aperture_jy"]["direct_thermal"] for s in seeds],
                         style, color=colors[geometry], label=geometry)
    first, second = _selection(rows, "coarse"), _selection(rows, "fine")
    seeds = sorted(set(first) & set(second))
    axes[2].plot([s - 42000 for s in seeds],
                 [100 * (second[s]["flux_jy"] - first[s]["flux_jy"]) / first[s]["flux_jy"] for s in seeds],
                 "o-", color="#525252", label="paired pixel-scale change")
    for ax, ylabel in zip(axes, ("Total-I aperture flux (Jy)", "Direct thermal aperture flux (Jy)",
                                "(Fine − coarse) / coarse (%)")):
        ax.set_xlabel("Image seed − 42000")
        ax.set_ylabel(ylabel)
        ax.set_xticks(range(1, 6))
        ax.grid(alpha=.18)
        ax.legend(fontsize=8)
    axes[2].axhline(0, color="grey", linewidth=.7)
    fig.suptitle(f"2.546 µm pixel scale at fixed T42004 and 6000-AU field: {summary['status']}\n"
                 f"{summary['validated_images']}/10 images; seed points, no outlier removal or production certification")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"crossover_diagnostics.{suffix}", dpi=170)
    plt.close(fig)


def analyze_bundle(bundle, output=None, make_plot=True):
    bundle = Path(bundle).resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mcfost-crossover-mpl"))
    output = Path(output).resolve() if output else bundle / "results"
    if output.is_relative_to(bundle) and not output.is_relative_to(bundle / "results"):
        raise ValueError("Reports inside the bundle must be under BUNDLE/results to protect frozen inputs")
    dispatcher, experiment = _frozen_dispatcher(bundle)
    tasks = experiment["tasks"]
    expected_tasks = {(TEMPERATURE_SEED, seed) for seed in IMAGE_SEEDS}
    if (len(tasks) != 5 or {task["index"] for task in tasks} != set(range(5))
            or {(task["temperature_seed"], task["image_seed"]) for task in tasks} != expected_tasks):
        raise ValueError("Expected exactly five distinct image-seed tasks at fixed temperature seed 42004")
    temperature_hashes = {task["temperature_sha256"] for task in tasks}
    if (len(temperature_hashes) != 1 or any(not isinstance(value, str) or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value) for value in temperature_hashes)):
        raise ValueError("All five tasks must pin the same valid temperature SHA-256")
    if ((output / "temperature_changes.csv").exists()
            or (output / "summary.json").exists()
            and read_json(output / "summary.json").get("diagnostic_id") != DIAGNOSTIC_ID):
        raise ValueError("Report destination contains another diagnostic; choose a new output directory")
    statuses, rows, errors, warnings = [], [], [], []
    for task in sorted(tasks, key=lambda value: value["index"]):
        try:
            status, task_rows = dispatcher.inspect_task(bundle, experiment, task, remeasure=True)
            if status.get("index") != task["index"]:
                raise ValueError("Task inspector returned the wrong task index")
            if status.get("state") not in {"complete", "partial", "failed", "not_started", "running", "invalid"}:
                raise ValueError("Task inspector returned an unrecognized state")
            for row in task_rows:
                if any(row.get(key) != task[key] for key in ("temperature_seed", "image_seed")):
                    raise ValueError("Inspected measurement has a different factor identity")
                if row.get("task_index") != task["index"] or row.get("temperature_sha256") != task["temperature_sha256"]:
                    raise ValueError("Inspected measurement has a different task or temperature hash")
                if row.get("threads") != 64 or not row.get("host") or row.get("quality_pass") is not True:
                    raise ValueError("Inspected measurement lacks valid thread/host/aperture evidence")
                if not isinstance(row.get("image_sha256"), str) or len(row["image_sha256"]) != 64:
                    raise ValueError("Inspected measurement lacks an image hash")
                if set(row.get("components_aperture_jy", {})) != set(COMPONENTS):
                    raise ValueError("Inspected measurement lacks the expected intensity components")
                if (row.get("image_npix") != {"coarse": 2401, "fine": 4801}.get(row.get("geometry"))
                        or row.get("image_size_au") != 6000):
                    raise ValueError("Inspected image geometry differs from the two-scale design")
            if status["state"] == "complete" and (len(task_rows) != 2 or {r["geometry"] for r in task_rows} != set(GEOMETRIES)):
                raise ValueError("A complete task must contain both unique geometries")
            if len({row["host"] for row in task_rows}) > 1:
                warnings.append(f"Task {task['index']}: coarse/fine images ran on different hosts, "
                                "possibly after resuming; this pair does not control the worker.")
            status["pair_host_matched"] = (len(task_rows) == 2 and len({row["host"] for row in task_rows}) == 1)
            statuses.append(status)
            rows.extend(task_rows)
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
            statuses.append({"index": task["index"], "temperature_seed": task["temperature_seed"],
                             "image_seed": task["image_seed"], "state": "invalid",
                             "error": f"{type(exc).__name__}: {exc}"})
        if statuses[-1]["state"] == "invalid":
            errors.append(f"Task {task['index']}: {statuses[-1].get('error') or statuses[-1].get('reason') or 'invalid'}")
    summary = summarize(rows)
    summary["errors"].extend(errors)
    completed = sum(status["state"] == "complete" for status in statuses)
    summary.update(schema_version=2, experiment_id=experiment.get("experiment_id"),
                   diagnostic_id=DIAGNOSTIC_ID, analysis="fixed_temperature_pixel_scale_diagnostic",
                   experiment_sha256=digest(bundle / "experiment.json"), analyzer_sha256=digest(__file__),
                   total_tasks=5, complete_tasks=completed, tasks=statuses,
                   fixed_temperature_seed=TEMPERATURE_SEED,
                   fixed_temperature_sha256=next(iter(temperature_hashes)), new_temperature_solves=0,
                   task_state_counts=dict(Counter(status["state"] for status in statuses)),
                   assumptions=ASSUMPTIONS, warnings=warnings, remeasurement_performed=True,
                   resources=experiment.get("resources"), production_authorized=False)
    if summary["errors"]:
        summary["status"] = "integrity_error"
    elif completed != 5 or not summary["design_complete"]:
        summary["status"] = "incomplete"
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    _write_csv(output / "seed_fluxes.csv", summary["rows"],
               ["task_index", "temperature_seed", "image_seed", "geometry", "flux_jy", "components_aperture_jy"])
    _write_csv(output / "scatter.csv", summary["scatter"], ["temperature_seed", "geometry", "n_seeds"])
    _write_csv(output / "spatial_changes.csv", summary["spatial_changes"], ["temperature_seed", "n_pairs"])
    (output / "REVIEW.md").write_text(_review(summary))
    if make_plot:
        plot_diagnostics(output, summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", nargs="?", type=Path, help="Frozen crossover run bundle")
    parser.add_argument("--bundle", dest="bundle_option", type=Path, help="Frozen crossover run bundle")
    parser.add_argument("--output", type=Path, help="Report directory (default: BUNDLE/results)")
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG and PDF generation")
    args = parser.parse_args(argv)
    if bool(args.bundle) == bool(args.bundle_option):
        parser.error("Supply the bundle once, as a positional path or with --bundle")
    try:
        summary = analyze_bundle(args.bundle or args.bundle_option, args.output, not args.no_plot)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ImportError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps({key: summary[key] for key in
                      ("status", "complete_tasks", "total_tasks", "validated_images", "production_convergence_certified")}, indent=2))
    return 0 if summary["status"] == "complete_diagnostic" else 2


if __name__ == "__main__":
    raise SystemExit(main())
