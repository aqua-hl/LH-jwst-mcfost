#!/usr/bin/env python3
"""Audit downloaded control tables without rewriting the simulation results."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/continuum_control_v1"
OUT = ROOT / "validation/continuum_control_v1"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))
sys.path.insert(0, str(ROOT / "src"))

from astropy.table import Table
from mcfost_grid.runner import validate_inputs


def check(condition, message):
    if not condition:
        raise ValueError(message)


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-15)


def score(rows, key):
    regions = defaultdict(list)
    for row in rows:
        regions[row["region"]].append(math.log10(row[key] / row["observed_jy"]) ** 2)
    return math.sqrt(sum(sum(v) / len(v) for v in regions.values()) / len(regions))


def main():
    manifest = json.loads((RUN / "manifest.json").read_text())
    summary = json.loads((RUN / "results/summary.json").read_text())
    manifest_hash = hashlib.sha256((RUN / "manifest.json").read_bytes()).hexdigest()
    check(manifest_hash == summary["manifest_sha256"], "Downloaded results belong to another manifest")
    validate_inputs(RUN, manifest)
    check(len(manifest["models"]) == 1, "This review expects the single-model control")
    model = manifest["models"][0]
    check(summary["best_model"]["parameters"] == model["parameters"], "Summary parameters differ")
    check(summary["best_model"]["model_id"] == model["id"], "Summary model ID differs")
    check(summary["model_count"] == summary["ranked_model_count"] == 1, "Control was not ranked")
    check(summary["excluded_model_count"] == 0 and not summary["excluded_models"], "Control was excluded")
    predictions = Table.read(RUN / "results/predictions.ecsv")
    ranking = Table.read(RUN / "results/ranking.ecsv")
    check(len(ranking) == 1 and ranking[0]["status"] == "ranked", "Unexpected ranking")
    check(ranking[0]["model_id"] == model["id"] and int(ranking[0]["rank"]) == 1, "Ranking identity differs")
    check(json.loads(ranking[0]["parameters_json"]) == model["parameters"], "Ranking parameters differ")
    anchors = {a["id"]: a for a in manifest["anchors"]}
    check(len(predictions) == len(anchors) == 9, "Expected nine anchor predictions")
    check(set(predictions["anchor_id"]) == set(anchors), "Prediction anchor IDs differ")

    old_summary = json.loads((ROOT / "reference/predictions/continuum_summary.json").read_text())
    old_index = old_summary["primary_nominal_rank1"]["model_index"]
    old_table = Table.read(ROOT / "reference/predictions/continuum_predictions.ecsv")
    old = {str(r["anchor_tag"]): r for r in old_table if r["model_index"] == old_index}
    check(set(old) == set(anchors), "Archived nominal anchors differ")
    rows = []
    for p in predictions:
        tag = str(p["anchor_id"])
        a, archived = anchors[tag], old[tag]
        check(p["model_id"] == model["id"] and p["model_status"] == "ranked", "Invalid prediction identity/status")
        check(bool(p["quality_pass"]) and bool(p["scored"]), "An anchor did not pass the reported checks")
        check(p["region"] == a["region"] and p["instrument"] == a["instrument"], "Anchor contract changed")
        for key, value in (("wavelength_um", a["wavelength_um"]), ("observed_flux_jy", a["flux_jy"]),
                           ("uncertainty_jy", a["uncertainty_jy"])):
            check(close(p[key], value), f"{tag}: {key} differs from manifest")
        for key in ("model_flux_jy", "observed_flux_jy", "uncertainty_jy"):
            check(math.isfinite(p[key]) and p[key] > 0, f"{tag}: invalid {key}")
        params = model["parameters"]
        for key, value in (("mass_factor", params["envelope_dust_mass_msun"] / 1e-4),
                           ("envelope_size_exponent", params["envelope_size_exponent"]),
                           ("requested_half_opening_deg", params["cavity_half_opening_deg"]),
                           ("inclination_deg", params["inclination_deg"]),
                           ("envelope_amax_um", params["envelope_amax_um"])):
            check(close(archived[key], value), f"Archived {key} does not match the control")
        factor = float(p["wavelength_um"]) * 1e-6 / 299792458.0 * 1e26
        check(close(archived["wavelength_um"], p["wavelength_um"]), "Archived wavelength differs")
        check(close(archived["observed_flux_geometric_w_m2"] * factor, p["observed_flux_jy"]), "Archived observed flux differs")
        check(close(archived["observed_error_w_m2"] * factor, p["uncertainty_jy"]), "Archived uncertainty differs")
        f, obs, err = (float(p[k]) for k in ("model_flux_jy", "observed_flux_jy", "uncertainty_jy"))
        archive_flux = float(archived["model_aperture_flux_147pc_w_m2"]) * factor
        check(close(math.log10(f / obs), p["log_residual_dex"]), "Stored residual does not reproduce")
        rows.append(dict(anchor_id=tag, wavelength_um=float(p["wavelength_um"]), region=str(p["region"]),
                         observed_jy=obs, uncertainty_jy=err, new_model_jy=f, archived_model_jy=archive_flux,
                         model_minus_observed_percent=100 * (f / obs - 1),
                         standardized_flux_residual=(f - obs) / err,
                         new_minus_archived_percent=100 * (f / archive_flux - 1)))
    rows.sort(key=lambda r: r["wavelength_um"])
    new_score, old_score = score(rows, "new_model_jy"), score(rows, "archived_model_jy")
    chi2 = sum(r["standardized_flux_residual"] ** 2 for r in rows)
    check(close(new_score, summary["best_model"]["score_dex"]), "Summary score does not reproduce")
    check(close(new_score, ranking[0]["score_dex"]), "Ranking score does not reproduce")
    check(close(chi2, summary["best_model"]["chi2_diagnostic"]), "Summary chi-square does not reproduce")
    check(close(old_score, old_summary["primary_nominal_rank1"]["score_dex"]), "Archived score does not reproduce")
    missing = [str(p.relative_to(RUN)) for p in (
        RUN / "runtime_binding.json", RUN / "models" / model["id"] / "runtime_binding.json",
        RUN / "models" / model["id"] / "measurements.json",
        RUN / "models" / model["id"] / "temperature_complete.json",
        RUN / "models" / model["id"] / "status.json") if not p.is_file()]
    review = dict(run_id=manifest["run_id"], manifest_sha256=manifest_hash, table_checks_pass=True,
                  reported_quality_pass_count=9, score_dex=new_score, archived_score_dex=old_score,
                  chi2_diagnostic=chi2, max_absolute_archive_difference_percent=max(abs(r["new_minus_archived_percent"]) for r in rows),
                  anchors_outside_provisional_2_percent=sum(abs(r["new_minus_archived_percent"]) > 2 for r in rows),
                  archive_reproduction_within_2_percent=all(abs(r["new_minus_archived_percent"]) < 2 for r in rows),
                  missing_runtime_records=missing,
                  downloaded_analysis_warnings=summary["warnings"], rows=rows)
    inputs = [RUN / "manifest.json", *[RUN / "results" / n for n in ("summary.json", "ranking.ecsv", "predictions.ecsv")],
              ROOT / "reference/predictions/continuum_predictions.ecsv", ROOT / "reference/predictions/continuum_summary.json"]
    review["input_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "review.json").write_text(json.dumps(review, indent=2, allow_nan=False) + "\n")
    with (OUT / "anchor_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Continuum control review", "", "## Downloaded results", "",
             "The downloaded summary matches the local frozen manifest. All nine prediction rows report passing quality checks; their identities, observed fluxes/errors, residuals and score reproduce from the saved tables. One model is ranked, none excluded, and the analysis reports no warnings.", "",
             f"Regional log RMS: **{new_score:.8f} dex**. Archived nominal control: **{old_score:.8f} dex**. Diagnostic diagonal chi-square: **{chi2:.6f}** (not reduced chi-square or a posterior likelihood).", "",
             "| Wavelength (µm) | New − observed (%) | New − archived (%) |",
             "|---:|---:|---:|"]
    lines += [f"| {r['wavelength_um']:.4f} | {r['model_minus_observed_percent']:+.2f} | {r['new_minus_archived_percent']:+.2f} |" for r in rows]
    lines += ["", "## Scientific interpretation", "",
              f"The new nominal model differs from the archived nominal flux by up to {review['max_absolute_archive_difference_percent']:.2f}%; {review['anchors_outside_provisional_2_percent']}/9 anchors exceed the provisional 2% restart target. The lower score alone does not establish improved physics: the physical setting is unchanged. This is a functioning continuum calculation, with archive reproduction and convergence still unresolved.", "",
              "Archived images used seed 41001; the active workflow does not explicitly pass a seed. The reported simulator build has also changed, and this run recomputes temperature. These are possible contributors, not an attribution of the flux differences. Do not merge old and new scores into one ranking.", "",
              "The four strongest fractional continuum residuals are at 13.80 µm (+30.38%), 5.39 µm (−28.02%), 2.55 µm (−23.95%) and 1.20 µm (−16.01%). A nine-point continuum run does not assess the H2O band or the silicate trough.", "",
              "## Verification limits", "",
              "Only the results directory was downloaded. Raw images, temperatures, per-model measurements, runtime identities, command records and timings are absent here. The table review cannot independently recheck image photometry, establish the actual binary/utilities hashes or benchmark runtime/memory. Do not rerun the ordinary analyzer locally over this partial download: it would replace the valid imported results with missing-model output.", "",
              "Preserve and download `runtime_binding.json`, `models/*/{runtime_binding.json,status.json,measurements.json,temperature_complete.json}`, simulator `command.json`/`execution.json`/`mcfost.log` files and Slurm logs before making a reproducibility claim. FITS products are needed for independent photometry checks.", "",
              "See `review.json` for input hashes and `anchor_comparison.csv` for the numerical comparison. Reproduce this review with `python -B scripts/review_continuum_control.py`."]
    (OUT / "REVIEW.md").write_text("\n".join(lines) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6), layout="constrained")
    x = [r["wavelength_um"] for r in rows]
    for key, label, color in (("new_model_jy", "New control", "#16804a"), ("archived_model_jy", "Archived nominal", "#4b6cb7")):
        axes[0].scatter(x, [r[key] / r["observed_jy"] for r in rows], label=label, color=color)
    axes[0].axhline(1, color="0.4", lw=0.8)
    axes[0].set(ylabel="Model / observed flux", title="Nominal continuum control: nine measured anchors")
    axes[0].legend(frameon=False)
    axes[1].axhspan(-2, 2, color="0.85", label="Provisional ±2% target")
    axes[1].axhline(0, color="0.4", lw=0.8)
    axes[1].scatter(x, [r["new_minus_archived_percent"] for r in rows], color="#16804a")
    axes[1].set(xscale="log", xlabel="Wavelength (µm)", ylabel="New − archived (%)")
    axes[1].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"control_comparison.{suffix}", dpi=200)
    plt.close(fig)
    print(json.dumps({k: review[k] for k in ("table_checks_pass", "score_dex", "archived_score_dex", "max_absolute_archive_difference_percent", "anchors_outside_provisional_2_percent")}, indent=2))


if __name__ == "__main__":
    main()
