#!/usr/bin/env python3
"""Check downloaded production tables without overwriting their results directory."""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation/continuum_production_v1"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))
sys.path.insert(0, str(ROOT / "src"))
from astropy.table import Table
from mcfost_grid.runner import validate_inputs


def check(ok, message):
    if not ok:
        raise ValueError(message)


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-15)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(run):
    manifest = json.loads((run / "manifest.json").read_text())
    summary = json.loads((run / "results/summary.json").read_text())
    ranking = Table.read(run / "results/ranking.ecsv")
    predictions = Table.read(run / "results/predictions.ecsv")
    validate_inputs(run, manifest)
    check(summary["manifest_sha256"] == digest(run / "manifest.json"), "Result manifest identity differs")
    check(summary["run_id"] == manifest["run_id"], "Run identity differs")
    models = {m["id"]: m for m in manifest["models"]}
    anchors = {a["id"]: a for a in manifest["anchors"]}
    check(len(anchors) == 9 and all(a["score"] for a in anchors.values()), "Expected nine scored anchors")
    check(summary["model_count"] == summary["ranked_model_count"] == len(models), "Some models were not ranked")
    check(summary["excluded_model_count"] == 0 and not summary["excluded_models"], "Some models were excluded")
    check(summary["status_counts"] == {"ranked": len(models)} and not summary["warnings"], "Unexpected status counts or warnings")
    check(summary["n_scored_anchors"] == 9 and summary["n_regions"] == 6, "Scoring contract differs")
    check(summary["expected_scored_anchor_ids"] == list(anchors), "Anchor order differs")
    check(len(ranking) == len(models) and set(ranking["model_id"]) == set(models), "Ranking model IDs differ")
    check(len(predictions) == 9 * len(models), "Prediction row count differs")
    check(set(predictions["model_id"]) == set(models), "Prediction model IDs differ")
    values, rows = {}, []
    for rank in ranking:
        mid = str(rank["model_id"])
        model = models[mid]
        check(rank["status"] == "ranked" and int(rank["model_index"]) == model["index"], "Ranking identity/status differs")
        check(json.loads(rank["parameters_json"]) == model["parameters"], "Ranking parameters differ")
        for key, value in model["parameters"].items():
            check(close(rank["param_" + key], value), "Ranking parameter column differs")
        group = predictions[predictions["model_id"] == mid]
        check(len(group) == 9 and set(group["anchor_id"]) == set(anchors), "Missing/duplicate anchor")
        regions, chi2, values[mid] = defaultdict(list), 0.0, {}
        for p in group:
            tag = str(p["anchor_id"])
            a = anchors[tag]
            check(bool(p["quality_pass"]) and bool(p["scored"]) and p["model_status"] == "ranked", "Prediction quality/status failed")
            check(p["instrument"] == a["instrument"] and p["region"] == a["region"], "Prediction contract differs")
            for key, field in (("wavelength_um", "wavelength_um"), ("observed_flux_jy", "flux_jy"), ("uncertainty_jy", "uncertainty_jy")):
                check(close(p[key], a[field]), f"{mid}/{tag}: {key} differs")
            flux = float(p["model_flux_jy"])
            check(math.isfinite(flux) and flux > 0, "Invalid model flux")
            residual = math.log10(flux / a["flux_jy"])
            check(close(residual, p["log_residual_dex"]), "Saved residual differs")
            regions[a["region"]].append(residual * residual)
            chi2 += ((flux - a["flux_jy"]) / a["uncertainty_jy"]) ** 2
            values[mid][tag] = flux
        score = math.sqrt(math.fsum(math.fsum(v) / len(v) for v in regions.values()) / len(regions))
        regional = json.loads(rank["region_rms_dex_json"])
        check(set(regional) == set(regions), "Saved region IDs differ")
        check(all(close(regional[k], math.sqrt(math.fsum(v) / len(v))) for k, v in regions.items()), "Saved region scores differ")
        check(close(score, rank["score_dex"]) and close(chi2, rank["chi2_diagnostic"]), "Saved model score differs")
        check(int(rank["n_scored"]) == 9, "Wrong scored count")
        rows.append(dict(rank=int(rank["rank"]), index=model["index"], model_id=mid,
                         score_dex=score, chi2_diagnostic=chi2, **model["parameters"]))
    rows.sort(key=lambda x: (x["score_dex"], x["index"]))
    check([r["rank"] for r in rows] == list(range(1, len(rows) + 1)), "Rank order does not reproduce")
    check(summary["best_model"]["model_id"] == rows[0]["model_id"], "Best identity differs")
    check(summary["best_model"]["parameters"] == models[rows[0]["model_id"]]["parameters"], "Best parameters differ")
    check(close(summary["best_model"]["score_dex"], rows[0]["score_dex"]), "Best score differs")
    check(close(summary["best_model"]["chi2_diagnostic"], rows[0]["chi2_diagnostic"]), "Best chi-square differs")
    return manifest, summary, rows, values


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    run = ROOT / "runs/continuum_production_v1"
    control = ROOT / "runs/continuum_control_v1"
    m, summary, rows, fluxes = audit(run)
    cm, cs, cr, cf = audit(control)
    check(len(rows) == 12 and len(cr) == 1, "Unexpected experiment sizes")
    check(m["anchors"] == cm["anchors"] and m["measurement"] == cm["measurement"], "Control observation contract differs")
    nominal = m["models"][0]
    check(nominal["parameters"] == cm["models"][0]["parameters"], "Control physical parameters differ")
    check(m["configuration"]["numerics"] == cm["configuration"]["numerics"], "Control numerics differ")
    mid, cid, best = nominal["id"], cr[0]["model_id"], rows[0]["model_id"]
    check(best == mid, "The best model changed; revise this experiment-specific interpretation")
    for model_fluxes in fluxes.values():
        for a in m["anchors"]:
            if a["id"] in ("w005", "w006", "w007"):
                check(model_fluxes[a["id"]] < a["flux_jy"], "The common underprediction pattern changed")
            if a["id"] in ("w008", "w009"):
                check(model_fluxes[a["id"]] > a["flux_jy"], "The common overprediction pattern changed")
    comparison = []
    for a in m["anchors"]:
        tag = a["id"]
        comparison.append(dict(anchor_id=tag, wavelength_um=a["wavelength_um"], observed_jy=a["flux_jy"],
            control_jy=cf[cid][tag], production_nominal_jy=fluxes[mid][tag], best_jy=fluxes[best][tag],
            nominal_repeat_difference_percent=100 * (fluxes[mid][tag] / cf[cid][tag] - 1),
            best_minus_observed_percent=100 * (fluxes[best][tag] / a["flux_jy"] - 1)))
    differences = [r["nominal_repeat_difference_percent"] for r in comparison]
    repeated_score = next(r["score_dex"] for r in rows if r["model_id"] == mid)
    max_difference = max(map(abs, differences))
    missing = [str(p.relative_to(run)) for p in [run / "runtime_binding.json", *[
        run / "models" / model["id"] / filename for model in m["models"]
        for filename in ("status.json", "measurements.json", "runtime_binding.json", "temperature_complete.json")]] if not p.is_file()]
    report = dict(run_id=m["run_id"], table_checks_pass=True, ranked_models=len(rows), prediction_rows=9 * len(rows),
        reported_quality_pass_count=9 * len(rows), warnings=summary["warnings"], best_model=rows[0],
        top_two_score_gap_dex=rows[1]["score_dex"] - rows[0]["score_dex"],
        nominal_score_shift_dex=repeated_score - cr[0]["score_dex"],
        max_absolute_nominal_repeat_difference_percent=max_difference,
        rms_nominal_repeat_difference_percent=math.sqrt(math.fsum(v*v for v in differences) / len(differences)),
        repeat_anchors_outside_2_percent=sum(abs(v) > 2 for v in differences),
        previously_failed_indices_now_ranked=[r["index"] for r in rows if r["index"] in (8, 9, 10)],
        missing_runtime_records=missing, ranking=rows, nominal_comparison=comparison)
    inputs = [r / name for r in (run, control) for name in
              ("manifest.json", "results/summary.json", "results/ranking.ecsv", "results/predictions.ecsv")]
    report["input_sha256"] = {str(p.relative_to(ROOT)): digest(p) for p in inputs}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results_review.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_csv(OUT / "reviewed_ranking.csv", rows)
    write_csv(OUT / "nominal_repeat.csv", comparison)
    lines = ["# Production continuum results review", "",
        "All 12 models are ranked, with 108 passing prediction rows, nine scored anchors per model, zero exclusions and no reported warnings. The downloaded manifest identity matches the prepared experiment; frozen input hashes, model parameters, anchor contracts, residuals, regional scores, diagnostic chi-square and rank order all pass independent table checks.", "",
        "Previously failed indices 8, 9 and 10 are now present with all nine passing predictions each. This establishes their inclusion in the completed analysis; scheduler/recovery logs were not downloaded, so the particular recovery method is not independently verified.", "",
        "## Ranking", "", "| Rank | Index | Inclination (°) | Envelope dust mass (M☉) | q | Cavity (°) | Score (dex) |",
        "|---:|---:|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['rank']} | {r['index']} | {r['inclination_deg']:g} | {r['envelope_dust_mass_msun']:.6g} | {r['envelope_size_exponent']:g} | {r['cavity_half_opening_deg']:g} | {r['score_dex']:.6f} |" for r in rows]
    lines += ["", "The nominal 70° model remains best. The lower-mass (2.0e-4 M☉) and lower-q (2.5) directions do not improve the score in this sampled set. The alternative mass2.5/q3.25 family changes two parameters together, so it does not isolate either effect.", "",
        "Lower inclinations improve some intermediate-wavelength fluxes but worsen the long-wavelength excess. Their poorer scores at fixed heating and dust do not exclude those inclinations after other parameters are refitted. This is a sparse continuum screen, not an inclination posterior or a test of ice/silicate absorption profiles.", "",
        "## Repeatability and remaining residuals", "",
        f"The identical nominal setting scores {repeated_score:.8f} dex here versus {cr[0]['score_dex']:.8f} dex in the preceding control. Their fluxes differ by up to **{max_difference:.2f}%** (RMS {report['rms_nominal_repeat_difference_percent']:.2f}%); {report['repeat_anchors_outside_2_percent']}/9 anchors exceed the provisional 2% target.", "",
        f"The score shift between nominal runs is {report['nominal_score_shift_dex']:.6f} dex, compared with the top-two gap of {report['top_two_score_gap_dex']:.6f} dex. The cavity17.5/cavity20 ordering is therefore fragile. One repeat does not estimate a statistical uncertainty; missing runtime records prevent attributing the differences solely to Monte Carlo noise.", "",
        "| Wavelength (µm) | Best − observed (%) | Production nominal − earlier control (%) |",
        "|---:|---:|---:|"]
    lines += [f"| {r['wavelength_um']:.4f} | {r['best_minus_observed_percent']:+.2f} | {r['nominal_repeat_difference_percent']:+.2f} |" for r in comparison]
    lines += ["", "![Repeated nominal continuum calculation](nominal_repeat.png)", "",
        "The principal residual pattern remains: every tested model underpredicts 2.55, 3.91 and 5.39 µm, while overpredicting 13.80 and 27.51 µm. The best model is about 28–29% low at 2.55/5.39 µm and 31% high at 13.80 µm.", "",
        "## Next decision and verification limits", "",
        "The evidence does not support simply extending the lower-mass or lower-q directions. Before using small score differences to refine the search, check repeatability with controlled seeds and photon counts and verify the actual runtime identities. The remaining residual pattern also motivates the deferred joint geometry/heating tests in GOAL.md; this batch held heating and dust composition fixed. No new simulations or submissions were made during this review.", "",
        "Only compact results are available here. Per-model measurements, temperature receipts, binary/utilities identities, Python environments, scheduler/recovery logs, timings and FITS products are absent. Raw photometry, resource usage and the exact retry environment cannot be independently verified. Preserve those records on the cluster. Do not rerun the ordinary analyzer locally on this results-only download.", "",
        "Reproduce this review with `python -B scripts/review_continuum_production.py`. Input hashes and numerical checks are in [results_review.json](results_review.json); CSV comparisons are [reviewed_ranking.csv](reviewed_ranking.csv) and [nominal_repeat.csv](nominal_repeat.csv). The original [spectrum plot](../../runs/continuum_production_v1/results/spectrum_comparison.png) is unchanged."]
    (OUT / "RESULTS_REVIEW.md").write_text("\n".join(lines) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6), layout="constrained")
    wave = [r["wavelength_um"] for r in comparison]
    for key, label, color in (("control_jy", "Earlier nominal control", "#4b6cb7"), ("production_nominal_jy", "Production nominal (rank 1)", "#16804a")):
        axes[0].scatter(wave, [r[key] / r["observed_jy"] for r in comparison], label=label, color=color)
    axes[0].axhline(1, color="0.5", lw=.8)
    axes[0].set(ylabel="Model / observed flux", title="Same physical model: continuum control versus production")
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].axhspan(-2, 2, color="0.9", label="Provisional ±2% target")
    axes[1].axhline(0, color="0.5", lw=.8)
    axes[1].scatter(wave, differences, color="#16804a")
    axes[1].set(xscale="log", xlabel="Wavelength (µm)", ylabel="Production − control (%)")
    axes[1].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"nominal_repeat.{suffix}", dpi=200)
    plt.close(fig)
    check(all(digest(ROOT / name) == h for name, h in report["input_sha256"].items()), "Review modified input results")
    print(json.dumps({k: report[k] for k in ("table_checks_pass", "ranked_models", "prediction_rows", "max_absolute_nominal_repeat_difference_percent", "top_two_score_gap_dex", "nominal_score_shift_dex")}, indent=2))


if __name__ == "__main__":
    main()
