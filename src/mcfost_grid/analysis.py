"""Rank compact measurements and plot an aperture-matched grid screen.

Only the run manifest, per-model measurement JSON and observational ECSV are
read. Image FITS are deliberately outside this analysis interface.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
from typing import Any


class AnalysisError(ValueError):
    """The shared observation contract cannot support a comparable ranking."""


def _positive(value: Any) -> bool:
    try:
        return not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0
    except (TypeError, ValueError):
        return False


def validate_anchors(anchors: list[dict]) -> list[dict]:
    """Validate the common scored set once, not independently for each model."""
    if not isinstance(anchors, list) or not anchors:
        raise AnalysisError("The manifest must contain a nonempty anchor list")
    ids = [str(a["id"]) for a in anchors]
    if len(set(ids)) != len(ids):
        raise AnalysisError("Duplicate anchor IDs in the manifest")
    for anchor in anchors:
        if not isinstance(anchor.get("score"), bool):
            raise AnalysisError(f"Anchor {anchor['id']} requires an explicit boolean score flag")
        if not _positive(anchor.get("wavelength_um")):
            raise AnalysisError(f"Invalid wavelength for anchor {anchor['id']}")
        if anchor["score"] and (not _positive(anchor.get("flux_jy")) or not anchor.get("region")):
            raise AnalysisError(f"Scored anchor {anchor['id']} requires positive flux and a region")
    scored = [a for a in anchors if a["score"]]
    if not scored:
        raise AnalysisError("No scientific anchors are marked for scoring")
    return scored


def region_balanced_log_rms(residuals_by_region: dict[str, list[float]]) -> float:
    """sqrt(mean_regions(mean_anchors(log10(model/observed)^2)))."""
    if not residuals_by_region or any(not values for values in residuals_by_region.values()):
        raise AnalysisError("Every scored region needs at least one residual")
    if any(not math.isfinite(x) for values in residuals_by_region.values() for x in values):
        raise AnalysisError("Cannot score nonfinite residuals")
    return math.sqrt(math.fsum(math.fsum(x*x for x in values) / len(values)
                               for values in residuals_by_region.values()) / len(residuals_by_region))


def evaluate_model(model: dict, anchors: list[dict], payload: dict | None) -> dict:
    """Return an explicit status, predictions and score for one manifest model.

    Missing/failed scored anchors exclude a model. Diagnostic anchors never
    alter its score, and no model can improve by dropping a difficult point.
    """
    scored = validate_anchors(anchors)
    result = {"model_id": str(model["id"]), "model_index": int(model["index"]),
              "parameters": model.get("parameters", {}), "status": "missing",
              "reason": "measurements.json is absent", "score_dex": None,
              "chi2_diagnostic": None, "n_scored": 0, "region_rms_dex": {},
              "predictions": [], "smoke_test": False}
    if payload is None:
        return result
    if not isinstance(payload, dict):
        return {**result, "status": "invalid", "reason": "Measurement payload must be an object"}
    result["smoke_test"] = payload.get("smoke_test") is True
    if (str(payload.get("model_id")) != result["model_id"]
            or payload.get("model_index") != result["model_index"]
            or payload.get("parameters") != result["parameters"]):
        return {**result, "status": "invalid", "reason": "Model identity/parameters disagree with manifest"}
    measurements = payload.get("measurements", [])
    if not isinstance(measurements, list) or any(not isinstance(m, dict) for m in measurements):
        return {**result, "status": "invalid", "reason": "Invalid measurements list"}
    ids = [str(m.get("anchor_id")) for m in measurements]
    if len(ids) != len(set(ids)):
        return {**result, "status": "invalid", "reason": "Duplicate measurement anchor IDs"}
    expected = {str(a["id"]) for a in anchors}
    if set(ids) - expected:
        return {**result, "status": "invalid", "reason": "Unexpected measurement anchor IDs"}
    by_id = dict(zip(ids, measurements))
    failures = []
    residuals: dict[str, list[float]] = defaultdict(list)
    chi2_terms = []
    errors_valid = True
    for anchor in anchors:
        measurement = by_id.get(str(anchor["id"]))
        valid = (measurement is not None and measurement.get("quality_pass") is True
                 and _positive(measurement.get("flux_jy")))
        if valid:
            try:
                valid = math.isclose(float(measurement["wavelength_um"]),
                                     float(anchor["wavelength_um"]), rel_tol=1e-7, abs_tol=1e-10)
            except (KeyError, TypeError, ValueError):
                valid = False
        flux = float(measurement["flux_jy"]) if valid else None
        observed = anchor.get("flux_jy")
        residual = math.log10(flux) - math.log10(float(observed)) if valid and _positive(observed) else None
        result["predictions"].append({"model_id": result["model_id"],
            "anchor_id": str(anchor["id"]), "wavelength_um": float(anchor["wavelength_um"]),
            "instrument": str(anchor.get("instrument", "")), "region": str(anchor.get("region", "")),
            "scored": anchor["score"], "quality_pass": valid,
            "model_flux_jy": flux, "observed_flux_jy": float(observed) if _positive(observed) else None,
            "uncertainty_jy": float(anchor["uncertainty_jy"]) if _positive(anchor.get("uncertainty_jy")) else None,
            "log_residual_dex": residual})
        if not anchor["score"]:
            continue
        if not valid:
            failures.append(str(anchor["id"]))
            continue
        residuals[str(anchor["region"])].append(residual)
        result["n_scored"] += 1
        if _positive(anchor.get("uncertainty_jy")):
            try:
                term = ((flux - float(observed)) / float(anchor["uncertainty_jy"])) ** 2
            except OverflowError:
                term = math.inf
            errors_valid = errors_valid and math.isfinite(term)
            chi2_terms.append(term)
        else:
            errors_valid = False
    if failures or payload.get("complete") is not True:
        reason = "Missing or failed scientific anchors: " + ", ".join(failures) if failures else "Completion flag is not true"
        return {**result, "status": "incomplete", "reason": reason}
    if result["n_scored"] != len(scored):
        raise AnalysisError("Internal error: inconsistent scored-anchor count")
    result.update(status="ranked", reason="", score_dex=region_balanced_log_rms(residuals),
                  chi2_diagnostic=math.fsum(chi2_terms) if errors_valid else None,
                  region_rms_dex={r: math.sqrt(math.fsum(x*x for x in v)/len(v)) for r, v in residuals.items()})
    return result


def _load_measurement(path: Path) -> tuple[dict | None, str | None]:
    try:
        return json.loads(path.read_text()), None
    except FileNotFoundError:
        return None, None
    except (OSError, ValueError) as exc:
        return None, f"Unreadable measurement JSON: {exc}"


def _sha256_string(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _verified_observation(manifest: dict, run_dir: Path) -> Path:
    """Reject edited shared context before creating potentially misleading plots."""
    relative = manifest.get("observation_spectrum", "inputs/observed_spectrum.ecsv")
    path = run_dir / relative
    if "input_hashes" not in manifest:
        return path
    hashes = manifest["input_hashes"]
    if (not isinstance(hashes, dict) or not _sha256_string(hashes.get(relative))
            or Path(relative).is_absolute() or not path.resolve().is_relative_to(run_dir)):
        raise AnalysisError("Pinned manifest requires a valid, in-run observation_spectrum input hash")
    try:
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        raise AnalysisError(f"Pinned observation spectrum is unavailable: {path}") from exc
    if actual_hash != hashes[relative]:
        raise AnalysisError("Observation spectrum changed from its pinned input hash; use an unmodified run snapshot")
    return path


def _measurement_fingerprint(payload: dict, manifest_hash: str, required: bool) -> tuple[tuple | None, str | None]:
    """Return the simulation identity, or a reason this result is not comparable."""
    if "fingerprint" not in payload and not required:
        return None, None  # Compact, unpinned analysis fixtures remain supported.
    fingerprint = payload.get("fingerprint")
    if not isinstance(fingerprint, dict) or fingerprint.get("manifest_sha256") != manifest_hash:
        return None, "Measurement fingerprint is missing or does not match this manifest"
    executable = fingerprint.get("executable_sha256")
    backend = fingerprint.get("backend")
    utilities = fingerprint.get("utilities_content_sha256")
    if (not _sha256_string(executable) or backend not in ("image_method2", "coeval_method2")
            or (utilities is not None and not _sha256_string(utilities))):
        return None, "Measurement fingerprint has an invalid executable, utilities hash or backend"
    for measurement in payload.get("measurements", []):
        if isinstance(measurement, dict) and "backend" in measurement and measurement["backend"] != backend:
            return None, "Measurement backend disagrees with its fingerprint"
    return (executable, backend, utilities), None


def _write_tables(results: list[dict], out: Path) -> None:
    import numpy as np
    from astropy.table import Table

    ranked = sorted((r for r in results if r["status"] == "ranked"), key=lambda r: (r["score_dex"], r["model_index"]))
    ranks = {r["model_id"]: i+1 for i, r in enumerate(ranked)}
    parameter_names = sorted({key for r in results for key in r["parameters"]})
    rows = []
    for result in sorted(results, key=lambda r: (ranks.get(r["model_id"], math.inf), r["model_index"])):
        row = {"rank": ranks.get(result["model_id"], 0), "model_id": result["model_id"],
               "model_index": result["model_index"], "status": result["status"], "reason": result["reason"],
               "score_dex": result["score_dex"] if result["score_dex"] is not None else np.nan,
               "chi2_diagnostic": result["chi2_diagnostic"] if result["chi2_diagnostic"] is not None else np.nan,
               "n_scored": result["n_scored"], "parameters_json": json.dumps(result["parameters"], sort_keys=True),
               "region_rms_dex_json": json.dumps(result["region_rms_dex"], sort_keys=True)}
        for name in parameter_names:
            row[f"param_{name}"] = result["parameters"].get(name)
        rows.append(row)
    ranking = Table(rows=rows) if rows else Table(names=["rank", "model_id", "score_dex"], dtype=[int, str, float])
    ranking.meta.update(metric="region-balanced log10 RMS", rank_zero="excluded", posterior=False,
                        chi2="Unreduced diagonal-error diagnostic; not the ranking objective")
    ranking.write(out / "ranking.ecsv", format="ascii.ecsv", overwrite=True)
    predictions = [{**p, "model_status": r["status"], "model_exclusion_reason": r["reason"]}
                   for r in results for p in r["predictions"]]
    numeric = {"model_flux_jy", "observed_flux_jy", "uncertainty_jy", "log_residual_dex"}
    for row in predictions:
        for key in numeric:
            if row[key] is None:
                row[key] = np.nan
    table = Table(rows=predictions) if predictions else Table(names=["model_id", "anchor_id", "model_flux_jy"], dtype=[str, str, float])
    table.meta["note"] = "Actual sampled wavelengths only; unscored diagnostics do not affect ranking"
    table.write(out / "predictions.ecsv", format="ascii.ecsv", overwrite=True)


def _plots(manifest: dict, results: list[dict], observed_path: Path, out: Path, smoke: bool) -> list[str]:
    import os
    os.environ.setdefault("MPLCONFIGDIR", str(out / ".mplconfig"))
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    from astropy.table import Table

    warnings = []
    ranked = sorted((r for r in results if r["status"] == "ranked"), key=lambda r: (r["score_dex"], r["model_index"]))
    title = "SMOKE TEST — not a scientific fit" if smoke else "TMC1A aperture-matched grid screen"
    style = {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
             "savefig.dpi": 300, "pdf.fonttype": 42}
    with plt.rc_context(style):
        fig, (ax, res) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1.3]}, layout="constrained")
        if observed_path.is_file():
            try:
                observed = Table.read(observed_path, format="ascii.ecsv")
                lam, flux = np.asarray(observed["wavelength_um"], float), np.asarray(observed["flux_jy"], float)
                valid = np.isfinite(lam) & np.isfinite(flux) & (lam > 0) & (flux > 0)
                ax.scatter(lam[valid], 2.99792458e-12*flux[valid]/lam[valid], s=5, c="0.65", alpha=.6,
                           linewidths=0, label="JWST spectrum (context)", rasterized=True)
            except (OSError, ValueError, KeyError) as exc:
                warnings.append(f"Observation context could not be plotted: {exc}")
        else:
            warnings.append("Observation context ECSV is absent; anchors remain plotted")
        for scored, marker, color, label in [(True, "o", "black", "Scored JWST anchors"),
                                             (False, "P", "#8750a1", "JWST diagnostics (unscored)")]:
            anchors = [a for a in manifest["anchors"] if a["score"] is scored and _positive(a.get("flux_jy"))]
            if not anchors:
                continue
            x = np.array([a["wavelength_um"] for a in anchors], float)
            y = np.array([a["flux_jy"] for a in anchors], float) * 2.99792458e-12/x
            ax.scatter(x, y, s=34, marker=marker, color=color, label=label, zorder=5)
            for a, xx, yy in zip(anchors, x, y):
                if _positive(a.get("uncertainty_jy")):
                    err = float(a["uncertainty_jy"])*2.99792458e-12/xx
                    ax.errorbar(xx, yy, yerr=[[min(err, .95*yy)], [err]], fmt="none", ecolor=color, lw=.8, capsize=2)
        if ranked:
            best = ranked[0]
            predictions = sorted((p for p in best["predictions"] if p["quality_pass"]), key=lambda p: p["wavelength_um"])
            for scored, marker, label in [(True, "s", "Best sampled model"),
                                           (False, "P", "Model diagnostics (unscored)")]:
                selected = [p for p in predictions if p["scored"] is scored]
                if not selected:
                    continue
                x = np.array([p["wavelength_um"] for p in selected])
                flux = np.array([p["model_flux_jy"] for p in selected])
                ax.scatter(x, flux*2.99792458e-12/x, marker=marker, s=33, facecolors="none", edgecolors="#14844c",
                           linewidths=1.4, label=label, zorder=6)
            residual_points = [p for p in predictions if p["log_residual_dex"] is not None]
            for scored, marker in [(True, "s"), (False, "P")]:
                rows = [p for p in residual_points if p["scored"] is scored]
                res.scatter([p["wavelength_um"] for p in rows], [p["log_residual_dex"] for p in rows],
                            marker=marker, s=33, color="#14844c", facecolors="#14844c" if scored else "none")
            ax.text(.02, .98, f"{best['model_id']}\nRegion-balanced RMS = {best['score_dex']:.4f} dex",
                    transform=ax.transAxes, va="top", fontsize=9)
        else:
            ax.text(.5, .5, "No complete, valid models available for ranking", transform=ax.transAxes, ha="center")
        ax.set(xscale="log", yscale="log", ylabel=r"$\nu F_\nu$ [W m$^{-2}$]", title=title)
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        res.axhline(0, color="0.3", lw=.8)
        res.set(xscale="log", xlabel=r"Wavelength [$\mu$m]", ylabel=r"$\log_{10}(F_{model}/F_{obs})$")
        res.xaxis.set_major_formatter(ScalarFormatter())
        res.grid(alpha=.15)
        fig.supxlabel("Model markers are actual evaluated wavelengths; no dense model spectrum is inferred.", fontsize=8)
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"spectrum_comparison.{suffix}")
        plt.close(fig)

        names = sorted({k for r in results for k in r["parameters"]})
        varying = [k for k in names if len({json.dumps(r["parameters"].get(k), sort_keys=True) for r in results}) > 1]
        if len(varying) > 6:
            warnings.append("Parameter-space figure shows the first six varying axes; all parameters are in ranking.ecsv")
            varying = varying[:6]
        if not varying:
            fig, axis = plt.subplots(figsize=(7, 4), layout="constrained")
            axis.scatter([r["model_index"] for r in ranked], [r["score_dex"] for r in ranked], c="#14844c")
            axis.set(xlabel="Model index", ylabel="Region-balanced RMS [dex]", title=title)
            if not ranked:
                axis.text(.5, .5, "No ranked models", transform=axis.transAxes, ha="center")
        else:
            dim = len(varying)
            fig, axes = plt.subplots(dim, dim, figsize=(max(6, 2.7*dim), max(4, 2.5*dim)), squeeze=False,
                                     layout="constrained")
            values, labels = {}, {}
            for name in varying:
                raw = [r["parameters"].get(name) for r in results]
                if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in raw):
                    values[name] = {r["model_id"]: float(r["parameters"][name]) for r in results}
                else:
                    categories = sorted({str(v) for v in raw})
                    labels[name] = categories
                    values[name] = {r["model_id"]: categories.index(str(r["parameters"].get(name))) for r in results}
            scores = [r["score_dex"] for r in ranked]
            scatter = None
            for i, yname in enumerate(varying):
                for j, xname in enumerate(varying):
                    axis = axes[i, j]
                    if j > i:
                        axis.set_visible(False)
                        continue
                    if i == j:
                        prof: dict[float, list[float]] = defaultdict(list)
                        for row in ranked:
                            prof[values[xname][row["model_id"]]].append(row["score_dex"])
                        xx = sorted(prof)
                        axis.plot(xx, [min(prof[x]) for x in xx], "o-", color="#14844c", ms=4)
                        axis.set_ylabel("Minimum RMS [dex]")
                    else:
                        axis.scatter([values[xname][r["model_id"]] for r in results],
                                     [values[yname][r["model_id"]] for r in results], c="0.8", marker="x", s=17)
                        if ranked:
                            scatter = axis.scatter([values[xname][r["model_id"]] for r in ranked],
                                                   [values[yname][r["model_id"]] for r in ranked], c=scores,
                                                   cmap="viridis_r", s=35, vmin=min(scores), vmax=max(scores) or 1)
                            axis.scatter(values[xname][ranked[0]["model_id"]], values[yname][ranked[0]["model_id"]],
                                         marker="*", s=150, facecolors="none", edgecolors="crimson", linewidths=1.3)
                        axis.set_ylabel(yname)
                    axis.set_xlabel(xname)
                    if xname in labels:
                        axis.set_xticks(range(len(labels[xname])), labels[xname], rotation=35, fontsize=8)
                    if yname in labels and i != j:
                        axis.set_yticks(range(len(labels[yname])), labels[yname], fontsize=8)
                    axis.grid(alpha=.12)
            if scatter is not None:
                fig.colorbar(scatter, ax=[a for a in axes.flat if a.get_visible()], label="Region-balanced RMS [dex]", shrink=.55)
            fig.suptitle(title)
            fig.supxlabel("Grid samples and profile minima, not a posterior. Red star: best; grey crosses: proposed grid.", fontsize=8)
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"parameter_space.{suffix}")
        plt.close(fig)
    return warnings


def analyze_run(run_dir: Path, workers: int = 1) -> dict:
    """Analyze completed compact results; report partial models without ranking them.

    Workers parallelize independent small-file reads only (threads). Scoring and
    plotting are small, serial operations; no image reads or BLAS pools occur.
    """
    run_dir = Path(run_dir).resolve()
    if isinstance(workers, bool) or int(workers) != workers or workers < 1:
        raise AnalysisError("workers must be a positive integer")
    manifest_bytes = (run_dir / "manifest.json").read_bytes()
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = json.loads(manifest_bytes)
    scored = validate_anchors(manifest["anchors"])
    observed = _verified_observation(manifest, run_dir)
    models = manifest["models"]
    if len({str(m["id"]) for m in models}) != len(models) or len({m["index"] for m in models}) != len(models):
        raise AnalysisError("Duplicate model IDs or indices in the manifest")
    template = manifest.get("model_directory", "models/{id}")
    paths = [run_dir / template.format(id=m["id"], index=m["index"]) / "measurements.json" for m in models]
    if workers > 1 and len(paths) > 1:
        with ThreadPoolExecutor(max_workers=min(int(workers), len(paths))) as pool:
            payloads = list(pool.map(_load_measurement, paths))
    else:
        payloads = [_load_measurement(p) for p in paths]
    results = []
    runtime_identities = {}
    for model, (payload, error) in zip(models, payloads):
        result = evaluate_model(model, manifest["anchors"], payload)
        if error:
            result.update(status="invalid", reason=error)
        elif isinstance(payload, dict) and result["status"] != "invalid":
            identity, fingerprint_error = _measurement_fingerprint(payload, manifest_hash, "input_hashes" in manifest)
            if fingerprint_error:
                result.update(status="invalid", reason=fingerprint_error,
                              score_dex=None, chi2_diagnostic=None)
            elif identity is not None:
                runtime_identities[result["model_id"]] = identity
        results.append(result)
    provenance_warnings = []
    if len(set(runtime_identities.values())) > 1:
        reason = "Mixed executable, utilities or backend fingerprints in one run; no cross-runtime ranking is valid"
        provenance_warnings.append(reason)
        for result in results:
            if result["model_id"] in runtime_identities:
                result.update(status="invalid", reason=reason, score_dex=None, chi2_diagnostic=None)
    config = manifest.get("configuration", {})
    smoke = manifest.get("smoke_test") is True or config.get("smoke_test") is True or any(r["smoke_test"] for r in results)
    out = run_dir / "results"
    out.mkdir(parents=True, exist_ok=True)
    _write_tables(results, out)
    warnings = provenance_warnings + _plots(manifest, results, observed, out, smoke)
    ranked = sorted((r for r in results if r["status"] == "ranked"), key=lambda r: (r["score_dex"], r["model_index"]))
    best = {k: v for k, v in ranked[0].items() if k != "predictions"} if ranked else None
    summary = {"schema_version": 1, "run_id": manifest.get("run_id", run_dir.name),
               "manifest_sha256": manifest_hash,
               "smoke_test": smoke, "interpretation": "SMOKE TEST, not scientific inference" if smoke else "Grid screen, not a posterior",
               "metric": "sqrt(mean_regions(mean_anchors(log10(model/observed)^2)))",
               "normalization_fitted": False, "expected_scored_anchor_ids": [str(a["id"]) for a in scored],
               "n_scored_anchors": len(scored), "n_regions": len({str(a["region"]) for a in scored}),
               "model_count": len(models), "ranked_model_count": len(ranked), "excluded_model_count": len(models)-len(ranked),
               "status_counts": dict(Counter(r["status"] for r in results)), "best_model": best,
               "chi2_note": "Unreduced diagonal-error diagnostic only; null if any scored uncertainty is invalid. No reduced chi-square or posterior inferred.",
               "excluded_models": [{"model_id": r["model_id"], "status": r["status"], "reason": r["reason"]} for r in results if r["status"] != "ranked"],
               "warnings": warnings}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary
