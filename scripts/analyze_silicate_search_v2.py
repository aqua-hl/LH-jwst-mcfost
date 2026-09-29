#!/usr/bin/env python3
"""All-model, predeclared composition diagnostics at matched NIR extinction.

The raw-flux success table and normalized feature diagnostics are deliberately
separate. Profiled scores never select comparisons or determine success.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_extinction_ice_production as production
from analyze_silicate_structure_production import feature_shapes

require = production.require
COUNTS = {"carbon": 18, "ice": 16, "material": 8, "profile": 6}
SUCCESS_RULE = {
    "silicate_9p7": "abs(log10(model core/continuum / observed core/continuum)) <= 0.10 dex",
    "miri_18_20": "Each raw s05 and s06 flux: abs(log10(model/observed)) <= 0.15 dex",
    "ice_2p95_3p11": "Each raw h02 and h03 flux: abs(log10(model/observed)) <= 0.15 dex",
    "nir_continuum": "Unweighted RMS log10(model/observed) across all nine c01-c09 bands < 0.20 dex",
    "joint": "One individual model meets all four conditions; no combining successes from different models",
    "ice_interpretation": "The design emphasizes normalized ice features but its success table specifies the 2.95/3.11-micron bands without saying normalized. This implementation interprets those two thresholds as RAW flux, and prominently reports normalized h02/h03 residuals and their separate diagnostic threshold. Normalized ice agreement is not silently substituted into the declared joint rule.",
    "observed_silicate_reference": "Computed from frozen observed band means; 0.463 in the design is a rounded summary, not a replacement datum",
    "h04": "The 3.36-micron red wing is reported in full but excluded from the declared success rule",
    "quadrature": "Retain and flag unresolved quadrature; any success with a warning remains provisional",
    "statistics": "Descriptive tolerances, not statistical acceptance regions, a posterior, or unique abundance constraints",
}


def _coordinate(record):
    return (record["material_id"], record["carbon_mass_fraction"],
            record["ice_mass_fraction"], record["density_exponent"], record["inclination_deg"])


def predeclared_contrasts(catalogue):
    """Validate the complete design and declare its 48 contrasts without fit data."""
    require(len(catalogue) == 34 and [r["index"] for r in catalogue] == list(range(34)),
            "Require the indexed 34-model silicate-search-v2 catalogue")
    expected = {( "draine", carbon, ice, -1.5, inc)
                for carbon in (0., .1, .2, .3) for ice in (.04, .08, .12) for inc in (60., 70.)}
    expected |= {(material, .2, ice, exponent, 60.)
                 for material in ("draine", "pyroxene_mg50", "olivine_mg50")
                 for exponent in (-1.5, -1.75) for ice in (.08, .12)}
    grid = {_coordinate(r): r for r in catalogue}
    require(set(grid) == expected and len(grid) == 34, "Catalogue coordinates differ from the predeclared design")
    for record in catalogue:
        p = record["parameters"]
        require(p["inclination_deg"] == record["inclination_deg"]
                and p["envelope_density_exponent"] == record["density_exponent"],
                "Catalogue coordinate differs from physical parameters")
        require(production.positive(p["envelope_dust_mass_msun"]), "Matched dust mass must be finite and positive")
        for key, value in (("envelope_carbon_mass_fraction", record["carbon_mass_fraction"]),
                           ("envelope_ice_mass_fraction", record["ice_mass_fraction"])):
            require(key not in p or p[key] == value, "Catalogue fraction differs from physical parameters")
    contrasts = []

    def add(kind, reference, changed):
        row = {"contrast_id": f"{kind}_{reference['index']:03d}_{changed['index']:03d}",
               "kind": kind, "reference_index": reference["index"], "changed_index": changed["index"],
               "inclination_deg": reference["inclination_deg"]}
        for prefix, record in (("reference", reference), ("changed", changed)):
            for key in ("material_id", "carbon_mass_fraction", "ice_mass_fraction", "density_exponent"):
                row[prefix + "_" + key] = record[key]
            row[prefix + "_dust_mass_msun"] = record["parameters"]["envelope_dust_mass_msun"]
        contrasts.append(row)

    for ice in (.04, .08, .12):
        for inc in (60., 70.):
            for carbon in (.1, .2, .3):
                add("carbon", grid["draine", 0., ice, -1.5, inc], grid["draine", carbon, ice, -1.5, inc])
    for carbon in (0., .1, .2, .3):
        for inc in (60., 70.):
            for ice in (.08, .12):
                add("ice", grid["draine", carbon, .04, -1.5, inc], grid["draine", carbon, ice, -1.5, inc])
    for exponent in (-1.5, -1.75):
        for ice in (.08, .12):
            for material in ("pyroxene_mg50", "olivine_mg50"):
                add("material", grid["draine", .2, ice, exponent, 60.], grid[material, .2, ice, exponent, 60.])
    for material in ("draine", "pyroxene_mg50", "olivine_mg50"):
        for ice in (.08, .12):
            add("profile", grid[material, .2, ice, -1.5, 60.], grid[material, .2, ice, -1.75, 60.])
    require(Counter(r["kind"] for r in contrasts) == COUNTS, "Predeclared contrast count differs")
    return contrasts


def _observed(contract):
    bands = contract["bands"]
    require(len(bands) == 19 and len({b["id"] for b in bands}) == 19, "Require 19 distinct broad bands")
    continuum = [b["id"] for b in bands if b["block"] == "nir_continuum"]
    require(set(continuum) == {f"c{i:02d}" for i in range(1, 10)}, "Require all nine NIR continuum bands")
    flux = np.array([b["observed_flux_jy"] for b in bands])
    return bands, flux, feature_shapes(flux, bands)


def _lookup(models):
    require(len({r["model_index"] for r in models}) == len(models), "Duplicate analyzed model indices")
    return {r["model_index"]: r for r in models}


def model_diagnostics(catalogue, models, contract):
    """Evaluate every model separately using raw flux, never nuisance-adjusted flux."""
    predeclared_contrasts(catalogue)
    bands, observed, observed_shapes = _observed(contract)
    lookup = _lookup(models)
    rows = []
    for record in catalogue:
        result = lookup.get(record["index"])
        metadata = {key: record[key] for key in
                    ("index", "material_id", "carbon_mass_fraction", "ice_mass_fraction", "density_exponent", "inclination_deg")}
        metadata["dust_mass_msun"] = record["parameters"]["envelope_dust_mass_msun"]
        metadata["model_id"] = result["model_id"] if result else None
        if not result or result["status"] != "screened":
            rows.append({**metadata, "status": result["status"] if result else "absent",
                         "reason": result.get("reason", "") if result else "No analyzed receipt",
                         "joint_success": None})
            continue
        require(result["parameters"] == record["parameters"], "Analyzed parameters differ from catalogue")
        flux = np.asarray(result["band_flux_jy"])
        shapes = feature_shapes(flux, bands)
        residuals = {b["id"]: float(np.log10(f / o)) for b, f, o in zip(bands, flux, observed)}
        shape_residuals = {key: value["log10_core_over_continuum"] - observed_shapes[key]["log10_core_over_continuum"]
                           for key, value in shapes.items()}
        rms = float(np.sqrt(np.mean([residuals[band["id"]] ** 2 for band in bands if band["block"] == "nir_continuum"])))
        checks = {"silicate_9p7": abs(shape_residuals["silicate_9p7"]) <= .10 + 1e-12,
                  "miri_18_20": all(abs(residuals[key]) <= .15 + 1e-12 for key in ("s05", "s06")),
                  "ice_2p95_3p11_raw": all(abs(residuals[key]) <= .15 + 1e-12 for key in ("h02", "h03")),
                  "nir_continuum": rms < .20}
        rows.append({**metadata, "status": "screened", "band_residuals_dex": residuals,
            "feature_shapes": shapes, "feature_log_residuals_dex": shape_residuals,
            "nir_continuum_log_rms_dex": rms, "success_conditions": checks,
            "joint_success": all(checks.values()),
            "normalized_ice_within_0p15_dex_diagnostic": all(abs(shape_residuals[f"ice_{key}"]) <= .15 + 1e-12 for key in ("h02", "h03")),
            "quadrature_unresolved": result["quadrature_unresolved"],
            "unresolved_quadrature_bands": result.get("unresolved_quadrature_bands", []),
            "success_provisional_due_to_quadrature": all(checks.values()) and result["quadrature_unresolved"],
            "baseline_profiled_score_diagnostic_only": result["baseline_score"]})
    return rows


def paired_responses(catalogue, models, contract):
    declarations = predeclared_contrasts(catalogue)
    bands, observed, observed_shapes = _observed(contract)
    lookup = _lookup(models)
    complete, missing, band_rows, feature_rows = [], [], [], []
    for declaration in declarations:
        indices = [declaration["reference_index"], declaration["changed_index"]]
        available = [lookup.get(index) for index in indices]
        if any(row is None or row["status"] != "screened" for row in available):
            missing.append({**declaration, "states": [r["status"] if r else "absent" for r in available],
                "reasons": [r.get("reason", "") if r else "No analyzed receipt" for r in available],
                "reason": "Missing or invalid matched model; no paired response computed"})
            continue
        for index, result in zip(indices, available):
            require(result["parameters"] == catalogue[index]["parameters"], "Analyzed parameters differ from catalogue")
        flux = np.asarray([r["band_flux_jy"] for r in available])
        shapes = [feature_shapes(vector, bands) for vector in flux]
        residual_dex = np.log10(flux / observed)
        metadata = {**declaration, "reference_model_id": available[0]["model_id"], "changed_model_id": available[1]["model_id"]}
        for j, band in enumerate(bands):
            band_rows.append({**metadata, "band_id": band["id"], "block": band["block"],
                "lower_um": band["lower_um"], "upper_um": band["upper_um"],
                "observed_flux_jy": float(observed[j]),
                "reference_flux_jy": float(flux[0, j]), "changed_flux_jy": float(flux[1, j]),
                "changed_over_reference_flux_dex": float(np.log10(flux[1, j] / flux[0, j])),
                "changed_over_reference_flux_ratio": float(flux[1, j] / flux[0, j]),
                "reference_residual_dex": float(residual_dex[0, j]), "changed_residual_dex": float(residual_dex[1, j]),
                "absolute_residual_improvement_dex": float(abs(residual_dex[0, j]) - abs(residual_dex[1, j])),
                "reference_quadrature_unresolved": band["id"] in available[0].get("unresolved_quadrature_bands", []),
                "changed_quadrature_unresolved": band["id"] in available[1].get("unresolved_quadrature_bands", [])})
        for feature_id, obs in observed_shapes.items():
            reference, changed = [s[feature_id] for s in shapes]
            residual = [s["log10_core_over_continuum"] - obs["log10_core_over_continuum"] for s in (reference, changed)]
            involved = {obs[key] for key in ("core_band", "left_band", "right_band")}
            feature_rows.append({**metadata, "feature_id": feature_id,
                **{key: obs[key] for key in ("core_band", "left_band", "right_band")},
                "observed_core_over_continuum": obs["core_over_continuum"],
                "reference_core_over_continuum": reference["core_over_continuum"],
                "changed_core_over_continuum": changed["core_over_continuum"],
                "reference_log_shape_residual_dex": residual[0], "changed_log_shape_residual_dex": residual[1],
                "changed_over_reference_shape_dex": residual[1] - residual[0],
                "absolute_shape_residual_improvement_dex": abs(residual[0]) - abs(residual[1]),
                "reference_quadrature_unresolved": bool(involved.intersection(available[0].get("unresolved_quadrature_bands", []))),
                "changed_quadrature_unresolved": bool(involved.intersection(available[1].get("unresolved_quadrature_bands", [])))})
        complete.append({**metadata, "quadrature_unresolved": any(r["quadrature_unresolved"] for r in available),
                         "baseline_profiled_scores_diagnostic_only": [r["baseline_score"] for r in available]})
    return {"status": "complete_campaign" if len(complete) == len(declarations) else "partial_campaign",
        "pairs_expected": 48, "pairs_complete": len(complete),
        "pairs_expected_by_kind": COUNTS, "pairs_complete_by_kind": dict(Counter(r["kind"] for r in complete)),
        "declarations": declarations, "pairs": complete, "missing_pairs": missing,
        "band_responses": band_rows, "feature_responses": feature_rows,
        "observed_feature_shapes": observed_shapes}


def _plots(output, summary, contract):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kinds = tuple(COUNTS)
    bands = [b["id"] for b in contract["bands"]]
    values = {(r["contrast_id"], r["band_id"]): r for r in summary["band_responses"]}
    vmax = max([abs(r["changed_over_reference_flux_dex"]) for r in values.values()] + [.05])
    fig, axes = plt.subplots(4, 1, figsize=(11, 16), constrained_layout=True,
                             gridspec_kw={"height_ratios": list(COUNTS.values())})
    for ax, kind in zip(axes, kinds):
        declarations = [r for r in summary["declarations"] if r["kind"] == kind]
        array = np.array([[values.get((r["contrast_id"], band), {}).get("changed_over_reference_flux_dex", np.nan)
                           for band in bands] for r in declarations])
        color = ax.imshow(np.ma.masked_invalid(array), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_yticks(range(len(declarations)), [f"{r['reference_index']:02d}→{r['changed_index']:02d}" for r in declarations], fontsize=7)
        ax.set_xticks(range(len(bands)), bands, fontsize=8)
        ax.set_title(kind.capitalize(), fontsize=11)
        for j, declaration in enumerate(declarations):
            for i, band in enumerate(bands):
                row = values.get((declaration["contrast_id"], band))
                if row and (row["reference_quadrature_unresolved"] or row["changed_quadrature_unresolved"]):
                    ax.plot(i, j, ".", color="black", markersize=2)
    fig.colorbar(color, ax=axes.tolist(), label="log₁₀(changed / reference raw band flux)", shrink=.5)
    fig.suptitle(f"All declared composition contrasts: {summary['pairs_complete']}/48\nDots: quadrature warnings; blank: unavailable pair")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"matched_response_map.{suffix}", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    available = [r for r in summary["model_diagnostics"] if r["status"] == "screened"]
    for row in available:
        x = row["index"]
        marker = "o" if not row["quadrature_unresolved"] else "x"
        axes[0, 0].scatter(x, row["feature_log_residuals_dex"]["silicate_9p7"], marker=marker, color="#0072B2", s=22)
        for key, color in (("s05", "#0072B2"), ("s06", "#D55E00")):
            axes[0, 1].scatter(x, row["band_residuals_dex"][key], marker=marker, color=color, s=22)
        for key, color in (("h02", "#0072B2"), ("h03", "#D55E00")):
            axes[1, 0].scatter(x, row["band_residuals_dex"][key], marker=marker, color=color, s=22)
            axes[1, 0].scatter(x, row["feature_log_residuals_dex"]["ice_" + key], marker="_", color=color, s=28)
        axes[1, 1].scatter(x, row["nir_continuum_log_rms_dex"], marker=marker, color="#0072B2", s=22)
    titles = ("9.7 µm normalized feature", "18.45 (blue) / 19.8 µm (orange): raw",
              "2.95 (blue) / 3.11 µm (orange): raw; bars: normalized", "Nine-band near-IR continuum RMS")
    for index, (ax, title, threshold) in enumerate(zip(axes.flat, titles, (.10, .15, .15, .20))):
        ax.axhline(threshold, color="grey", ls="--", lw=.8)
        if index != 3:
            ax.axhline(-threshold, color="grey", ls="--", lw=.8)
            ax.axhline(0., color="grey", lw=.5)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Predeclared model index")
        ax.set_ylabel("Log residual / RMS (dex)")
        ax.set_xlim(-1, 34)
    fig.suptitle("Individual-model criteria; all available models shown\nCrosses: quadrature warning; tolerances are descriptive")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"all_model_criteria.{suffix}", dpi=180)
    plt.close(fig)


def analyze_campaign(run, output=None, make_plot=True):
    run = Path(run).resolve()
    output = Path(output).resolve() if output else run / "results"
    experiment = production.read_json(run / "production_experiment.json")
    declarations = predeclared_contrasts(experiment["catalogue"])
    require(experiment.get("predeclared_contrasts") == declarations, "Frozen predeclared contrasts changed or missing")
    generic = production.analyze_run(run, output, make_plot=False)
    contract = production.read_json(run / "inputs/production_observation_contract.json")
    require(len(contract["probes"]) == 98, "Campaign requires all 98 wavelength probes")
    comparison = paired_responses(experiment["catalogue"], generic["models"], contract)
    diagnostics = model_diagnostics(experiment["catalogue"], generic["models"], contract)
    successful = [r["index"] for r in diagnostics if r["joint_success"] is True]
    note = "Supplied H2O30K, carbon and silicate are separate DHS species; all 48 predeclared composition/profile contrasts enter independently of fit. Stage 0 sets masses from opacities before model construction."
    generic["limitations"] = [v for v in generic["limitations"] if v != "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used."] + [note]
    generic["controlled_contrasts"] = {"selection": "All 48 predeclared comparisons in campaign_summary.json",
        "ice_responses": [], "ice_mass_svd": [], "missing_contrasts": comparison["missing_pairs"],
        "conditional_only": True, "number_of_physical_constraints": None,
        "limitations": "No legacy generic-ice or bare-silicate ladder is interpreted; use campaign contrasts."}
    generic["ranking_role"] = "All nine weighting/calibration scenarios are descriptive diagnostics only; no model, pair, success or discussion selection by fit."
    generic["campaign_context"] = {key: experiment[key] for key in
        ("dust_prescription", "ice_input_check", "material_input_check", "campaign", "stage0_opacity", "opacity_matching", "stage0") if key in experiment}
    (output / "summary.json").write_text(json.dumps(generic, indent=2, sort_keys=True, allow_nan=False) + "\n")
    summary = {"schema_version": 1, "analysis": "predeclared_silicate_search_v2",
        "manifest_sha256": generic["manifest_sha256"], "analyzer_sha256": production.digest(__file__),
        "production_summary_sha256": production.digest(output / "summary.json"),
        "observation_contract_sha256": generic["observation_contract_sha256"],
        "models_total": generic["models_total"], "models_screened": generic["models_screened"],
        "quadrature_flagged_models": generic["quadrature_flagged_models"], **comparison,
        "model_diagnostics": diagnostics, "success_rule": SUCCESS_RULE,
        "joint_success_count": len(successful), "joint_success_indices": successful,
        "campaign_context": generic["campaign_context"], "catalogue_with_opacity_matching": experiment["catalogue"],
        "posterior_intervals_reported": False, "production_convergence_certified": False,
        "physical_cause_uniquely_identified": False, "abundances_uniquely_identified": False,
        "discrete_geometry_column_validated": False,
        "limitations": [
            "Joint success describes a tested configuration only. It does not uniquely determine abundances, causal mechanisms, or the number of independent physical constraints.",
            "Joint failure cannot by itself establish envelope temperature structure, carbon size, or optical constants as the cause; those remain hypotheses conditional on this finite grid.",
            "DHS opacity-matched masses and analytic radial-column normalization do not certify identical sightline optical depths on the discrete, cavity-bearing grid.",
            "Normalized features use power-law interpolation between broad-band means; they are descriptive emergent flux ratios, not measured optical depths.",
            "The success table is implemented as normalized silicate, raw ice and 18–20-micron residuals, plus raw NIR continuum RMS. Normalized ice residuals remain prominent diagnostics; the document's ambiguity is recorded rather than silently changing the rule.",
            "The h04 red ice wing remains in every band/shape output and the descriptive score, but lies outside the success rule.",
            "A single image seed gives no numerical covariance; unresolved five-versus-three-node quadrature warnings remain visible and do not select away models.",
            "All nine nuisance-profiled scenarios remain in summary.json/ranking.csv as diagnostics only. Success uses raw physical flux, without a fitted amplitude or segment gain.",
            "Carbon uses the silicate size distribution by assumption; native optical-table coverage/tails are accepted inputs rather than physically validated extrapolations.",
            "Compact receipts and frozen inputs are checked; raw FITS contents are not reopened or rehashed by this analysis."]}
    (output / "campaign_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    for filename, rows, fields in (
        ("matched_band_responses.csv", summary["band_responses"], ["contrast_id", "kind", "reference_index", "changed_index", "band_id"]),
        ("matched_feature_shapes.csv", summary["feature_responses"], ["contrast_id", "kind", "reference_index", "changed_index", "feature_id"]),
        ("missing_pairs.csv", summary["missing_pairs"], ["contrast_id", "kind", "reference_index", "changed_index", "reason"]),
        ("model_success.csv", diagnostics, ["index", "model_id", "status", "joint_success"]),
        ("opacity_matched_catalogue.csv", experiment["catalogue"], ["index", "material_id", "carbon_mass_fraction", "ice_mass_fraction"])):
        production._write_csv(output / filename, rows, fields)
    lines = ["# Silicate composition search v2", "",
        f"**{summary['status']}**: {summary['models_screened']}/34 models and {summary['pairs_complete']}/48 predeclared comparisons.", "",
        f"Joint descriptive success: {len(successful)} models; indices {successful}. Quadrature-flagged cases remain provisional.", "",
        "Every model and declared pair is reported regardless of fit. All nine weighting/calibration scores are diagnostics only; no score selects discussion or success.", "",
        "## Declared success and ice interpretation", "",
        *[f"- **{key}**: {value}" for key, value in SUCCESS_RULE.items()], "",
        "## All individual models", "",
        "Raw h02/h03 ice residuals and normalized ice residuals are shown together. h04 is retained in the machine-readable diagnostics and excluded only from success.", "",
        "| Index | Status | 9.7 shape dex | Ice raw h02/h03 dex | Ice shape h02/h03 dex | 18–20 raw s05/s06 dex | NIR RMS dex | Joint | Quadrature |",
        "|---:|---|---:|---|---|---|---:|---|---|"]
    for row in diagnostics:
        if row["status"] != "screened":
            lines.append(f"| {row['index']} | {row['status']} | — | — | — | — | — | unavailable | — |")
            continue
        raw, shape = row["band_residuals_dex"], row["feature_log_residuals_dex"]
        lines.append(f"| {row['index']} | screened | {shape['silicate_9p7']:.4f} | {raw['h02']:.4f} / {raw['h03']:.4f} | {shape['ice_h02']:.4f} / {shape['ice_h03']:.4f} | {raw['s05']:.4f} / {raw['s06']:.4f} | {row['nir_continuum_log_rms_dex']:.4f} | {row['joint_success']} | {row['quadrature_unresolved']} |")
    lines += ["", "## Comparisons", "", "| Kind | Complete | Expected |", "|---|---:|---:|"]
    lines += [f"| {kind} | {summary['pairs_complete_by_kind'].get(kind, 0)} | {count} |" for kind, count in COUNTS.items()]
    lines += ["", "## Limits", "", *[f"- {v}" for v in summary["limitations"]], "",
              "Stage 0 species opacities, mixture S=κabs(9.7)/κext(2.2), the 18/9.7 absorption ratio and frozen mass matching are retained in the campaign context and opacity_matched_catalogue.csv. They precede all image fits.", ""]
    review = "\n".join(lines)
    (output / "CAMPAIGN_REVIEW.md").write_text(review)
    (output / "REVIEW.md").write_text(review)
    if make_plot:
        _plots(output, summary, contract)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-plots", "--no-plot", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = analyze_campaign(args.run, args.output, not args.no_plots)
    except (OSError, ValueError, TypeError, KeyError, ImportError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps({key: result[key] for key in
        ("status", "models_screened", "models_total", "pairs_complete", "pairs_expected", "joint_success_count")}, indent=2))
    return 0 if result["status"] == "complete_campaign" else 2


if __name__ == "__main__":
    raise SystemExit(main())
