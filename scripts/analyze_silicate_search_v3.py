#!/usr/bin/env python3
"""Predeclared profile/column search with seed and cross-campaign checks.

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
COUNTS = {"profile": 18, "column": 6, "carbon": 16, "material": 12, "seed": 2}
CONFIGURATIONS = {"anchor": (-1.5, 1.), "shallow": (-1.25, .86),
                  "flat_high": (-1., .83), "flat_low": (-1., .66)}
MATERIALS = ("draine", "pyroxene_mg50")
LIMITS = {"silicate_shape_dex": .10, "s05_raw_dex": .15, "s06_raw_dex": .15,
          "h02_raw_dex": .15, "h03_raw_dex": .15, "nir_continuum_rms_dex": .20}
BRIDGE_METRICS = ("silicate_shape_dex", "s05_raw_dex", "s06_raw_dex", "nir_continuum_rms_dex")
SUCCESS_RULE = {
    "silicate_9p7": "Absolute log10 residual of s02 / power-law continuum(s01,s03) <= 0.10 dex",
    "miri_18_20": "Each RAW s05 and s06 log10(model/observed) residual <= 0.15 dex in absolute value",
    "ice_2p95_3p11_raw": "Each RAW h02 and h03 log10(model/observed) residual <= 0.15 dex in absolute value",
    "nir_continuum": "Unweighted RMS of nine RAW c01-c09 log10 residuals < 0.20 dex",
    "joint": "One run meets all four; primary physical configurations and seed replicas are counted separately",
    "marginal": "A passing run is marginal when any criterion margin is no larger than the conservative paired-seed absolute difference assigned to that diagnostic",
    "ice_and_h04": "Normalized ice shapes and h04 remain reported; neither replaces the raw h02/h03 success criterion, and h04 is outside success",
    "scope": "Descriptive thresholds, not statistical acceptance regions or posterior intervals; unresolved quadrature remains flagged",
}


def ice_fraction(record):
    return record.get("ice_mass_fraction", record.get("derived_ice_mass_fraction"))


def predeclared_contrasts(catalogue):
    """Fix all 54 comparisons without using any flux or fitted score."""
    require(len(catalogue) == 26 and [r["index"] for r in catalogue] == list(range(26)),
            "Require the indexed 26-run v3 catalogue")
    expected = [(config, material, carbon) for config in CONFIGURATIONS
                for material in MATERIALS for carbon in (.15, .2, .25)]
    require([(r["configuration_id"], r["material_id"], r["carbon_mass_fraction"])
             for r in catalogue[:24]] == expected, "Primary catalogue coordinates or ordering changed")
    for row in catalogue:
        p = row["parameters"]
        require((row["density_exponent"], row["column_factor"]) == CONFIGURATIONS[row["configuration_id"]],
                "Profile/column configuration changed")
        require(row["inclination_deg"] == p["inclination_deg"] == 70.
                and p["envelope_density_exponent"] == row["density_exponent"]
                and p["envelope_carbon_mass_fraction"] == row["carbon_mass_fraction"]
                and p["envelope_ice_mass_fraction"] == ice_fraction(row), "Catalogue and physical parameters differ")
        require(production.positive(p["envelope_dust_mass_msun"]), "Invalid matched envelope mass")
        require(row["random_seed"] == (43001 if row["index"] < 24 else 43002), "Wrong primary or replica seed")
        require(row.get("replica_of") == (None if row["index"] < 24 else {24: 7, 25: 10}[row["index"]]),
                "Replica mapping changed")
    for index, parent in ((24, 7), (25, 10)):
        require(catalogue[index]["parameters"] == catalogue[parent]["parameters"], "Seed pair physical parameters differ")
        require(all(catalogue[index][key] == catalogue[parent][key] for key in
                    ("configuration_id", "material_id", "density_exponent", "column_factor", "carbon_mass_fraction")),
                "Seed pair metadata differ")
    grid = {(r["configuration_id"], r["material_id"], r["carbon_mass_fraction"]): r for r in catalogue[:24]}
    declarations = []

    def add(kind, reference, changed):
        item = {"contrast_id": f"{kind}_{reference['index']:03d}_{changed['index']:03d}",
                "kind": kind, "reference_index": reference["index"], "changed_index": changed["index"], "inclination_deg": 70.}
        for label, record in (("reference", reference), ("changed", changed)):
            item.update({label + "_" + key: record[key] for key in
                         ("configuration_id", "material_id", "density_exponent", "column_factor", "carbon_mass_fraction", "random_seed")})
            item[label + "_ice_mass_fraction"] = ice_fraction(record)
            item[label + "_dust_mass_msun"] = record["parameters"]["envelope_dust_mass_msun"]
        declarations.append(item)

    for config in tuple(CONFIGURATIONS)[1:]:
        for material in MATERIALS:
            for carbon in (.15, .2, .25):
                add("profile", grid["anchor", material, carbon], grid[config, material, carbon])
    for material in MATERIALS:
        for carbon in (.15, .2, .25):
            add("column", grid["flat_high", material, carbon], grid["flat_low", material, carbon])
    for config in CONFIGURATIONS:
        for material in MATERIALS:
            for carbon in (.15, .25):
                add("carbon", grid[config, material, .2], grid[config, material, carbon])
        for carbon in (.15, .2, .25):
            add("material", grid[config, "draine", carbon], grid[config, "pyroxene_mg50", carbon])
    for index in (24, 25):
        add("seed", catalogue[catalogue[index]["replica_of"]], catalogue[index])
    require(Counter(r["kind"] for r in declarations) == COUNTS, "Wrong predeclared comparison counts")
    return declarations


def _flux_diagnostics(flux, contract):
    bands, observed, observed_shapes = _observed(contract)
    flux = np.asarray(flux, dtype=float)
    shapes = feature_shapes(flux, bands)
    raw = {b["id"]: float(np.log10(f/o)) for b, f, o in zip(bands, flux, observed)}
    shape_residuals = {key: value["log10_core_over_continuum"] - observed_shapes[key]["log10_core_over_continuum"]
                       for key, value in shapes.items()}
    rms = float(np.sqrt(np.mean([raw[f"c{i:02d}"]**2 for i in range(1, 10)])))
    metrics = {"silicate_shape_dex": shape_residuals["silicate_9p7"], "nir_continuum_rms_dex": rms,
               **{key + "_raw_dex": raw[key] for key in ("s05", "s06", "h02", "h03")}}
    margins = {key: limit - abs(metrics[key]) for key, limit in LIMITS.items()}
    checks = {"silicate_9p7": margins["silicate_shape_dex"] >= -1e-12,
              "miri_18_20": all(margins[k + "_raw_dex"] >= -1e-12 for k in ("s05", "s06")),
              "ice_2p95_3p11_raw": all(margins[k + "_raw_dex"] >= -1e-12 for k in ("h02", "h03")),
              "nir_continuum": rms < .20}
    return dict(band_residuals_dex=raw, feature_shapes=shapes, feature_log_residuals_dex=shape_residuals,
                nir_continuum_log_rms_dex=rms, criterion_metrics=metrics, criterion_margins_dex=margins,
                success_conditions=checks, joint_success=all(checks.values()),
                normalized_ice_within_0p15_dex_diagnostic=all(abs(shape_residuals["ice_"+key]) <= .15+1e-12 for key in ("h02", "h03")))


def model_diagnostics(catalogue, models, contract):
    predeclared_contrasts(catalogue)
    lookup = _lookup(models)
    rows = []
    for record in catalogue:
        result = lookup.get(record["index"])
        metadata = {key: record[key] for key in ("index", "configuration_id", "material_id", "carbon_mass_fraction",
                    "density_exponent", "column_factor", "inclination_deg", "random_seed", "replica_of")}
        metadata.update(ice_mass_fraction=ice_fraction(record), is_primary=record["replica_of"] is None,
                        dust_mass_msun=record["parameters"]["envelope_dust_mass_msun"],
                        model_id=result["model_id"] if result else None)
        if not result or result["status"] != "screened":
            rows.append({**metadata, "status": result["status"] if result else "absent", "joint_success": None,
                         "reason": result.get("reason", "") if result else "No analyzed receipt"})
            continue
        require(result["parameters"] == record["parameters"], "Analyzed parameters differ from frozen catalogue")
        rows.append({**metadata, **_flux_diagnostics(result["band_flux_jy"], contract), "status": "screened",
            "quadrature_unresolved": result["quadrature_unresolved"],
            "unresolved_quadrature_bands": result.get("unresolved_quadrature_bands", []),
            "baseline_profiled_score_diagnostic_only": result["baseline_score"]})
    return rows


def seed_scatter(diagnostics):
    lookup = {r["index"]: r for r in diagnostics}
    pairs = []
    for parent, replica in ((7, 24), (10, 25)):
        a, b = lookup[parent], lookup[replica]
        if a["status"] != "screened" or b["status"] != "screened":
            continue
        pairs.append(dict(material_id=a["material_id"], primary_index=parent, replica_index=replica,
            criterion_absolute_paired_differences_dex={key: abs(b["criterion_metrics"][key]-a["criterion_metrics"][key]) for key in LIMITS},
            band_absolute_paired_log_differences_dex={key: abs(b["band_residuals_dex"][key]-a["band_residuals_dex"][key]) for key in a["band_residuals_dex"]},
            feature_absolute_paired_log_differences_dex={key: abs(b["feature_log_residuals_dex"][key]-a["feature_log_residuals_dex"][key]) for key in a["feature_log_residuals_dex"]},
            quadrature_unresolved=a["quadrature_unresolved"] or b["quadrature_unresolved"]))
    complete = len(pairs) == 2
    conservative = {key: max(p["criterion_absolute_paired_differences_dex"][key] for p in pairs) for key in LIMITS} if complete else None
    criterion_bands = [f"c{i:02d}" for i in range(1, 10)] + ["h02", "h03", "s01", "s02", "s03", "s05", "s06"]
    exceedances = [{"material_id": p["material_id"], "band_id": key, "absolute_paired_difference_dex": p["band_absolute_paired_log_differences_dex"][key]}
                   for p in pairs for key in criterion_bands if p["band_absolute_paired_log_differences_dex"][key] > .03]
    return dict(status="complete" if complete else "incomplete", pairs=pairs,
        conservative_criterion_scatter_dex=conservative, criterion_band_exceeds_0p03_dex=exceedances,
        small_responses_unresolved_warning=bool(exceedances),
        definition="Absolute paired log differences between seeds 43001 and 43002; no division by sqrt(2), standard deviation, covariance or confidence interval is inferred",
        transfer_assumption="For marginality and the bridge, each diagnostic uses the larger of the Draine and pyroxene differences measured only at shallow/C20%. Applying it across materials, carbon levels, columns and profiles is an explicit conservative transfer assumption, not a demonstrated noise bound.")


def apply_marginality(diagnostics, scatter):
    spread = scatter["conservative_criterion_scatter_dex"]
    for row in diagnostics:
        if row["status"] != "screened":
            continue
        row["seed_scatter_available"] = spread is not None
        row["criterion_margins_within_seed_scatter"] = ([key for key in LIMITS
            if abs(row["criterion_margins_dex"][key]) <= spread[key]+1e-12] if spread is not None else None)
        row["marginal_success"] = (bool(row["criterion_margins_within_seed_scatter"]) if row["joint_success"] and spread is not None
                                   else None if row["joint_success"] else False)
        row["success_provisional_due_to_quadrature"] = row["joint_success"] and row["quadrature_unresolved"]


def bridge_check(diagnostics, reference, contract, scatter):
    """Require the new anchor to fall inside the two v2 controls plus seed scatter."""
    controls = reference["bridge_models"]
    require([r["model_index"] for r in controls] == [15, 17], "Require v2 bridge models 15 and 17 in order")
    for model, ice in zip(controls, (.08, .12)):
        p = model["parameters"]
        require(model["status"] == "screened" and p["inclination_deg"] == 70.
                and p["envelope_density_exponent"] == -1.5 and p["envelope_carbon_mass_fraction"] == .2
                and p["envelope_ice_mass_fraction"] == ice and p["envelope_silicate_file"] == "Draine_Si_sUV.dat",
                "V2 bridge control identity differs")
    values = [_flux_diagnostics(model["band_flux_jy"], contract)["criterion_metrics"] for model in controls]
    anchor = diagnostics[1]
    require(anchor["configuration_id"] == "anchor" and anchor["material_id"] == "draine"
            and anchor["carbon_mass_fraction"] == .2 and .08 <= anchor["ice_mass_fraction"] <= .12,
            "V3 bridge anchor does not interpolate the v2 ice bracket")
    spread = scatter["conservative_criterion_scatter_dex"]
    if anchor["status"] != "screened" or spread is None:
        return dict(status="unavailable", passed=None, anchor_index=1, reference_indices=[15, 17],
                    reason="The anchor and both seed pairs are required before cross-campaign interpretation", metrics=[])
    rows = []
    for key in BRIDGE_METRICS:
        low, high = min(r[key] for r in values), max(r[key] for r in values)
        actual = anchor["criterion_metrics"][key]
        rows.append(dict(metric=key, v2_range=[low, high], seed_extension_dex=spread[key],
            allowed_range=[low-spread[key], high+spread[key]], actual=actual,
            passed=low-spread[key]-1e-12 <= actual <= high+spread[key]+1e-12))
    passed = all(r["passed"] for r in rows)
    return dict(status="passed" if passed else "failed", passed=passed, anchor_index=1,
        reference_indices=[15, 17], metrics=rows, interpretation="A failed bridge blocks interpretive outcome conclusions; raw fluxes, descriptive criteria and all pairs remain reported")


def declared_predictions(diagnostics):
    primary = [r for r in diagnostics if r["is_primary"]]
    ice = []
    for config in CONFIGURATIONS:
        expected = [r for r in primary if r["configuration_id"] == config]
        available = [r for r in expected if r["status"] == "screened"]
        ice.append(dict(configuration_id=config, expected_models=6, available_models=len(available),
            raw_ice_pass_indices=[r["index"] for r in available if r["success_conditions"]["ice_2p95_3p11_raw"]],
            raw_ice_fail_indices=[r["index"] for r in available if not r["success_conditions"]["ice_2p95_3p11_raw"]],
            residual_ranges_dex={key: [min(r["band_residuals_dex"][key] for r in available), max(r["band_residuals_dex"][key] for r in available)] for key in ("h02", "h03")} if available else {},
            all_pass=all(r["success_conditions"]["ice_2p95_3p11_raw"] for r in available) if len(available) == 6 else None))
    lookup = {(r["configuration_id"], r["material_id"], r["carbon_mass_fraction"]): r for r in primary}
    nir = []
    for material in MATERIALS:
        for carbon in (.15, .2, .25):
            rows = [lookup[config, material, carbon] for config in ("shallow", "flat_high", "flat_low")]
            valid = all(r["status"] == "screened" for r in rows)
            nir.append(dict(material_id=material, carbon_mass_fraction=carbon, indices=[r["index"] for r in rows],
                configuration_pass={r["configuration_id"]: r["success_conditions"]["nir_continuum"] if r["status"] == "screened" else None for r in rows},
                shallow_and_at_least_one_flat_pass=(rows[0]["success_conditions"]["nir_continuum"] and any(r["success_conditions"]["nir_continuum"] for r in rows[1:])) if valid else None))
    return dict(ice_column_prediction_by_configuration=ice, nir_rematch_by_composition=nir,
                ice_prediction_all_24_pass=all(r["all_pass"] for r in ice) if all(r["all_pass"] is not None for r in ice) else None,
                nir_rematch_all_compositions_pass=all(r["shallow_and_at_least_one_flat_pass"] for r in nir) if all(r["shallow_and_at_least_one_flat_pass"] is not None for r in nir) else None)


def campaign_outcomes(diagnostics, comparison, bridge, scatter, predictions):
    """Return bounded discussion signals; never launch or prescribe a new grid."""
    by_index = {r["index"]: r for r in diagnostics}
    primary = diagnostics[:24]
    profile_rows = []
    for config in tuple(CONFIGURATIONS)[1:]:
        pairs = [r for r in comparison["pairs"] if r["kind"] == "profile" and r["changed_configuration_id"] == config]
        changes = [float(np.mean([by_index[r["changed_index"]]["band_residuals_dex"][key]-by_index[r["reference_index"]]["band_residuals_dex"][key] for key in ("s05", "s06")])) for r in pairs]
        median = float(np.median(changes)) if len(changes) == 6 else None
        profile_rows.append(dict(configuration_id=config, complete_pairs=len(changes), expected_pairs=6,
            individual_pair_mean_changes_dex=changes, median_long_miri_change_dex=median,
            decrease_dex=-median if median is not None else None,
            response_class=("at_least_0p10_dex" if -median >= .10-1e-12 else "0p05_to_0p10_dex" if -median >= .05-1e-12 else "below_0p05_dex") if median is not None else "incomplete"))
    complete = all(r["status"] == "screened" for r in diagnostics)
    gate = bridge["passed"] is True and complete
    raw_success = [r["index"] for r in primary if r.get("joint_success") is True]
    replica_success = [r["index"] for r in diagnostics[24:] if r.get("joint_success") is True]
    any_success = raw_success + replica_success
    replica_only = []
    for index in replica_success:
        row = by_index[index]
        parent = by_index[row["replica_of"]]
        if parent.get("joint_success") is True:
            continue
        # Keep the original primary flux in every profile comparison. A
        # boundary crossed by the alternate seed is still a passing run, but
        # supports only a seed-dependent, marginal configuration claim.
        seed_dependent = parent.get("joint_success") is False
        replica_only.append(dict(replica_index=index, primary_index=parent["index"],
            configuration_id=row["configuration_id"], material_id=row["material_id"],
            carbon_mass_fraction=row["carbon_mass_fraction"], primary_status=parent["status"],
            seed_dependent=seed_dependent, marginal=True if seed_dependent else row.get("marginal_success"),
            interpretation="Replica-only success is seed-dependent and marginal; it prevents a categorical null/stop result" if seed_dependent
                else "Replica passes but its primary is unavailable; the configuration comparison remains incomplete"))
    weak = all(r["decrease_dex"] is not None and r["decrease_dex"] < .05-1e-12 for r in profile_rows)
    contingency = []
    for row in profile_rows:
        relevant = [r for r in primary if r["configuration_id"] == row["configuration_id"] and r["status"] == "screened"
                    and r["success_conditions"]["silicate_9p7"] and r["success_conditions"]["miri_18_20"]]
        if row["decrease_dex"] is not None and row["decrease_dex"] >= .10-1e-12 and relevant and all(
                not (r["success_conditions"]["nir_continuum"] and r["success_conditions"]["ice_2p95_3p11_raw"]) for r in relevant):
            contingency.append(row["configuration_id"])
    flat_limits = [dict(index=r["index"], configuration_id=r["configuration_id"], material_id=r["material_id"],
                        too_deep=r["criterion_metrics"]["silicate_shape_dex"] < -.10-1e-12)
                   for r in primary if r["density_exponent"] == -1. and r["carbon_mass_fraction"] == .25 and r["status"] == "screened"]
    return dict(interpretive_gate_passed=gate,
        interpretation_blocked_reason=None if gate else "Bridge failed/unavailable or required runs are incomplete; no interpretive outcome conclusion is assigned",
        profile_response_by_configuration=profile_rows, raw_primary_joint_success_indices=raw_success,
        raw_replica_joint_success_indices=replica_success, raw_any_run_joint_success_indices=any_success,
        replica_only_success_configurations=replica_only,
        success_statement_supported=bool(any_success) if gate else None,
        success_statement_seed_dependent=bool(replica_only) and not raw_success if gate else None,
        stop_global_profile_search_signal=(weak and not any_success) if gate else None,
        inner_flattened_contingency_discussion_configurations=contingency if gate and not any_success else [],
        partial_response_discussion_configurations=[r["configuration_id"] for r in profile_rows if r["response_class"] == "0p05_to_0p10_dex"] if gate and not any_success else [],
        ice_prediction_failure_signal=(not predictions["ice_prediction_all_24_pass"]) if gate else None,
        flat_high_carbon_depth_limits=flat_limits, seed_small_response_warning=scatter["small_responses_unresolved_warning"],
        actions_launched=False, posterior_or_unique_solution=False,
        rule="Profile response is the median across six declared primary pairs of the mean change in RAW s05/s06 log residuals. Replica fluxes never replace primaries. Either a primary or replica passing the raw criteria prevents a categorical null/stop result; a replica-only pass with a failing primary is a seed-dependent, marginal configuration. Signals are descriptive discussion aids, conditional on the bridge and complete runs; no causal proof or automatic contingency launch.")

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
        "pairs_expected": 54, "pairs_complete": len(complete),
        "pairs_expected_by_kind": COUNTS, "pairs_complete_by_kind": dict(Counter(r["kind"] for r in complete)),
        "declarations": declarations, "pairs": complete, "missing_pairs": missing,
        "band_responses": band_rows, "feature_responses": feature_rows,
        "observed_feature_shapes": observed_shapes}


def _plots(output, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    colors = {"anchor": "#666666", "shallow": "#0072B2", "flat_high": "#D55E00", "flat_low": "#009E73"}
    for row in summary["model_diagnostics"]:
        if row["status"] != "screened":
            continue
        color = colors[row["configuration_id"]]
        marker = "x" if not row["is_primary"] else "o" if row["material_id"] == "draine" else "s"
        index = row["index"]
        series = ([row["criterion_metrics"]["silicate_shape_dex"]],
                  [row["criterion_metrics"][key+"_raw_dex"] for key in ("s05", "s06")],
                  [row["criterion_metrics"][key+"_raw_dex"] for key in ("h02", "h03")],
                  [row["criterion_metrics"]["nir_continuum_rms_dex"]])
        for ax, values in zip(axes.flat, series):
            for offset, value in enumerate(values):
                ax.scatter(index+(offset-.5)*.15, value, color=color, marker=marker, s=24,
                           alpha=.6 if row["quadrature_unresolved"] else 1.)
    for i, (ax, title, threshold) in enumerate(zip(axes.flat,
            ("9.7 µm normalized feature", "18.45 and 19.8 µm raw flux", "2.95 and 3.11 µm raw flux", "Nine-band NIR continuum RMS"),
            (.10, .15, .15, .20))):
        ax.axhline(threshold, color="grey", ls="--", lw=.8)
        if i != 3:
            ax.axhline(-threshold, color="grey", ls="--", lw=.8)
            ax.axhline(0, color="grey", lw=.5)
        ax.axvline(23.5, color="grey", lw=.5)
        ax.set(title=title, xlabel="Declared run index; 24–25 are seed replicas", ylabel="Residual / RMS (dex)", xlim=(-1, 26))
    fig.suptitle(f"All v3 runs; bridge: {summary['bridge_to_v2']['status']}\n"
        "Grey: anchor; blue: shallow; orange/green: flat high/low column. Circles: Draine; squares: pyroxene; crosses: replicas.", fontsize=10)
    for suffix in ("png", "pdf"):
        fig.savefig(output/f"all_model_criteria.{suffix}", dpi=180)
    plt.close(fig)


def analyze_campaign(run, output=None, make_plot=True):
    run = Path(run).resolve()
    output = Path(output).resolve() if output else run/"results"
    experiment = production.read_json(run/"production_experiment.json")
    catalogue = experiment["catalogue"]
    declarations = predeclared_contrasts(catalogue)
    require(experiment.get("predeclared_contrasts") == declarations, "Frozen v3 contrasts changed or missing")
    generic = production.analyze_run(run, output, make_plot=False)
    contract_path = run/"inputs/production_observation_contract.json"
    contract = production.read_json(contract_path)
    require(len(contract["probes"]) == 98, "Require all 98 production probes")
    reference = production.read_json(run/"inputs/v2_reference.json")
    require(reference["observation_contract_sha256"] == production.digest(contract_path),
            "Bridge observation contract differs from v2")
    comparison = paired_responses(catalogue, generic["models"], contract)
    diagnostics = model_diagnostics(catalogue, generic["models"], contract)
    scatter = seed_scatter(diagnostics)
    apply_marginality(diagnostics, scatter)
    bridge = bridge_check(diagnostics, reference, contract, scatter)
    predictions = declared_predictions(diagnostics)
    outcomes = campaign_outcomes(diagnostics, comparison, bridge, scatter, predictions)
    primary = diagnostics[:24]
    replicas = diagnostics[24:]
    screened = [r for r in primary if r["status"] == "screened"]
    primary_success = [r["index"] for r in primary if r.get("joint_success") is True]
    replica_success = [r["index"] for r in replicas if r.get("joint_success") is True]
    material_note = "V3 tests profiles/columns with separate DHS silicate, carbon and supplied H2O30K; ice fraction is derived from a material-specific ice-column target. All 24 primary configurations and two seed replicas are reported separately, regardless of score."
    generic["limitations"] = [value for value in generic["limitations"] if value != "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used."]+[material_note]
    generic["controlled_contrasts"] = dict(selection="All 54 predeclared v3 pairs in campaign_summary.json",
        ice_responses=[], ice_mass_svd=[], missing_contrasts=comparison["missing_pairs"],
        conditional_only=True, number_of_physical_constraints=None,
        limitations="No obsolete generic-ice or bare-silicate ladder is interpreted")
    generic["ranking_role"] = "All nine scenario scores are descriptive only; no fit-selected models, comparisons or outcome rule"
    generic["campaign_context"] = {key: experiment[key] for key in
        ("dust_prescription", "campaign", "stage0_opacity", "opacity_matching", "design_provenance", "ice_column_targets_msun", "observation_policy") if key in experiment}
    (output/"summary.json").write_text(json.dumps(generic, indent=2, sort_keys=True, allow_nan=False)+"\n")
    summary = dict(schema_version=1, analysis="predeclared_silicate_search_v3",
        manifest_sha256=generic["manifest_sha256"], analyzer_sha256=production.digest(__file__),
        production_summary_sha256=production.digest(output/"summary.json"),
        observation_contract_sha256=generic["observation_contract_sha256"], v2_reference_sha256=production.digest(run/"inputs/v2_reference.json"),
        models_total=26, models_screened=generic["models_screened"], primary_models_total=24,
        primary_models_screened=len(screened), seed_replicas_total=2,
        seed_replicas_screened=sum(r["status"] == "screened" for r in replicas),
        quadrature_flagged_models=generic["quadrature_flagged_models"],
        primary_criterion_pass_counts={key: sum(r["success_conditions"][key] for r in screened) for key in
            ("silicate_9p7", "miri_18_20", "ice_2p95_3p11_raw", "nir_continuum")},
        primary_joint_success_indices=primary_success, primary_joint_success_count=len(primary_success),
        replica_joint_success_indices=replica_success, joint_success_indices=primary_success+replica_success,
        joint_success_count=len(primary_success)+len(replica_success),
        model_diagnostics=diagnostics, success_rule=SUCCESS_RULE, seed_scatter=scatter,
        bridge_to_v2=bridge, declared_predictions=predictions, outcomes=outcomes,
        campaign_context=generic["campaign_context"], catalogue_with_opacity_matching=catalogue,
        posterior_intervals_reported=False, physical_cause_uniquely_identified=False,
        abundances_uniquely_identified=False, production_convergence_certified=False,
        discrete_geometry_column_validated=False, **comparison,
        limitations=[
            "The raw four-part criterion is descriptive. Joint success is a tested configuration, not a unique abundance solution, posterior or causal identification.",
            "Bridge failure or unavailable seed comparisons suppress interpretive outcome conclusions, while all raw criteria and comparisons remain visible.",
            "Only two seed differences at shallow/C20% are measured. Their maximum per diagnostic is transferred across configurations and materials as an explicitly conservative heuristic, not a validated upper bound, covariance or standard error.",
            "Quadrature warnings remain in every report; seed marginality and integration sensitivity are different checks, and neither certifies convergence.",
            "The current frozen MIRI reduction is retained; no artifact-corrected replacement was supplied. If its broad 12–20-micron level changes by more than about 0.1 dex, revisit these criteria and projections.",
            "Derived mass, ice fractions and radial-column factors are analytic/interpolated and do not certify equal sightline or scattered-light optical depths on MCFOST's discrete cavity-bearing grid.",
            "The ice targets are inherited from v2; pyroxene's two bracketed crossings were at 60 degrees. All v3 runs are conditional on 70 degrees.",
            "A flat-profile success needs comparison of its envelope mass with independent TMC1A estimates before physical interpretation; no independent mass bound is imposed here.",
            "Normalized ice features and the h04 wing are retained; success uses raw h02/h03 and excludes h04. No ice survival or optical-tail validation is inferred.",
            "Nine nuisance-profiled score scenarios are diagnostics only. Raw success, bridge and profile responses never apply fitted amplitudes or calibration gains.",
            "The inner-flattened contingency is not implemented or launched; discussion signals do not authorize a further campaign.",
            "Compact receipts and frozen inputs are verified; this analysis does not reopen or rehash raw FITS images."])
    (output/"campaign_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)+"\n")
    exports = [
        ("matched_band_responses.csv", summary["band_responses"], ["contrast_id", "kind", "reference_index", "changed_index", "band_id"]),
        ("matched_feature_shapes.csv", summary["feature_responses"], ["contrast_id", "kind", "reference_index", "changed_index", "feature_id"]),
        ("missing_pairs.csv", summary["missing_pairs"], ["contrast_id", "kind", "reference_index", "changed_index", "reason"]),
        ("model_success.csv", diagnostics, ["index", "model_id", "is_primary", "status", "joint_success"]),
        ("opacity_matched_catalogue.csv", catalogue, ["index", "configuration_id", "material_id", "carbon_mass_fraction", "ice_mass_fraction"]),
        ("seed_scatter.csv", scatter["pairs"], ["material_id", "primary_index", "replica_index"]),
        ("bridge_to_v2.csv", bridge["metrics"], ["metric", "actual", "passed"]),
        ("ice_column_prediction.csv", predictions["ice_column_prediction_by_configuration"], ["configuration_id", "all_pass"]),
        ("nir_rematch.csv", predictions["nir_rematch_by_composition"], ["material_id", "carbon_mass_fraction", "shallow_and_at_least_one_flat_pass"]),
        ("profile_outcomes.csv", outcomes["profile_response_by_configuration"], ["configuration_id", "median_long_miri_change_dex", "response_class"])]
    for name, rows, fields in exports:
        production._write_csv(output/name, rows, fields)
    lines = ["# Silicate profile/column search v3", "",
        f"**{summary['status']}**: {len(screened)}/24 primary models, {summary['seed_replicas_screened']}/2 seed replicas, {summary['pairs_complete']}/54 declared comparisons.", "",
        "## Bridge before interpretation", "",
        f"Bridge to v2: **{bridge['status']}**. Interpretive outcome gate: **{outcomes['interpretive_gate_passed']}**.",
        outcomes["interpretation_blocked_reason"] or "The bridge and required runs are available; the following signals remain descriptive.", "",
        f"Primary raw joint successes: {primary_success}. Replica raw joint successes: {replica_success}. These counts remain separate.", "",
        f"Replica-only configurations: {outcomes['replica_only_success_configurations']}. A replica-only pass with a failing primary is seed-dependent and marginal, and prevents a categorical null/stop interpretation.", "",
        "## Declared prediction checks", "",
        f"All 24 ice-column predictions pass: {predictions['ice_prediction_all_24_pass']}.",
        f"Shallow and at least one flat column pass NIR for every composition: {predictions['nir_rematch_all_compositions_pass']}.",
        "Per-configuration ice ranges and per-composition NIR outcomes are retained in the companion CSVs, including every miss.", "",
        "## Seed comparison and marginality", "", scatter["definition"], "", scatter["transfer_assumption"], "",
        f"Criterion-band seed differences above 0.03 dex: {scatter['criterion_band_exceeds_0p03_dex']}", "",
        "## Individual runs", "",
        "| Index | Primary | Status | Silicate shape | Ice raw h02/h03 | Ice shape h02/h03 | Long MIRI s05/s06 | NIR RMS | Joint | Marginal |",
        "|---:|---|---|---:|---|---|---|---:|---|---|"]
    for row in diagnostics:
        if row["status"] != "screened":
            lines.append(f"| {row['index']} | {row['is_primary']} | {row['status']} | — | — | — | — | — | unavailable | unavailable |")
            continue
        raw, shape = row["band_residuals_dex"], row["feature_log_residuals_dex"]
        lines.append(f"| {row['index']} | {row['is_primary']} | screened | {shape['silicate_9p7']:.4f} | {raw['h02']:.4f}/{raw['h03']:.4f} | {shape['ice_h02']:.4f}/{shape['ice_h03']:.4f} | {raw['s05']:.4f}/{raw['s06']:.4f} | {row['nir_continuum_log_rms_dex']:.4f} | {row['joint_success']} | {row['marginal_success']} |")
    lines += ["", "## Declared rules", "", *[f"- {key}: {value}" for key, value in SUCCESS_RULE.items()], "",
        "## Conditional discussion signals", "", outcomes["rule"], "",
        f"Stop-global-profile-search signal: {outcomes['stop_global_profile_search_signal']}.",
        f"Inner-flattened contingency discussion configurations: {outcomes['inner_flattened_contingency_discussion_configurations']}.",
        f"Partial-response discussion configurations: {outcomes['partial_response_discussion_configurations']}.", "",
        "## Limits", "", *[f"- {value}" for value in summary["limitations"]], ""]
    for name in ("REVIEW.md", "CAMPAIGN_REVIEW.md"):
        (output/name).write_text("\n".join(lines))
    if make_plot:
        _plots(output, summary)
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
    print(json.dumps({key: result[key] for key in ("status", "models_screened", "models_total", "pairs_complete",
        "pairs_expected", "primary_joint_success_count")}, indent=2))
    return 0 if result["status"] == "complete_campaign" else 2


if __name__ == "__main__":
    raise SystemExit(main())
