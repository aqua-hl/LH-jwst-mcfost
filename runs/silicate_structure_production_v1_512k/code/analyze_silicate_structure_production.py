#!/usr/bin/env python3
"""Receipt-validated, predeclared mass/size/material/envelope comparisons.

All declared pairs are reported, independently of fit. Feature/continuum ratios
are descriptive emergent-flux diagnostics, not optical depths or significance.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
import silicate_structure_design as design

_SPEC = importlib.util.spec_from_file_location(
    "_structure_production_analysis", Path(__file__).with_name("analyze_extinction_ice_production.py"))
production = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(production)
require = production.require


def predeclared_contrasts(catalogue):
    """Declare 75 comparisons using physical coordinates, never model scores."""
    require(catalogue == design.catalogue(), "Campaign catalogue differs from the frozen 60-model design")
    grid = {(r["parameters"]["inclination_deg"], r["reference_mass_msun"],
             r["parameters"]["envelope_amax_um"]): r
            for r in catalogue if r["role"] == "mass_size"}
    masses = sorted({key[1] for key in grid})
    inclinations = sorted({key[0] for key in grid})
    contrasts = []

    def add(kind, reference, changed):
        contrasts.append({"contrast_id": f"{kind}_{reference['index']:03d}_{changed['index']:03d}",
            "kind": kind, "reference_index": reference["index"], "changed_index": changed["index"],
            "inclination_deg": changed["parameters"]["inclination_deg"],
            "reference_column_mass_msun": reference["reference_mass_msun"],
            "changed_column_mass_msun": changed["reference_mass_msun"],
            "reference_actual_mass_msun": reference["parameters"]["envelope_dust_mass_msun"],
            "changed_actual_mass_msun": changed["parameters"]["envelope_dust_mass_msun"],
            "reference_material_id": reference["material_id"], "changed_material_id": changed["material_id"],
            "reference_amax_um": reference["parameters"]["envelope_amax_um"],
            "changed_amax_um": changed["parameters"]["envelope_amax_um"],
            "reference_density_exponent": reference["parameters"]["envelope_density_exponent"],
            "changed_density_exponent": changed["parameters"]["envelope_density_exponent"]})

    for inclination in inclinations:
        for mass in masses:
            add("size", grid[inclination, mass, .4], grid[inclination, mass, 1.])
        for amax in (.4, 1.):
            for low, high in zip(masses, masses[1:]):
                add("mass", grid[inclination, high, amax], grid[inclination, low, amax])
        for mass in masses[:-1]:
            add("mass_size_compensation", grid[inclination, masses[-1], .4], grid[inclination, mass, 1.])
    for changed in catalogue:
        if changed["role"] not in {"material", "density_profile"}:
            continue
        reference = catalogue[changed["parent_index"]]
        require(reference["role"] == "mass_size" and reference["material_id"] == "draine"
                and reference["parameters"]["envelope_amax_um"] == .4,
                "Matched reference must be the Draine small-grain control")
        add(changed["role"], reference, changed)
    require(Counter(r["kind"] for r in contrasts) == {
        "size": 12, "mass": 18, "mass_size_compensation": 9, "material": 18, "density_profile": 18},
        "Predeclared comparison counts differ")
    return contrasts


def feature_shapes(flux, bands):
    """Raw core/shoulder ratios with a power-law interpolated local reference.

Band midpoints are only coordinates for this descriptive construction. These
wide-band emergent flux ratios do not measure an absorption optical depth.
"""
    require(len(flux) == len(bands) and np.all(np.isfinite(flux)) and np.all(np.array(flux) > 0),
            "Feature fluxes must be finite and positive")
    lookup = {band["id"]: (float(value), .5 * (band["lower_um"] + band["upper_um"]))
              for value, band in zip(flux, bands)}
    require(len(lookup) == len(bands), "Duplicate broad-band identifiers")
    result = {}
    for feature, core, left, right in (("silicate_9p7", "s02", "s01", "s03"),
            *((f"ice_{band}", band, "c07", "c08") for band in ("h01", "h02", "h03", "h04"))):
        fcore, xcore = lookup[core]
        fleft, xleft = lookup[left]
        fright, xright = lookup[right]
        require(xleft < xcore < xright, f"Invalid feature continuum bracket: {feature}")
        fraction = float(np.log(xcore / xleft) / np.log(xright / xleft))
        continuum = float(np.exp((1 - fraction) * np.log(fleft) + fraction * np.log(fright)))
        result[feature] = {"core_band": core, "left_band": left, "right_band": right,
            "interpolation_fraction_log_wavelength": fraction,
            "interpolated_continuum_flux_jy": continuum, "core_flux_jy": fcore,
            "core_over_continuum": fcore / continuum,
            "log10_core_over_continuum": float(np.log10(fcore / continuum))}
    return result


def paired_responses(catalogue, models, contract):
    declarations = predeclared_contrasts(catalogue)
    require(len({r["model_index"] for r in models}) == len(models), "Duplicate analyzed model indices")
    lookup = {r["model_index"]: r for r in models}
    bands = contract["bands"]
    require(len(bands) == 19 and len({b["id"] for b in bands}) == 19, "Require 19 distinct broad bands")
    observed = np.array([b["observed_flux_jy"] for b in bands])
    observed_shapes = feature_shapes(observed, bands)
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
        flux = np.array([r["band_flux_jy"] for r in available])
        require(flux.shape == (2, len(bands)) and np.all(np.isfinite(flux)) and np.all(flux > 0),
                "Paired model flux must be finite and positive")
        residual = flux / observed - 1
        residual_dex = np.log10(flux / observed)
        metadata = {**declaration, "reference_model_id": available[0]["model_id"],
                    "changed_model_id": available[1]["model_id"]}
        pair_bands = []
        for j, band in enumerate(bands):
            row = {**metadata, "band_id": band["id"], "block": band["block"],
                "lower_um": band["lower_um"], "upper_um": band["upper_um"],
                "observed_flux_jy": float(observed[j]),
                "reference_flux_jy": float(flux[0, j]), "changed_flux_jy": float(flux[1, j]),
                "changed_over_reference_flux_ratio": float(flux[1, j] / flux[0, j]),
                "changed_over_reference_flux_dex": float(np.log10(flux[1, j] / flux[0, j])),
                "reference_residual_percent": float(100 * residual[0, j]),
                "changed_residual_percent": float(100 * residual[1, j]),
                "reference_residual_dex": float(residual_dex[0, j]),
                "changed_residual_dex": float(residual_dex[1, j]),
                "absolute_residual_improvement_percentage_points": float(100 * (abs(residual[0, j]) - abs(residual[1, j]))),
                "reference_quadrature_unresolved": band["id"] in available[0].get("unresolved_quadrature_bands", []),
                "changed_quadrature_unresolved": band["id"] in available[1].get("unresolved_quadrature_bands", [])}
            band_rows.append(row)
            pair_bands.append(row)
        shapes = [feature_shapes(vector, bands) for vector in flux]
        pair_features = []
        for feature_id, obs in observed_shapes.items():
            reference, changed = [s[feature_id] for s in shapes]
            ref_residual = reference["log10_core_over_continuum"] - obs["log10_core_over_continuum"]
            new_residual = changed["log10_core_over_continuum"] - obs["log10_core_over_continuum"]
            involved = {obs[key] for key in ("core_band", "left_band", "right_band")}
            row = {**metadata, "feature_id": feature_id,
                "core_band": obs["core_band"], "left_band": obs["left_band"], "right_band": obs["right_band"],
                "observed_core_over_continuum": obs["core_over_continuum"],
                "reference_core_over_continuum": reference["core_over_continuum"],
                "changed_core_over_continuum": changed["core_over_continuum"],
                "reference_log_shape_residual_dex": ref_residual,
                "changed_log_shape_residual_dex": new_residual,
                "absolute_shape_residual_improvement_dex": abs(ref_residual) - abs(new_residual),
                "reference_quadrature_unresolved": bool(involved.intersection(available[0].get("unresolved_quadrature_bands", []))),
                "changed_quadrature_unresolved": bool(involved.intersection(available[1].get("unresolved_quadrature_bands", [])))}
            feature_rows.append(row)
            pair_features.append(row)
        regions = {}
        for region, blocks in (("nir_all", {"nir_continuum", "nir_ice"}),
                ("nir_continuum", {"nir_continuum"}), ("nir_ice", {"nir_ice"}), ("miri", {"miri"})):
            selected = np.array([b["block"] in blocks for b in bands])
            require(np.any(selected), f"Missing bands in {region}")
            regions[region] = {"number_of_bands": int(selected.sum())}
            for index, label in enumerate(("reference", "changed")):
                regions[region][label + "_fractional_residual_rms_percent"] = float(100 * np.sqrt(np.mean(residual[index, selected]**2)))
                regions[region][label + "_log_residual_rms_dex"] = float(np.sqrt(np.mean(residual_dex[index, selected]**2)))
        complete.append({**metadata, "region_residuals": regions,
            "s02": next(r for r in pair_bands if r["band_id"] == "s02"),
            "feature_shapes": pair_features,
            "baseline_profiled_scores": [r["baseline_score"] for r in available],
            "quadrature_unresolved": any(r["quadrature_unresolved"] for r in available)})
    return {"status": "complete_campaign" if len(complete) == len(declarations) else "partial_campaign",
        "pairs_expected": len(declarations), "pairs_complete": len(complete),
        "pairs_expected_by_kind": dict(Counter(r["kind"] for r in declarations)),
        "pairs_complete_by_kind": dict(Counter(r["kind"] for r in complete)),
        "declarations": declarations, "pairs": complete, "missing_pairs": missing,
        "band_responses": band_rows, "feature_responses": feature_rows,
        "observed_feature_shapes": observed_shapes,
        "flux_convention": "Raw physical signed aperture Fnu at the target distance; no fitted amplitude or segment gain applied.",
        "feature_definition": "Core divided by continuum interpolated as a power law between flanking band means, using band-midpoint wavelengths. Silicate: s02 with s01/s03; each ice band h01-h04 with c07/c08. These descriptive emergent-flux shapes are not optical depths or independent ice-abundance measurements.",
        "residual_definition": "Flux: 100*(model/observed-1). Shape: log10(model core/continuum) minus log10(observed core/continuum). Positive absolute-residual improvement means closer to observed.",
        "rms_definition": "Unweighted RMS of the declared band fractional or log10 residuals within each region; descriptive, without statistical significance."}


def _plots(output, summary, contract):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kinds = ("size", "mass", "mass_size_compensation", "material", "density_profile")
    fig, axes = plt.subplots(5, 1, figsize=(11, 19), constrained_layout=True,
                             gridspec_kw={"height_ratios": [12, 18, 9, 18, 18]})
    bands = [b["id"] for b in contract["bands"]]
    values = {(row["contrast_id"], row["band_id"]): row for row in summary["band_responses"]}
    vmax = max([abs(row["changed_over_reference_flux_dex"]) for row in summary["band_responses"]] + [.05])
    for ax, kind in zip(axes, kinds):
        declarations = [r for r in summary["declarations"] if r["kind"] == kind]
        array = np.array([[values.get((r["contrast_id"], band), {}).get("changed_over_reference_flux_dex", np.nan)
                           for band in bands] for r in declarations])
        color = ax.imshow(np.ma.masked_invalid(array), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_yticks(range(len(declarations)), [f"{r['reference_index']:02d}→{r['changed_index']:02d} / {r['inclination_deg']:g}°" for r in declarations], fontsize=7)
        ax.set_xticks(range(len(bands)), bands, fontsize=8)
        ax.set_title(kind.replace("_", " "), fontsize=11)
        for j, declaration in enumerate(declarations):
            for i, band in enumerate(bands):
                row = values.get((declaration["contrast_id"], band))
                if row and (row["reference_quadrature_unresolved"] or row["changed_quadrature_unresolved"]):
                    ax.plot(i, j, ".", color="black", markersize=2)
        for boundary in (8.5, 12.5):
            ax.axvline(boundary, color="black", lw=.6)
    fig.colorbar(color, ax=axes.tolist(), label="log₁₀(changed / reference raw band flux)", shrink=.45)
    fig.suptitle(f"All predeclared contrasts: {summary['pairs_complete']}/{summary['pairs_expected']} complete\n"
                 "White cells: unavailable pair; dots: quadrature warning; model indices identify comparisons", fontsize=12)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"matched_response_map.{suffix}", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)
    colors = {50.: "#0072B2", 60.: "#009E73", 70.: "#D55E00"}
    for ax, kind in zip(axes.flat, kinds):
        for pair in summary["pairs"]:
            if pair["kind"] != kind:
                continue
            nir = pair["region_residuals"]["nir_all"]
            x = nir["reference_log_residual_rms_dex"] - nir["changed_log_residual_rms_dex"]
            feature = next(r for r in pair["feature_shapes"] if r["feature_id"] == "silicate_9p7")
            y = feature["absolute_shape_residual_improvement_dex"]
            color = colors[pair["inclination_deg"]]
            ax.scatter(x, y, s=30, edgecolors=color, facecolors="none" if pair["quadrature_unresolved"] else color)
            ax.annotate(f"{pair['reference_index']}→{pair['changed_index']}", (x, y), xytext=(3, 2), textcoords="offset points", fontsize=5)
        ax.axhline(0, color="grey", lw=.7)
        ax.axvline(0, color="grey", lw=.7)
        ax.set_title(kind.replace("_", " "))
        ax.set_xlabel("NIR log-residual RMS improvement (dex)", fontsize=8)
        ax.set_ylabel("Silicate shape residual improvement (dex)", fontsize=8)
        ax.grid(alpha=.15)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0, .9, "Positive = closer to observations\nBoth axes descriptive; no significance\n\nBlue: 50°; green: 60°; orange: 70°\nOpen: quadrature warning in either model\nLabels: reference → changed model index\n\nEvery available predeclared pair shown.\nSilicate shape: s02 / interpolated(s01, s03).\nA better shape alone is not a successful fit.", va="top", fontsize=10)
    fig.suptitle("Near-IR agreement and 9.7-µm core/shoulder shape — raw physical fluxes")
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"shape_tradeoffs.{suffix}", dpi=200)
    plt.close(fig)


def analyze_campaign(run, output=None, make_plot=True):
    run = Path(run).resolve()
    output = Path(output).resolve() if output else run / "results"
    experiment = production.read_json(run / "production_experiment.json")
    declarations = predeclared_contrasts(experiment["catalogue"])
    if "predeclared_contrasts" in experiment:
        require(experiment["predeclared_contrasts"] == declarations, "Frozen predeclared contrasts changed")
    # Generic analysis verifies frozen inputs, parameters, runtime identity,
    # signed total-I/aperture receipts and the 19-band integration contract.
    generic = production.analyze_run(run, output, make_plot=False)
    contract = production.read_json(run / "inputs/production_observation_contract.json")
    require(len(contract["probes"]) == 98, "Campaign requires all 98 wavelength probes")
    comparison = paired_responses(experiment["catalogue"], generic["models"], contract)
    material_note = "This campaign fixes supplied H2O30K at 4% by mass in separate DHS grains; silicate material, size, envelope mass and density profile are tested in declared matched arms. There is no bare/ice abundance ladder."
    generic["limitations"] = [value for value in generic["limitations"]
        if value != "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used."] + [material_note]
    generic["controlled_contrasts"] = {"selection": "All 75 predeclared comparisons are reported in campaign_summary.json",
        "ice_responses": [], "ice_mass_svd": [], "missing_contrasts": comparison["missing_pairs"],
        "conditional_only": True, "number_of_physical_constraints": None,
        "limitations": "No independent ice-abundance contrast; use the fixed-ice campaign diagnostics."}
    generic["campaign_context"] = {key: experiment.get(key, {}) for key in
        ("dust_prescription", "ice_input_check", "material_input_check", "campaign")}
    (output / "summary.json").write_text(json.dumps(generic, indent=2, sort_keys=True, allow_nan=False) + "\n")
    generic_review = output / "REVIEW.md"
    review = generic_review.read_text().replace(
        "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used.", material_note)
    review = review.replace("All available predeclared bare/ice ladders enter the controlled contrasts, regardless of fit. Full conditional covariance is whitened and a free common amplitude removed before ice/mass response diagnostics. Small or collinear responses are valid outcomes; no count of physical constraints is assigned.",
        "All 75 predeclared mass, size, compensation, material and density-profile comparisons are reported in campaign_summary.json, independently of fit. Ice abundance is fixed; no bare/ice ladder is interpreted.")
    review = "\n".join(
        f"Missing/invalid matched contrasts: {len(comparison['missing_pairs'])}."
        if line.startswith("Missing/invalid matched contrasts:") else line
        for line in review.split("\n"))
    generic_review.write_text(review)
    summary = {"schema_version": 1, "analysis": "predeclared_silicate_structure_production",
        "manifest_sha256": generic["manifest_sha256"], "analyzer_sha256": production.digest(__file__),
        "production_summary_sha256": production.digest(output / "summary.json"),
        "observation_contract_sha256": generic["observation_contract_sha256"],
        "models_total": generic["models_total"], "models_screened": generic["models_screened"],
        "quadrature_flagged_models": generic["quadrature_flagged_models"],
        "dust_prescription": experiment.get("dust_prescription", {}),
        "ice_input_check": experiment.get("ice_input_check", {}),
        "material_input_check": experiment.get("material_input_check", {}),
        "campaign": experiment.get("campaign", {}), **comparison,
        "posterior_intervals_reported": False, "production_convergence_certified": False,
        "physical_cause_uniquely_identified": False, "discrete_geometry_column_validated": False,
        "limitations": [
            "Finite exploratory grid and matched comparisons only: no posterior, significance, uniquely identified physical cause or count of independent constraints.",
            "Density-profile masses are matched by an analytic radial-column prescription. This is not validation of the actual discretized, cavity-bearing MCFOST geometry or its observer sightline optical depth.",
            "The H2O mass fraction and size distribution are fixed. Changing total dust mass also changes absolute ice mass; these runs do not independently constrain ice abundance.",
            "Silicate composition comparisons include the material density and full optical-constant response recorded in the frozen prescription. They are not changes only to the 9.7-micron opacity.",
            "Power-law-interpolated feature/continuum ratios are descriptive emergent-flux shapes. Scattering, emission, broad-band integration and an imperfect continuum can affect them; they are not measured optical depths.",
            "A single image seed supplies no numerical covariance or formal significance. Fresh temperatures are recomputed for every physical model.",
            "Five-versus-three-node integration warnings remain visible without removing bands or relaxing mandatory positive signed total-I and aperture checks.",
            "Historical H2O and supplied silicate coverage/tail policies are recorded, not physically validated by this experiment.",
            "All nine weighting/calibration scenario scores remain in summary.json and the scenario_scores column of ranking.csv. CSV rank/order follows the baseline score only; pair membership does not depend on fit.",
            "Compact receipts and frozen inputs are checked; this analysis does not reopen or rehash raw FITS images."]}
    (output / "campaign_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    production._write_csv(output / "matched_band_responses.csv", summary["band_responses"],
                          ["contrast_id", "kind", "reference_index", "changed_index", "band_id"])
    production._write_csv(output / "matched_feature_shapes.csv", summary["feature_responses"],
                          ["contrast_id", "kind", "reference_index", "changed_index", "feature_id"])
    production._write_csv(output / "missing_pairs.csv", summary["missing_pairs"],
                          ["contrast_id", "kind", "reference_index", "changed_index", "reason"])
    lines = ["# Silicate mass, size, material and envelope campaign", "",
        f"**{summary['status']}**: {summary['models_screened']}/60 models; {summary['pairs_complete']}/75 complete predeclared comparisons.", "",
        "512,000 photons; fresh temperatures; approximately 1 AU/pixel. Fixed supplied H2O30K at 4% mass in separate DHS grains. Every available declared comparison enters, independently of fit.", "",
        "| Comparison | Complete | Expected |", "|---|---:|---:|"]
    for kind, count in summary["pairs_expected_by_kind"].items():
        lines.append(f"| {kind} | {summary['pairs_complete_by_kind'].get(kind, 0)} | {count} |")
    lines += ["", summary["flux_convention"], "", summary["feature_definition"], "",
        summary["residual_definition"], "", summary["rms_definition"], "",
        "Assess NIR continuum, all four ice bands, the silicate core and its shoulders together, with 18–20 µm as a further check. Neither a reduced combined score nor a better normalized feature alone establishes success.", "",
        "The mass-size compensation pairs compare each reduced-mass 1.0-µm Draine case with the 2.25e-4-solar-mass 0.4-µm case at fixed inclination. Adjacent-mass pairs hold size fixed; size pairs hold mass fixed. Material and density-profile pairs use their declared Draine small-grain parents.", "",
        "## Missing comparisons", ""]
    lines += ([f"- {row['contrast_id']}: {row['reason']} ({', '.join(row['states'])})." for row in summary["missing_pairs"]]
              or ["None."])
    lines += ["", "## Limits", "", *[f"- {value}" for value in summary["limitations"]], "",
        "`summary.json` retains every nuisance fit and receipt exclusion. `ranking.csv` is ordered by baseline score, with all nine weighting/calibration scores in its `scenario_scores` column. `matched_band_responses.csv` reports every band; `matched_feature_shapes.csv` records continuum-normalized shapes. The response maps display every declared pair, leaving missing pairs blank. No numerical error bars are inferred.", ""]
    (output / "CAMPAIGN_REVIEW.md").write_text("\n".join(lines))
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
        ("status", "models_screened", "models_total", "pairs_complete", "pairs_expected")}, indent=2))
    return 0 if result["status"] == "complete_campaign" else 2


if __name__ == "__main__":
    raise SystemExit(main())
