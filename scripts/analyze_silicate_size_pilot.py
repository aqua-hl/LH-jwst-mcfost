#!/usr/bin/env python3
"""Analyze all four fixed-mass grain-size models as matched inclination pairs.

The inherited production analyzer validates frozen inputs and signed-aperture
receipts before any response is calculated. No opacity/structure verdict or
parameter confidence interval is inferred from this conditional experiment.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


_SPEC = importlib.util.spec_from_file_location(
    "_silicate_production_analysis", Path(__file__).with_name("analyze_extinction_ice_production.py"))
production = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(production)
require = production.require


def validate_catalogue(catalogue):
    """Require the declared 2 x 2 experiment; drift is an error, not a result."""
    require(len(catalogue) == 4, "Silicate pilot requires exactly four models")
    require([r["index"] for r in catalogue] == list(range(4)), "Invalid pilot model indices")
    groups = defaultdict(list)
    fixed = []
    for record in catalogue:
        require(record["role"] == "silicate_size_pair", "Unexpected pilot catalogue role")
        params = record["parameters"]
        require(production.positive(params["inclination_deg"]), "Invalid inclination")
        require(math.isclose(params["envelope_dust_mass_msun"], 2.25e-4, rel_tol=1e-12),
                "Pilot mass must be 2.25e-4 solar masses")
        fraction = record["ice_mass_fraction"]
        require(isinstance(fraction, (int, float)) and not isinstance(fraction, bool)
                and math.isfinite(fraction) and 0 <= fraction < 1, "Invalid ice mass fraction")
        require(isinstance(record["contrast_group_id"], str) and bool(record["contrast_group_id"]),
                "Missing pair identifier")
        groups[record["contrast_group_id"]].append(record)
        fixed.append({k: v for k, v in params.items() if k not in {"inclination_deg", "envelope_amax_um"}})
    require(len(groups) == 2, "Require two predeclared inclination pairs")
    require(all(p == fixed[0] for p in fixed), "Other physical parameters drift across the four models")
    require(all(r["ice_mass_fraction"] == catalogue[0]["ice_mass_fraction"] for r in catalogue),
            "Ice mass fraction must remain fixed")
    inclinations = []
    pairs = []
    for group_id, records in sorted(groups.items()):
        require(len(records) == 2 and {r["parameters"]["envelope_amax_um"] for r in records} == {.4, 1.},
                "Each inclination requires amax 0.4 and 1.0 microns")
        incl = {r["parameters"]["inclination_deg"] for r in records}
        require(len(incl) == 1, "Pair members must have the same inclination")
        inclinations.extend(incl)
        pairs.append((group_id, sorted(records, key=lambda r: r["parameters"]["envelope_amax_um"])))
    require(len(set(inclinations)) == 2, "Inclination pairs must be distinct")
    return pairs


def opacity_screen_comparison(screen, ice_mass_fraction):
    """Keep fixed-column absorption separate from fixed-NIR-normalized opacity."""
    baseline = screen["baseline"]
    require(baseline["amax_um"] == .4 and baseline["porosity"] == 0,
            "Screen reference must be compact amax=0.4 grains")
    selected = [r for r in screen["variations"] if r["amax_um"] == 1. and r["porosity"] == 0]
    require(len(selected) == 1, "Missing/duplicate compact amax=1.0 screen row")
    changed = selected[0]
    for row in (baseline, changed):
        require(all(production.positive(row[key]) for key in ("kappa_abs_9p7", "kappa_ext_2p2", "S_9p7_per_nir")),
                "Nonpositive/nonfinite opacity screen values")
        require(math.isclose(row["S_9p7_per_nir"], row["kappa_abs_9p7"] / row["kappa_ext_2p2"], rel_tol=1e-9),
                "Screen opacity-ratio identity fails")
    absorption = changed["kappa_abs_9p7"] / baseline["kappa_abs_9p7"]
    extinction = changed["kappa_ext_2p2"] / baseline["kappa_ext_2p2"]
    return {"material": screen["material"], "material_sha256": screen["sha256"],
            "bare_grain_kappa_abs_9p7_ratio_at_fixed_mass": absorption,
            "bare_grain_kappa_abs_9p7_percent_change_at_fixed_mass": 100 * (absorption - 1),
            "bare_grain_kappa_ext_2p2_ratio_at_fixed_mass": extinction,
            "bare_grain_S_9p7_per_nir_ratio_at_fixed_nir_extinction": absorption / extinction,
            "pilot_ice_mass_fraction": ice_mass_fraction,
            "bare_screen_vs_icy_pilot_mismatch": ice_mass_fraction > 0,
            "interpretation": "The screen is bare-Mie opacity, while the pilot has separate DHS silicate and ice populations. Its S ratio holds NIR extinction fixed by changing column; this pilot holds total dust mass fixed. Neither opacity ratio predicts aperture flux after scattering and fresh radiative equilibrium.",
            "screen_flux_gain_dex_used_as_prediction": False}


def paired_responses(catalogue, models, contract):
    """Use every predeclared pair, without choosing by fit or rescaling flux."""
    pairs = validate_catalogue(catalogue)
    require(len({r["model_index"] for r in models}) == len(models), "Duplicate analyzed models")
    lookup = {r["model_index"]: r for r in models}
    bands = contract["bands"]
    ids = [b["id"] for b in bands]
    require(len(set(ids)) == len(ids) and "s02" in ids, "Unique bands including s02 are required")
    observed = np.array([b["observed_flux_jy"] for b in bands])
    require(np.all(np.isfinite(observed)) and np.all(observed > 0), "Invalid observed band flux")
    comparisons, missing, rows = [], [], []
    for group_id, records in pairs:
        available = [lookup.get(r["index"]) for r in records]
        if any(r is None or r["status"] != "screened" for r in available):
            missing.append({"group_id": group_id, "model_indices": [r["index"] for r in records],
                            "reason": "Missing or invalid matched model; no paired response computed",
                            "states": [r["status"] if r else "absent" for r in available]})
            continue
        for result, record in zip(available, records):
            require(result["parameters"] == record["parameters"], "Screened parameters differ from catalogue")
        flux = np.array([r["band_flux_jy"] for r in available])
        require(flux.shape == (2, len(bands)) and np.all(np.isfinite(flux)) and np.all(flux > 0),
                "Paired model flux must be finite and positive")
        residual = flux / observed - 1
        residual_dex = np.log10(flux / observed)
        pair_rows = []
        for j, band in enumerate(bands):
            flags = [band["id"] in r.get("unresolved_quadrature_bands", []) for r in available]
            row = {"group_id": group_id, "inclination_deg": records[0]["parameters"]["inclination_deg"],
                   "band_id": band["id"], "block": band["block"],
                   "lower_um": band["lower_um"], "upper_um": band["upper_um"],
                   "observed_flux_jy": float(observed[j]),
                   "reference_model_id": available[0]["model_id"], "changed_model_id": available[1]["model_id"],
                   "reference_amax_um": .4, "changed_amax_um": 1.,
                   "reference_flux_jy": float(flux[0, j]), "changed_flux_jy": float(flux[1, j]),
                   "changed_over_reference_flux_ratio": float(flux[1, j] / flux[0, j]),
                   "changed_over_reference_flux_dex": float(np.log10(flux[1, j] / flux[0, j])),
                   "reference_residual_percent": float(100 * residual[0, j]),
                   "changed_residual_percent": float(100 * residual[1, j]),
                   "reference_residual_dex": float(residual_dex[0, j]),
                   "changed_residual_dex": float(residual_dex[1, j]),
                   "absolute_residual_improvement_percentage_points": float(100 * (abs(residual[0, j]) - abs(residual[1, j]))),
                   "absolute_residual_improvement_dex": float(abs(residual_dex[0, j]) - abs(residual_dex[1, j])),
                   "reference_quadrature_unresolved": flags[0], "changed_quadrature_unresolved": flags[1]}
            rows.append(row)
            pair_rows.append(row)
        regions = {}
        for region, blocks in (("nir_all", {"nir_continuum", "nir_ice"}), ("nir_continuum", {"nir_continuum"}),
                               ("nir_ice", {"nir_ice"}), ("miri", {"miri"})):
            selected = np.array([b["block"] in blocks for b in bands])
            require(np.any(selected), f"Missing bands for {region}")
            values = {}
            for k, label in enumerate(("reference", "changed")):
                values[label + "_fractional_residual_rms_percent"] = float(100 * np.sqrt(np.mean(residual[k, selected]**2)))
                values[label + "_log_residual_rms_dex"] = float(np.sqrt(np.mean(residual_dex[k, selected]**2)))
            values["number_of_bands"] = int(selected.sum())
            regions[region] = values
        comparisons.append({"group_id": group_id, "inclination_deg": records[0]["parameters"]["inclination_deg"],
            "model_indices": [r["index"] for r in records], "region_residuals": regions,
            "s02": next(r for r in pair_rows if r["band_id"] == "s02"),
            "miri_bands": [r for r in pair_rows if r["block"] == "miri"],
            "baseline_profiled_scores": [r["baseline_score"] for r in available],
            "quadrature_unresolved": any(r["quadrature_unresolved"] for r in available)})
    return {"status": "complete_pilot" if len(comparisons) == 2 else "incomplete_pilot",
            "pairs_expected": 2, "pairs_complete": len(comparisons), "pairs": comparisons,
            "missing_pairs": missing, "band_responses": rows,
            "flux_convention": "Raw physical aperture Fnu at the target distance, without common-scale or segment-gain adjustments.",
            "residual_definition": "100*(model/observed-1); dex=log10(model/observed). Positive absolute-residual improvement means closer to observed flux.",
            "rms_definition": "Unweighted RMS across declared broad bands within each region, calculated separately in fractional flux and log10 flux; no statistical significance implied."}


def _plot(output, summary, contract):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex="col", constrained_layout=True)
    colors = ("#0072B2", "#D55E00")
    for col, region in enumerate(("NIR", "MIRI")):
        selected = [b for b in contract["bands"] if (b["block"] == "miri") == (col == 1)]
        for band in selected:
            x = .5 * (band["lower_um"] + band["upper_um"])
            axes[0, col].errorbar(x, band["observed_flux_jy"], xerr=.5 * (band["upper_um"] - band["lower_um"]),
                                  fmt="s", color="black", markersize=3, lw=.7)
        for color, pair in zip(colors, summary["pairs"]):
            rows = [r for r in summary["band_responses"] if r["group_id"] == pair["group_id"]
                    and (r["block"] == "miri") == (col == 1)]
            x = [(.5 * (r["lower_um"] + r["upper_um"])) for r in rows]
            for key, marker, size in (("reference", "o", .4), ("changed", "^", 1.)):
                label = f"{pair['inclination_deg']:g}°, silicate amax={size:g} µm"
                axes[0, col].scatter(x, [r[key + "_flux_jy"] for r in rows], color=color, marker=marker, s=24, label=label)
                axes[1, col].scatter(x, [r[key + "_residual_dex"] for r in rows], color=color, marker=marker, s=24)
            axes[2, col].scatter(x, [r["changed_over_reference_flux_dex"] for r in rows], color=color, s=24,
                                 label=f"{pair['inclination_deg']:g}°")
        axes[0, col].set(title=region + " broad bands", yscale="log")
        for row in (1, 2):
            axes[row, col].axhline(0, color="grey", lw=.7)
        axes[2, col].set_xlabel("Wavelength (µm)")
        if not summary["pairs"]:
            axes[1, col].text(.5, .5, "No complete matched pair", ha="center", transform=axes[1, col].transAxes)
        for ax in axes[:, col]:
            ax.grid(alpha=.15)
    axes[0, 0].set_ylabel("Aperture band mean Fν (Jy)")
    axes[1, 0].set_ylabel("log₁₀(model / observed)")
    axes[2, 0].set_ylabel("log₁₀(Fν[1.0 µm] / Fν[0.4 µm])")
    if summary["pairs"]:
        axes[0, 0].legend(fontsize=8)
        axes[2, 0].legend(fontsize=8)
    fig.suptitle(f"Fixed-mass grain-size pilot: {summary['pairs_complete']}/2 matched pairs\n"
                 "Silicate size varies; raw fluxes; black horizontal bars are band widths; no numerical error bars", fontsize=11)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"pilot_plot.{suffix}", dpi=220)
    plt.close(fig)


def analyze_pilot(run, output=None, make_plot=True):
    run = Path(run).resolve()
    output = Path(output).resolve() if output else run / "results"
    # This validates every pinned input and every mandatory receipt check before
    # any paired response. Invalid/negative total-I receipts cannot enter pairs.
    generic = production.analyze_run(run, output, make_plot=False)
    manifest = production.read_json(run / "manifest.json")
    experiment = production.read_json(run / "production_experiment.json")
    contract = production.read_json(run / "inputs/production_observation_contract.json")
    # Generic production statistics remain useful, but its hardcoded material
    # prose describes the older production campaign rather than this pilot.
    material_note = "The actual frozen material and grain prescription is recorded in pilot_summary.json; no material constants are varied within this size-only pilot."
    generic["limitations"] = [value for value in generic["limitations"]
        if value != "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used."] + [material_note]
    generic["controlled_contrasts"]["limitations"] = "No ice-abundance or mass contrasts are defined in this size-only pilot. Use the predeclared pairs in pilot_summary.json."
    generic["pilot_context"] = {"dust_prescription": experiment.get("dust_prescription", {}),
                                "ice_input_check": experiment.get("ice_input_check", {}),
                                "pilot": experiment.get("pilot", {})}
    (output / "summary.json").write_text(json.dumps(generic, indent=2, sort_keys=True, allow_nan=False) + "\n")
    generic_review = output / "REVIEW.md"
    generic_review.write_text(generic_review.read_text().replace(
        "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used.", material_note))
    require(len(contract["bands"]) == 19 and len(contract["probes"]) == 98,
            "Pilot must retain the complete 19-band/98-probe production contract")
    screen_relative = "inputs/silicate_screen_draine.json"
    require(screen_relative in manifest["input_hashes"], "Opacity screen must be a frozen hashed input")
    comparison = paired_responses(experiment["catalogue"], generic["models"], contract)
    screen = opacity_screen_comparison(production.read_json(run / screen_relative),
                                       experiment["catalogue"][0]["ice_mass_fraction"])
    summary = {"schema_version": 1, "analysis": "fixed_mass_silicate_size_pilot",
               "manifest_sha256": generic["manifest_sha256"], "analyzer_sha256": production.digest(__file__),
               "production_summary_sha256": production.digest(output / "summary.json"),
               "opacity_screen_sha256": manifest["input_hashes"][screen_relative],
               "dust_prescription": experiment.get("dust_prescription", {}),
               "ice_input_check": experiment.get("ice_input_check", {}),
               "pilot": experiment.get("pilot", {}),
               "models_screened": generic["models_screened"], "models_total": generic["models_total"],
               **comparison, "opacity_screen_comparison": screen,
               "posterior_intervals_reported": False, "production_convergence_certified": False,
               "controlling_variable_identified": False,
               "limitations": ["This is a conditional grain-size experiment at fixed mass, ice prescription and background parameters; it does not isolate material constants from thermal or scattering responses.",
                   "Improvement motivates further controlled comparisons; a weak response does not establish envelope structure as the unique cause.",
                   "Neither the bare-grain opacity screen nor its fixed-NIR-extinction normalization predicts the fixed-mass emergent aperture-flux change. Grain geometry, mixture and ice-table differences from the screen must also be considered using the frozen dust prescription.",
                   "A single production image seed supplies no numerical covariance or significance; fresh temperatures are part of the coupled physical response.",
                   "Five-versus-three-node integration warnings remain visible and do not remove any band or relax mandatory signed total-I checks.",
                   "The inherited weighted screening scores are secondary descriptive comparisons, not likelihoods, parameter intervals or hypothesis-test significance.",
                   "The frozen ice_input_check records table coverage and the accepted extrapolation policy. Agreement within this pilot does not independently validate unmeasured ice optical constants or their thermal effect."]}
    (output / "pilot_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    production._write_csv(output / "paired_band_responses.csv", summary["band_responses"],
                          ["group_id", "inclination_deg", "band_id", "block", "reference_flux_jy", "changed_flux_jy"])
    lines = ["# Fixed-mass silicate grain-size pilot", "", f"**{summary['status']}**: {summary['pairs_complete']}/2 matched pairs, {summary['models_screened']}/4 valid models.", "",
             "All predeclared pairs enter the comparison, independent of fit. Envelope dust mass is 2.25e-4 solar masses; silicate amax changes from 0.4 to 1.0 µm at each inclination. All other physical parameters are fixed; see the frozen dust prescription for the separate ice size distribution.", "",
             f"H2O coverage policy: `{summary['ice_input_check'].get('coverage_policy', 'not recorded')}`. Native extrapolation required: {summary['ice_input_check'].get('native_extrapolation_required', 'not recorded')}. The input-policy check is distinct from physical validation of the tails.", "",
             summary["flux_convention"], summary["residual_definition"], summary["rms_definition"], "",
             "## Paired results", ""]
    for pair in summary["pairs"]:
        s02 = pair["s02"]
        nir = pair["region_residuals"]["nir_all"]
        lines += [f"- {pair['inclination_deg']:g}°: s02 (9.3–10.1 µm) residual {s02['reference_residual_percent']:+.3f}% → {s02['changed_residual_percent']:+.3f}%; flux response {s02['changed_over_reference_flux_dex']:+.4f} dex. Absolute residual improvement {s02['absolute_residual_improvement_percentage_points']:+.3f} percentage points. NIR fractional RMS {nir['reference_fractional_residual_rms_percent']:.3f}% → {nir['changed_fractional_residual_rms_percent']:.3f}%. Quadrature warning in either model: {pair['quadrature_unresolved']}."]
    for missing in summary["missing_pairs"]:
        lines += [f"- {missing['group_id']}: {missing['reason']} ({', '.join(missing['states'])})."]
    lines += ["", "The CSV reports all 19 bands, including every MIRI shoulder and both NIR regions. Positive absolute-residual improvement means closer to the observation; no improvement threshold or physical-cause classification is imposed.", "",
              "## Opacity-screen context", "",
              f"The bare Draine screen changes 9.7-µm absorption opacity at fixed mass by {screen['bare_grain_kappa_abs_9p7_percent_change_at_fixed_mass']:+.4f}% and 2.2-µm extinction by a factor {screen['bare_grain_kappa_ext_2p2_ratio_at_fixed_mass']:.5f}. Its fixed-NIR-extinction S ratio is {screen['bare_grain_S_9p7_per_nir_ratio_at_fixed_nir_extinction']:.5f}. These are different comparisons; the S ratio is not a fixed-mass flux prediction.", "",
              f"Pilot ice mass fraction: {screen['pilot_ice_mass_fraction']:.6g}. Bare-screen/icy-pilot mismatch: {screen['bare_screen_vs_icy_pilot_mismatch']}.", "",
              "The frozen `dust_prescription` and `pilot` metadata in `pilot_summary.json` identify the material files, grain geometry and exactly which species changes size. The bare-grain screen is contextual rather than a radiative-transfer prediction for this prescription.", "",
              "## Limits", "", *[f"- {value}" for value in summary["limitations"]], "",
              "`summary.json` retains receipt validation, exclusions and all nuisance-profiled weighting/covariance scenarios. `pilot_summary.json` and `paired_band_responses.csv` contain the raw-flux experiment. Horizontal figure bars are observed band widths, not uncertainty bars.", ""]
    (output / "PILOT_REVIEW.md").write_text("\n".join(lines))
    if make_plot:
        _plot(output, summary, contract)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-plots", "--no-plot", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = analyze_pilot(args.run, args.output, not args.no_plots)
    except (OSError, ValueError, TypeError, KeyError, ImportError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps({key: result[key] for key in ("status", "models_screened", "models_total", "pairs_complete", "pairs_expected")}, indent=2))
    return 0 if result["status"] == "complete_pilot" else 2


if __name__ == "__main__":
    raise SystemExit(main())
