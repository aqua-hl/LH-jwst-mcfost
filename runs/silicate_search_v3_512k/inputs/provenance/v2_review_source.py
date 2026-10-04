#!/usr/bin/env python3
"""Read-only review of the returned silicate composition search v2.

Reads only the compact downloaded results of runs/silicate_search_v2_512k and,
for comparison, the earlier campaigns' band predictions. It recomputes nothing
with MCFOST, removes no model, and changes neither the declared success rule
nor the analyser's tables. It extracts findings the campaign review does not
make explicit:

  F1  the declared v2 criteria applied to every earlier campaign's models;
  F2  raw MIRI and check-band residuals along the matched carbon ladder;
  F3  the Stage 0 mixture opacity ratios behind the 9.7/18-um coupling;
  F4  material, profile and inclination responses, with the v1 density-profile
      and grain-size responses at each column for comparison;
  F5  the ice mass, and its column equivalent, at which the mean raw
      2.95/3.11-um residual crosses zero on every ice ladder;
  F6  an additive projection for the untested shallower-profile/carbon
      combination at a near-IR-matched column, under two stated assumptions.

Residuals are log10(model/observed) of broad-band F_nu means. They are
descriptive, not a likelihood. F6 is arithmetic on measured responses, not a
radiative-transfer prediction.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

import numpy as np

from silicate_structure_design import column_normalization

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "runs/silicate_search_v2_512k/results"
OUT = ROOT / "validation/silicate_search_v2_512k"
EARLIER = {
    "extinction_ice_production_1au_v1": ROOT / "runs/extinction_ice_production_1au_v1_review/results",
    "silicate_size_pilot_v1_512k": ROOT / "runs/silicate_size_pilot_v1_512k/results",
    "silicate_structure_production_v1_512k": ROOT / "runs/silicate_structure_production_v1_512k/results",
}
V1 = EARLIER["silicate_structure_production_v1_512k"]

NIR_CONTINUUM = tuple(f"c{i:02d}" for i in range(1, 10))
MIRI = tuple(f"s{i:02d}" for i in range(1, 7))
CHECKS = ("w007", "w008", "w009")
LEVER_KEYS = ("shape", "s01", "s02", "s03", "s04", "s05", "s06", "h02", "h03", "nir_rms")
# Declared v2 tolerances (docs/SILICATE_SEARCH_V2.md) as implemented: raw ice flux.
SHAPE_TOL, ICE_TOL, MIRI_TOL, NIR_TOL = 0.10, 0.15, 0.15, 0.20
REFERENCE_COLUMN = 2.25e-4  # M_sun, the near-IR-preferred column at p = -1.5

LIMITATIONS = [
    "One photon seed per model and unresolved five-versus-three-node quadrature "
    "warnings on every v2 model; differences of a few hundredths of a dex are "
    "not numerically certified.",
    "Earlier campaigns used other ice prescriptions (the 1-AU rerun used the "
    "generic MCFOST ice file), so F1 compares regions jointly fitted, not "
    "identical model families.",
    "Ice crossings interpolate log-linearly between two or three sampled ice "
    "fractions; flagged rows extrapolate from the nearest segment. Column "
    "equivalents use the analytic smooth-envelope factor, not a certified "
    "sightline column on the cavity-bearing grid.",
    "An ice mass here is the envelope H2O mass of this DHS template, conditional "
    "on the supplied constants, grain sizes, geometry and inclination. It is not "
    "a measured abundance or an H2O column density.",
    "F6 adds responses measured in different models and assumes they do not "
    "interact. Its profile step comes from v1 (no carbon) at each profile's "
    "near-IR-optimal column, interpolated between sampled columns; at p = -1.5 "
    "that optimum lies at the edge of v1's sampled range. The shallower step was "
    "never run with carbon; the weak variant infers its carbon scaling from the "
    "steeper step.",
]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_bands(results: Path):
    """Model and observed band fluxes and effective wavelengths."""
    flux, observed, wavelength = defaultdict(dict), {}, {}
    for row in read_csv(results / "band_predictions.csv"):
        flux[row["model_id"]][row["band_id"]] = float(row["model_flux_jy"])
        observed[row["band_id"]] = float(row["observed_flux_jy"])
        wavelength[row["band_id"]] = float(row["effective_wavelength_um"])
    complete = {m: f for m, f in flux.items() if all(b in f for b in observed)}
    return complete, observed, wavelength


def core_over_continuum(f: dict, wavelength: dict) -> float:
    """9.7-um core over a power law through the 8.4- and 11.0-um bands."""
    t = math.log(wavelength["s02"] / wavelength["s01"]) / math.log(wavelength["s03"] / wavelength["s01"])
    return f["s02"] / math.exp((1 - t) * math.log(f["s01"]) + t * math.log(f["s03"]))


def passes(r: dict) -> dict:
    return {"shape": abs(r["shape"]) <= SHAPE_TOL,
            "ice": max(abs(r["h02"]), abs(r["h03"])) <= ICE_TOL,
            "miri": max(abs(r["s05"]), abs(r["s06"])) <= MIRI_TOL,
            "nir": r["nir_rms"] < NIR_TOL}


def criteria(f: dict, observed: dict, wavelength: dict) -> dict:
    r = {b: math.log10(f[b] / observed[b]) for b in observed}
    r["shape"] = math.log10(core_over_continuum(f, wavelength) / core_over_continuum(observed, wavelength))
    r["nir_rms"] = math.sqrt(st.fmean(r[b] ** 2 for b in NIR_CONTINUUM))
    return r


def all_criteria(results: Path) -> dict:
    flux, observed, wavelength = load_bands(results)
    return {m: criteria(f, observed, wavelength) for m, f in flux.items()}


def cross_campaign(crit_v2: dict) -> list[dict]:
    rows = []
    for name, crit in {"silicate_search_v2_512k": crit_v2,
                       **{k: all_criteria(v) for k, v in EARLIER.items()}}.items():
        flags = [passes(r) for r in crit.values()]
        best = min(crit.values(), key=lambda r: abs(r["shape"]))
        rows.append({
            "campaign": name, "models": len(flags),
            "shape_9p7": sum(p["shape"] for p in flags),
            "shape_and_nir": sum(p["shape"] and p["nir"] for p in flags),
            "shape_nir_raw_ice": sum(p["shape"] and p["nir"] and p["ice"] for p in flags),
            "three_of_four": sum(sum(p.values()) >= 3 for p in flags),
            "joint": sum(all(p.values()) for p in flags),
            "smallest_abs_shape_dex": round(best["shape"], 4),
            "nir_rms_of_that_model_dex": round(best["nir_rms"], 4),
            "lowest_nir_rms_with_shape_dex": round(min(
                (r["nir_rms"] for r, p in zip(crit.values(), flags) if p["shape"]),
                default=float("nan")), 4),
        })
    return rows


def catalogue() -> dict[int, dict]:
    ids = {int(r["index"]): r["model_id"] for r in read_csv(RESULTS / "model_success.csv")}
    out = {}
    for row in read_csv(RESULTS / "opacity_matched_catalogue.csv"):
        i = int(row["index"])
        out[i] = {"index": i, "model_id": ids[i], "material": row["material_id"],
                  "carbon": round(float(row["carbon_mass_fraction"]), 6),
                  "ice": round(float(row["ice_mass_fraction"]), 6),
                  "p": float(row["density_exponent"]), "inc": float(row["inclination_deg"]),
                  **json.loads(row["opacity_matching"])}
    return out


def find(cat: dict, **attrs) -> dict:
    hits = [m for m in cat.values() if all(m[k] == v for k, v in attrs.items())]
    if len(hits) != 1:
        raise SystemExit(f"expected one model for {attrs}, found {len(hits)}")
    return hits[0]


def carbon_ladder(cat: dict, crit: dict) -> list[dict]:
    checks = defaultdict(dict)
    for row in read_csv(RESULTS / "consistency_checks.csv"):
        checks[row["model_id"]][row["id"]] = math.log10(
            float(row["model_flux_jy"]) / float(row["observed_flux_jy"]))
    rows = []
    for m in sorted((m for m in cat.values() if m["material"] == "draine" and m["p"] == -1.5
                     and m["ice"] == 0.08 and m["inc"] == 60.0), key=lambda m: m["carbon"]):
        r = crit[m["model_id"]]
        row = {"index": m["index"], "carbon_percent": round(100 * m["carbon"]),
               "dust_mass_msun": m["actual_mass_msun"], "shape_9p7_dex": round(r["shape"], 4)}
        row.update({f"{b}_dex": round(r[b], 4) for b in MIRI})
        row.update({f"{b}_dex": round(checks[m["model_id"]][b], 4) for b in CHECKS})
        rows.append(row)
    return rows


def opacity_ratios(cat: dict) -> list[dict]:
    seen, rows = set(), []
    for m in sorted(cat.values(), key=lambda m: (m["material"] != "draine", m["material"],
                                                  m["carbon"], m["ice"])):
        if (m["material"], m["carbon"], m["ice"]) in seen:
            continue
        seen.add((m["material"], m["carbon"], m["ice"]))
        rows.append({
            "material": m["material"], "carbon_percent": round(100 * m["carbon"]),
            "ice_percent": round(100 * m["ice"]), "S": round(m["S"], 4),
            "S_over_reference": round(m["S_over_reference"], 4),
            "absorption_18_over_9p7": round(m["absorption_18_over_9p7"], 4),
            "absorption_18_per_extinction_2p2": round(
                m["mixture_absorption_18_cm2_g"] / m["mixture_extinction_cm2_g"], 4),
            "composition_mass_factor": round(m["composition_mass_factor"], 4)})
    return rows


def difference(a: dict, b: dict) -> dict:
    return {f"d_{k}": round(b[k] - a[k], 4) for k in LEVER_KEYS}


def declared_contrasts(source: str, results: Path, crit: dict, kinds: set[str]) -> list[dict]:
    pairs = {}
    for row in read_csv(results / "matched_feature_shapes.csv"):
        if row["kind"] in kinds:
            pairs.setdefault(row["contrast_id"], row)
    rows = []
    for cid, row in pairs.items():
        if row["kind"] == "material":
            change = f"{row['reference_material_id']}->{row['changed_material_id']}"
        elif row["kind"] == "size":
            change = f"amax {row['reference_amax_um']}->{row['changed_amax_um']} um"
        else:
            change = f"p {row['reference_density_exponent']}->{row['changed_density_exponent']}"
        rows.append({
            "source": source, "contrast_id": cid, "kind": row["kind"], "change": change,
            "reference_model_id": row["reference_model_id"], "changed_model_id": row["changed_model_id"],
            "carbon_percent": round(100 * float(row.get("changed_carbon_mass_fraction") or 0)),
            "inclination_deg": float(row["inclination_deg"]),
            "reference_column_mass_msun": float(row.get("reference_column_mass_msun")
                                                or row.get("reference_dust_mass_msun")),
            **difference(crit[row["reference_model_id"]], crit[row["changed_model_id"]])})
    return rows


def inclination_contrasts(cat: dict, crit: dict) -> list[dict]:
    groups = defaultdict(dict)
    for m in cat.values():
        groups[(m["material"], m["carbon"], m["ice"], m["p"])][m["inc"]] = m
    rows = []
    for (material, carbon, ice, p), g in sorted(groups.items()):
        if 60.0 in g and 70.0 in g:
            a, b = g[60.0], g[70.0]
            rows.append({
                "source": "v2", "contrast_id": f"inclination_{a['index']:03d}_{b['index']:03d}",
                "kind": "inclination", "change": "60->70 deg",
                "reference_model_id": a["model_id"], "changed_model_id": b["model_id"],
                "carbon_percent": round(100 * carbon), "inclination_deg": 70.0,
                "reference_column_mass_msun": a["reference_mass_msun"],
                **difference(crit[a["model_id"]], crit[b["model_id"]])})
    return rows


def ice_crossings(cat: dict, crit: dict) -> list[dict]:
    groups = defaultdict(list)
    for m in cat.values():
        r = crit[m["model_id"]]
        groups[(m["material"], m["carbon"], m["p"], m["inc"])].append(
            (m["actual_mass_msun"] * m["ice"], 0.5 * (r["h02"] + r["h03"]), m))
    rows = []
    for (material, carbon, p, inc), g in sorted(groups.items()):
        g.sort(key=lambda x: x[0])
        segment, kind = None, "interpolated"
        for lo, hi in zip(g, g[1:]):
            if lo[1] * hi[1] <= 0:
                segment = (lo, hi)
                break
        if segment is None:  # every residual has one sign: extend the end nearest zero
            segment, kind = ((g[-2], g[-1]) if g[-1][1] > 0 else (g[0], g[1])), "extrapolated"
        (m1, r1, a), (m2, r2, b) = segment
        t = r1 / (r1 - r2)
        slope = (r2 - r1) / math.log10(m2 / m1)
        crossing = 10 ** (math.log10(m1) + t * math.log10(m2 / m1))
        rows.append({
            "material": material, "carbon_percent": round(100 * carbon), "density_exponent": p,
            "inclination_deg": inc,
            "ladder": " ".join(f"{100 * m['ice']:.0f}%:{r:+.3f}" for _, r, m in g),
            "crossing": kind,
            "crossing_ice_fraction": round(a["ice"] + t * (b["ice"] - a["ice"]), 4),
            "crossing_ice_mass_msun": crossing,
            "segment_slope_dex_per_dex": round(slope, 4),
            "column_equivalent_ice_mass_msun": crossing / a["profile_mass_factor"],
            "total_dust_mass_at_8pct_ice_msun": next(
                (m["actual_mass_msun"] for _, _, m in g if m["ice"] == 0.08), float("nan"))})
    return rows


def v1_column_ladders(v1_crit: dict) -> dict:
    """v1 Draine/0.4-um models by (exponent, inclination), as (reference column, criteria)."""
    ladders, seen = defaultdict(list), set()
    for row in read_csv(V1 / "ranking.csv"):
        mid = row["model_id"]
        if mid in seen or mid not in v1_crit:
            continue
        seen.add(mid)
        p = json.loads(row["parameters"])
        if p["envelope_amax_um"] != 0.4 or not p["envelope_silicate_file"].startswith("Draine"):
            continue
        exponent = p["envelope_density_exponent"]
        factor = column_normalization(exponent, REFERENCE_COLUMN)["mass_ratio"]
        ladders[(exponent, p["inclination_deg"])].append((p["envelope_dust_mass_msun"] / factor, v1_crit[mid]))
    return ladders


def nir_optimum(ladder: list) -> dict:
    """Interpolate a column ladder, log-linearly per band, to its near-IR RMS minimum."""
    ladder = sorted(ladder, key=lambda t: t[0])
    x = np.log10([column for column, _ in ladder])
    fine = np.linspace(x[0], x[-1], 2001)
    values = {k: np.interp(fine, x, [r[k] for _, r in ladder]) for k in ladder[0][1] if k != "nir_rms"}
    rms = np.sqrt(np.mean([values[b] ** 2 for b in NIR_CONTINUUM], axis=0))
    i = int(np.argmin(rms))
    return {"reference_column_msun": float(10 ** fine[i]), "at_sampled_edge": i in (0, len(fine) - 1),
            "nir_rms": float(rms[i]), **{k: float(v[i]) for k, v in values.items()}}


def projection(cat: dict, crit: dict, v1_crit: dict, v1_levers: list[dict], v2_profile: list[dict],
               crossings: list[dict]) -> tuple[list[dict], dict]:
    """Additive arithmetic: 20% carbon, 8% ice, 70 deg, p=-1.25 at a near-IR-matched column."""
    c = lambda **kw: crit[find(cat, **kw)["model_id"]]
    base60 = c(material="draine", carbon=0.2, ice=0.08, p=-1.5, inc=60.0)
    base70 = c(material="draine", carbon=0.2, ice=0.08, p=-1.5, inc=70.0)
    pyro60 = c(material="pyroxene_mg50", carbon=0.2, ice=0.08, p=-1.5, inc=60.0)
    ladders = v1_column_ladders(v1_crit)
    a, b = nir_optimum(ladders[(-1.5, 70.0)]), nir_optimum(ladders[(-1.25, 70.0)])
    keys = ("shape", "s05", "s06") + NIR_CONTINUUM
    step = {q: b[q] - a[q] for q in keys}
    delta = lambda rows, table, q: st.median(table[r["changed_model_id"]][q] - table[r["reference_model_id"]][q]
                                             for r in rows)
    matched = [r for r in v1_levers if r["reference_column_mass_msun"] == REFERENCE_COLUMN]
    steep = [r for r in matched if r["change"].endswith("-1.75")]
    shallow = [r for r in matched if r["change"].endswith("-1.25")]
    # Carbon weakens the MIRI profile response; measure that on the steeper step,
    # the only profile step run both with and without carbon.
    ratio = {q: delta(v2_profile, crit, q) / delta(steep, v1_crit, q) for q in ("shape", "s05", "s06", "c01", "c02")}
    assumptions = {"strong_v1_step": step,
                   "weak_carbon_scaled": {q: step[q] * ratio.get(q, 1.0) if q in ("shape", "s05", "s06")
                                          else step[q] for q in keys}}
    material = {q: pyro60[q] - base60[q] for q in keys}
    rows = []
    for name, profile in [("measured_p-1.5", {q: 0.0 for q in keys}), *assumptions.items()]:
        for silicate, offset in (("draine", {q: 0.0 for q in keys}), ("pyroxene_mg50", material)):
            r = {q: base70[q] + offset[q] + profile[q] for q in keys}
            r["nir_rms"] = math.sqrt(st.fmean(r[band] ** 2 for band in NIR_CONTINUUM))
            flags = {"shape": abs(r["shape"]) <= SHAPE_TOL,
                     "miri": max(abs(r["s05"]), abs(r["s06"])) <= MIRI_TOL, "nir": r["nir_rms"] < NIR_TOL}
            rows.append({"profile_response": name, "silicate": silicate,
                         **{q: round(r[q], 3) for q in ("shape", "s05", "s06", "nir_rms")},
                         **{f"pass_{f}": v for f, v in flags.items()},
                         "criteria_met_of_three": sum(flags.values())})
    factor = b["reference_column_msun"] / REFERENCE_COLUMN
    ice_change = 0.5 * (b["h02"] + b["h03"] - a["h02"] - a["h03"])
    slope = st.median(r["segment_slope_dex_per_dex"] for r in crossings if r["crossing"] == "interpolated")
    meta = {"optimum_p-1.5": a, "optimum_p-1.25": b, "column_factor_p-1.25": factor,
            "step": step, "carbon_ratio": ratio, "ice_change_fixed_fraction": ice_change,
            "ice_slope": slope, "ice_offset_column_held": slope * math.log10(1 / factor),
            "analytic_shallow": {q: delta(shallow, v1_crit, q) for q in ("c01", "c02")}}
    return rows, meta


def markdown_table(rows: list[dict], columns: list[tuple[str, str, str]]) -> str:
    head = "| " + " | ".join(label for _, label, _ in columns) + " |"
    rule = "|" + "|".join("---:" if fmt else "---" for _, _, fmt in columns) + "|"
    body = ["| " + " | ".join(format(r[key], fmt) if fmt else str(r[key]) for key, _, fmt in columns) + " |"
            for r in rows]
    return "\n".join([head, rule, *body])


def lever_summary(rows: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for r in rows:
        if r["source"] == "v2":
            label = {"material": f"v2 {r['change']} (20% C)", "profile": "v2 p -1.5 -> -1.75 (20% C)",
                     "inclination": "v2 60 -> 70 deg (all carbon levels)"}[r["kind"]]
        elif r["kind"] == "size":
            label = f"v1 {r['change']}, all columns (0% C)"
        else:
            label = f"v1 {r['change']}, column {r['reference_column_mass_msun']:.2e} (0% C)"
        groups[label].append(r)
    return [{"lever": label, "n": len(g),
             **{k: st.median(r[f"d_{k}"] for r in g) for k in ("shape", "s04", "s05", "s06", "h02", "nir_rms")}}
            for label, g in groups.items()]


def findings(cross, ladder, ratios, levers, crossings, proj, meta) -> str:
    v2 = cross[0]
    earlier = cross[1:]
    n_earlier = sum(r["models"] for r in earlier)
    joint_earlier = sum(r["shape_and_nir"] for r in earlier)
    lowest = min(r["lowest_nir_rms_with_shape_dex"] for r in earlier if r["shape_9p7"])
    ten = next(r for r in ladder if r["carbon_percent"] == 10)
    inner = [ten[f"{b}_dex"] for b in ("s01", "s02", "s03")]
    outer = [ten[f"{b}_dex"] for b in ("s04", "w008", "s05", "s06")]
    ratio = lambda m: next(r["absorption_18_over_9p7"] for r in ratios if r["material"] == m
                           and r["carbon_percent"] == 20 and r["ice_percent"] == 8)
    raised = [ratio(m) / ratio("draine") - 1 for m in ("pyroxene_mg50", "olivine_mg50")]
    drop = [r[f"d_{b}"] for r in levers if r["kind"] == "material" for b in ("s05", "s06")]
    size = [r for r in levers if r["kind"] == "size"]
    size_long = [r[f"d_{b}"] for r in size for b in ("s05", "s06")]
    bracketed = [r["column_equivalent_ice_mass_msun"] for r in crossings if r["crossing"] == "interpolated"]
    extended = [r["column_equivalent_ice_mass_msun"] for r in crossings if r["crossing"] != "interpolated"]
    fractions = [r["crossing_ice_fraction"] for r in crossings]
    draine = [r for r in ratios if r["material"] == "draine"]
    dust = [r["total_dust_mass_at_8pct_ice_msun"] for r in crossings if r["density_exponent"] == -1.5]
    parts = [
        "# Silicate search v2: review findings",
        "",
        "Generated by `scripts/review_silicate_search_v2.py` from the returned compact results. "
        "Read-only: the declared success rule, the analyser's tables and `CAMPAIGN_REVIEW.md` are "
        "unchanged. Residuals are log10(model/observed) of broad-band means and are descriptive.",
        "",
        "## F1 The first joint fit of the near-IR continuum, 3-µm ice and 9.7-µm silicate",
        "",
        f"In v2, {v2['shape_and_nir']} of {v2['models']} models meet the 9.7-µm shape and near-IR "
        f"criteria together, and {v2['shape_nir_raw_ice']} also meet the raw ice criterion. "
        f"Across the {n_earlier} models of the three earlier campaigns, {joint_earlier} meet the "
        f"first two together: all {sum(r['shape_9p7'] for r in earlier)} earlier models that "
        f"reproduced the 9.7-µm shape had a near-IR RMS of at least {lowest:.2f} dex. No model in "
        "any campaign meets all four criteria.",
        "",
        markdown_table(cross, [("campaign", "Campaign", ""), ("models", "Models", "d"),
                               ("shape_9p7", "9.7 shape", "d"), ("shape_and_nir", "+ NIR", "d"),
                               ("shape_nir_raw_ice", "+ raw ice", "d"), ("three_of_four", "≥3 of 4", "d"),
                               ("joint", "All 4", "d"), ("lowest_nir_rms_with_shape_dex",
                                                        "Lowest NIR RMS with 9.7 shape", ".3f")]),
        "",
        "## F2 The remaining tension is a 12–20 µm excess",
        "",
        "Matched carbon ladder (Draine, 8% ice, p = −1.5, 60°). The 13.8- and 27.5-µm points are the "
        "unweighted check bands.",
        "",
        markdown_table(ladder, [("carbon_percent", "C %", "d"), ("shape_9p7_dex", "9.7 shape", "+.3f"),
                                ("s01_dex", "8.4", "+.3f"), ("s02_dex", "9.7", "+.3f"),
                                ("s03_dex", "11.0", "+.3f"), ("s04_dex", "12.7", "+.3f"),
                                ("w008_dex", "13.8", "+.3f"), ("s05_dex", "18.45", "+.3f"),
                                ("s06_dex", "19.8", "+.3f"), ("w009_dex", "27.5", "+.3f")]),
        "",
        f"With 10% carbon, the 8.4–11.0-µm bands lie within {max(map(abs, inner)):.2f} dex of the "
        f"data while 12.7–19.8 µm is {min(outer):+.2f} to {max(outer):+.2f} dex bright. From 0% to "
        "30% carbon the 27.5-µm check runs "
        + ", ".join(f"{r['w009_dex']:+.2f}" for r in ladder) + ": it does not brighten.",
        "",
        "## F3 Carbon lowers 9.7- and 18-µm absorption together",
        "",
        f"At matched 2.2-µm extinction, the Draine mixtures' 18/9.7 absorption ratio stays within "
        f"{min(r['absorption_18_over_9p7'] for r in draine):.3f}–"
        f"{max(r['absorption_18_over_9p7'] for r in draine):.3f} from 0% to 30% carbon, while S falls "
        "fourfold. The carbon fraction therefore cannot shallow the 9.7-µm feature without "
        "reducing 18-µm absorption per unit near-IR extinction in nearly the same proportion.",
        "",
        markdown_table(ratios, [("material", "Silicate", ""), ("carbon_percent", "C %", "d"),
                                ("ice_percent", "Ice %", "d"), ("S", "S", ".3f"),
                                ("absorption_18_over_9p7", "κ18/κ9.7", ".3f"),
                                ("absorption_18_per_extinction_2p2", "κ18/κext(2.2)", ".3f")]),
        "",
        "## F4 Lever responses (median changed − reference, dex)",
        "",
        markdown_table(lever_summary(levers), [("lever", "Lever", ""), ("n", "n", "d"),
                                               ("shape", "9.7 shape", "+.3f"), ("s04", "12.7", "+.3f"),
                                               ("s05", "18.45", "+.3f"), ("s06", "19.8", "+.3f"),
                                               ("h02", "2.95", "+.3f"), ("nir_rms", "NIR RMS", "+.3f")]),
        "",
        f"Raising κ18/κ9.7 by {100 * min(raised):.0f}–{100 * max(raised):.0f}% (pyroxene, olivine) "
        f"changes 18.45/19.8 µm by only {min(drop):+.3f} to {max(drop):+.3f} dex. The "
        "steeper profile moves 18–20 µm the wrong way, and so do larger silicate grains: in v1, "
        f"amax 0.4 → 1.0 µm changed the 9.7-µm shape by {min(r['d_shape'] for r in size):+.3f} to "
        f"{max(r['d_shape'] for r in size):+.3f} dex and 18.45/19.8 µm by {min(size_long):+.3f} to "
        f"{max(size_long):+.3f} dex. In v1 the profile response grew with "
        "column; at 20% carbon the steeper-profile response resembles v1's lowest column, "
        "consistent with carbon cutting the silicate optical depth per unit near-IR extinction.",
        "",
        "## F5 The 3-µm band fixes the ice column, not the ice fraction",
        "",
        f"Across {len(crossings)} ice ladders, the mean raw 2.95/3.11-µm residual crosses zero at a "
        f"column-equivalent envelope ice mass of {min(bracketed):.2e}–{max(bracketed):.2e} M☉ where "
        f"the crossing is bracketed ({len(bracketed)} ladders) and {min(extended):.2e}–"
        f"{max(extended):.2e} M☉ where it is extrapolated ({len(extended)}). Over the same ladders the "
        f"ice fraction at the crossing runs {min(fractions):.1%}–{max(fractions):.1%}, and the total "
        f"dust mass at p = −1.5 spans {min(dust):.2e}–{max(dust):.2e} M☉. The fitting ice fraction "
        "rises with carbon because carbon lowers the total dust mass at matched near-IR extinction.",
        "",
        markdown_table(crossings, [("material", "Silicate", ""), ("carbon_percent", "C %", "d"),
                                   ("density_exponent", "p", "+.2f"), ("inclination_deg", "Inc.", ".0f"),
                                   ("ladder", "Ice %: residual", ""), ("crossing", "Crossing", ""),
                                   ("crossing_ice_fraction", "Ice fraction", ".3f"),
                                   ("column_equivalent_ice_mass_msun", "Column-equiv. ice M☉", ".2e")]),
        "",
        "## F6 Untested: shallower profile with carbon at a re-matched column (arithmetic projection)",
        "",
        "At the matched central column, v1's step to p = −1.25 dimmed the 1.1- and 1.25-µm bands by "
        f"{meta['analytic_shallow']['c01']:+.2f}/{meta['analytic_shallow']['c02']:+.2f} dex, so a shallower "
        "profile needs its column re-matched to the near-IR. Interpolating v1's column ladders "
        "(Draine, no carbon, 70°), the near-IR RMS is lowest at "
        f"{meta['optimum_p-1.5']['reference_column_msun']:.2e} M☉ for p = −1.5"
        + (" (edge of the sampled range)" if meta['optimum_p-1.5']['at_sampled_edge'] else "")
        + f" and at {meta['optimum_p-1.25']['reference_column_msun']:.2e} M☉ for p = −1.25, "
        f"{meta['column_factor_p-1.25']:.2f} of the analytic match, with RMS "
        f"{meta['optimum_p-1.5']['nir_rms']:.3f} and {meta['optimum_p-1.25']['nir_rms']:.3f}. Between the two "
        f"optima the 9.7-µm shape changes by {meta['step']['shape']:+.3f} dex and 18.45/19.8 µm by "
        f"{meta['step']['s05']:+.3f}/{meta['step']['s06']:+.3f} dex. At a fixed ice fraction the mean "
        f"2.95/3.11-µm residual rises by {meta['ice_change_fixed_fraction']:+.3f} dex; holding the radial ice "
        f"column instead offsets {meta['ice_offset_column_held']:+.3f} dex of that, at F5's median ladder "
        f"slope of {meta['ice_slope']:.2f} dex per dex.",
        "",
        "The projection adds that step to the measured Draine 20% C, 8% ice, 70° model; pyroxene adds "
        "the measured Draine→pyroxene change at 60°. The strong variant keeps the v1 step. The weak "
        "variant scales its 9.7-µm and 18–20-µm parts by the carbon weakening measured on the steeper "
        "step (v2 at 20% carbon against v1 at the matched column: "
        f"{meta['carbon_ratio']['shape']:.2f}, {meta['carbon_ratio']['s05']:.2f} and "
        f"{meta['carbon_ratio']['s06']:.2f}). The near-IR part is not scaled: on the steeper step the "
        f"1.1/1.25-µm responses with carbon were {meta['carbon_ratio']['c01']:.2f}/"
        f"{meta['carbon_ratio']['c02']:.2f} of v1's. Ice is not projected.",
        "",
        markdown_table(proj, [("profile_response", "Profile response", ""), ("silicate", "Silicate", ""),
                              ("shape", "9.7 shape", "+.3f"), ("s05", "18.45", "+.3f"), ("s06", "19.8", "+.3f"),
                              ("nir_rms", "NIR RMS", ".3f"),
                              ("criteria_met_of_three", "Met (9.7, 18–20, NIR)", "d")]),
        "",
        "These rows motivate a test; they do not predict its outcome.",
        "",
        "## Limits",
        "",
        *[f"- {item}" for item in LIMITATIONS],
        "",
    ]
    return "\n".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    cat = catalogue()
    crit = all_criteria(RESULTS)
    cross = cross_campaign(crit)
    ladder = carbon_ladder(cat, crit)
    ratios = opacity_ratios(cat)
    v1_crit = all_criteria(V1)
    v1_levers = declared_contrasts("v1", V1, v1_crit, {"density_profile"})
    v1_size = declared_contrasts("v1", V1, v1_crit, {"size"})
    v2_declared = declared_contrasts("v2", RESULTS, crit, {"material", "profile"})
    levers = v2_declared + inclination_contrasts(cat, crit) + v1_levers + v1_size
    crossings = ice_crossings(cat, crit)
    proj, meta = projection(cat, crit, v1_crit, v1_levers,
                            [r for r in v2_declared if r["kind"] == "profile"], crossings)

    outputs = {"cross_campaign_criteria.csv": cross, "carbon_ladder_miri.csv": ladder,
               "opacity_ratios.csv": ratios, "lever_responses.csv": levers,
               "ice_column_crossings.csv": crossings, "crossed_profile_projection.csv": proj}
    for name, rows in outputs.items():
        write_csv(args.output / name, rows)
    (args.output / "FINDINGS.md").write_text(findings(cross, ladder, ratios, levers, crossings, proj, meta))
    try:
        shown = args.output.resolve().relative_to(ROOT)
    except ValueError:
        shown = args.output
    print(f"wrote {len(outputs) + 1} files to {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
