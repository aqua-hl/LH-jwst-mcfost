#!/usr/bin/env python3
"""Planning numbers for docs/SILICATE_SEARCH_V3.md. Builds no run package.

Every composition in v2 is a mass-weighted sum of separate DHS species, so the
pure-species opacities at 2.2, 9.7 and 18 um follow exactly from the frozen
v2 catalogue by least squares. This script recovers them and derives:

  - the column factor at each shallower profile, from the column at which v1's
    near-IR RMS is lowest (Draine, no carbon, 70 deg), relative to the analytic
    central-column match;
  - each silicate's target ice column, from the v2 review's bracketed 3-um
    crossings, and the ice fraction that holds it at every column factor;
  - the opacity ratios and matched envelope masses of every v3 model.

The production builder must recompute masses and ice fractions from its own
Stage 0 DHS opacities. These planning values never become run inputs.

    python -B scripts/silicate_search_v3_plan.py [--json PATH]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
from pathlib import Path

import numpy as np

from review_silicate_search_v2 import REFERENCE_COLUMN, V1, all_criteria, nir_optimum, v1_column_ladders
from silicate_structure_design import column_normalization

ROOT = Path(__file__).resolve().parents[1]
V2_CATALOGUE = ROOT / "runs/silicate_search_v2_512k/results/opacity_matched_catalogue.csv"
V2_ICE_CROSSINGS = ROOT / "validation/silicate_search_v2_512k/ice_column_crossings.csv"

REFERENCE_MASS = REFERENCE_COLUMN  # 96% Draine / 4% ice at p = -1.5, as in v2
MATERIALS = ("draine", "pyroxene_mg50")
CARBON = (0.15, 0.20, 0.25)
INCLINATION_DEG = 70.0
P_MINUS_1_BRACKET_DEX = 0.05  # half-width of the p = -1.0 column bracket
FIELDS = {"ext_2p2": "mixture_extinction_cm2_g", "abs_9p7": "mixture_absorption_9p7_cm2_g",
          "abs_18": "mixture_absorption_18_cm2_g"}


def species_opacities() -> tuple[dict, float]:
    """Pure-species DHS opacities recovered from the v2 mixtures."""
    rows = []
    with V2_CATALOGUE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            o = json.loads(row["opacity_matching"])
            rows.append((row["material_id"], float(row["carbon_mass_fraction"]),
                         float(row["ice_mass_fraction"]), {k: o[v] for k, v in FIELDS.items()}))
    out, worst = {s: {} for s in ("draine", "pyroxene_mg50", "carbon", "ice")}, 0.0
    for field in FIELDS:
        draine = [(1 - c - i, c, i, k[field]) for m, c, i, k in rows if m == "draine"]
        design = np.array([r[:3] for r in draine])
        target = np.array([r[3] for r in draine])
        solution = np.linalg.lstsq(design, target, rcond=None)[0]
        worst = max(worst, float(np.max(np.abs(design @ solution / target - 1))))
        out["draine"][field], out["carbon"][field], out["ice"][field] = map(float, solution)
        pyroxene = [(k[field] - c * solution[1] - i * solution[2]) / (1 - c - i)
                    for m, c, i, k in rows if m == "pyroxene_mg50"]
        worst = max(worst, (max(pyroxene) - min(pyroxene)) / st.fmean(pyroxene))
        out["pyroxene_mg50"][field] = st.fmean(pyroxene)
    return out, worst


def ice_targets() -> dict[str, tuple[float, int]]:
    """Per-silicate median column-equivalent ice mass of the bracketed v2 crossings."""
    with V2_ICE_CROSSINGS.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["crossing"] == "interpolated"]
    out = {}
    for material in MATERIALS:
        values = [float(r["column_equivalent_ice_mass_msun"]) for r in rows if r["material"] == material]
        out[material] = (st.median(values), len(values))
    return out


def configurations() -> tuple[list[tuple[float, float]], dict]:
    """(exponent, column factor) pairs; factors relative to the analytic central-column match."""
    ladders = v1_column_ladders(all_criteria(V1))
    optimum = nir_optimum(ladders[(-1.25, INCLINATION_DEG)])
    measured = optimum["reference_column_msun"] / REFERENCE_COLUMN
    shallow = round(measured, 2)
    # Extend the measured log-factor linearly to p = -1.0 and bracket it.
    centre = shallow ** 2
    pair = [round(centre * 10 ** (sign * P_MINUS_1_BRACKET_DEX), 2) for sign in (1, -1)]
    configs = [(-1.5, 1.0), (-1.25, shallow), (-1.0, pair[0]), (-1.0, pair[1])]
    return configs, {"v1_optimum_column_msun": optimum["reference_column_msun"],
                     "v1_optimum_nir_rms": optimum["nir_rms"], "measured_factor": measured,
                     "p-1.0_centre": centre}


def plan() -> dict:
    kappa, worst = species_opacities()
    targets = ice_targets()
    configs, derivation = configurations()
    reference = {f: 0.96 * kappa["draine"][f] + 0.04 * kappa["ice"][f] for f in FIELDS}
    s_reference = reference["abs_9p7"] / reference["ext_2p2"]
    rows = []
    for exponent, factor in configs:
        profile = column_normalization(exponent, REFERENCE_MASS)["mass_ratio"]
        for material in MATERIALS:
            target = targets[material][0]
            for carbon in CARBON:
                # Radial ice column f*M_ref*c*k_ref/k_mix(f) = target, with k_mix linear
                # in the ice fraction f: closed form, no iteration.
                base = (1 - carbon) * kappa[material]["ext_2p2"] + carbon * kappa["carbon"]["ext_2p2"]
                scaled = REFERENCE_MASS * factor
                ice = target * base / (scaled * reference["ext_2p2"]
                                       - target * (kappa["ice"]["ext_2p2"] - kappa[material]["ext_2p2"]))
                fractions = {material: 1 - carbon - ice, "carbon": carbon, "ice": ice}
                mix = {f: sum(x * kappa[s][f] for s, x in fractions.items()) for f in FIELDS}
                composition = reference["ext_2p2"] / mix["ext_2p2"]
                rows.append({
                    "exponent": exponent, "column_factor": factor, "material": material,
                    "carbon": carbon, "ice": ice, "silicate": 1 - carbon - ice,
                    "S_over_reference": mix["abs_9p7"] / mix["ext_2p2"] / s_reference,
                    "absorption_18_over_9p7": mix["abs_18"] / mix["abs_9p7"],
                    "radial_ice_column_msun": ice * scaled * composition,
                    "dust_mass_msun": scaled * composition * profile})
    replicas = [(material, 0.20, configs[1]) for material in MATERIALS]
    return {"species_opacities_cm2_g": kappa, "max_relative_fit_residual": worst,
            "reference_mass_msun": REFERENCE_MASS, "reference_S": s_reference,
            "ice_column_targets_msun": {m: t for m, (t, _) in targets.items()},
            "ice_column_target_crossings": {m: n for m, (_, n) in targets.items()},
            "configurations": configs, "column_factor_derivation": derivation,
            "profile_mass_factors": {str(p): column_normalization(p, REFERENCE_MASS)["mass_ratio"]
                                     for p in sorted({p for p, _ in configs})},
            "inclination_deg": INCLINATION_DEG, "seed_replicas": replicas,
            "physical_models": len(rows), "runs": len(rows) + len(replicas), "rows": rows,
            "note": "Planning values only; the builder recomputes them from its own Stage 0."}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()
    p = plan()
    print(f"Recovered DHS species opacities (max relative residual {p['max_relative_fit_residual']:.1e}):")
    for s, v in p["species_opacities_cm2_g"].items():
        ratio = f"{v['abs_18'] / v['abs_9p7']:.3f}" if v["abs_9p7"] > 1.0 else "--"
        print(f"  {s:>14}: ext(2.2) {v['ext_2p2']:9.1f}  abs(9.7) {v['abs_9p7']:8.1f}  "
              f"abs(18) {v['abs_18']:8.1f}  abs18/abs9.7 {ratio}")
    for m, t in p["ice_column_targets_msun"].items():
        print(f"Ice column target, {m}: {t:.3e} Msun "
              f"(median of {p['ice_column_target_crossings'][m]} bracketed v2 crossings)")
    d = p["column_factor_derivation"]
    print(f"v1 near-IR optimum at p=-1.25, {INCLINATION_DEG:.0f} deg: {d['v1_optimum_column_msun']:.3e} Msun "
          f"(factor {d['measured_factor']:.3f}, RMS {d['v1_optimum_nir_rms']:.3f}); "
          f"p=-1.0 centre {d['p-1.0_centre']:.3f}")
    print("Configurations (p, column factor): " + ", ".join(f"({e}, {c:.2f})" for e, c in p["configurations"]))
    print("Profile mass factors: " + ", ".join(f"p={e}: {f:.11g}" for e, f in p["profile_mass_factors"].items()))
    print("\n| p | Column factor | Silicate | Carbon | Ice | Silicate frac | S/S_ref | k18/k9.7 | Dust mass (Msun) |")
    print("|---:|---:|---|---:|---:|---:|---:|---:|---:|")
    for r in p["rows"]:
        print(f"| {r['exponent']} | {r['column_factor']:.2f} | {r['material']} | {100 * r['carbon']:.0f}% | "
              f"{100 * r['ice']:.2f}% | {100 * r['silicate']:.2f}% | {r['S_over_reference']:.3f} | "
              f"{r['absorption_18_over_9p7']:.3f} | {r['dust_mass_msun']:.3e} |")
    print(f"\n{p['physical_models']} physical models + {len(p['seed_replicas'])} seed replicas = {p['runs']} runs")
    if args.json:
        args.json.write_text(json.dumps(p, indent=2) + "\n")
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
