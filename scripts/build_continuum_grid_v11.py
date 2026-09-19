#!/usr/bin/env python3
"""Build the bounded 96-model continuum v1.1 design; never prepare or submit it.

The physical catalogue and its design provenance use only the standard library.
Matplotlib is needed only for the optional coverage figure.  Once the run exists,
normal generation is refused; --check verifies the already-written catalogue.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import tempfile


ROOT = Path(__file__).resolve().parents[1]
RUN_NAME = "continuum_production_v1.1"
CONFIG_PATH = ROOT / "config/grid.continuum-production-v1.1.json"
OUT = ROOT / "validation" / RUN_NAME
DESIGN_VERSION = "matched-inclinations-maximin-latin-v1"
DESIGN_SEED = 20260919
MCFOST_SEED = 41001
N_CANDIDATES = 2048
N_LATIN = 13
INCLINATIONS = [70, 45, 50, 55, 60, 65]

# Explicit SI constants, matching Astropy 7.1.1's nominal solar constants,
# CODATA G/sigma and the Julian year. These constants are frozen design inputs;
# the generator does not depend on the installed Astropy version.
CONSTANTS = {
    "G_m3_kg_s2": 6.6743e-11,
    "solar_mass_kg": 1.988409870698051e30,
    "solar_radius_m": 695700000.0,
    "solar_luminosity_w": 3.828e26,
    "stefan_boltzmann_w_m2_k4": 5.6703744191844314e-8,
    "julian_year_seconds": 31557600.0,
}
FIXED = {
    "envelope_dust_mass_msun": 0.000225,
    "envelope_amax_um": 0.4,
    "envelope_size_exponent": 2.75,
    "cavity_half_opening_deg": 17.5,
    "envelope_ice_volume_fraction": 0.05,
    "distance_pc": 140.0,
    "stellar_temperature_k": 4000.0,
    "stellar_radius_rsun": 2.5,
    "stellar_mass_msun": 0.5,
    "accretion_rate_msun_per_year": 7.3011756855e-8,
}
AXES = [
    ("envelope_dust_mass_msun", 7.5e-5, 4e-4, "log", "Msun"),
    ("envelope_amax_um", 0.1, 3.0, "log", "micrometre"),
    ("envelope_size_exponent", 2.5, 4.0, "linear", "dimensionless"),
    ("cavity_half_opening_deg", 10.0, 30.0, "linear", "degree"),
    ("source_heating_target_lsun", 1.5, 3.0, "linear", "Lsun"),
]


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def digest(value):
    return hashlib.sha256(value).hexdigest()


def photospheric_luminosity():
    radius_m = FIXED["stellar_radius_rsun"] * CONSTANTS["solar_radius_m"]
    return (4 * math.pi * radius_m**2 * CONSTANTS["stefan_boltzmann_w_m2_k4"]
            * FIXED["stellar_temperature_k"]**4 / CONSTANTS["solar_luminosity_w"])


def accretion_luminosity(mdot):
    c = CONSTANTS
    return (c["G_m3_kg_s2"] * FIXED["stellar_mass_msun"] * c["solar_mass_kg"]
            * mdot * c["solar_mass_kg"] / c["julian_year_seconds"]
            / (FIXED["stellar_radius_rsun"] * c["solar_radius_m"])
            / c["solar_luminosity_w"])


def mdot_for_source_heating(target):
    """Source target = 4 pi R^2 sigma T^4 + G M Mdot / R, not observed Lbol."""
    extra_lsun = target - photospheric_luminosity()
    if extra_lsun <= 0:
        raise ValueError("Heating target must exceed the fixed photosphere")
    return extra_lsun / accretion_luminosity(1.0)


def maximin_latin():
    """Select the first best of deterministic random rank permutations.

    Each dimension visits ranks 0..12 exactly once. Ranks are divided by 12,
    so both bounds are sampled. Distance is Euclidean in the unit cube after
    the stated log/linear transformations; no physics weighting is inferred.
    """
    rng = random.Random(DESIGN_SEED)
    best_ranks, best_distance2, best_candidate = None, -1, None
    for candidate in range(N_CANDIDATES):
        columns = []
        for _ in AXES:
            column = list(range(N_LATIN))
            rng.shuffle(column)
            columns.append(column)
        rows = list(zip(*columns))
        distance2 = min(sum((a - b)**2 for a, b in zip(left, right))
                        for i, left in enumerate(rows) for right in rows[i + 1:])
        if distance2 > best_distance2:
            best_ranks, best_distance2, best_candidate = rows, distance2, candidate
    return [list(row) for row in best_ranks], math.sqrt(best_distance2) / (N_LATIN - 1), best_candidate


def scale(rank, low, high, transform):
    # Return boundary tokens exactly, without exp/log roundoff.
    if rank == 0:
        return low
    if rank == N_LATIN - 1:
        return high
    fraction = rank / (N_LATIN - 1)
    if transform == "log":
        return math.exp(math.log(low) + fraction * math.log(high / low))
    return low + fraction * (high - low)


def make_outputs():
    ranks, separation, candidate = maximin_latin()
    nominal_source = photospheric_luminosity() + accretion_luminosity(FIXED["accretion_rate_msun_per_year"])
    families = []
    controls = [
        ("nominal_control", {}),
        ("cavity20_control", {"cavity_half_opening_deg": 20.0}),
        ("mass2p5_q3p25_control", {"envelope_dust_mass_msun": 0.00025,
                                 "envelope_size_exponent": 3.25}),
    ]
    for label, override in controls:
        families.append({"family_index": len(families), "family": label,
                         "kind": "archived_control", "parameters": {**FIXED, **override},
                         "source_heating_target_lsun": 1.9, "latin_ranks": None})
    for row in ranks:
        sampled = {key: scale(rank, low, high, transform)
                   for rank, (key, low, high, transform, _) in zip(row, AXES)}
        target = sampled.pop("source_heating_target_lsun")
        sampled["accretion_rate_msun_per_year"] = mdot_for_source_heating(target)
        families.append({"family_index": len(families), "family": f"latin_{len(families) - 3:02d}",
                         "kind": "space_filling", "parameters": {**FIXED, **sampled},
                         "source_heating_target_lsun": target, "latin_ranks": row})

    models, catalogue, identities = [], [], set()
    for family in families:
        for inclination in INCLINATIONS:
            parameters = {key: float(value) for key, value in family["parameters"].items()}
            parameters["inclination_deg"] = float(inclination)
            index = len(models)
            # Match configuration.expand_models' physical model identity exactly.
            identity = digest(json.dumps(parameters, sort_keys=True, allow_nan=False).encode())
            if identity in identities:
                raise ValueError("Duplicate physical model")
            identities.add(identity)
            models.append(parameters)
            catalogue.append({
                "index": index, "model_id": f"m{index:04d}_{identity[:8]}",
                "family_index": family["family_index"], "family": family["family"],
                "kind": family["kind"], "source_heating_target_lsun": family["source_heating_target_lsun"],
                "photospheric_luminosity_lsun": photospheric_luminosity(),
                "accretion_luminosity_lsun": accretion_luminosity(parameters["accretion_rate_msun_per_year"]),
                **parameters,
            })

    config = {
        "schema_version": 1, "run_name": RUN_NAME, "output_dir": "../runs",
        "smoke_test": False, "max_models": 96,
        "notes": (
            "User-selected 96-model joint continuum search: 16 dust/heating settings each at "
            "inclinations [70,45,50,55,60,65] degrees. Three exact archived structural controls "
            "plus 13 endpoint-inclusive maximin Latin-hypercube settings vary envelope mass, "
            "amax, grain exponent, cavity half-opening, and stellar accretion heating together. "
            "These are exploratory bounds, not inferred or literature confidence limits. "
            "The fixed 4000 K, 2.5 Rsun, 0.5 Msun photosphere has approximately 1.4415 Lsun; "
            "Mdot converts approximate source-heating targets 1.5–3 Lsun using Lacc=G*M*Mdot/R. "
            "Source heating is distinct from observed bolometric or aperture luminosity. "
            "Coated-Mie envelope dust with 5% ice mantle volume and the nine strict continuum "
            "anchors are retained; no ice/silicate-band inference. Fresh temperature for every "
            "model. All packet counts rise fourfold to 512000; explicit MCFOST seed 41001 is "
            "new relative to the previous unseeded runs, so control differences cannot isolate "
            "photon count alone. Numerical convergence and unique parameter inference remain "
            "unestablished. See validation/continuum_production_v1.1/design.json for design "
            "seed, constants, model order, and catalogue hashes."
        ),
        "fixed": FIXED.copy(), "models": models,
        "numerics": {"photons_temperature": 512000, "photons_image": 512000,
                     "photons_sed": 512000, "grid_nr": 100, "grid_ntheta": 70,
                     "grid_n_inner": 20, "grains": 50, "random_seed": MCFOST_SEED},
        "observations": {"preset": "strict9", "aperture_radius_arcsec": 1.0,
                         "target_distance_pc": 147.0, "aperture_subpixels": 64},
    }
    if len(models) != 96 or models[0] != {**{k: float(v) for k, v in FIXED.items()}, "inclination_deg": 70.0}:
        raise ValueError("Unexpected model count or nominal control")
    csv_stream = io.StringIO(newline="")
    writer = csv.DictWriter(csv_stream, fieldnames=list(catalogue[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(catalogue)
    config_data = json_bytes(config)
    catalogue_data = csv_stream.getvalue().encode()
    design = {
        "schema_version": 1, "design_version": DESIGN_VERSION, "run_name": RUN_NAME,
        "status": "specified_catalogue_not_a_completed_simulation_or_posterior",
        "model_count": 96, "structural_heating_families": 16,
        "continuum_anchors_per_model": 9, "image_calculations": 864, "temperature_calculations": 96,
        "inclination_order_deg": INCLINATIONS,
        "model_order": "Family-major; index = 6 * family_index + position in inclination_order_deg. Index 0 is nominal at 70 degrees.",
        "inclination_priority_range_deg": [50, 65], "models_in_priority_range": 64,
        "sampling": {
            "method": "Endpoint-inclusive Latin hypercube; best minimum Euclidean distance among random permutations in transformed unit cube",
            "design_seed": DESIGN_SEED, "generator": "Python standard-library random.Random / shuffle (MT19937)",
            "candidate_count": N_CANDIDATES, "selected_candidate_index_zero_based": candidate,
            "space_filling_family_count": N_LATIN, "rank_denominator": N_LATIN - 1,
            "minimum_normalized_pair_distance": separation, "selected_latin_rank_rows": ranks,
            "ties": "First candidate wins; distance comparisons use squared integer rank distances",
            "axes": [{"name": key, "bounds": [low, high], "transform": transform, "unit": unit}
                     for key, low, high, transform, unit in AXES],
        },
        "heating_conversion": {
            "photospheric_formula": "Lstar = 4*pi*Rstar^2*sigma*Teff^4",
            "accretion_formula": "Lacc = G*Mstar*Mdot/Rstar; Mdot = (Lsource_target-Lstar)*Rstar/(G*Mstar)",
            "constants_si": CONSTANTS,
            "constants_provenance": "Explicit values matching Astropy 7.1.1 solar/CODATA constants and Julian year; no runtime Astropy dependency",
            "fixed_photospheric_luminosity_lsun": photospheric_luminosity(),
            "nominal_accretion_luminosity_lsun": accretion_luminosity(FIXED["accretion_rate_msun_per_year"]),
            "nominal_source_luminosity_lsun": nominal_source,
            "scope": "Approximate input source-heating targets, not observed Lbol or JWST-aperture luminosity. Confirm emitted luminosities in the actual MCFOST logs.",
            "source": "https://github.com/cpinte/mcfost/blob/main/src/stars.f90",
            "source_check": "Current official main source reviewed; exact af0dec17 cluster source could not be retrieved in this session",
        },
        "numerics": config["numerics"], "observations": config["observations"],
        "families": families,
        "limits": [
            "Only 16 shared dust/heating settings are sampled; six inclinations per setting permit matched geometry comparisons, not exhaustive refitting.",
            "No parameter bounds are claimed as literature confidence limits; disk structure, dust representation, ice fraction and photosphere stay fixed.",
            "MCFOST random seed 41001 is distinct from design seed 20260919. Older runs had no explicit simulation seed.",
            "Fourfold packets may reduce Monte Carlo scatter but do not establish convergence; fixed-seed controls are not independent noise replicates.",
            "Nine aperture continuum measurements cannot uniquely constrain all searched axes; report best tested settings and structured residuals.",
        ],
        "source_documents": ["GOAL.md", "config/restart_proposal.json",
                             "reports/restart_review_2026-09-18/literature_audit.md",
                             "validation/continuum_production_v1/RESULTS_REVIEW.md"],
        "file_sha256": {
            str(CONFIG_PATH.relative_to(ROOT)): digest(config_data),
            str((OUT / "model_catalogue.csv").relative_to(ROOT)): digest(catalogue_data),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).read_bytes()),
        },
    }
    return config_data, catalogue_data, json_bytes(design), families


def coverage_figure(families):
    with tempfile.TemporaryDirectory(prefix="mcfost-v11-plot-") as cache:
        os.environ.setdefault("MPLCONFIGDIR", cache)
        os.environ.setdefault("XDG_CACHE_HOME", cache)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(2, 3, figsize=(12.5, 7.5), constrained_layout=True)
        labels = ["Envelope dust mass (M☉)", "Envelope maximum grain size (µm)",
                  "Grain exponent q", "Cavity half-opening (°)", "Source heating target (L☉)"]
        for ax, axis, label in zip(axes.flat, AXES, labels):
            key, low, high, transform, _ = axis
            values = [family[key] if key in family else family["parameters"][key] for family in families]
            ax.scatter(range(3, 16), values[3:], color="#246a9b", s=38, label="Joint search settings")
            ax.scatter(range(3), values[:3], color="#a84324", marker="*", s=110, label="Archived controls")
            if transform == "log":
                ax.set_yscale("log")
            ax.axhline(low, color="0.7", lw=0.6, ls=":")
            ax.axhline(high, color="0.7", lw=0.6, ls=":")
            ax.set(xlabel="Shared setting index", ylabel=label, xticks=[0, 3, 6, 9, 12, 15])
            ax.grid(alpha=0.2)
        ax = axes.flat[-1]
        ordered = sorted(INCLINATIONS)
        ax.bar(ordered, [16] * 6, width=3, color=["#b8cbd7" if value in (45, 70) else "#246a9b" for value in ordered])
        ax.set(xlabel="Inclination from polar axis (°)", ylabel="Models per inclination", xticks=ordered, ylim=(0, 19))
        ax.text(0.5, 0.95, "The same 16 settings at every inclination", transform=ax.transAxes,
                ha="center", va="top", fontsize=9)
        axes.flat[0].legend(fontsize=8, loc="best")
        figure.suptitle("Continuum v1.1: 16 shared settings × 6 inclinations = 96 models\nExploratory design bounds; nine continuum anchors; 512,000 packets", fontsize=13)
        figure.savefig(OUT / "coverageplot.png", dpi=180)
        figure.savefig(OUT / "coverageplot.pdf", metadata={"CreationDate": None, "ModDate": None})
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify existing config/catalogue/provenance without writing")
    parser.add_argument("--skip-plot", action="store_true", help="Write only standard-library catalogue artifacts")
    args = parser.parse_args()
    run = ROOT / "runs" / RUN_NAME
    if not args.check and run.exists():
        raise SystemExit(f"Refusing to rewrite the design because the run exists: {run}; use --check")
    config, catalogue, design, families = make_outputs()
    files = {CONFIG_PATH: config, OUT / "model_catalogue.csv": catalogue, OUT / "design.json": design}
    if args.check:
        for path, expected in files.items():
            if not path.is_file() or path.read_bytes() != expected:
                raise SystemExit(f"Generated content differs or is missing: {path}")
        print("Verified 96-model config, catalogue and design provenance; no files changed.")
        return
    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(content)
        temporary.replace(path)
    if not args.skip_plot:
        coverage_figure(families)
    print(f"Wrote {CONFIG_PATH.relative_to(ROOT)}")
    print(f"Wrote {OUT.relative_to(ROOT)}: 96 unique models, 16 settings, 864 images, 96 temperature solves")
    print("This command did not prepare or submit a run.")


if __name__ == "__main__":
    main()
