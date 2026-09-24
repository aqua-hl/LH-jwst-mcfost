#!/usr/bin/env python3
"""Compact broad-band production screen and predeclared ice/mass contrasts.

Read frozen inputs and measurement receipts, never raw FITS. Scores describe a
finite model screen; response singular values are conditional diagnostics, not
posterior constraints or a count of identifiable physical parameters.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np


BLOCKS = ("nir_continuum", "nir_ice", "miri")
WEIGHTING = {"baseline": {"nir_continuum": .45, "nir_ice": .30, "miri": .25},
             "nir_only": {"nir_continuum": .6, "nir_ice": .4, "miri": 0.},
             "nir_miri_equal": {"nir_continuum": .3, "nir_ice": .2, "miri": .5}}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    def pairs(values):
        result = {}
        for key, value in values:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def invalid(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")

    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def positive(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value) and value > 0


def hash_string(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def contained(root, relative):
    relative = Path(relative)
    path = root / relative
    require(not relative.is_absolute() and ".." not in relative.parts and path.resolve().is_relative_to(root),
            f"Unsafe bundle path: {relative}")
    return path


def _matrix(value, count, name, positive_definite=False):
    array = np.asarray(value, dtype=float)
    require(array.shape == (count, count) and np.all(np.isfinite(array)), f"Invalid {name} dimensions/values")
    require(np.allclose(array, array.T, rtol=1e-10, atol=1e-20), f"Nonsymmetric {name}")
    if positive_definite:
        try:
            np.linalg.cholesky(array)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"Nonpositive-definite {name}") from exc
    else:
        eigen = np.linalg.eigvalsh(array)
        require(eigen.min() >= -1e-10 * max(abs(eigen).max(), 1e-30), f"Nonpositive-semidefinite {name}")
    return array


def segment_key(band):
    return band.get("gain_group", band.get("segment_id", f"{band['instrument']}/{band['segment']}"))


def validate_contract(contract):
    probes, bands = contract["probes"], contract["bands"]
    require(bool(probes) and bool(bands), "Empty production observation contract")
    probe_map = {p["id"]: p for p in probes}
    require(len(probe_map) == len(probes), "Duplicate contract probes")
    require(len({b["id"] for b in bands}) == len(bands), "Duplicate band IDs")
    require(contract["block_weights"] == WEIGHTING["baseline"], "Changed declared block weighting")
    require({b["block"] for b in bands} == set(BLOCKS), "All three observation blocks are required")
    for probe in probes:
        require(positive(probe["wavelength_um"]) and probe.get("score") is False,
                "All monochromatic probes must be positive-wavelength unscored samples")
    for check in contract.get("consistency_checks", []):
        require(check["probe_id"] in probe_map and check["primary_score_weight"] == 0
                and positive(check["observed_flux_jy"]), "Legacy checks must have zero primary weight and valid observations")
    for band in bands:
        require(positive(band["observed_flux_jy"]) and positive(band["formal_sigma_jy"])
                and positive(band["lower_um"]) and band["upper_um"] > band["lower_um"], "Invalid band observation/support")
        require(len(band["probe_ids"]) == 5 and len(set(band["probe_ids"])) == 5
                and all(p in probe_map for p in band["probe_ids"]), "A band needs five unique frozen probes")
        require(np.allclose(band["weights"], np.array([1, 4, 2, 4, 1]) / 12, rtol=0, atol=1e-14),
                "Five-node band weights must be composite Simpson means")
        require(band["coarse_probe_ids"] == band["probe_ids"][::2]
                and np.allclose(band["coarse_weights"], np.array([1, 4, 1]) / 6, rtol=0, atol=1e-14),
                "Three-node quadrature must be nested Simpson means")
        expected = np.linspace(band["lower_um"], band["upper_um"], 5)
        require(np.allclose([probe_map[p]["wavelength_um"] for p in band["probe_ids"]], expected,
                            rtol=1e-10, atol=1e-10), "Band nodes do not span their frozen wavelength interval")
    order = [band["id"] for band in bands]
    baseline_count = 0
    for name, scenario in contract["covariance_scenarios"].items():
        require(scenario["band_order"] == order, f"Covariance band order differs: {name}")
        base = _matrix(scenario["formal_plus_structure_covariance_jy2"], len(bands), name + " base", True)
        gain = _matrix(scenario["relative_gain_covariance_jy2"], len(bands), name + " gain")
        total = _matrix(scenario["conditional_covariance_jy2"], len(bands), name + " total", True)
        require(np.allclose(base + gain, total, rtol=1e-10, atol=1e-20), "Covariance components do not sum")
        keys = {segment_key(b) for b in bands}
        require(keys <= set(scenario["relative_gain_fraction_by_segment"]), "Missing relative-gain prior")
        require(all(positive(scenario["relative_gain_fraction_by_segment"][key]) for key in keys), "Invalid gain prior")
        z, sigmas, _ = gain_design(bands, scenario)
        require(np.allclose((z * sigmas) @ (z * sigmas).T, gain, rtol=1e-9, atol=1e-20),
                "Segment prior design disagrees with gain covariance")
        baseline_count += scenario.get("baseline") is True
    require(baseline_count == 1, "Require exactly one baseline covariance scenario")
    common = contract["common_scale_nuisance"]
    require(common["kind"] == "gaussian_multiplicative" and common["mean"] == 1
            and positive(common["sigma_fraction"]), "Invalid common-scale nuisance prior")
    return contract


def gain_design(bands, scenario):
    segments = sorted({segment_key(b) for b in bands})
    design = np.array([[band["observed_flux_jy"] if segment_key(band) == segment else 0.
                        for segment in segments] for band in bands])
    sigmas = np.array([scenario["relative_gain_fraction_by_segment"][key] for key in segments])
    return design, sigmas, segments


def block_precision(covariance, blocks, weights):
    """D C_base^-1 D, with D_ii=sqrt(w_block / number_of_bands_in_block).

    This is a declared weighted screening metric, not an unweighted Gaussian
    likelihood. Zero-weight blocks are removed before inversion, so the NIR-only
    alternative cannot depend on MIRI observations through covariance terms.
    """
    counts = Counter(blocks)
    require(set(weights) == set(BLOCKS) and all(v >= 0 for v in weights.values())
            and math.isclose(sum(weights.values()), 1), "Invalid normalized block weights")
    active = np.array([weights[block] > 0 for block in blocks])
    selected = np.flatnonzero(active)
    d = np.array([math.sqrt(weights[blocks[i]] / counts[blocks[i]]) for i in selected])
    precision = np.zeros_like(covariance)
    precision[np.ix_(selected, selected)] = d[:, None] * np.linalg.inv(covariance[np.ix_(selected, selected)]) * d[None, :]
    return precision


def profile_screen(model_flux, contract, scenario, weights=None):
    """Joint common-scale and shared-segment gain fit without double counting.

    Minimize (m-o + a*m - Z*g)'P(m-o + a*m - Z*g) plus Gaussian
    nuisance penalties. P uses only formal+structure covariance. Gain covariance
    is represented by g's priors, and is not also inserted into P.
    """
    bands = contract["bands"]
    model = np.asarray(model_flux, dtype=float)
    observed = np.array([b["observed_flux_jy"] for b in bands])
    require(model.shape == observed.shape and np.all(np.isfinite(model)) and np.all(model > 0), "Invalid model band vector")
    base = np.asarray(scenario["formal_plus_structure_covariance_jy2"], dtype=float)
    precision = block_precision(base, [b["block"] for b in bands], weights or WEIGHTING["baseline"])
    z, gain_sigmas, segments = gain_design(bands, scenario)
    design = np.column_stack((model, -z))
    priors = np.r_[contract["common_scale_nuisance"]["sigma_fraction"], gain_sigmas]
    prior_precision = np.diag(1 / priors**2)
    residual = model - observed
    nuisance = np.linalg.solve(design.T @ precision @ design + prior_precision,
                               -design.T @ precision @ residual)
    adjusted = residual + design @ nuisance
    data_term = float(adjusted @ precision @ adjusted)
    penalty = float(nuisance @ prior_precision @ nuisance)
    block_rms, standardized_rms, contributions = {}, {}, {}
    precision_residual = precision @ adjusted
    for block in BLOCKS:
        selected = np.array([band["block"] == block for band in bands])
        block_rms[block] = float(np.sqrt(np.mean(adjusted[selected]**2)))
        standardized_rms[block] = float(np.sqrt(np.mean(adjusted[selected]**2 / np.diag(base)[selected])))
        contributions[block] = float(np.sum(adjusted[selected] * precision_residual[selected]))
    return {"score": data_term + penalty, "weighted_data_term": data_term, "nuisance_penalty": penalty,
            "block_adjusted_residual_rms_jy": block_rms, "block_standardized_residual_rms": standardized_rms,
            "block_weighted_data_contribution": contributions,
            "block_contribution_note": "Sum of r_i*(P*r)_i within each block; sums to weighted_data_term. Cross-block base covariance could make individual contributions signed.",
            "common_scale": float(1 + nuisance[0]),
            "common_scale_prior_sigma_fraction": float(priors[0]),
            "segment_gain_fraction": {key: float(value) for key, value in zip(segments, nuisance[1:])},
            "segment_gain_prior_sigma_fraction": dict(zip(segments, map(float, gain_sigmas))),
            "adjusted_residual_jy": adjusted.tolist(), "block_weights": weights or WEIGHTING["baseline"],
            "interpretation": "regularized block-normalized screening objective; not a posterior or formal chi-square"}


def integrate_bands(flux_by_probe, contract):
    rows = []
    for band in contract["bands"]:
        fine = float(np.dot(band["weights"], [flux_by_probe[p] for p in band["probe_ids"]]))
        coarse = float(np.dot(band["coarse_weights"], [flux_by_probe[p] for p in band["coarse_probe_ids"]]))
        require(positive(fine) and positive(coarse), "Nonpositive integrated band flux")
        change = (fine - coarse) / fine
        rows.append({"band_id": band["id"], "block": band["block"], "model_flux_jy": fine,
                     "three_node_flux_jy": coarse, "five_minus_three_over_five": change,
                     "quadrature_unresolved": abs(change) > .01,
                     "observed_flux_jy": band["observed_flux_jy"],
                     "effective_wavelength_um": .5 * (band["lower_um"] + band["upper_um"])})
    return rows


def project_common_scale(responses, reference_flux, covariance):
    """Whiten then remove a free common-amplitude direction conservatively."""
    response = np.asarray(responses, dtype=float)
    if response.ndim == 1:
        response = response[:, None]
    reference = np.asarray(reference_flux, dtype=float)
    lower = np.linalg.cholesky(np.asarray(covariance, dtype=float))
    white = np.linalg.solve(lower, response)
    common = np.linalg.solve(lower, reference)
    require(np.all(np.isfinite(white)) and float(common @ common) > 0, "Invalid response or common-scale direction")
    projected = white - common[:, None] * ((common @ white) / (common @ common))[None, :]
    return projected


def response_svd(ice_response, mass_response, reference_flux, covariance):
    projected = project_common_scale(np.column_stack((ice_response, mass_response)), reference_flux, covariance)
    singular = np.linalg.svd(projected, compute_uv=False)
    norms = np.linalg.norm(projected, axis=0)
    cosine = float(projected[:, 0] @ projected[:, 1] / np.prod(norms)) if np.all(norms > 1e-12) else None
    ratio = float(singular[1] / singular[0]) if len(singular) >= 2 and singular[0] > 1e-12 else None
    return {"singular_values": singular.tolist(), "small_to_large_singular_value_ratio": ratio,
            "ice_response_norm": float(norms[0]), "mass_response_norm": float(norms[1]),
            "response_cosine": cosine, "near_collinear_diagnostic": ratio is not None and ratio < .05,
            "response_below_one_covariance_unit": bool(max(norms) < 1),
            "diagnostic_cutoff_note": "Ratio <0.05 flags near-collinearity for inspection; it is not an identifiability threshold.",
            "number_of_physical_constraints": None}


def controlled_contrasts(catalogue, results, contract):
    """Use every predeclared matched group, irrespective of observational score."""
    by_index = {r["model_index"]: r for r in results if r["status"] == "screened"}
    groups = defaultdict(list)
    for record in catalogue:
        if record["role"] == "ice_ladder":
            groups[record["contrast_group_id"]].append(record)
    ice_responses, missing, svds = [], [], []
    for group_id, records in sorted(groups.items()):
        records = sorted(records, key=lambda r: r["ice_mass_fraction"])
        bare = next((r for r in records if r["ice_mass_fraction"] == 0), None)
        if bare is None:
            missing.append({"group_id": group_id, "reason": "No predeclared true bare-silicate reference"})
            continue
        for record in records:
            if record is bare:
                continue
            stripped = [{k: v for k, v in r["parameters"].items() if k != "envelope_ice_volume_fraction"}
                        for r in (bare, record)]
            if stripped[0] != stripped[1]:
                missing.append({"group_id": group_id, "model_index": record["index"], "reason": "Other physical parameters differ from bare reference"})
                continue
            if bare["index"] not in by_index or record["index"] not in by_index:
                missing.append({"group_id": group_id, "model_index": record["index"], "reason": "Missing or invalid matched model"})
                continue
            reference = np.array(by_index[bare["index"]]["band_flux_jy"])
            delta = np.array(by_index[record["index"]]["band_flux_jy"]) - reference
            for scenario_id, scenario in contract["covariance_scenarios"].items():
                response = project_common_scale(delta, reference, scenario["conditional_covariance_jy2"])[:, 0]
                ice_responses.append({"group_id": group_id, "model_index": record["index"], "bare_model_index": bare["index"],
                    "inclination_deg": record["parameters"]["inclination_deg"],
                    "envelope_dust_mass_msun": record["parameters"]["envelope_dust_mass_msun"],
                    "ice_mass_fraction": record["ice_mass_fraction"], "covariance_scenario": scenario_id,
                    "whitened_common_scale_removed_response": response.tolist(),
                    "response_norm": float(np.linalg.norm(response)),
                    "quadrature_unresolved": bool(by_index[bare["index"]]["quadrature_unresolved"] or by_index[record["index"]]["quadrature_unresolved"])})
    ladders = [r for r in catalogue if r["role"] == "ice_ladder"]
    inclinations = sorted({r["parameters"]["inclination_deg"] for r in ladders})
    for inclination in inclinations:
        selected = [r for r in ladders if r["parameters"]["inclination_deg"] == inclination]
        masses = sorted({r["parameters"]["envelope_dust_mass_msun"] for r in selected})
        if len(masses) != 2:
            missing.append({"inclination_deg": inclination, "reason": "Two predeclared envelope masses required for matched SVD"})
            continue
        lookup = {(r["parameters"]["envelope_dust_mass_msun"], r["ice_mass_fraction"]): r for r in selected}
        for fraction in sorted({r["ice_mass_fraction"] for r in selected} - {0.}):
            keys = [(masses[0], 0.), (masses[0], fraction), (masses[1], fraction)]
            if not all(k in lookup and lookup[k]["index"] in by_index for k in keys):
                missing.append({"inclination_deg": inclination, "ice_mass_fraction": fraction, "reason": "Incomplete matched ice/mass rectangle"})
                continue
            records = [lookup[k] for k in keys]
            allowed = {"envelope_dust_mass_msun", "envelope_ice_volume_fraction"}
            stripped = [{k: v for k, v in r["parameters"].items() if k not in allowed} for r in records]
            if not all(value == stripped[0] for value in stripped):
                missing.append({"inclination_deg": inclination, "ice_mass_fraction": fraction, "reason": "Other physical parameters differ in matched rectangle"})
                continue
            fluxes = [np.array(by_index[r["index"]]["band_flux_jy"]) for r in records]
            ice = (fluxes[1] - fluxes[0]) * (.04 / fraction)
            mass = (fluxes[2] - fluxes[1]) * (.000075 / (masses[1] - masses[0]))
            for scenario_id, scenario in contract["covariance_scenarios"].items():
                svds.append({"inclination_deg": inclination, "ice_mass_fraction": fraction,
                    "model_indices": [r["index"] for r in records], "covariance_scenario": scenario_id,
                    "ice_coordinate_step": .04, "mass_coordinate_step_msun": .000075,
                    "definition": "Ice secant from bare at lower mass, scaled to delta ice mass fraction 0.04; mass secant at this ice fraction, scaled to 7.5e-5 solar masses.",
                    "quadrature_unresolved": any(by_index[r["index"]]["quadrature_unresolved"] for r in records),
                    **response_svd(ice, mass, fluxes[0], scenario["conditional_covariance_jy2"])})
    return {"selection": "all predeclared ice ladders, without selecting on fit", "ice_responses": ice_responses,
            "ice_mass_svd": svds, "missing_contrasts": missing,
            "conditional_only": True, "number_of_physical_constraints": None,
            "common_scale_policy": "Project out a free common amplitude after full conditional-covariance whitening; more conservative than the 3% scoring prior.",
            "limitations": "Fixed inclination/background comparisons in one generic ice material family; no posterior, global nuisance marginalization, numerical-noise covariance, material selection or H2O30K inference."}


def _load_quality(run, manifest):
    relative = "code/src/mcfost_grid/quality.py"
    require(relative in manifest["input_hashes"], "Frozen quality policy is not pinned")
    name = "_production_quality_" + hashlib.sha256(str(run).encode()).hexdigest()[:16]
    spec = importlib.util.spec_from_file_location(name, run / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate_model(model, payload, manifest, manifest_hash, quality, bound_fingerprint):
    result = {"model_id": model["id"], "model_index": model["index"], "parameters": model["parameters"],
              "status": "missing", "reason": "measurements.json is absent"}
    if payload is None:
        return result
    try:
        require(isinstance(payload, dict), "Measurement payload must be an object")
        require(payload.get("model_id") == model["id"] and payload.get("model_index") == model["index"]
                and payload.get("parameters") == model["parameters"], "Model identity/parameters differ from manifest")
        fingerprint = payload.get("fingerprint", {})
        require(fingerprint.get("manifest_sha256") == manifest_hash
                and fingerprint.get("backend") == "image_method2"
                and all(hash_string(fingerprint.get(k)) for k in ("executable_sha256", "utilities_content_sha256")),
                "Invalid manifest/executable/utilities/backend fingerprint")
        require(bound_fingerprint is not None and fingerprint == bound_fingerprint, "Fingerprint differs from the shared runtime binding")
        require(payload.get("complete") is True, "Model is incomplete")
        require(payload.get("smoke_test") is not True, "Smoke-test products cannot enter a production screen")
        measurements = payload.get("measurements", [])
        require(isinstance(measurements, list) and all(isinstance(r, dict) for r in measurements), "Invalid measurement list")
        lookup = {r.get("anchor_id"): r for r in measurements}
        anchors = {a["id"]: a for a in manifest["anchors"]}
        require(len(lookup) == len(measurements) and set(lookup) == set(anchors), "Missing, duplicate or unexpected probes")
        for key, row in lookup.items():
            anchor = anchors[key]
            require(row.get("quality_pass") is True and row.get("valid") is True and positive(row.get("flux_jy")),
                    f"Invalid aperture flux/quality: {key}")
            error = quality.cached_policy_error(row, manifest["measurement"]["quality_policy"])
            require(not error, f"{key}: {error}")
            require(row.get("backend") == "image_method2" and hash_string(row.get("image_sha256")), "Missing image provenance")
            require(math.isclose(row["wavelength_um"], anchor["wavelength_um"], rel_tol=1e-10)
                    and row["instrument"] == anchor["instrument"], f"Probe wavelength/instrument mismatch: {key}")
            for field, expected in (("model_distance_pc", model["parameters"]["distance_pc"]),
                                    ("target_distance_pc", manifest["measurement"]["target_distance_pc"]),
                                    ("aperture_radius_arcsec", manifest["measurement"]["aperture_radius_arcsec"])):
                require(math.isclose(row[field], expected, rel_tol=1e-10), f"Measurement {field} differs: {key}")
            diagnostics = row["diagnostics"]
            expected_npix = anchor.get("image_npix", manifest["configuration"]["numerics"]["image_npix"])
            require(diagnostics.get("nx") == diagnostics.get("ny") == expected_npix, f"Wrong image dimensions: {key}")
            size = anchor.get("image_size_au", manifest["configuration"]["numerics"]["image_size_au"])
            pixel = size / expected_npix / model["parameters"]["distance_pc"]
            require(math.isclose(diagnostics["model_pixel_arcsec"], pixel, rel_tol=5e-6), f"Wrong pixel scale: {key}")
            ratio = model["parameters"]["distance_pc"] / manifest["measurement"]["target_distance_pc"]
            require(math.isclose(diagnostics["distance_flux_factor"], ratio**2, rel_tol=1e-10)
                    and math.isclose(diagnostics["target_pixel_arcsec"], diagnostics["model_pixel_arcsec"] * ratio, rel_tol=1e-10)
                    and diagnostics.get("zero_clipping_used_for_scoring") is False, f"Changed signed-distance estimator: {key}")
            require(row.get("canonical_estimator") == "signed_direct_method2_polarized_total_i_plane0",
                    f"Changed total-I estimator: {key}")
        return {**result, "status": "valid", "reason": "", "flux_by_probe": {key: row["flux_jy"] for key, row in lookup.items()},
                "quality_warnings": {key: row.get("quality_warnings", []) for key, row in lookup.items() if row.get("quality_warnings")}}
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as exc:
        return {**result, "status": "excluded", "reason": str(exc)}


def _write_csv(path, rows, fields):
    fields = list(fields) + sorted({key for row in rows for key in row} - set(fields))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, sort_keys=True, allow_nan=False) if isinstance(value, (dict, list)) else value
                             for key, value in row.items()})


def _plots(out, summary, contract):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bands = contract["bands"]
    x = np.array([.5 * (b["lower_um"] + b["upper_um"]) for b in bands])
    observed = np.array([b["observed_flux_jy"] for b in bands])
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    axes[0].scatter(x, observed, color="black", label="Observed broad bands", zorder=5)
    screened = sorted((r for r in summary["models"] if r["status"] == "screened"), key=lambda r: r["baseline_score"])
    for result in screened[:5]:
        scale = result["baseline_profile"]["common_scale"]
        flux = np.array(result["band_flux_jy"]) * scale
        axes[0].plot(x, flux, "o-", alpha=.65, label=result["model_id"] + (" [quadrature flag]" if result["quadrature_unresolved"] else ""))
        axes[1].plot(x, np.log10(flux / observed), "o-", alpha=.65)
    if not screened:
        axes[0].text(.5, .5, "No complete valid models", ha="center", transform=axes[0].transAxes)
    axes[0].set(ylabel="Band mean Fν (Jy)", yscale="log", title="Broad-band production screen; five smallest screening scores for display")
    axes[1].set(xlabel="Band midpoint wavelength (µm)", ylabel="log10(scaled model / observed)", xscale="log")
    axes[1].axhline(0, color="grey", lw=.7)
    axes[0].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.15)
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"band_comparison.{suffix}", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    responses = summary["controlled_contrasts"]["ice_responses"]
    selected = [r for r in responses if r["covariance_scenario"] == summary["baseline_covariance_scenario"]]
    for group in sorted({r["group_id"] for r in selected}):
        rows = sorted((r for r in selected if r["group_id"] == group), key=lambda r: r["ice_mass_fraction"])
        ax.plot([r["ice_mass_fraction"] for r in rows], [r["response_norm"] for r in rows], "o-", label=group)
    if selected:
        ax.legend(fontsize=7, ncol=3)
    else:
        ax.text(.5, .5, "Matched ice/bare responses unavailable", ha="center", transform=ax.transAxes)
    ax.set(xlabel="Ice mass fraction", ylabel="Covariance-whitened response norm after removing common scale",
           title="Every available predeclared ice ladder; conditional response, not significance or posterior")
    ax.grid(alpha=.15)
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"ice_response.{suffix}", dpi=180)
    plt.close(fig)


def analyze_run(run, output=None, make_plot=True):
    run = Path(run).resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mcfost-production-mpl"))
    output = Path(output).resolve() if output else run / "results"
    require(not output.is_relative_to(run) or output.is_relative_to(run / "results"), "Reports within a run must be under results/")
    manifest = read_json(run / "manifest.json")
    sidecar = run / "manifest.sha256"
    if sidecar.is_file():
        require(sidecar.read_text().split() == [digest(run / "manifest.json"), "manifest.json"], "Manifest checksum differs")
    hashes = manifest.get("input_hashes", {})
    require(isinstance(hashes, dict) and bool(hashes), "Missing frozen input hashes")
    for relative, expected in hashes.items():
        path = contained(run, relative)
        require(not path.is_symlink() and path.is_file() and digest(path) == expected, f"Frozen input changed or missing: {relative}")
    for required in ("production_experiment.json", "inputs/production_observation_contract.json"):
        require(required in hashes, f"Unpinned production input: {required}")
    experiment = read_json(run / "production_experiment.json")
    contract = validate_contract(read_json(run / "inputs/production_observation_contract.json"))
    require(manifest["configuration"].get("numerical_only") is True, "Production runner must use explicit unscored numerical probes")
    anchors = manifest["anchors"]
    require(len(anchors) == len(contract["probes"]) and {a["id"]: a["wavelength_um"] for a in anchors}
            == {p["id"]: p["wavelength_um"] for p in contract["probes"]}
            and all(a.get("score") is False for a in anchors), "Manifest probes differ from observation contract")
    catalogue = experiment["catalogue"]
    require(len(catalogue) == len(manifest["models"]) and [r["index"] for r in catalogue] == list(range(len(catalogue))),
            "Production catalogue differs from manifest model indices")
    for model, record in zip(manifest["models"], catalogue):
        require(model["index"] == record["index"] and model["parameters"] == record["parameters"], "Catalogue parameters differ from manifest")
    quality = _load_quality(run, manifest)
    manifest_hash = digest(run / "manifest.json")
    binding_path = run / "runtime_binding.json"
    binding = read_json(binding_path).get("fingerprint") if binding_path.is_file() else None
    production_runtime_path = run / "production_runtime_binding.json"
    production_runtime = ({"sha256": digest(production_runtime_path), "record": read_json(production_runtime_path)}
                          if production_runtime_path.is_file() else None)
    baseline_id = next(key for key, scenario in contract["covariance_scenarios"].items() if scenario.get("baseline") is True)
    results = []
    for model in manifest["models"]:
        path = contained(run, f"models/{model['id']}/measurements.json")
        try:
            payload = read_json(path) if path.is_file() else None
            result = evaluate_model(model, payload, manifest, manifest_hash, quality, binding)
        except (OSError, ValueError, TypeError) as exc:
            result = {"model_id": model["id"], "model_index": model["index"], "parameters": model["parameters"],
                      "status": "excluded", "reason": str(exc)}
        result["catalogue"] = catalogue[model["index"]]
        if result["status"] == "valid":
            predictions = integrate_bands(result["flux_by_probe"], contract)
            vector = [row["model_flux_jy"] for row in predictions]
            profiles = {scenario_id: {weight_id: profile_screen(vector, contract, scenario, weights)
                                     for weight_id, weights in WEIGHTING.items()}
                        for scenario_id, scenario in contract["covariance_scenarios"].items()}
            result.update(status="screened", band_predictions=predictions, band_flux_jy=vector,
                          quadrature_unresolved=any(row["quadrature_unresolved"] for row in predictions),
                          unresolved_quadrature_bands=[row["band_id"] for row in predictions if row["quadrature_unresolved"]],
                          profiles=profiles, baseline_profile=profiles[baseline_id]["baseline"],
                          baseline_score=profiles[baseline_id]["baseline"]["score"],
                          interpretation="Finite-grid screening only; unresolved quadrature retained and flagged")
            result["consistency_checks"] = [{**check, "model_flux_jy": result["flux_by_probe"][check["probe_id"]],
                "common_scale_adjusted_model_flux_jy": result["flux_by_probe"][check["probe_id"]]
                    * result["baseline_profile"]["common_scale"]} for check in contract.get("consistency_checks", [])]
        results.append(result)
    contrasts = controlled_contrasts(catalogue, results, contract)
    screened = sorted((r for r in results if r["status"] == "screened"), key=lambda r: (r["baseline_score"], r["model_index"]))
    summary = {"schema_version": 1, "analysis": "conditional_extinction_ice_production_screen",
               "run": str(run), "manifest_sha256": manifest_hash, "analyzer_sha256": digest(__file__),
               "observation_contract_sha256": hashes["inputs/production_observation_contract.json"],
               "status": "complete_screen" if len(screened) == len(results) else "partial_screen",
               "models_total": len(results), "models_screened": len(screened),
               "model_state_counts": dict(Counter(r["status"] for r in results)), "models": results,
               "baseline_covariance_scenario": baseline_id, "block_weighting": WEIGHTING,
               "score_definition": "min (m-o+a*m-Z*g)' D C_formal+structure^-1 D (m-o+a*m-Z*g) + (a/sigma_scale)^2 + sum(g/sigma_gain)^2; D_ii=sqrt(block_weight/number_of_bands_in_block). Gains shared jointly by instrument/segment, including segments shared across NIR blocks. No gain-covariance double counting.",
               "controlled_contrasts": contrasts, "raw_images_required": False, "raw_fits_revalidated": False,
               "runtime_fingerprint": binding, "production_convergence_certified": False,
               "production_runtime_binding": production_runtime,
               "posterior_intervals_reported": False, "number_of_physical_constraints": None,
               "quadrature_flagged_models": sum(r["quadrature_unresolved"] for r in screened),
               "limitations": ["This is a conditional screen of a finite catalogue, not a likelihood/posterior or proof of identified ice properties.",
                   "MIRI artifacts are not mapped; structure covariance is a provisional descriptive model, not measured calibration truth.",
                   "Three-versus-five-node quadrature differences flag unresolved bands without dropping observations or models.",
                   "A single image seed per model does not estimate numerical covariance; no numerical-convergence certification is made.",
                   "Source-centred 1-arcsecond aperture and both model/target distances are retained; no 0.35-arcsecond H2O operator is inferred.",
                   "Generic ice and true bare-silicate controls only; H2O30K laboratory constants are not used.",
                   "Compact receipts and all frozen inputs are verified; raw FITS contents are not rehashed by this analysis."]}
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    ranking = [{"rank": rank, "model_id": row["model_id"], "model_index": row["model_index"], "score": row["baseline_score"],
                "quadrature_unresolved": row["quadrature_unresolved"], "unresolved_bands": row["unresolved_quadrature_bands"],
                "common_scale": row["baseline_profile"]["common_scale"], "role": row["catalogue"]["role"],
                "scenario_scores": {scenario_id: {weight_id: profile["score"] for weight_id, profile in profiles.items()}
                                    for scenario_id, profiles in row["profiles"].items()},
                "parameters": row["parameters"]} for rank, row in enumerate(screened, 1)]
    excluded = [r for r in results if r["status"] != "screened"]
    _write_csv(output / "ranking.csv", ranking, ["rank", "model_id", "model_index", "score", "quadrature_unresolved"])
    _write_csv(output / "excluded_models.csv", excluded, ["model_id", "model_index", "status", "reason"])
    _write_csv(output / "band_predictions.csv", [{"model_id": r["model_id"], **b} for r in screened for b in r["band_predictions"]],
               ["model_id", "band_id", "block", "model_flux_jy", "observed_flux_jy"])
    _write_csv(output / "ice_responses.csv", contrasts["ice_responses"], ["group_id", "model_index", "ice_mass_fraction", "covariance_scenario", "response_norm"])
    _write_csv(output / "ice_mass_svd.csv", contrasts["ice_mass_svd"], ["inclination_deg", "ice_mass_fraction", "covariance_scenario"])
    _write_csv(output / "consistency_checks.csv", [{"model_id": r["model_id"], **check}
               for r in screened for check in r["consistency_checks"]], ["model_id", "id", "model_flux_jy", "observed_flux_jy", "primary_score_weight"])
    lines = ["# Extinction/ice production screen", "", f"**{summary['status']}**: {len(screened)}/{len(results)} models screened.", "",
             "The ranking is a block-normalized, nuisance-profiled finite-grid screen. It is not a posterior or a numerical-convergence certificate.", "",
             f"{summary['quadrature_flagged_models']} screened models have at least one band with >1% three-versus-five-node quadrature change. These models and bands remain in all reported comparisons.", "",
             "Baseline block weights are 45% NIR continuum, 30% NIR ice and 25% MIRI, divided by each block's band count. NIR-only and equal NIR/MIRI alternatives are included. Legacy monochromatic checks do not enter these scores.", "",
             "All available predeclared bare/ice ladders enter the controlled contrasts, regardless of fit. Full conditional covariance is whitened and a free common amplitude removed before ice/mass response diagnostics. Small or collinear responses are valid outcomes; no count of physical constraints is assigned.", "",
             f"Missing/invalid matched contrasts: {len(contrasts['missing_contrasts'])}.", "", "## Limitations", "",
             *[f"- {item}" for item in summary["limitations"]], "", "See `summary.json` for every nuisance fit, covariance scenario, model exclusion and matched contrast; CSV files and figures provide compact views.", ""]
    (output / "REVIEW.md").write_text("\n".join(lines))
    if make_plot:
        _plots(output, summary, contract)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = analyze_run(args.run, args.output, not args.no_plot)
    except (OSError, ValueError, TypeError, KeyError, ImportError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(json.dumps({key: result[key] for key in ("status", "models_total", "models_screened", "quadrature_flagged_models")}, indent=2))
    return 0 if result["status"] == "complete_screen" else 2


if __name__ == "__main__":
    raise SystemExit(main())
