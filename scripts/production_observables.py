#!/usr/bin/env python3
"""Immutable one-arcsecond broadband contract for the 1 AU production run.

The data operator is an exact wavelength-width average of the existing R100
piecewise-constant F_nu estimates, using one named detector segment per band.
Model means use nested five/three-node composite Simpson quadrature. The
three-node estimate is a diagnostic, never an additional fitted observation.
No empirical calibration covariance is available: all gain priors and the
diagonal treatment of measured within-bin structure are declared assumptions.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table


BLOCK_WEIGHTS = {"nir_continuum": .45, "nir_ice": .30, "miri": .25}
# Fixed before any new model result is inspected. Detector gaps and the
# continuum-mask intervals constrain the clean continuum choices.
WINDOWS = (
    ("c01", "nir_continuum", "NIRSpec", 1.10, 1.16, "G140H_F100LP"),
    ("c02", "nir_continuum", "NIRSpec", 1.18, 1.24, "G140H_F100LP"),
    ("c03", "nir_continuum", "NIRSpec", 1.492, 1.55, "G140H_F100LP"),
    ("c04", "nir_continuum", "NIRSpec", 1.70, 1.78, "G140H_F100LP"),
    ("c05", "nir_continuum", "NIRSpec", 2.16, 2.24, "G235H_F170LP"),
    ("c06", "nir_continuum", "NIRSpec", 2.26, 2.34, "G235H_F170LP"),
    ("c07", "nir_continuum", "NIRSpec", 2.52, 2.58, "G235H_F170LP"),
    ("c08", "nir_continuum", "NIRSpec", 3.86, 3.96, "G395H_F290LP"),
    ("c09", "nir_continuum", "NIRSpec", 4.505, 4.545, "G395H_F290LP"),
    ("h01", "nir_ice", "NIRSpec", 2.70, 2.84, "G235H_F170LP"),
    ("h02", "nir_ice", "NIRSpec", 2.88, 3.02, "G395H_F290LP"),
    ("h03", "nir_ice", "NIRSpec", 3.04, 3.18, "G395H_F290LP"),
    ("h04", "nir_ice", "NIRSpec", 3.22, 3.50, "G395H_F290LP"),
    ("s01", "miri", "MIRI", 8.2, 8.6, "CH2_short"),
    ("s02", "miri", "MIRI", 9.3, 10.1, "CH2_medium"),
    ("s03", "miri", "MIRI", 10.7, 11.3, "CH2_long"),
    ("s04", "miri", "MIRI", 12.4, 13.0, "CH3_short"),
    ("s05", "miri", "MIRI", 18.1, 18.8, "CH4_short"),
    ("s06", "miri", "MIRI", 19.3, 20.3, "CH4_short"),
)
CHECKS = (("w007", 5.390628315339, "CH1_short"),
          ("w008", 13.7999, "CH3_medium"),
          ("w009", 27.5131, "CH4_long"))
FINE_WEIGHTS = (1 / 12, 4 / 12, 2 / 12, 4 / 12, 1 / 12)
COARSE_WEIGHTS = (1 / 6, 4 / 6, 1 / 6)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_window(table: Table, native: Table, window: tuple) -> dict:
    """Select exact partial-bin widths; reject gaps, duplicate rows and stitching.

    Native support is checked separately from R100 geometric support. Original
    masks, clipping and retention counts remain visible in the returned rows.
    Partial bins inherit their existing full-bin flux estimate and uncertainty.
    """
    name, block, instrument, low, high, segment = window
    if not np.isfinite([low, high]).all() or not 0 < low < high:
        raise ValueError(f"{name}: invalid window bounds")
    nmask = ((np.asarray(native["instrument"], str) == instrument) &
             (np.asarray(native["segment"], str) == segment))
    nlam = np.asarray(native["wavelength_um"], float)[nmask]
    nflux = np.asarray(native["flux_jy"], float)[nmask]
    nerr = np.asarray(native["error_jy"], float)[nmask]
    nvalid = np.isfinite(nlam) & np.isfinite(nflux) & np.isfinite(nerr) & (nerr > 0)
    if "coverage_fraction" in native.colnames:
        nvalid &= np.asarray(native["coverage_fraction"], float)[nmask] >= .9
    native_in_window = np.isfinite(nlam) & (nlam >= low) & (nlam <= high)
    native_samples = int(np.count_nonzero(native_in_window))
    native_quality_samples = int(np.count_nonzero(native_in_window & nvalid))
    nlam = nlam[nvalid]
    if len(nlam) < 2 or low < nlam.min() or high > nlam.max():
        raise ValueError(f"{name}: window exceeds native support of {segment}; no stitching allowed")
    left = np.asarray(table["bin_left_um"], float)
    right = np.asarray(table["bin_right_um"], float)
    selected = np.flatnonzero(
        (np.asarray(table["instrument"], str) == instrument) &
        (np.asarray(table["segment"], str) == segment) & (left < high) & (right > low))
    selected = selected[np.argsort(left[selected], kind="stable")]
    if not len(selected):
        raise ValueError(f"{name}: no R100 bins selected")
    start, stop = np.maximum(left[selected], low), np.minimum(right[selected], high)
    tol = 1e-10 * max(1., high)
    if (abs(start[0] - low) > tol or abs(stop[-1] - high) > tol or
            np.any(start[1:] - stop[:-1] > tol)):
        raise ValueError(f"{name}: incomplete R100 bin coverage")
    if np.any(start[1:] - stop[:-1] < -tol):
        raise ValueError(f"{name}: duplicate or overlapping source bins")
    widths = stop - start
    if np.any(widths <= 0) or not np.isclose(widths.sum(), high - low, rtol=0, atol=tol):
        raise ValueError(f"{name}: invalid overlap widths")
    weights = widths / (high - low)
    flux = np.asarray(table["flux_jy"], float)[selected]
    formal = np.asarray(table["formal_error_jy"], float)[selected]
    scatter = np.asarray(table["within_bin_scatter_jy"], float)[selected]
    if (not np.isfinite(np.concatenate([flux, formal, scatter])).all() or
            np.any(flux <= 0) or np.any(formal <= 0) or np.any(scatter < 0)):
        raise ValueError(f"{name}: selected bins require positive finite flux/formal error and nonnegative scatter")
    if "n_used" in table.colnames and np.any(np.asarray(table["n_used"])[selected] <= 0):
        raise ValueError(f"{name}: selected bin has no retained samples")
    mean = float(weights @ flux)
    # No independent-bin or sqrt(N) assumption is available from the source.
    formal_sigma = float(weights @ formal)
    structure_fraction = float(np.sqrt(weights @ (scatter / flux) ** 2))
    rows = []
    for index, lo, hi, width, weight in zip(selected, start, stop, widths, weights):
        row = table[int(index)]
        entry = {
            "source_row_zero_based": int(index), "instrument": instrument, "segment": segment,
            "bin_index": int(row["bin_index"]), "wavelength_um": float(row["wavelength_um"]),
            "bin_left_um": float(left[index]), "bin_right_um": float(right[index]),
            "selected_lower_um": float(lo), "selected_upper_um": float(hi),
            "overlap_width_um": float(width), "normalized_weight": float(weight),
            "flux_jy": float(row["flux_jy"]), "formal_error_jy": float(row["formal_error_jy"]),
            "within_bin_scatter_jy": float(row["within_bin_scatter_jy"]),
        }
        for column in ("n_input", "n_quality_valid", "n_manual_masked", "n_auto_clipped", "n_used"):
            if column in table.colnames:
                entry[column] = int(row[column])
        for column in ("used_fraction_of_quality_samples", "median_coverage_fraction"):
            if column in table.colnames:
                entry[column] = float(row[column])
        rows.append(entry)
    return {
        "id": name, "block": block, "instrument": instrument, "lower_um": low, "upper_um": high,
        "segment": segment, "gain_group": f"{instrument}/{segment}",
        "native_support_um": [float(nlam.min()), float(nlam.max())],
        "geometric_bin_coverage_fraction": float(widths.sum() / (high - low)),
        "sum_observation_weights": float(weights.sum()), "observed_flux_jy": mean,
        "formal_sigma_jy": formal_sigma, "structure_fraction": structure_fraction,
        "structure_sigma_jy": mean * structure_fraction,
        "r100_bins": len(selected), "selected_rows": rows,
        "native_sample_quality_coverage_threshold": .9,
        "native_samples_in_window": native_samples,
        "native_quality_valid_samples_in_window": native_quality_samples,
        "native_quality_valid_fraction_in_window": (native_quality_samples / native_samples
                                                     if native_samples else None),
        "continuum_feature_mask_applied": block == "nir_continuum",
    }


def covariance_scenario(bands: list[dict], miri_gain: float, nir_gain: float = .03) -> dict:
    """Separate nuisance covariance from formal/structure terms to avoid double counting."""
    if not np.isfinite([miri_gain, nir_gain]).all() or min(miri_gain, nir_gain) <= 0:
        raise ValueError("Gain fractions must be positive and finite")
    flux = np.array([b["observed_flux_jy"] for b in bands])
    formal = np.diag([b["formal_sigma_jy"] ** 2 for b in bands])
    structure = np.diag([b["structure_sigma_jy"] ** 2 for b in bands])
    groups = np.array([b["gain_group"] for b in bands])
    fractions = np.array([miri_gain if b["instrument"] == "MIRI" else nir_gain for b in bands])
    gain = np.outer(flux * fractions, flux * fractions) * (groups[:, None] == groups[None, :])
    total = formal + structure + gain
    np.linalg.cholesky(total)
    sigma = np.sqrt(np.diag(total))
    return {
        "band_order": [b["id"] for b in bands],
        "block_indices": {block: [i for i, band in enumerate(bands) if band["block"] == block]
                          for block in BLOCK_WEIGHTS},
        "relative_gain_fraction_by_segment": dict(zip(groups.tolist(), fractions.tolist())),
        "formal_covariance_jy2": formal.tolist(),
        "structure_covariance_jy2": structure.tolist(),
        "formal_plus_structure_covariance_jy2": (formal + structure).tolist(),
        "relative_gain_covariance_jy2": gain.tolist(),
        "conditional_covariance_jy2": total.tolist(),
        "conditional_correlation": (total / np.outer(sigma, sigma)).tolist(),
        "baseline": miri_gain == .02 and nir_gain == .03,
        "common_scale_in_covariance": False,
        "treatment": "Use conditional covariance OR profile segment gains with formal_plus_structure; never both.",
    }


def _probe(identifier: str, wavelength: float, instrument: str, region: str) -> dict:
    return {
        "id": identifier, "wavelength_um": float(wavelength), "instrument": instrument,
        "region": region, "score": False, "image_npix": 6001, "image_size_au": 6000.,
        "psf_fwhm_arcsec": .033 * wavelength + (.106 if instrument == "MIRI" else 0.),
    }


def build_contract(root: Path) -> dict:
    """Build from the three tracked immutable reference inputs; write nothing."""
    root = Path(root).resolve()
    r100_path = root / "reference/observations/continuum_sed_R100.ecsv"
    native_path = root / "reference/observations/tmc1a_sed_unstitched.ecsv"
    mask_path = root / "reference/observations/continuum_fit_mask_v1_manifest.json"
    table, native = Table.read(r100_path), Table.read(native_path)
    mask = json.loads(mask_path.read_text())
    if table.meta.get("segment_rescaling_applied") is not False or native.meta.get("segment_rescaling_applied") is not False:
        raise ValueError("Require explicitly unstitched, unrescaled source tables")
    if float(native.meta.get("aperture_radius_arcsec", -1)) != 1.0:
        raise ValueError("Require the native one-arcsecond aperture extraction")
    if mask["source"]["sha256"] != digest(r100_path):
        raise ValueError("Continuum-mask source hash does not match the R100 input")
    intervals = mask["scenarios"]["defer_major_ice_and_both_silicates"]["intervals_um"]
    bands, probes = [], []
    for window in WINDOWS:
        band = select_window(table, native, window)
        if band["block"] == "nir_continuum":
            for interval in intervals:
                lo, hi = interval["lower_um"], interval["upper_um"]
                if band["lower_um"] <= hi and band["upper_um"] >= lo:
                    raise ValueError(f"{band['id']}: continuum window intersects {interval['reason']}")
                if any(lo <= row["wavelength_um"] <= hi for row in band["selected_rows"]):
                    raise ValueError(f"{band['id']}: selected source-bin centre is feature-masked")
        nodes = np.linspace(band["lower_um"], band["upper_um"], 5)
        ids = [f"{band['id']}_q{i}" for i in range(5)]
        probes.extend(_probe(identifier, float(wavelength), band["instrument"], band["block"])
                      for identifier, wavelength in zip(ids, nodes))
        band.update({"probe_ids": ids, "weights": list(FINE_WEIGHTS),
                     "coarse_probe_ids": [ids[0], ids[2], ids[4]],
                     "coarse_weights": list(COARSE_WEIGHTS)})
        bands.append(band)
    # Reusing an R100 estimate in two fitted means creates additional covariance.
    # This fixed design avoids that ambiguity entirely.
    selected_keys = [(r["instrument"], r["segment"], r["source_row_zero_based"])
                     for band in bands for r in band["selected_rows"]]
    if len(selected_keys) != len(set(selected_keys)):
        raise ValueError("A source R100 row contributes to multiple primary bands")
    checks = []
    for identifier, wavelength, segment in CHECKS:
        selected = np.flatnonzero((np.asarray(table["instrument"], str) == "MIRI") &
                                  (np.asarray(table["segment"], str) == segment))
        nearest = int(selected[np.argmin(abs(np.asarray(table["wavelength_um"], float)[selected] - wavelength))])
        # Exact original observed R100 wavelength is used for these check images.
        actual = float(table["wavelength_um"][nearest])
        probes.append(_probe(identifier, actual, "MIRI", "consistency_check"))
        checks.append({"id": identifier, "probe_id": identifier, "wavelength_um": actual,
                       "segment": segment, "observed_flux_jy": float(table["flux_jy"][nearest]),
                       "formal_sigma_jy": float(table["formal_error_jy"][nearest]),
                       "source_row_zero_based": nearest, "primary_score_weight": 0.})
    paths = [r100_path, native_path, mask_path]
    return {
        "schema_version": 1, "contract_id": "extinction_ice_production_1au_observables_v1",
        "observable": "wavelength-width mean Fnu in Jy",
        "observation_operator": "exact R100 overlap-width average, one explicit detector segment per band",
        "model_operator": "five-node composite Simpson wavelength-width mean Fnu",
        "quadrature_diagnostic": "nested three-node Simpson comparison; no extra score term or automatic convergence certification",
        "quadrature_fraction_tolerance": .01,
        "probe_count": len(probes), "primary_band_count": len(bands),
        "probes": probes, "bands": bands, "consistency_checks": checks,
        "block_weights": dict(BLOCK_WEIGHTS),
        "region_weighting_rule": "Weighted mean of each block's nuisance-profiled squared residuals: 0.45 NIR continuum, 0.30 NIR ice, 0.25 MIRI; divide each block sum by its number of bands. Priors are counted once globally.",
        "covariance_scenarios": {f"miri_gain_{gain:.2f}": covariance_scenario(bands, gain)
                                 for gain in (.01, .02, .05)},
        "common_scale_nuisance": {"kind": "gaussian_multiplicative", "mean": 1., "sigma_fraction": .03,
                                  "assumption": "Sensitivity prior, not an empirically measured absolute calibration uncertainty"},
        "aperture_radius_arcsec": 1., "model_distance_pc": 140., "target_distance_pc": 147.,
        "image_npix": 6001, "image_size_au": 6000., "pixel_size_au": 6000. / 6001.,
        "psf_model": "historical Gaussian FWHM=0.033 wavelength_um arcsec; MIRI adds 0.106 arcsec",
        "continuum_mask_scenario": "defer_major_ice_and_both_silicates",
        "feature_bands_unmasked": True, "silicate_9p7_um_retained": True,
        "h2o_30k_constants_required": False,
        "sources": {str(path.relative_to(root)): digest(path) for path in paths},
        "builder_sha256": digest(Path(__file__).resolve()),
        "limitations": [
            "Gain priors and diagonal measured-structure treatment are descriptive sensitivity assumptions, not calibrated likelihoods.",
            "Within-bin structure includes astrophysical slope, residual lines and artifacts; it is not a clean detector-noise estimate.",
            "Formal errors are averaged without sqrt(N) reduction because source inter-bin covariance is unavailable.",
            "Partial R100 bins inherit the full-bin flux estimate; selected-bin masking and clipping diagnostics are retained.",
            "Nine NIR continuum windows stop at 4.545 micrometres because the original feature mask excludes 4.55–5.25 micrometres.",
            "The four broad 3-micrometre bands are one-arcsecond flux proxies, not the separate 0.35-arcsecond R2400 H2O optical-depth contract.",
            "The nested 3/5-node difference is an integration diagnostic, not a bound on unresolved spectral structure.",
            "The 6000 AU map has 6001 pixels for an odd centred grid: 0.999833 AU per pixel, approximately 1 AU.",
            "MIRI legacy anchors w007/w008/w009 have zero primary weight; the six broad feature-inclusive bands carry 25 percent of the regional score.",
            "This contract tests conditional constraints within the supplied dust/ice family and does not establish material identifiability across untested optical constants.",
        ],
    }
