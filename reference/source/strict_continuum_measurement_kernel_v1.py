#!/usr/bin/env python3
"""Versioned Method-2 signed-aperture measurement kernel.

This module is the reusable, rank-blind measurement layer for future
strict-continuum RT tasks.  The numerical operations and scientific constants
are a direct extraction of ``analyze_strict_continuum_refinement_v1.py``:

* polarized Method-2 total-I plane 0 is the canonical observable;
* the aperture is centred on FITS CRPIX and has a nominal 1 arcsec radius;
* boundary pixels use the frozen 64x64 midpoint rule;
* the MCFOST 140 pc product is evaluated at 147 pc for both angular sampling
  and flux scaling;
* the wavelength-dependent Gaussian PSF and all closure, component, edge,
  convolution, and signed-versus-clipped gates are unchanged.

``measure_method2_signed_aperture_v1`` is pure: it receives arrays, scalar
metadata, and an already-decoded closure object, and has no filesystem side
effects.  The FITS/JSON functions below it are deliberately thin, fail-closed
adapters suitable for one Slurm RT task at a time.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter


KERNEL_ID = "strict_continuum_method2_signed_aperture_v1"
SCHEMA_VERSION = 1

# Frozen scientific constants.  Do not make these runtime options: changing
# one requires a new versioned kernel and a new parity exercise.
MODEL_DISTANCE_PC = 140.0
TARGET_DISTANCE_PC = 147.0
DISTANCE_FACTOR = (MODEL_DISTANCE_PC / TARGET_DISTANCE_PC) ** 2
APERTURE_RADIUS_ARCSEC = 1.0
APERTURE_SUBPIXELS = 64
PSF_TRUNCATION_SIGMA = 6.0
COEVAL_TOLERANCE = 1.0e-5
NEGATIVE_FRACTION_MAX = 1.0e-3
EDGE_FRACTION_MAX = 1.0e-3
CONVOLUTION_LOSS_MAX = 1.0e-4
APERTURE_CLIP_MAX = 1.0e-3
COMPONENT_SIGNIFICANT_FRACTION = 1.0e-8
COMPONENT_ACTIVE_FRACTION = 0.01
COMPONENT_APERTURE_FRACTION_MIN = -1.0e-6
COMPONENT_APERTURE_FRACTION_MAX = 1.0001

COMPONENTS = {
    4: "direct_star",
    5: "scattered_star",
    6: "direct_thermal",
    7: "scattered_thermal",
}
CLOSURE_LABELS = {
    0: "total_I",
    4: "direct_star",
    5: "scattered_star",
    6: "thermal_dust",
    7: "scattered_thermal",
}


class MeasurementContractError(ValueError):
    """An input or serialized record violates the v1 measurement contract."""


class MeasurementGateError(MeasurementContractError):
    """A numerically valid product fails one or more frozen quality gates."""

    def __init__(self, message: str, *, failed_checks: tuple[str, ...] = ()) -> None:
        super().__init__(message)
        self.failed_checks = failed_checks


def fixed_policy_v1() -> dict[str, Any]:
    """Return a fresh JSON-safe copy of the immutable v1 policy."""
    return {
        "model_distance_pc": MODEL_DISTANCE_PC,
        "target_distance_pc": TARGET_DISTANCE_PC,
        "distance_flux_factor_140pc_to_147pc": DISTANCE_FACTOR,
        "aperture_radius_arcsec": APERTURE_RADIUS_ARCSEC,
        "aperture_subpixels_per_axis": APERTURE_SUBPIXELS,
        "psf_truncation_sigma": PSF_TRUNCATION_SIGMA,
        "significant_component_coeval_relative_difference_max": COEVAL_TOLERANCE,
        "component_significant_fraction": COMPONENT_SIGNIFICANT_FRACTION,
        "component_active_fraction": COMPONENT_ACTIVE_FRACTION,
        "active_component_negative_fraction_max_exclusive": NEGATIVE_FRACTION_MAX,
        "positive_edge_fraction_max_exclusive": EDGE_FRACTION_MAX,
        "psf_convolution_flux_loss_fraction_max_exclusive": CONVOLUTION_LOSS_MAX,
        "component_aperture_fraction_min_inclusive": COMPONENT_APERTURE_FRACTION_MIN,
        "component_aperture_fraction_max_inclusive": COMPONENT_APERTURE_FRACTION_MAX,
        "aperture_clip_symmetric_fraction_max_inclusive": APERTURE_CLIP_MAX,
        "canonical_plane_zero_based": 0,
        "signed_flux_used_for_scoring": True,
        "zero_clipped_flux_used_for_scoring": False,
        "full_plane_negative_fraction_role": "diagnostic_only",
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "intrinsic_separated_component_product_integrity_gates_applied": True,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_regular(path: Path, label: str) -> str:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise MeasurementContractError(
            f"Missing regular non-symlink {label}: {path}"
        )
    return sha256(path)


def psf_fwhm_arcsec(instrument: str, wavelength_um: float) -> float:
    if instrument == "NIRSpec":
        return 0.033 * wavelength_um
    if instrument == "MIRI":
        return 0.033 * wavelength_um + 0.106
    raise MeasurementContractError(f"Unsupported instrument: {instrument}")


def symmetric_fraction(first: float, second: float) -> float:
    denominator = 0.5 * (abs(first) + abs(second))
    if denominator == 0.0:
        return 0.0 if first == second else math.inf
    return abs(first - second) / denominator


def relative_difference(value: float, reference: float) -> float:
    if reference == 0.0:
        return 0.0 if value == 0.0 else math.inf
    return (value - reference) / reference


def fractional_circle_weights(
    shape: tuple[int, int], centre_x: float, centre_y: float, radius: float
) -> np.ndarray:
    """Frozen 64x64 midpoint boundary algorithm, copied without alteration."""
    ny, nx = shape
    weights = np.zeros(shape, dtype=np.float64)
    x_min = max(0, int(np.floor(centre_x - radius - 1.0)))
    x_max = min(nx - 1, int(np.ceil(centre_x + radius + 1.0)))
    y_min = max(0, int(np.floor(centre_y - radius - 1.0)))
    y_max = min(ny - 1, int(np.ceil(centre_y + radius + 1.0)))
    x_index = np.arange(x_min, x_max + 1)
    y_index = np.arange(y_min, y_max + 1)
    offsets = (
        (np.arange(APERTURE_SUBPIXELS, dtype=np.float64) + 0.5)
        / APERTURE_SUBPIXELS
        - 0.5
    )
    dx = x_index[:, None] + offsets[None, :] - centre_x
    for start in range(0, y_index.size, 8):
        selected = y_index[start : start + 8]
        dy = selected[:, None] + offsets[None, :] - centre_y
        inside = dy[:, None, :, None] ** 2 + dx[None, :, None, :] ** 2 <= radius**2
        weights[np.ix_(selected, x_index)] = inside.mean(axis=(2, 3))
    return weights


def _closure_rows(closure: Mapping[str, Any]) -> dict[int, Mapping[str, Any]]:
    rows = closure.get("planes")
    if not isinstance(rows, list) or len(rows) != len(CLOSURE_LABELS):
        raise MeasurementContractError("Closure report plane table changed")
    selected: dict[int, Mapping[str, Any]] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            raise MeasurementContractError("Closure report contains a non-object plane")
        plane = row.get("plane")
        if isinstance(plane, bool) or not isinstance(plane, int):
            raise MeasurementContractError("Closure report plane index is not an integer")
        if plane in CLOSURE_LABELS:
            if plane in selected:
                raise MeasurementContractError(f"Duplicate closure plane {plane}")
            selected[plane] = row
    if set(selected) != set(CLOSURE_LABELS):
        raise MeasurementContractError(
            "Closure report is not the exact plane 0,4,5,6,7 set"
        )
    return selected


def _validate_reported_closure_values(
    closure: Mapping[str, Any],
    *,
    sed: np.ndarray,
    image: np.ndarray,
    wavelength_um: float,
) -> tuple[dict[int, Mapping[str, Any]], dict[int, float], dict[int, float]]:
    if not (
        closure.get("schema_version") == 1
        and closure.get("passes_requested_closure_tolerance") is True
        and math.isclose(
            float(closure.get("tolerance", math.nan)),
            COEVAL_TOLERANCE,
            abs_tol=0.0,
        )
        and math.isclose(
            float(closure.get("significant_fraction", math.nan)),
            COMPONENT_SIGNIFICANT_FRACTION,
            abs_tol=0.0,
        )
        and math.isfinite(
            float(
                closure.get(
                    "max_significant_custom_image_relative_difference", math.inf
                )
            )
        )
        and float(
            closure.get("max_significant_custom_image_relative_difference", math.inf)
        )
        <= COEVAL_TOLERANCE
        and closure.get("stock_component_closure_relative") is None
        and math.isclose(
            float(closure.get("wavelength_um", math.nan)),
            wavelength_um,
            abs_tol=1.0e-6,
        )
    ):
        raise MeasurementContractError("Coeval SED/map closure report failed")

    rows = _closure_rows(closure)
    numpy_totals = {
        plane: float(np.sum(image[plane], dtype=np.float64))
        for plane in CLOSURE_LABELS
    }
    total_scale = max(abs(float(sed[0])), abs(numpy_totals[0]))
    significance_cut = COMPONENT_SIGNIFICANT_FRACTION * total_scale
    significant_errors: dict[int, float] = {}
    for plane, label in CLOSURE_LABELS.items():
        sed_value = float(sed[plane])
        map_value = numpy_totals[plane]
        significant = plane == 0 or max(abs(sed_value), abs(map_value)) > significance_cut
        coeval = relative_difference(sed_value, map_value)
        reported = rows[plane]
        if not (
            reported.get("label") == label
            and reported.get("significant") is significant
            and math.isclose(
                float(reported.get("custom_sed_flux", math.nan)),
                sed_value,
                rel_tol=1.0e-13,
                abs_tol=0.0,
            )
            and math.isclose(
                float(reported.get("image_flux", math.nan)),
                map_value,
                rel_tol=1.0e-13,
                abs_tol=0.0,
            )
            and math.isclose(
                float(reported.get("custom_relative_difference", math.nan)),
                coeval,
                rel_tol=1.0e-10,
                abs_tol=1.0e-14,
            )
        ):
            raise MeasurementContractError(
                f"Independent closure report binding failed for plane {plane}"
            )
        if significant:
            error = abs(coeval)
            if not math.isfinite(error) or error > COEVAL_TOLERANCE:
                raise MeasurementGateError(
                    f"Independent significant-plane closure failed for plane {plane}: {error}",
                    failed_checks=(f"coeval_plane_{plane}",),
                )
            significant_errors[plane] = error
    recomputed_max = max(significant_errors.values(), default=math.inf)
    if not math.isclose(
        recomputed_max,
        float(closure["max_significant_custom_image_relative_difference"]),
        rel_tol=1.0e-10,
        abs_tol=1.0e-14,
    ):
        raise MeasurementContractError(
            "Closure report maximum does not match independently reopened planes"
        )
    return rows, numpy_totals, significant_errors


def _measure_component_integrity(
    *,
    sed: np.ndarray,
    image: np.ndarray,
    closure: Mapping[str, Any],
    sigma_pixels: float,
    weights: np.ndarray,
    border_mask: np.ndarray,
) -> dict[str, Any]:
    _, numpy_totals, significant_errors = _validate_reported_closure_values(
        closure,
        sed=sed,
        image=image,
        wavelength_um=float(closure["wavelength_um"]),
    )
    recomputed_max = max(significant_errors.values(), default=math.inf)
    component_l1 = float(np.sum(np.abs(sed[4:8]), dtype=np.float64))
    component_rows: dict[str, dict[str, Any]] = {}

    for plane, name in COMPONENTS.items():
        values = np.asarray(image[plane], dtype=np.float64)
        map_total = math.fsum(float(value) for value in values.ravel(order="C"))
        map_total_numpy = numpy_totals[plane]
        sed_component = float(sed[plane])
        positive = float(np.sum(values[values > 0.0], dtype=np.float64))
        negative = float(np.sum(np.abs(values[values < 0.0]), dtype=np.float64))
        negative_fraction = (
            negative / positive
            if positive > 0.0
            else (0.0 if negative == 0.0 else float(np.finfo(float).max))
        )
        edge_fraction = (
            float(np.sum(np.clip(values, 0.0, None)[border_mask], dtype=np.float64))
            / positive
            if positive > 0.0
            else 0.0
        )
        active = (
            abs(sed_component) >= COMPONENT_ACTIVE_FRACTION * component_l1
            if component_l1 > 0.0
            else False
        )
        significant = plane in significant_errors
        coeval = relative_difference(sed_component, map_total_numpy)
        if positive > 0.0:
            convolved = gaussian_filter(
                values,
                sigma_pixels,
                mode="constant",
                cval=0.0,
                truncate=PSF_TRUNCATION_SIGMA,
            )
            convolved_total = float(np.sum(convolved, dtype=np.float64))
            map_aperture = float(np.sum(convolved * weights, dtype=np.float64))
            convolution_loss = abs(convolved_total - map_total) / max(
                abs(map_total), np.finfo(float).tiny
            )
        else:
            map_aperture = 0.0
            convolution_loss = (
                0.0 if map_total == 0.0 else float(np.finfo(float).max)
            )
        zero_map_contract = map_total != 0.0 or sed_component == 0.0
        aperture_fraction = 0.0 if map_total == 0.0 else map_aperture / map_total
        finite_map_contract = all(
            math.isfinite(value)
            for value in (
                sed_component,
                map_total,
                map_total_numpy,
                positive,
                negative,
                negative_fraction,
                edge_fraction,
                map_aperture,
                aperture_fraction,
                convolution_loss,
                coeval,
            )
        )
        aperture_fraction_contract = finite_map_contract and (
            COMPONENT_APERTURE_FRACTION_MIN
            <= aperture_fraction
            <= COMPONENT_APERTURE_FRACTION_MAX
        )
        active_morphology_contract = (not active) or (
            0.0 <= negative_fraction < NEGATIVE_FRACTION_MAX
            and 0.0 <= edge_fraction < EDGE_FRACTION_MAX
            and 0.0 <= convolution_loss < CONVOLUTION_LOSS_MAX
        )
        quality = (
            zero_map_contract
            and finite_map_contract
            and aperture_fraction_contract
            and active_morphology_contract
            and ((not significant) or abs(coeval) <= COEVAL_TOLERANCE)
        )
        if not quality:
            raise MeasurementGateError(
                f"Intrinsic separated-component gate failed: {name}",
                failed_checks=(f"component_{name}",),
            )
        component_rows[name] = {
            "plane_zero_based": plane,
            "significant_ge_1e-8_total_scale": significant,
            "active_ge_1pct_component_l1": active,
            "sed_flux_140pc_w_m2": sed_component,
            "map_total_fixed_fsum_140pc_w_m2": map_total,
            "map_total_numpy_140pc_w_m2": map_total_numpy,
            "map_min_pixel_flux_w_m2": float(np.min(values)),
            "map_max_pixel_flux_w_m2": float(np.max(values)),
            "map_positive_flux_140pc_w_m2": positive,
            "map_negative_magnitude_140pc_w_m2": negative,
            "coeval_sed_minus_map_over_map": coeval,
            "negative_flux_fraction": negative_fraction,
            "edge_positive_flux_fraction": edge_fraction,
            "convolution_flux_loss_fraction": convolution_loss,
            "psf_map_aperture_flux_140pc_w_m2": map_aperture,
            "psf_map_aperture_fraction": aperture_fraction,
            "map_pixels_finite": True,
            "zero_map_requires_zero_sed_component_pass": zero_map_contract,
            "finite_map_contract_pass": finite_map_contract,
            "aperture_fraction_contract_pass": aperture_fraction_contract,
            "active_morphology_contract_pass": active_morphology_contract,
            "quality_pass": quality,
        }

    significant_components = [
        item
        for item in component_rows.values()
        if item["significant_ge_1e-8_total_scale"]
    ]
    active_components = [
        item for item in component_rows.values() if item["active_ge_1pct_component_l1"]
    ]
    return {
        "evidence_kind": "direct_recomputation_from_coeval_sed_and_image_products",
        "intrinsic_component_product_integrity_gates_applied": True,
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "significant_fraction_of_total_scale": COMPONENT_SIGNIFICANT_FRACTION,
        "active_fraction_of_component_l1": COMPONENT_ACTIVE_FRACTION,
        "component_l1_sed_flux_140pc_w_m2": component_l1,
        "significant_component_count": len(significant_components),
        "active_component_count": len(active_components),
        "max_abs_significant_all_plane_coeval_relative_difference": recomputed_max,
        "max_abs_significant_component_coeval_relative_difference": max(
            (
                abs(float(item["coeval_sed_minus_map_over_map"]))
                for item in significant_components
            ),
            default=0.0,
        ),
        "max_active_component_negative_flux_fraction": max(
            (float(item["negative_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_edge_positive_flux_fraction": max(
            (float(item["edge_positive_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_convolution_flux_loss_fraction": max(
            (
                float(item["convolution_flux_loss_fraction"])
                for item in active_components
            ),
            default=0.0,
        ),
        "components": component_rows,
        "quality_pass": True,
    }


def measure_method2_signed_aperture_v1(
    *,
    sed_flux_140pc_w_m2: np.ndarray,
    image_planes_140pc_w_m2_pixel: np.ndarray,
    closure: Mapping[str, Any],
    wavelength_um: float,
    instrument: str,
    cdelt1_deg: float,
    cdelt2_deg: float,
    crpix1_fits: float,
    crpix2_fits: float,
    expected_cdelt_deg: float,
    expected_psf_fwhm_arcsec: float,
) -> dict[str, Any]:
    """Measure one coeval single-wavelength Method-2 SED/map realization.

    The function either returns a fully valid record or raises; it never emits
    a partially passing measurement.  Arrays are treated as read-only.
    """
    sed = np.asarray(sed_flux_140pc_w_m2, dtype=np.float64)
    image = np.asarray(image_planes_140pc_w_m2_pixel)
    if sed.shape != (8,):
        raise MeasurementContractError(f"Expected eight SED planes, got {sed.shape}")
    if image.ndim != 3 or image.shape[0] != 8:
        raise MeasurementContractError(
            f"Expected image planes shaped (8,ny,nx), got {image.shape}"
        )
    if not np.all(np.isfinite(sed)) or not np.all(np.isfinite(image)):
        raise MeasurementContractError("SED/map contains non-finite values")
    if not math.isfinite(wavelength_um) or wavelength_um <= 0.0:
        raise MeasurementContractError("Wavelength must be finite and positive")
    if not isinstance(closure, Mapping):
        raise MeasurementContractError("Closure report root is not an object")

    ny, nx = int(image.shape[1]), int(image.shape[2])
    centre_x = float(crpix1_fits) - 1.0
    centre_y = float(crpix2_fits) - 1.0
    if not (
        math.isclose(abs(float(cdelt1_deg)), expected_cdelt_deg, rel_tol=5.0e-6)
        and math.isclose(abs(float(cdelt2_deg)), expected_cdelt_deg, rel_tol=5.0e-6)
        and math.isclose(centre_x, (nx - 1) / 2.0, abs_tol=1.0e-6)
        and math.isclose(centre_y, (ny - 1) / 2.0, abs_tol=1.0e-6)
    ):
        raise MeasurementContractError("Map scale/CRPIX centre changed")

    target_pixel_arcsec = (
        abs(float(cdelt1_deg))
        * 3600.0
        * MODEL_DISTANCE_PC
        / TARGET_DISTANCE_PC
    )
    expected_target_pixel = (
        expected_cdelt_deg
        * 3600.0
        * MODEL_DISTANCE_PC
        / TARGET_DISTANCE_PC
    )
    radius_pixels = APERTURE_RADIUS_ARCSEC / target_pixel_arcsec
    fwhm = psf_fwhm_arcsec(instrument, wavelength_um)
    sigma_pixels = fwhm / 2.354820045 / target_pixel_arcsec
    half_field = min(
        (centre_x + 0.5) * target_pixel_arcsec,
        (nx - 0.5 - centre_x) * target_pixel_arcsec,
        (centre_y + 0.5) * target_pixel_arcsec,
        (ny - 0.5 - centre_y) * target_pixel_arcsec,
    )
    support_minimum = (
        APERTURE_RADIUS_ARCSEC
        + PSF_TRUNCATION_SIGMA * fwhm / 2.354820045
    )
    if not (
        math.isclose(target_pixel_arcsec, expected_target_pixel, rel_tol=5.0e-6)
        and math.isclose(expected_psf_fwhm_arcsec, fwhm, abs_tol=1.0e-14)
        and half_field >= support_minimum
    ):
        raise MeasurementGateError(
            "PSF/aperture support geometry failed",
            failed_checks=("aperture_plus_6sigma_support",),
        )

    closure_rows, _, _ = _validate_reported_closure_values(
        closure,
        sed=sed,
        image=image,
        wavelength_um=wavelength_um,
    )
    plane0 = np.asarray(image[0], dtype=np.float64)
    sed_plane0 = float(sed[0])
    map_total = math.fsum(float(value) for value in plane0.ravel(order="C"))
    map_total_numpy = float(np.sum(plane0, dtype=np.float64))
    coeval_relative = relative_difference(sed_plane0, map_total)
    reported_relative = float(closure_rows[0]["custom_relative_difference"])
    if not (
        abs(coeval_relative) <= COEVAL_TOLERANCE
        and math.isclose(
            float(closure_rows[0]["custom_sed_flux"]),
            sed_plane0,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        )
        and math.isclose(
            float(closure_rows[0]["image_flux"]),
            map_total_numpy,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        )
        and math.isclose(
            reported_relative,
            relative_difference(sed_plane0, map_total_numpy),
            rel_tol=1.0e-10,
            abs_tol=1.0e-14,
        )
    ):
        raise MeasurementGateError(
            "Independently recomputed plane-0 closure failed",
            failed_checks=("coeval_plane0_abs_relative_difference_le_1e-5",),
        )

    positive = float(np.sum(plane0[plane0 > 0.0], dtype=np.float64))
    negative = float(np.sum(np.abs(plane0[plane0 < 0.0]), dtype=np.float64))
    negative_fraction = negative / positive if positive > 0.0 else math.inf
    border = max(1, min(nx, ny) // 100)
    border_mask = np.ones((ny, nx), dtype=bool)
    border_mask[border:-border, border:-border] = False
    weights = fractional_circle_weights((ny, nx), centre_x, centre_y, radius_pixels)
    component_integrity = _measure_component_integrity(
        sed=sed,
        image=image,
        closure=closure,
        sigma_pixels=sigma_pixels,
        weights=weights,
        border_mask=border_mask,
    )
    edge_fraction = (
        float(np.sum(np.clip(plane0, 0.0, None)[border_mask], dtype=np.float64))
        / positive
        if positive > 0.0
        else math.inf
    )
    signed_convolved = gaussian_filter(
        plane0,
        sigma_pixels,
        mode="constant",
        cval=0.0,
        truncate=PSF_TRUNCATION_SIGMA,
    )
    signed_convolved_total = float(np.sum(signed_convolved, dtype=np.float64))
    signed_aperture_140 = float(np.sum(signed_convolved * weights, dtype=np.float64))
    convolution_loss = abs(signed_convolved_total - map_total) / max(
        abs(map_total), np.finfo(float).tiny
    )
    if negative == 0.0:
        clipped_aperture_140 = signed_aperture_140
    else:
        clipped_convolved = gaussian_filter(
            np.clip(plane0, 0.0, None),
            sigma_pixels,
            mode="constant",
            cval=0.0,
            truncate=PSF_TRUNCATION_SIGMA,
        )
        clipped_aperture_140 = float(
            np.sum(clipped_convolved * weights, dtype=np.float64)
        )
    dclip = symmetric_fraction(signed_aperture_140, clipped_aperture_140)
    checks = {
        "finite_positive_plane0_flux": (
            math.isfinite(map_total) and map_total > 0.0 and positive > 0.0
        ),
        "coeval_plane0_abs_relative_difference_le_1e-5": (
            math.isfinite(coeval_relative)
            and abs(coeval_relative) <= COEVAL_TOLERANCE
        ),
        "edge_positive_fraction_lt_1e-3": (
            math.isfinite(edge_fraction)
            and 0.0 <= edge_fraction < EDGE_FRACTION_MAX
        ),
        "convolution_flux_loss_lt_1e-4": (
            math.isfinite(convolution_loss)
            and 0.0 <= convolution_loss < CONVOLUTION_LOSS_MAX
        ),
        "finite_positive_signed_aperture": (
            math.isfinite(signed_aperture_140) and signed_aperture_140 > 0.0
        ),
        "direct_aperture_clip_symmetric_fraction_le_1e-3": (
            math.isfinite(dclip) and 0.0 <= dclip <= APERTURE_CLIP_MAX
        ),
        "aperture_plus_6sigma_support": half_field >= support_minimum,
    }
    if not all(checks.values()):
        failed = tuple(sorted(key for key, value in checks.items() if not value))
        raise MeasurementGateError(
            f"Canonical quality gate failed: {list(failed)}",
            failed_checks=failed,
        )

    prediction_147 = signed_aperture_140 * DISTANCE_FACTOR
    return {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "valid": True,
        "canonical_estimator": "signed_direct_coeval_method2_polarized_total_i_plane0",
        "policy": fixed_policy_v1(),
        "wavelength_um": float(wavelength_um),
        "instrument": instrument,
        "signed_aperture_flux_140pc_w_m2": signed_aperture_140,
        "signed_aperture_flux_147pc_w_m2": prediction_147,
        "zero_clipped_aperture_flux_140pc_w_m2": clipped_aperture_140,
        "aperture_clip_symmetric_fraction": dclip,
        "aperture_clip_evidence_kind": (
            "direct_post_psf_nominal_1arcsec_symmetric_fraction"
        ),
        "full_plane_negative_fraction_diagnostic": negative_fraction,
        "plane0_edge_positive_flux_fraction": edge_fraction,
        "plane0_convolution_flux_loss_fraction": convolution_loss,
        "coeval_plane0_relative_difference": coeval_relative,
        "coeval_max_significant_relative_difference": float(
            closure["max_significant_custom_image_relative_difference"]
        ),
        "intrinsic_component_integrity": component_integrity,
        "geometry": {
            "nx": nx,
            "ny": ny,
            "cdelt1_deg": float(cdelt1_deg),
            "cdelt2_deg": float(cdelt2_deg),
            "crpix1_fits": float(crpix1_fits),
            "crpix2_fits": float(crpix2_fits),
            "centre_x_zero_based": centre_x,
            "centre_y_zero_based": centre_y,
            "target_pixel_arcsec": target_pixel_arcsec,
            "aperture_radius_pixels": radius_pixels,
            "psf_fwhm_arcsec": fwhm,
            "psf_sigma_pixels": sigma_pixels,
            "available_half_field_arcsec": half_field,
            "required_half_field_arcsec": support_minimum,
        },
        "quality_checks": checks,
    }


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementContractError(f"Cannot read {label}: {path}: {error}") from error
    if not isinstance(value, dict):
        raise MeasurementContractError(f"{label} root is not an object")
    return value


def _validate_external_closure_identity(
    closure: Mapping[str, Any],
    *,
    sed_path: Path,
    image_path: Path,
    sed_sha256: str,
    image_sha256: str,
) -> None:
    inputs = closure.get("inputs")
    if not isinstance(inputs, Mapping):
        raise MeasurementContractError("Closure report input table is absent")
    try:
        reported_sed = Path(str(inputs.get("custom_sed"))).resolve()
        reported_image = Path(str(inputs.get("image"))).resolve()
    except (OSError, RuntimeError) as error:
        raise MeasurementContractError("Closure input path cannot be resolved") from error
    if not (
        reported_sed == sed_path.resolve()
        and reported_image == image_path.resolve()
        and inputs.get("custom_sed_sha256") == sed_sha256
        and inputs.get("image_sha256") == image_sha256
        and inputs.get("stock_sed") is None
        and inputs.get("stock_sed_sha256") is None
    ):
        raise MeasurementContractError("Closure report product identity changed")


def measure_method2_fits_task_v1(
    *,
    sed_path: Path,
    image_path: Path,
    closure_path: Path,
    wavelength_um: float,
    instrument: str,
    expected_nx: int,
    expected_ny: int,
    expected_cdelt_deg: float,
    expected_psf_fwhm_arcsec: float,
    task_identity: Mapping[str, Any],
) -> dict[str, Any]:
    """Read, bind, measure, and return one fail-closed per-RT-task record."""
    sed_path = Path(sed_path)
    image_path = Path(image_path)
    closure_path = Path(closure_path)
    sed_sha = require_regular(sed_path, "coeval Method-2 SED")
    image_sha = require_regular(image_path, "coeval Method-2 image")
    closure_sha = require_regular(closure_path, "coeval closure report")
    closure = _load_json_object(closure_path, "coeval closure report")
    _validate_external_closure_identity(
        closure,
        sed_path=sed_path,
        image_path=image_path,
        sed_sha256=sed_sha,
        image_sha256=image_sha,
    )

    with fits.open(sed_path, memmap=False) as hdus:
        sed_raw = np.array(hdus[0].data, dtype=np.float64, copy=True)
        sed_shape = tuple(sed_raw.shape)
        sed_header = hdus[0].header.copy()
        sed_bitpix = int(sed_header.get("BITPIX", 0))
        sed_unit = str(sed_header.get("BUNIT", "")).lower().replace(" ", "")
        wave = np.asarray(hdus[1].data, dtype=np.float64).ravel()
    with fits.open(image_path, memmap=False) as hdus:
        image_raw = np.array(hdus[0].data, copy=True)
        image_shape = tuple(image_raw.shape)
        image_bitpix = int(hdus[0].header.get("BITPIX", 0))
        image_header = hdus[0].header.copy()
        image_unit = str(image_header.get("BUNIT", "")).lower().replace(" ", "")

    expected_component_tokens = {
        5: ("direct", "star"),
        6: ("scattered", "star"),
        7: ("direct", "thermal"),
        8: ("scattered", "thermal"),
    }
    component_semantics_pass = all(
        all(
            token in str(image_header.get(f"FLUX_{header_plane}", "")).lower()
            for token in tokens
        )
        for header_plane, tokens in expected_component_tokens.items()
    )
    if not (
        sed_shape == (8, 1, 1, 1)
        and sed_bitpix == -64
        and wave.shape == (1,)
        and np.all(np.isfinite(sed_raw))
        and image_shape == (8, 1, 1, expected_ny, expected_nx)
        and image_bitpix == -32
        and np.all(np.isfinite(image_raw))
        and sed_unit in {"w.m-2", "w/m2"}
        and image_unit == "w.m-2.pixel-1"
        and "total flux" in str(image_header.get("FLUX_1", "")).lower()
        and component_semantics_pass
        and math.isclose(float(wave[0]), wavelength_um, rel_tol=1.0e-6, abs_tol=1.0e-8)
        and math.isclose(
            float(image_header.get("WAVE", math.nan)),
            float(np.float32(wavelength_um)),
            rel_tol=1.0e-6,
            abs_tol=1.0e-8,
        )
    ):
        raise MeasurementContractError(
            "FITS shape/precision/unit/wavelength gate failed"
        )
    for header, label in ((sed_header, "SED"), (image_header, "map")):
        if any(key in header for key in ("BSCALE", "BZERO")):
            raise MeasurementContractError(f"{label} uses forbidden FITS scaling")
    forbidden_wcs = (
        "PC1_1",
        "PC1_2",
        "PC2_1",
        "PC2_2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
        "CROTA1",
        "CROTA2",
    )
    if any(key in image_header for key in forbidden_wcs):
        raise MeasurementContractError("Map WCS matrix/rotation overrides are forbidden")

    measurement = measure_method2_signed_aperture_v1(
        sed_flux_140pc_w_m2=sed_raw[:, 0, 0, 0],
        image_planes_140pc_w_m2_pixel=image_raw[:, 0, 0],
        closure=closure,
        wavelength_um=wavelength_um,
        instrument=instrument,
        cdelt1_deg=float(image_header.get("CDELT1", math.nan)),
        cdelt2_deg=float(image_header.get("CDELT2", math.nan)),
        crpix1_fits=float(image_header.get("CRPIX1", math.nan)),
        crpix2_fits=float(image_header.get("CRPIX2", math.nan)),
        expected_cdelt_deg=expected_cdelt_deg,
        expected_psf_fwhm_arcsec=expected_psf_fwhm_arcsec,
    )
    # Bind again after all reads and expensive calculations.
    if not (
        require_regular(sed_path, "post-read coeval Method-2 SED") == sed_sha
        and require_regular(image_path, "post-read coeval Method-2 image") == image_sha
        and require_regular(closure_path, "post-read coeval closure report") == closure_sha
    ):
        raise MeasurementContractError("A measurement input changed while being read")
    record = {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "valid": True,
        "task_identity": dict(task_identity),
        "inputs": {
            "sed_rt": {"path": str(sed_path.resolve()), "sha256": sed_sha},
            "image_rt": {"path": str(image_path.resolve()), "sha256": image_sha},
            "closure": {"path": str(closure_path.resolve()), "sha256": closure_sha},
        },
        "measurement": measurement,
    }
    validate_task_measurement_record_v1(record)
    return record


def _require_finite_json_tree(value: Any, path: str = "$") -> None:
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise MeasurementContractError(f"Non-finite JSON number at {path}")
        return
    if isinstance(value, list):
        for index, child in enumerate(value):
            _require_finite_json_tree(child, f"{path}[{index}]")
        return
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                raise MeasurementContractError(f"Non-string JSON key at {path}")
            _require_finite_json_tree(child, f"{path}.{key}")
        return
    raise MeasurementContractError(
        f"Non-JSON value at {path}: {type(value).__name__}"
    )


def validate_measurement_v1(measurement: Mapping[str, Any]) -> None:
    """Fail closed if a decoded pure-kernel measurement is incomplete/tampered."""
    if not isinstance(measurement, Mapping):
        raise MeasurementContractError("Measurement root is not an object")
    if not (
        measurement.get("schema_version") == SCHEMA_VERSION
        and measurement.get("kernel_id") == KERNEL_ID
        and measurement.get("valid") is True
        and measurement.get("canonical_estimator")
        == "signed_direct_coeval_method2_polarized_total_i_plane0"
        and measurement.get("policy") == fixed_policy_v1()
    ):
        raise MeasurementContractError("Measurement identity/policy changed")
    checks = measurement.get("quality_checks")
    if not isinstance(checks, Mapping) or set(checks) != {
        "finite_positive_plane0_flux",
        "coeval_plane0_abs_relative_difference_le_1e-5",
        "edge_positive_fraction_lt_1e-3",
        "convolution_flux_loss_lt_1e-4",
        "finite_positive_signed_aperture",
        "direct_aperture_clip_symmetric_fraction_le_1e-3",
        "aperture_plus_6sigma_support",
    } or any(value is not True for value in checks.values()):
        raise MeasurementContractError("Measurement quality-check table is not all-pass")
    signed140 = float(measurement.get("signed_aperture_flux_140pc_w_m2", math.nan))
    signed147 = float(measurement.get("signed_aperture_flux_147pc_w_m2", math.nan))
    clipped = float(
        measurement.get("zero_clipped_aperture_flux_140pc_w_m2", math.nan)
    )
    dclip = float(measurement.get("aperture_clip_symmetric_fraction", math.nan))
    if not (
        math.isfinite(signed140)
        and signed140 > 0.0
        and math.isclose(
            signed147, signed140 * DISTANCE_FACTOR, rel_tol=1.0e-14, abs_tol=0.0
        )
        and math.isfinite(clipped)
        and math.isclose(
            dclip, symmetric_fraction(signed140, clipped), rel_tol=1.0e-14, abs_tol=0.0
        )
        and 0.0 <= dclip <= APERTURE_CLIP_MAX
        and abs(float(measurement.get("coeval_plane0_relative_difference", math.inf)))
        <= COEVAL_TOLERANCE
        and 0.0
        <= float(measurement.get("plane0_edge_positive_flux_fraction", math.inf))
        < EDGE_FRACTION_MAX
        and 0.0
        <= float(
            measurement.get("plane0_convolution_flux_loss_fraction", math.inf)
        )
        < CONVOLUTION_LOSS_MAX
    ):
        raise MeasurementContractError("Measurement scalar gates or scaling changed")
    integrity = measurement.get("intrinsic_component_integrity")
    if not (
        isinstance(integrity, Mapping)
        and integrity.get("quality_pass") is True
        and integrity.get("intrinsic_component_product_integrity_gates_applied")
        is True
        and integrity.get("plane0_vs_sum_planes4to7_gate_applied") is False
        and integrity.get("aggregate_componentwise_flux_gate_applied") is False
        and isinstance(integrity.get("components"), Mapping)
        and set(integrity["components"]) == set(COMPONENTS.values())
        and all(
            isinstance(row, Mapping) and row.get("quality_pass") is True
            for row in integrity["components"].values()
        )
    ):
        raise MeasurementContractError("Intrinsic component integrity changed")
    geometry = measurement.get("geometry")
    if not (
        isinstance(geometry, Mapping)
        and int(geometry.get("nx", 0)) > 0
        and int(geometry.get("ny", 0)) > 0
        and math.isclose(
            float(geometry.get("crpix1_fits", math.nan)) - 1.0,
            float(geometry.get("centre_x_zero_based", math.nan)),
            abs_tol=0.0,
        )
        and math.isclose(
            float(geometry.get("crpix2_fits", math.nan)) - 1.0,
            float(geometry.get("centre_y_zero_based", math.nan)),
            abs_tol=0.0,
        )
        and float(geometry.get("available_half_field_arcsec", -math.inf))
        >= float(geometry.get("required_half_field_arcsec", math.inf))
    ):
        raise MeasurementContractError("Measurement geometry changed")
    _require_finite_json_tree(measurement)


def validate_task_measurement_record_v1(record: Mapping[str, Any]) -> None:
    """Validate the envelope used as one RT task's durable measurement."""
    if not isinstance(record, Mapping):
        raise MeasurementContractError("Task measurement root is not an object")
    if not (
        record.get("schema_version") == SCHEMA_VERSION
        and record.get("kernel_id") == KERNEL_ID
        and record.get("valid") is True
        and isinstance(record.get("task_identity"), Mapping)
        and bool(record["task_identity"])
    ):
        raise MeasurementContractError("Task measurement identity changed")
    inputs = record.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != {
        "sed_rt",
        "image_rt",
        "closure",
    }:
        raise MeasurementContractError("Task measurement input table changed")
    for name, item in inputs.items():
        if not (
            isinstance(item, Mapping)
            and isinstance(item.get("path"), str)
            and item.get("path")
            and isinstance(item.get("sha256"), str)
            and len(item["sha256"]) == 64
            and all(character in "0123456789abcdef" for character in item["sha256"])
        ):
            raise MeasurementContractError(f"Malformed task measurement input: {name}")
    validate_measurement_v1(record.get("measurement"))
    _require_finite_json_tree(record)


def validate_task_measurement_input_files_v1(record: Mapping[str, Any]) -> None:
    """Reopen and hash-bind every product referenced by a task record."""
    validate_task_measurement_record_v1(record)
    inputs = record["inputs"]
    for name in ("sed_rt", "image_rt", "closure"):
        item = inputs[name]
        path = Path(str(item["path"]))
        actual = require_regular(path, f"task measurement input {name}")
        if actual != item["sha256"]:
            raise MeasurementContractError(
                f"Task measurement input {name} SHA-256 changed: "
                f"{actual}; expected {item['sha256']}"
            )
    closure = _load_json_object(
        Path(str(inputs["closure"]["path"])), "task measurement closure report"
    )
    _validate_external_closure_identity(
        closure,
        sed_path=Path(str(inputs["sed_rt"]["path"])),
        image_path=Path(str(inputs["image_rt"]["path"])),
        sed_sha256=str(inputs["sed_rt"]["sha256"]),
        image_sha256=str(inputs["image_rt"]["sha256"]),
    )


def atomic_write_task_measurement_json_v1(
    path: Path,
    record: Mapping[str, Any],
    *,
    replace_existing: bool = False,
) -> str:
    """Validate, serialize strictly, and atomically publish one task record."""
    validate_task_measurement_input_files_v1(record)
    path = Path(path)
    if not replace_existing and (path.exists() or path.is_symlink()):
        raise MeasurementContractError(f"Refusing pre-existing measurement JSON: {path}")
    try:
        payload = json.dumps(
            record,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        ) + "\n"
        decoded = json.loads(payload)
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise MeasurementContractError(f"Measurement is not strict JSON: {error}") from error
    validate_task_measurement_record_v1(decoded)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        try:
            directory_fd = os.open(path.parent, os.O_RDONLY)
        except OSError:
            directory_fd = None
        if directory_fd is not None:
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return sha256(path)


def load_task_measurement_json_v1(
    path: Path,
    *,
    expected_sha256: str | None = None,
    verify_input_files: bool = True,
) -> dict[str, Any]:
    """Read a durable record, optionally hash-pin it, and validate fail closed."""
    path = Path(path)
    actual_sha = require_regular(path, "task measurement JSON")
    if expected_sha256 is not None and actual_sha != expected_sha256:
        raise MeasurementContractError(
            f"Task measurement SHA-256 changed: {actual_sha}; expected {expected_sha256}"
        )
    record = _load_json_object(path, "task measurement JSON")
    validate_task_measurement_record_v1(record)
    if verify_input_files:
        validate_task_measurement_input_files_v1(record)
    if require_regular(path, "post-read task measurement JSON") != actual_sha:
        raise MeasurementContractError("Task measurement JSON changed while being read")
    return record
