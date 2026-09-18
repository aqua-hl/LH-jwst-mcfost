#!/usr/bin/env python3
"""Compare the v1 measurement kernel with the frozen legacy implementation.

The utility uses real local MCFOST image morphology where products are
available.  Those older convergence products do not contain the paired strict
SED/closure bundle, so a self-closing SED and closure report are constructed
from each *unchanged real image* solely for numerical old/new parity.  This is
not a scientific remeasurement of those convergence runs.

The strict support-expanded and task-024 recovery products are cluster-only in
the current checkout.  Their absence is recorded explicitly rather than being
silently replaced by synthetic evidence.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

import analyze_strict_continuum_refinement_v1 as legacy
import strict_continuum_measurement_kernel_v1 as kernel


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "output/strict_continuum_measurement_kernel_v1_parity.json"

LOCAL_CASES = (
    {
        "case": "full_n2401",
        "path": HERE
        / "runs/aperture_spatial_convergence_stage4b/tasks/"
        "t001_I_w000_seed41001_z1_n2401/image/seed=41001/"
        "data_1.0046728324804879/RT.fits.gz",
    },
    {
        "case": "zoom_n1201",
        "path": HERE
        / "runs/aperture_spatial_convergence_stage4/tasks/"
        "t007_F_w002_seed41001_z3p9933471945199667_n1201/image/seed=41001/"
        "data_12.486673413753987/RT.fits.gz",
    },
    {
        "case": "support_n2401",
        "path": HERE
        / "runs/strict_continuum_refinement_v1_support_v1/representative/"
        "paired/seed=41001/data_th/RT.fits.gz",
    },
    {
        "case": "recovered_negative_task024_seed41002",
        "path": HERE
        / "runs/strict_continuum_3d_q_v1_task024_recovery_v1/seed41002_128k/"
        "paired/seed=41002/data_th/RT.fits.gz",
    },
)


def closure_for_arrays(
    sed: np.ndarray, image: np.ndarray, wavelength_um: float
) -> dict[str, Any]:
    totals = {
        plane: float(np.sum(image[plane], dtype=np.float64))
        for plane in kernel.CLOSURE_LABELS
    }
    total_scale = max(abs(float(sed[0])), abs(totals[0]))
    significant_cut = kernel.COMPONENT_SIGNIFICANT_FRACTION * total_scale
    rows: list[dict[str, Any]] = []
    errors: list[float] = []
    for plane, label in kernel.CLOSURE_LABELS.items():
        sed_value = float(sed[plane])
        map_value = totals[plane]
        significant = (
            plane == 0
            or max(abs(sed_value), abs(map_value)) > significant_cut
        )
        relative = kernel.relative_difference(sed_value, map_value)
        if significant:
            errors.append(abs(relative))
        rows.append(
            {
                "plane": plane,
                "label": label,
                "significant": significant,
                "image_flux": map_value,
                "custom_sed_flux": sed_value,
                "custom_relative_difference": relative,
            }
        )
    maximum = max(errors, default=math.inf)
    return {
        "schema_version": 1,
        "wavelength_um": float(wavelength_um),
        "tolerance": kernel.COEVAL_TOLERANCE,
        "significant_fraction": kernel.COMPONENT_SIGNIFICANT_FRACTION,
        "stock_component_closure_relative": None,
        "max_significant_custom_image_relative_difference": maximum,
        "passes_requested_closure_tolerance": maximum <= kernel.COEVAL_TOLERANCE,
        "planes": rows,
    }


def legacy_measure(
    *,
    sed: np.ndarray,
    image: np.ndarray,
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
    """Invoke legacy helpers and reproduce its plane-0 block verbatim."""
    ny, nx = image.shape[1:]
    centre_x = crpix1_fits - 1.0
    centre_y = crpix2_fits - 1.0
    if not (
        math.isclose(abs(cdelt1_deg), expected_cdelt_deg, rel_tol=5.0e-6)
        and math.isclose(abs(cdelt2_deg), expected_cdelt_deg, rel_tol=5.0e-6)
        and math.isclose(centre_x, (nx - 1) / 2.0, abs_tol=1.0e-6)
        and math.isclose(centre_y, (ny - 1) / 2.0, abs_tol=1.0e-6)
    ):
        raise ValueError("legacy geometry gate failed")
    target_pixel_arcsec = (
        abs(cdelt1_deg)
        * 3600.0
        * legacy.MODEL_DISTANCE_PC
        / legacy.TARGET_DISTANCE_PC
    )
    expected_target_pixel = (
        expected_cdelt_deg
        * 3600.0
        * legacy.MODEL_DISTANCE_PC
        / legacy.TARGET_DISTANCE_PC
    )
    radius_pixels = legacy.APERTURE_RADIUS_ARCSEC / target_pixel_arcsec
    fwhm = legacy.psf_fwhm_arcsec(instrument, wavelength_um)
    sigma_pixels = fwhm / 2.354820045 / target_pixel_arcsec
    half_field = min(
        (centre_x + 0.5) * target_pixel_arcsec,
        (nx - 0.5 - centre_x) * target_pixel_arcsec,
        (centre_y + 0.5) * target_pixel_arcsec,
        (ny - 0.5 - centre_y) * target_pixel_arcsec,
    )
    support_minimum = (
        legacy.APERTURE_RADIUS_ARCSEC
        + legacy.PSF_TRUNCATION_SIGMA * fwhm / 2.354820045
    )
    if not (
        math.isclose(target_pixel_arcsec, expected_target_pixel, rel_tol=5.0e-6)
        and math.isclose(expected_psf_fwhm_arcsec, fwhm, abs_tol=1.0e-14)
        and half_field >= support_minimum
    ):
        raise ValueError("legacy support gate failed")

    closure_rows = legacy._closure_rows(dict(closure))
    plane0 = np.asarray(image[0], dtype=np.float64)
    sed_plane0 = float(sed[0])
    map_total = math.fsum(float(value) for value in plane0.ravel(order="C"))
    map_total_numpy = float(np.sum(plane0, dtype=np.float64))
    coeval_relative = legacy.relative_difference(sed_plane0, map_total)
    reported_relative = float(closure_rows[0]["custom_relative_difference"])
    if not (
        abs(coeval_relative) <= legacy.COEVAL_TOLERANCE
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
            legacy.relative_difference(sed_plane0, map_total_numpy),
            rel_tol=1.0e-10,
            abs_tol=1.0e-14,
        )
    ):
        raise ValueError("legacy plane-0 closure gate failed")

    positive = float(np.sum(plane0[plane0 > 0.0], dtype=np.float64))
    negative = float(np.sum(np.abs(plane0[plane0 < 0.0]), dtype=np.float64))
    negative_fraction = negative / positive if positive > 0.0 else math.inf
    border = max(1, min(nx, ny) // 100)
    border_mask = np.ones((ny, nx), dtype=bool)
    border_mask[border:-border, border:-border] = False
    weights = legacy.fractional_circle_weights(
        (ny, nx), centre_x, centre_y, radius_pixels
    )
    component_integrity = legacy.measure_fresh_component_integrity(
        sed=sed,
        image=image[:, None, None],
        closure=dict(closure),
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
        truncate=legacy.PSF_TRUNCATION_SIGMA,
    )
    signed_convolved_total = float(np.sum(signed_convolved, dtype=np.float64))
    signed_aperture_140 = float(
        np.sum(signed_convolved * weights, dtype=np.float64)
    )
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
            truncate=legacy.PSF_TRUNCATION_SIGMA,
        )
        clipped_aperture_140 = float(
            np.sum(clipped_convolved * weights, dtype=np.float64)
        )
    dclip = legacy.symmetric_fraction(signed_aperture_140, clipped_aperture_140)
    checks = {
        "finite_positive_plane0_flux": (
            math.isfinite(map_total) and map_total > 0.0 and positive > 0.0
        ),
        "coeval_plane0_abs_relative_difference_le_1e-5": (
            math.isfinite(coeval_relative)
            and abs(coeval_relative) <= legacy.COEVAL_TOLERANCE
        ),
        "edge_positive_fraction_lt_1e-3": (
            math.isfinite(edge_fraction)
            and 0.0 <= edge_fraction < legacy.EDGE_FRACTION_MAX
        ),
        "convolution_flux_loss_lt_1e-4": (
            math.isfinite(convolution_loss)
            and 0.0 <= convolution_loss < legacy.CONVOLUTION_LOSS_MAX
        ),
        "finite_positive_signed_aperture": (
            math.isfinite(signed_aperture_140) and signed_aperture_140 > 0.0
        ),
        "direct_aperture_clip_symmetric_fraction_le_1e-3": (
            math.isfinite(dclip) and 0.0 <= dclip <= legacy.APERTURE_CLIP_MAX
        ),
        "aperture_plus_6sigma_support": half_field >= support_minimum,
    }
    if not all(checks.values()):
        failed = sorted(key for key, value in checks.items() if not value)
        raise ValueError(f"legacy canonical quality gate failed: {failed}")
    return {
        "signed_aperture_flux_140pc_w_m2": signed_aperture_140,
        "signed_aperture_flux_147pc_w_m2": (
            signed_aperture_140 * legacy.DISTANCE_FACTOR
        ),
        "zero_clipped_aperture_flux_140pc_w_m2": clipped_aperture_140,
        "aperture_clip_symmetric_fraction": dclip,
        "full_plane_negative_fraction_diagnostic": negative_fraction,
        "plane0_edge_positive_flux_fraction": edge_fraction,
        "plane0_convolution_flux_loss_fraction": convolution_loss,
        "coeval_plane0_relative_difference": coeval_relative,
        "intrinsic_component_integrity": component_integrity,
        "geometry": {
            "target_pixel_arcsec": target_pixel_arcsec,
            "aperture_radius_pixels": radius_pixels,
            "psf_fwhm_arcsec": fwhm,
            "available_half_field_arcsec": half_field,
        },
        "quality_checks": checks,
    }


def _capture(action: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    try:
        return {"outcome": "pass", "value": action()}
    except (ValueError, FloatingPointError) as error:
        return {
            "outcome": "reject",
            "exception": type(error).__name__,
            "message": str(error),
        }


def _compare_numbers(old: Any, new: Any, path: str, differences: list[str]) -> None:
    if isinstance(old, bool) or isinstance(new, bool):
        if old is not new:
            differences.append(f"{path}: {old!r} != {new!r}")
        return
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        if not math.isclose(float(old), float(new), rel_tol=1.0e-13, abs_tol=1.0e-15):
            differences.append(f"{path}: {old!r} != {new!r}")
        return
    if isinstance(old, Mapping) and isinstance(new, Mapping):
        for key in old:
            if key not in new:
                differences.append(f"{path}.{key}: absent from new")
            else:
                _compare_numbers(old[key], new[key], f"{path}.{key}", differences)
        return
    if old != new:
        differences.append(f"{path}: {old!r} != {new!r}")


def compare_case(case: str, path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {
            "case": case,
            "status": "not_available_locally",
            "path": str(path),
            "cluster_parity_remaining": True,
        }
    with fits.open(path, memmap=False) as hdus:
        raw = np.array(hdus[0].data, copy=True)
        header = hdus[0].header.copy()
    if raw.ndim != 5 or raw.shape[:3] != (8, 1, 1):
        raise ValueError(f"Unexpected representative map shape: {raw.shape}")
    image = raw[:, 0, 0]
    sed = np.sum(image, axis=(-2, -1), dtype=np.float64)
    wavelength = float(header["WAVE"])
    instrument = "NIRSpec" if wavelength < 5.5 else "MIRI"
    expected_cdelt = abs(float(header["CDELT1"]))
    expected_fwhm = kernel.psf_fwhm_arcsec(instrument, wavelength)
    closure = closure_for_arrays(sed, image, wavelength)
    arguments = {
        "sed": sed,
        "image": image,
        "closure": closure,
        "wavelength_um": wavelength,
        "instrument": instrument,
        "cdelt1_deg": float(header["CDELT1"]),
        "cdelt2_deg": float(header["CDELT2"]),
        "crpix1_fits": float(header["CRPIX1"]),
        "crpix2_fits": float(header["CRPIX2"]),
        "expected_cdelt_deg": expected_cdelt,
        "expected_psf_fwhm_arcsec": expected_fwhm,
    }
    old = _capture(lambda: legacy_measure(**arguments))
    new_arguments = dict(arguments)
    new_arguments["sed_flux_140pc_w_m2"] = new_arguments.pop("sed")
    new_arguments["image_planes_140pc_w_m2_pixel"] = new_arguments.pop("image")
    new = _capture(
        lambda: kernel.measure_method2_signed_aperture_v1(**new_arguments)
    )
    differences: list[str] = []
    if old["outcome"] != new["outcome"]:
        differences.append(
            f"outcome: legacy={old['outcome']} new={new['outcome']}"
        )
    elif old["outcome"] == "pass":
        _compare_numbers(old["value"], new["value"], "$", differences)
    return {
        "case": case,
        "status": "parity_pass" if not differences else "parity_fail",
        "path": str(path),
        "image_sha256": kernel.sha256(path),
        "image_shape_yx": [int(image.shape[1]), int(image.shape[2])],
        "wavelength_um": wavelength,
        "evidence_kind": "unchanged_real_image_morphology_with_self_closing_sed_for_numerical_parity_only",
        "legacy_outcome": old["outcome"],
        "new_outcome": new["outcome"],
        "differences": differences,
        "cluster_parity_remaining": False,
    }


def _synthetic_arguments(*, plane0_negative: float = 0.0) -> dict[str, Any]:
    ny = nx = 161
    image = np.zeros((8, ny, nx), dtype=np.float32)
    centre = nx // 2
    image[0, centre, centre] = 1.0
    if plane0_negative:
        image[0, centre, centre + 20] = -plane0_negative
    image[4:8, centre, centre] = 0.25
    sed = np.sum(image, axis=(-2, -1), dtype=np.float64)
    wavelength = 1.5
    target_pixel_arcsec = 0.02
    cdelt = (
        target_pixel_arcsec
        / 3600.0
        * kernel.TARGET_DISTANCE_PC
        / kernel.MODEL_DISTANCE_PC
    )
    closure = closure_for_arrays(sed, image, wavelength)
    return {
        "sed": sed,
        "image": image,
        "closure": closure,
        "wavelength_um": wavelength,
        "instrument": "NIRSpec",
        "cdelt1_deg": -cdelt,
        "cdelt2_deg": cdelt,
        "crpix1_fits": 81.0,
        "crpix2_fits": 81.0,
        "expected_cdelt_deg": cdelt,
        "expected_psf_fwhm_arcsec": 0.033 * wavelength,
    }


def synthetic_case() -> dict[str, Any]:
    arguments = _synthetic_arguments()
    old = _capture(lambda: legacy_measure(**arguments))
    new_arguments = dict(arguments)
    new_arguments["sed_flux_140pc_w_m2"] = new_arguments.pop("sed")
    new_arguments["image_planes_140pc_w_m2_pixel"] = new_arguments.pop("image")
    new = _capture(
        lambda: kernel.measure_method2_signed_aperture_v1(**new_arguments)
    )
    differences: list[str] = []
    if old["outcome"] != new["outcome"]:
        differences.append("synthetic outcome mismatch")
    elif old["outcome"] == "pass":
        _compare_numbers(old["value"], new["value"], "$", differences)
    return {
        "case": "deterministic_synthetic_structure",
        "status": "parity_pass" if not differences else "parity_fail",
        "legacy_outcome": old["outcome"],
        "new_outcome": new["outcome"],
        "differences": differences,
    }


def synthetic_negative_rejection_case() -> dict[str, Any]:
    """Exercise the missing cluster recovery's decisive signed/clip rejection."""
    arguments = _synthetic_arguments(plane0_negative=0.01)
    old = _capture(lambda: legacy_measure(**arguments))
    new_arguments = dict(arguments)
    new_arguments["sed_flux_140pc_w_m2"] = new_arguments.pop("sed")
    new_arguments["image_planes_140pc_w_m2_pixel"] = new_arguments.pop("image")
    new = _capture(
        lambda: kernel.measure_method2_signed_aperture_v1(**new_arguments)
    )
    gate = "direct_aperture_clip_symmetric_fraction_le_1e-3"
    differences: list[str] = []
    if old["outcome"] != "reject" or new["outcome"] != "reject":
        differences.append(
            f"expected reject/reject; got {old['outcome']}/{new['outcome']}"
        )
    if gate not in str(old.get("message", "")):
        differences.append("legacy rejection did not identify the clip gate")
    if gate not in str(new.get("message", "")):
        differences.append("new rejection did not identify the clip gate")
    return {
        "case": "deterministic_negative_plane0_clip_rejection",
        "status": "parity_pass" if not differences else "parity_fail",
        "legacy_outcome": old["outcome"],
        "new_outcome": new["outcome"],
        "expected_failed_gate": gate,
        "differences": differences,
        "scope_note": (
            "structural fallback only; cluster task024 recovery product parity remains"
        ),
    }


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def run(*, include_real_products: bool) -> dict[str, Any]:
    cases = [synthetic_case(), synthetic_negative_rejection_case()]
    if include_real_products:
        cases.extend(compare_case(str(item["case"]), Path(item["path"])) for item in LOCAL_CASES)
    failed = [item["case"] for item in cases if item["status"] == "parity_fail"]
    unavailable = [
        item["case"]
        for item in cases
        if item["status"] == "not_available_locally"
    ]
    exercised = [
        item["case"]
        for item in cases
        if item["status"] == "parity_pass"
    ]
    return {
        "schema_version": 1,
        "parity_id": "strict_continuum_measurement_kernel_v1_old_vs_new",
        "kernel_id": kernel.KERNEL_ID,
        "legacy_module": str(Path(legacy.__file__).resolve()),
        "legacy_module_sha256": kernel.sha256(Path(legacy.__file__).resolve()),
        "kernel_module": str(Path(kernel.__file__).resolve()),
        "kernel_module_sha256": kernel.sha256(Path(kernel.__file__).resolve()),
        "parity_utility_sha256": kernel.sha256(Path(__file__).resolve()),
        "constants_exactly_equal": all(
            getattr(kernel, name) == getattr(legacy, name)
            for name in (
                "MODEL_DISTANCE_PC",
                "TARGET_DISTANCE_PC",
                "DISTANCE_FACTOR",
                "APERTURE_RADIUS_ARCSEC",
                "APERTURE_SUBPIXELS",
                "PSF_TRUNCATION_SIGMA",
                "COEVAL_TOLERANCE",
                "NEGATIVE_FRACTION_MAX",
                "EDGE_FRACTION_MAX",
                "CONVOLUTION_LOSS_MAX",
                "APERTURE_CLIP_MAX",
                "COMPONENT_SIGNIFICANT_FRACTION",
                "COMPONENT_ACTIVE_FRACTION",
                "COMPONENT_APERTURE_FRACTION_MIN",
                "COMPONENT_APERTURE_FRACTION_MAX",
            )
        ),
        "local_parity_passed": not failed,
        "complete_representative_product_coverage": not unavailable,
        "exercised_cases": exercised,
        "failed_cases": failed,
        "cluster_parity_remaining": unavailable,
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--synthetic-only",
        action="store_true",
        help="skip real local FITS maps (useful for a quick CI smoke test)",
    )
    args = parser.parse_args()
    report = run(include_real_products=not args.synthetic_only)
    atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["constants_exactly_equal"] or not report["local_parity_passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
