"""Fast signed Method-2 image photometry for the restart workflow.

The measurement is a direct aperture integral of total-I plane zero; an SED
is never used to rescale it.  For a symmetric, zero-padded Gaussian operator
G, <G image, aperture> = <image, G aperture>.  Cached adjoint weights thus
avoid convolving every new image, and all separated components, repeatedly.

This module implements lightweight input/geometry/flux checks, not the full
historical component/provenance audit or a convergence experiment.  Method-2
execution and the provenance of an optionally supplied coeval SED must be
established by the runner.  The Gaussian PSF is an approximation inherited
from the reference kernel, not a wavelength-dependent empirical JWST PSF.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.io import fits
from scipy.ndimage import gaussian_filter, gaussian_filter1d


KERNEL_ID = "restart_signed_method2_adjoint_aperture_v1"
PSF_TRUNCATION_SIGMA = 6.0
FWHM_TO_SIGMA = 2.354820045  # Frozen historical convention.
SPEED_OF_LIGHT_M_S = 299792458.0


class PhotometryError(ValueError):
    """Invalid input or a failed numerical gate; no score should be produced."""

    def __init__(self, message: str, *, failed_checks: tuple[str, ...] = (),
                 diagnostics: dict[str, Any] | None = None):
        super().__init__(message)
        self.failed_checks = failed_checks
        self.diagnostics = diagnostics or {}


def _positive(value: float, name: str) -> float:
    try:
        value = float(value)
    except (ValueError, TypeError) as exc:
        raise PhotometryError(f"{name} must be a finite positive number") from exc
    if not math.isfinite(value) or value <= 0:
        raise PhotometryError(f"{name} must be a finite positive number")
    return value


def gaussian_psf_fwhm_arcsec(instrument: str, wavelength_um: float) -> float:
    """Historical Gaussian FWHM: 0.033 lambda (+0.106 for MIRI), arcsec."""
    wavelength_um = _positive(wavelength_um, "wavelength_um")
    if not isinstance(instrument, str):
        raise PhotometryError("instrument must be NIRSpec or MIRI")
    if instrument.upper() == "NIRSPEC":
        return 0.033 * wavelength_um
    if instrument.upper() == "MIRI":
        return 0.033 * wavelength_um + 0.106
    raise PhotometryError(f"Unsupported instrument: {instrument!r}")


def fractional_circle_weights(shape: tuple[int, int], centre_x: float,
                              centre_y: float, radius_pixels: float,
                              subpixels: int = 16) -> np.ndarray:
    """Midpoint fractional circle, evaluating subpixels only at its boundary.

    Setting subpixels=64 reproduces the frozen kernel's midpoint convention.
    Fully interior/exterior pixels are classified from their square's farthest
    and nearest distances, rather than sampling the whole aperture N^2 times.
    """
    if (not isinstance(subpixels, (int, np.integer)) or isinstance(subpixels, bool)
            or not 1 <= subpixels <= 512):
        raise PhotometryError("aperture_subpixels must be an integer from 1 to 512")
    radius_pixels = _positive(radius_pixels, "radius_pixels")
    ny, nx = shape
    if ny < 1 or nx < 1 or not np.all(np.isfinite([centre_x, centre_y])):
        raise PhotometryError("Circle geometry is not finite or has an empty image")
    dx = np.abs(np.arange(nx, dtype=np.float64) - centre_x)[None, :]
    dy = np.abs(np.arange(ny, dtype=np.float64) - centre_y)[:, None]
    fully_inside = (dx + 0.5)**2 + (dy + 0.5)**2 <= radius_pixels**2
    touches = np.maximum(dx - 0.5, 0)**2 + np.maximum(dy - 0.5, 0)**2 <= radius_pixels**2
    weights = fully_inside.astype(np.float64)
    by, bx = np.nonzero(touches & ~fully_inside)
    offsets = (np.arange(subpixels, dtype=np.float64) + 0.5) / subpixels - 0.5
    # The chunk bounds intermediate memory even at the reference 64x64 rule.
    chunk = max(1, 524288 // (subpixels * subpixels))
    for first in range(0, len(bx), chunk):
        x, y = bx[first:first + chunk], by[first:first + chunk]
        xsq = (x[:, None] + offsets - centre_x)**2
        ysq = (y[:, None] + offsets - centre_y)**2
        weights[y, x] = np.mean(ysq[:, :, None] + xsq[:, None, :] <= radius_pixels**2,
                                axis=(1, 2))
    return weights


@lru_cache(maxsize=4)
def _adjoint_weights(shape: tuple[int, int], centre_x: float, centre_y: float,
                     radius_pixels: float, sigma_pixels: float,
                     subpixels: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cache aperture response and separable whole-field survival weights."""
    weights = fractional_circle_weights(shape, centre_x, centre_y, radius_pixels, subpixels)
    if sigma_pixels > 0:
        weights = gaussian_filter(weights, sigma_pixels, mode="constant", cval=0.0,
                                  truncate=PSF_TRUNCATION_SIGMA)
        sy = gaussian_filter1d(np.ones(shape[0]), sigma_pixels, mode="constant", cval=0.0,
                               truncate=PSF_TRUNCATION_SIGMA)
        sx = gaussian_filter1d(np.ones(shape[1]), sigma_pixels, mode="constant", cval=0.0,
                               truncate=PSF_TRUNCATION_SIGMA)
    else:
        sy, sx = np.ones(shape[0]), np.ones(shape[1])
    for array in (weights, sy, sx):
        array.flags.writeable = False
    return weights, sy, sx


def _unit_matches(header: fits.Header, expected: u.UnitBase) -> bool:
    try:
        return u.Unit(str(header.get("BUNIT", ""))) == expected
    except (ValueError, TypeError):
        return False


def _read_image(path: Path, wavelength_um: float) -> tuple[np.ndarray, fits.Header]:
    try:
        with fits.open(path, memmap=False) as hdus:
            header = hdus[0].header.copy()
            data = hdus[0].data
            if (data is None or data.ndim != 5 or data.shape[:3] != (8, 1, 1)
                    or min(data.shape[-2:]) < 3):
                raise PhotometryError("Expected polarized MCFOST image shape (8,1,1,ny,nx)")
            if not np.all(np.isfinite(data)):
                raise PhotometryError("Image contains non-finite pixels")
            plane = np.array(data[0, 0, 0], dtype=np.float64, copy=True)
    except (OSError, TypeError) as exc:
        raise PhotometryError(f"Could not read MCFOST image {path}: {exc}") from exc
    if not _unit_matches(header, u.W / u.m**2 / u.pix):
        raise PhotometryError("Image BUNIT must be W.m-2.pixel-1 (pixel-integrated lambda F_lambda)")
    if any(key in header for key in ("BSCALE", "BZERO")):
        raise PhotometryError("Unexpected FITS BSCALE/BZERO scaling")
    if "total flux" not in str(header.get("FLUX_1", "")).lower():
        raise PhotometryError("FLUX_1 does not identify the total-intensity flux plane")
    components = {5: ("direct", "star"), 6: ("scattered", "star"),
                  7: ("direct", "thermal"), 8: ("scattered", "thermal")}
    if not all(all(token in str(header.get(f"FLUX_{i}", "")).lower() for token in tokens)
               for i, tokens in components.items()):
        raise PhotometryError("Expected polarized image component labels FLUX_5 through FLUX_8")
    try:
        wave = float(header.get("WAVE", math.nan))
    except (ValueError, TypeError) as exc:
        raise PhotometryError("Image WAVE is invalid") from exc
    if not math.isclose(wave, wavelength_um, rel_tol=1e-6, abs_tol=1e-8):
        raise PhotometryError("Image WAVE does not match requested wavelength")
    return plane, header


def _geometry(header: fits.Header, shape: tuple[int, int]) -> tuple[float, float, float]:
    forbidden = ("PC1_1", "PC1_2", "PC2_1", "PC2_2", "CD1_1", "CD1_2", "CD2_1",
                 "CD2_2", "CROTA1", "CROTA2")
    if any(key in header for key in forbidden):
        raise PhotometryError("Rotated/matrix WCS is unsupported; require axis-aligned CDELT")
    for key in ("CUNIT1", "CUNIT2"):
        if key in header and str(header[key]).strip().lower() not in ("deg", "degree", "degrees"):
            raise PhotometryError("Spatial WCS CUNIT must be degrees")
    try:
        values = np.asarray([header.get(key, math.nan) for key in
                             ("CDELT1", "CDELT2", "CRPIX1", "CRPIX2")], dtype=float)
    except (TypeError, ValueError) as exc:
        raise PhotometryError("Invalid spatial WCS numeric values") from exc
    dx, dy, cx, cy = values
    cx, cy = cx - 1, cy - 1
    ny, nx = shape
    if (not np.all(np.isfinite(values)) or dx == 0 or dy == 0
            or not math.isclose(abs(dx), abs(dy), rel_tol=5e-6)
            or not math.isclose(cx, (nx - 1) / 2, abs_tol=1e-6, rel_tol=0)
            or not math.isclose(cy, (ny - 1) / 2, abs_tol=1e-6, rel_tol=0)):
        raise PhotometryError("Require finite square pixels and source-centred midpoint CRPIX")
    return float(cx), float(cy), float(abs(dx) * 3600)


def _coeval_total(path: Path, wavelength_um: float) -> float:
    try:
        with fits.open(path, memmap=False) as hdus:
            data = hdus[0].data
            header = hdus[0].header
            if data is None or data.shape != (8, 1, 1, 1) or len(hdus) < 2:
                raise PhotometryError("Coeval SED must have shape (8,1,1,1) and wavelength HDU")
            wave = np.asarray(hdus[1].data, dtype=float).ravel()
            if (not np.all(np.isfinite(data)) or wave.shape != (1,)
                    or not math.isclose(float(wave[0]), wavelength_um, rel_tol=1e-6, abs_tol=1e-8)
                    or not _unit_matches(header, u.W / u.m**2)):
                raise PhotometryError("Invalid coeval SED units, wavelength, or finite-flux contract")
            if any(key in header for key in ("BSCALE", "BZERO")):
                raise PhotometryError("Unexpected coeval SED FITS scaling")
            return float(data[0, 0, 0, 0])
    except (OSError, TypeError) as exc:
        raise PhotometryError(f"Could not read coeval SED {path}: {exc}") from exc


def measure_image(path: str | Path, wavelength_um: float, instrument: str,
                  aperture_radius_arcsec: float, model_distance_pc: float,
                  target_distance_pc: float, aperture_subpixels: int = 16,
                  psf_fwhm_arcsec: float | None = None,
                  coeval_sed_path: str | Path | None = None) -> dict[str, Any]:
    """Return JSON-serializable signed plane-0 aperture flux, or raise.

    ``flux_w_m2`` means lambda F_lambda at the target distance. ``flux_jy``
    is its F_nu conversion at wavelength_um. Image totals are also at target
    distance, with a separately named model-distance value in diagnostics.
    The optional SED validates plane-0 closure only; it cannot normalize flux.
    Explicit PSF FWHM zero disables convolution (useful for analytical tests).
    """
    wavelength_um = _positive(wavelength_um, "wavelength_um")
    aperture_radius_arcsec = _positive(aperture_radius_arcsec, "aperture_radius_arcsec")
    model_distance_pc = _positive(model_distance_pc, "model_distance_pc")
    target_distance_pc = _positive(target_distance_pc, "target_distance_pc")
    # Validate before the cache lookup: 16.0 and 16 compare equal as cache keys.
    if (not isinstance(aperture_subpixels, (int, np.integer))
            or isinstance(aperture_subpixels, bool) or not 1 <= aperture_subpixels <= 512):
        raise PhotometryError("aperture_subpixels must be an integer from 1 to 512")
    default_fwhm = gaussian_psf_fwhm_arcsec(instrument, wavelength_um)
    try:
        fwhm = default_fwhm if psf_fwhm_arcsec is None else float(psf_fwhm_arcsec)
    except (ValueError, TypeError) as exc:
        raise PhotometryError("psf_fwhm_arcsec must be finite and non-negative") from exc
    if not math.isfinite(fwhm) or fwhm < 0:
        raise PhotometryError("psf_fwhm_arcsec must be finite and non-negative")
    path = Path(path)
    plane, header = _read_image(path, wavelength_um)
    ny, nx = plane.shape
    cx, cy, model_pixel_arcsec = _geometry(header, plane.shape)
    distance_ratio = model_distance_pc / target_distance_pc
    try:
        flux_factor = distance_ratio**2
    except OverflowError as exc:
        raise PhotometryError("Distance ratio cannot be represented safely") from exc
    if not math.isfinite(flux_factor) or flux_factor <= 0:
        raise PhotometryError("Distance ratio cannot be represented safely")
    pixel_arcsec = model_pixel_arcsec * distance_ratio
    radius_pixels = aperture_radius_arcsec / pixel_arcsec
    sigma_pixels = fwhm / FWHM_TO_SIGMA / pixel_arcsec
    half_field = min((cx + 0.5), (nx - 0.5 - cx), (cy + 0.5), (ny - 0.5 - cy)) * pixel_arcsec
    required_field = aperture_radius_arcsec + PSF_TRUNCATION_SIGMA * fwhm / FWHM_TO_SIGMA
    if half_field < required_field:
        raise PhotometryError("Image does not contain aperture plus six-sigma PSF support",
                              failed_checks=("aperture_plus_6sigma_support",),
                              diagnostics={"available_half_field_arcsec": float(half_field),
                                           "required_half_field_arcsec": float(required_field)})
    weights, survival_y, survival_x = _adjoint_weights(plane.shape, cx, cy, radius_pixels,
                                                       sigma_pixels, aperture_subpixels)
    total = float(np.sum(plane, dtype=np.float64))
    positive = float(np.sum(plane, where=plane > 0, dtype=np.float64))
    negative = -float(np.sum(plane, where=plane < 0, dtype=np.float64))
    if total <= 0 or positive <= 0 or not math.isfinite(total):
        raise PhotometryError("Image total-I signed flux must be finite and positive",
                              failed_checks=("finite_positive_plane0_flux",))
    aperture = float(np.einsum("ij,ij->", plane, weights, optimize=False))
    # Optional packet-noise diagnostic, computed with the same adjoint weights.
    # It does not clip the scored image and does not introduce a second PSF
    # convolution. Historical v2/v3 also treated clipping sensitivity as a
    # diagnostic rather than rejecting signed realizations by this threshold.
    clipped_aperture = aperture
    if negative > 0:
        negative_pixels = plane < 0
        negative_aperture = float(np.einsum("i,i->", plane[negative_pixels],
                                          weights[negative_pixels], optimize=False))
        clipped_aperture -= negative_aperture
    clip_denominator = 0.5 * (abs(aperture) + abs(clipped_aperture))
    clipping_sensitivity = (abs(aperture - clipped_aperture) / clip_denominator
                            if clip_denominator > 0 else 0.0)
    convolved_total = float(np.einsum("ij,i,j->", plane, survival_y, survival_x, optimize=False))
    border = max(1, min(nx, ny) // 100)
    edge_positive = sum(float(np.sum(np.maximum(a, 0), dtype=np.float64)) for a in
                        (plane[:border], plane[-border:], plane[border:-border, :border],
                         plane[border:-border, -border:]))
    edge_fraction = edge_positive / positive
    convolution_loss = abs(convolved_total - total) / total
    closure_relative = None
    if coeval_sed_path is not None:
        closure_relative = (_coeval_total(Path(coeval_sed_path), wavelength_um) - total) / total
    checks = {"finite_positive_plane0_flux": True,
              "finite_positive_signed_aperture": bool(math.isfinite(aperture) and aperture > 0),
              "edge_positive_fraction_lt_1e-3": bool(edge_fraction < 1e-3),
              "convolution_flux_loss_lt_1e-4": bool(convolution_loss < 1e-4),
              "aperture_plus_6sigma_support": True}
    if closure_relative is not None:
        checks["coeval_plane0_abs_relative_difference_le_1e-5"] = bool(abs(closure_relative) <= 1e-5)
    diagnostics = {
        "quality_checks": checks,
        "plane0_edge_positive_flux_fraction": edge_fraction,
        "plane0_convolution_flux_loss_fraction": convolution_loss,
        "full_plane_negative_fraction_diagnostic": negative / positive,
        "aperture_clip_symmetric_fraction_diagnostic": clipping_sensitivity,
        "zero_clipped_aperture_flux_at_model_distance_w_m2_diagnostic": clipped_aperture,
        "zero_clipping_used_for_scoring": False,
        "coeval_plane0_relative_difference": closure_relative,
        "coeval_closure_checked": coeval_sed_path is not None,
        "image_total_at_model_distance_w_m2": total,
        "aperture_flux_at_model_distance_w_m2": aperture,
        "distance_flux_factor": flux_factor,
        "model_pixel_arcsec": model_pixel_arcsec,
        "target_pixel_arcsec": pixel_arcsec,
        "aperture_radius_pixels": radius_pixels,
        "centre_x_zero_based": cx, "centre_y_zero_based": cy,
        "available_half_field_arcsec": float(half_field),
        "required_half_field_arcsec": float(required_field),
        "nx": nx, "ny": ny,
        "psf_fwhm_arcsec": fwhm,
        "psf_sigma_pixels": sigma_pixels,
        "psf_truncation_sigma": PSF_TRUNCATION_SIGMA,
        "aperture_subpixels": int(aperture_subpixels),
        "component_integrity_audit_performed": False,
        "temperature_consistency_verified": False,
        "convergence_verified": False,
    }
    failed = tuple(key for key, passed in checks.items() if not passed)
    if failed:
        raise PhotometryError("Image photometry failed checks: " + ", ".join(failed),
                              failed_checks=failed, diagnostics=diagnostics)
    flux = aperture * flux_factor
    flux_jy = flux * wavelength_um * 1e-6 / SPEED_OF_LIGHT_M_S * 1e26
    image_total = total * flux_factor
    if not all(math.isfinite(x) and x > 0 for x in (flux, flux_jy, image_total)):
        raise PhotometryError("Distance/unit conversion produced non-finite or nonpositive flux",
                              failed_checks=("finite_positive_converted_flux",),
                              diagnostics=diagnostics)
    return {
        "schema_version": 1,
        "kernel_id": KERNEL_ID,
        "valid": True,
        "quality_pass": True,
        "quality_status": "passed" if coeval_sed_path is not None else "passed_without_coeval_closure",
        "canonical_estimator": "signed_direct_method2_polarized_total_i_plane0",
        "image_path": str(path.resolve()),
        "coeval_sed_path": None if coeval_sed_path is None else str(Path(coeval_sed_path).resolve()),
        "wavelength_um": wavelength_um,
        "instrument": "NIRSpec" if instrument.upper() == "NIRSPEC" else "MIRI",
        "aperture_radius_arcsec": aperture_radius_arcsec,
        "model_distance_pc": model_distance_pc,
        "target_distance_pc": target_distance_pc,
        "flux_w_m2": flux,
        "flux_jy": flux_jy,
        "image_total_w_m2": image_total,
        "diagnostics": diagnostics,
        "limitations": ["Gaussian PSF approximation; source-centred circular aperture",
                        "Method-2 and coeval execution provenance must be recorded by runner",
                        "No full historical component audit or convergence validation"],
    }
