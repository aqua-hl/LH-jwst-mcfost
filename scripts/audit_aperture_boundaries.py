#!/usr/bin/env python3
"""Read-only aperture/boundary audit of failed FITS; never accepts or reruns models."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"


def compare_crop(kernel, plane, diagnostics):
    """Compare the frozen adjoint sum with direct convolution of a central crop."""
    import numpy as np
    from scipy.ndimage import gaussian_filter

    d = diagnostics
    cx, cy = d["centre_x_zero_based"], d["centre_y_zero_based"]
    radius, sigma = d["aperture_radius_pixels"], d["psf_sigma_pixels"]
    truncate, subpixels = d["psf_truncation_sigma"], d["aperture_subpixels"]
    # Match scipy.ndimage's finite kernel radius, with two extra pixels beyond
    # the aperture boundary for fractional pixels and integer rounding.
    kernel_radius = int(truncate * sigma + 0.5)
    extent = radius + kernel_radius + 2
    x0, x1 = math.floor(cx - extent), math.ceil(cx + extent) + 1
    y0, y1 = math.floor(cy - extent), math.ceil(cy + extent) + 1
    ny, nx = plane.shape
    if not (0 < x0 < x1 < nx and 0 < y0 < y1 < ny):
        raise ValueError("No strictly interior crop contains this aperture and PSF support")
    weights, _, _ = kernel._adjoint_weights(plane.shape, cx, cy, radius, sigma, subpixels)
    outside_max = max(float(np.max(np.abs(a))) for a in
                      (weights[:y0], weights[y1:], weights[y0:y1, :x0], weights[y0:y1, x1:]))
    full = float(np.einsum("ij,ij->", plane, weights, optimize=False))
    crop = plane[y0:y1, x0:x1]
    aperture = kernel.fractional_circle_weights(crop.shape, cx-x0, cy-y0, radius, subpixels)
    blurred = gaussian_filter(crop, sigma, mode="constant", cval=0.0, truncate=truncate)
    direct = float(np.sum(blurred * aperture, dtype=np.float64))
    if not (math.isfinite(full) and full > 0 and math.isfinite(direct) and direct > 0):
        raise ValueError("Aperture flux is nonfinite or nonpositive")
    relative = abs(direct-full) / abs(full)
    return {
        "full_image_adjoint_aperture_w_m2": full,
        "central_crop_direct_aperture_w_m2": direct,
        "relative_difference": relative,
        "relative_tolerance": 1e-10,
        "crop_bounds_xy_exclusive": [x0, x1, y0, y1],
        "maximum_absolute_aperture_weight_outside_crop": outside_max,
        "aperture_crop_invariance_verified": bool(outside_max == 0 and relative <= 1e-10),
    }


def audit_model(run, manifest, index, kernel, sha256):
    model = manifest["models"][index]
    directory = run / "models" / model["id"]
    status = json.loads((directory / "status.json").read_text())
    permitted = {"edge_positive_fraction_lt_1e-3", "convolution_flux_loss_lt_1e-4"}
    failures = set(status.get("failed_checks", []))
    if status.get("state") != "failed" or not failures or not failures <= permitted:
        raise ValueError(f"Model {index}: expected only recorded whole-image boundary failures")
    saved = status["diagnostics"]
    anchors = [a for a in manifest["anchors"]
               if math.isclose(a["psf_fwhm_arcsec"], saved["psf_fwhm_arcsec"], rel_tol=1e-12)]
    if len(anchors) != 1:
        raise ValueError(f"Model {index}: cannot uniquely identify failed anchor from PSF width")
    anchor = anchors[0]
    attempts = sorted((directory / "anchors" / anchor["id"] / "attempts").glob("attempt_*"),
                      key=lambda p: int(p.name.split("_")[-1]))
    images = list((attempts[-1] / "output").rglob("RT.fits.gz")) if attempts else []
    if len(images) != 1:
        raise ValueError(f"Model {index}: latest failed attempt must contain exactly one RT.fits.gz")
    path = images[0]
    digest = sha256(path)
    try:
        kernel.measure_image(path, wavelength_um=anchor["wavelength_um"], instrument=anchor["instrument"],
                             model_distance_pc=model["parameters"]["distance_pc"],
                             psf_fwhm_arcsec=anchor["psf_fwhm_arcsec"], **manifest["measurement"])
    except kernel.PhotometryError as exc:
        if set(exc.failed_checks) != failures:
            raise ValueError(f"Model {index}: remeasurement failure checks differ from saved status") from exc
        measured = exc.diagnostics
    else:
        raise ValueError(f"Model {index}: latest FITS no longer reproduces the recorded failure")
    for key in ("nx", "ny", "aperture_flux_at_model_distance_w_m2",
                "plane0_edge_positive_flux_fraction", "plane0_convolution_flux_loss_fraction"):
        if not math.isclose(measured[key], saved[key], rel_tol=1e-9, abs_tol=0):
            raise ValueError(f"Model {index}: latest FITS and saved status disagree for {key}")
    plane, _ = kernel._read_image(path, anchor["wavelength_um"])
    comparison = compare_crop(kernel, plane, measured)
    if sha256(path) != digest:
        raise ValueError(f"Model {index}: FITS changed while it was audited")
    return {"index": index, "model_id": model["id"], "anchor_id": anchor["id"],
            "wavelength_um": anchor["wavelength_um"], "image": str(path.relative_to(run)),
            "image_sha256": digest, "original_failed_checks": sorted(failures),
            "original_diagnostics_reproduced": True, "original_diagnostics": measured,
            "crop_comparison": comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--indices", required=True, help="Comma-separated indices, e.g. 48,49")
    args = parser.parse_args()
    run = args.run.resolve()
    manifest_path = run / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    indices = sorted({int(x) for x in args.indices.split(",")})
    if not indices or min(indices) < 0 or max(indices) >= len(manifest["models"]):
        parser.error("Indices are outside this catalogue")
    # Use the run's measurement implementation; never replace its quality policy.
    sys.path.insert(0, str(run / "code/src"))
    from mcfost_grid import photometry as kernel
    from mcfost_grid.configuration import sha256
    from mcfost_grid.runner import validate_inputs
    validate_inputs(run, manifest)
    import numpy, scipy, astropy
    models = [audit_model(run, manifest, i, kernel, sha256) for i in indices]
    report = {
        "audit": "read_only_aperture_crop_invariance_v1", "run_id": manifest["run_id"],
        "manifest_sha256": sha256(manifest_path),
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "packages": {p.__name__: p.__version__ for p in (numpy, scipy, astropy)},
        "models": models, "original_acceptance_and_results_unchanged": True,
        "scope": "Same-FITS aperture operator check; no new MCFOST simulation, field expansion, photon convergence or physical-model validation",
    }
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if all(m["crop_comparison"]["aperture_crop_invariance_verified"] for m in models) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
