"""Versioned acceptance policy; the signed aperture estimator is unchanged."""
from copy import deepcopy
import math

STRICT_POLICY = "whole_image_v1"
APERTURE_POLICY = "aperture_v2"
BOUNDARY_CHECKS = frozenset({"edge_positive_fraction_lt_1e-3", "convolution_flux_loss_lt_1e-4"})
APERTURE_CHECKS = frozenset({"finite_positive_plane0_flux", "finite_positive_signed_aperture",
                             "aperture_plus_6sigma_support"})


def discrete_support_clear(diagnostics):
    """Conservative pixel bound including fractional pixels and kernel rounding."""
    try:
        d = diagnostics
        radius, sigma, truncate = d["aperture_radius_pixels"], d["psf_sigma_pixels"], d["psf_truncation_sigma"]
        cx, cy, nx, ny = (d[k] for k in ("centre_x_zero_based", "centre_y_zero_based", "nx", "ny"))
        if not all(math.isfinite(v) for v in (radius, sigma, truncate, cx, cy, nx, ny)):
            return False
        extent = math.ceil(radius + 1) + int(truncate*sigma + .5)
        return bool(radius > 0 and sigma >= 0 and truncate > 0 and min(cx, cy, nx-1-cx, ny-1-cy) > extent)
    except (KeyError, TypeError, ValueError, OverflowError):
        return False


def validate_policy(policy):
    if policy not in (STRICT_POLICY, APERTURE_POLICY):
        raise ValueError(f"Unknown quality_policy: {policy!r}")
    return policy


def classify_checks(checks, policy):
    validate_policy(policy)
    blocking = {k: v for k, v in checks.items() if policy == STRICT_POLICY or k not in BOUNDARY_CHECKS}
    warnings = sorted(k for k, v in checks.items() if policy == APERTURE_POLICY and k in BOUNDARY_CHECKS and not v)
    return blocking, warnings


def migrate_strict_measurement(item):
    """A previously passing stricter record implies aperture_v2 acceptance.

    This is metadata migration, not a claim of remeasurement. Its caller must
    verify the source record, simulator identity and raw product hashes.
    """
    result = deepcopy(item)
    checks = result.get("diagnostics", {}).get("quality_checks", {})
    if (result.get("quality_policy", STRICT_POLICY) != STRICT_POLICY
            or result.get("quality_pass") is not True or result.get("valid") is not True
            or not (APERTURE_CHECKS | BOUNDARY_CHECKS) <= set(checks)
            or any(v is not True for v in checks.values())):
        raise ValueError("Only a fully passing strict measurement can be migrated")
    checks["finite_kernel_aperture_support"] = discrete_support_clear(result["diagnostics"])
    if not checks["finite_kernel_aperture_support"]:
        raise ValueError("Source measurement lacks a conservative discrete aperture/PSF support margin")
    blocking, warnings = classify_checks(checks, APERTURE_POLICY)
    result.update(quality_policy=APERTURE_POLICY, quality_warnings=warnings)
    result["diagnostics"]["blocking_quality_checks"] = blocking
    return result


def cached_policy_error(item, expected):
    validate_policy(expected)
    if item.get("quality_policy", STRICT_POLICY) != expected:
        return "Measurement quality policy differs from the run"
    if expected == APERTURE_POLICY:
        diagnostics = item.get("diagnostics", {})
        if not isinstance(diagnostics, dict) or not isinstance(diagnostics.get("quality_checks"), dict):
            return "Aperture policy record is missing diagnostic checks"
        checks = diagnostics["quality_checks"]
        blocking, warnings = classify_checks(checks, expected)
        if (not (APERTURE_CHECKS | BOUNDARY_CHECKS | {"finite_kernel_aperture_support"}) <= set(checks)
                or any(not isinstance(v, bool) for v in checks.values())
                or any(v is not True for v in blocking.values())
                or not discrete_support_clear(diagnostics)
                or item.get("valid") is not True
                or item.get("quality_warnings") != warnings
                or item.get("diagnostics", {}).get("blocking_quality_checks") != blocking):
            return "Aperture policy record is missing or fails mandatory checks"
    return None
