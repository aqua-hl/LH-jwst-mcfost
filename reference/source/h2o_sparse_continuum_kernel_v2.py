#!/usr/bin/env python3
"""Scale-invariant six-support Huber continuum for H2O screening v2.

The v1 helper used absolute convergence thresholds in the native flux units.
That is unsafe when the same operator is applied to Jy observations and to
opacity values with a very different numerical scale.  This successor first
normalizes the six positive support values by their median, runs Huber IRLS in
that dimensionless space, uses a relative coefficient-convergence criterion,
and finally restores the input scale.  Multiplying all support values by any
positive constant therefore cannot change the fitted continuum shape.

This module is numerical only.  It does not read or write workflow products.
It also exposes the explicit ``lambda F_lambda`` to F_nu-proportional shape
conversion required before a future Method-2 RT continuum fit.  The opacity-
only Stage-1 rescore does not use that flux conversion.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np


KERNEL_ID = "six_support_scale_invariant_huber_continuum_v2"
SCHEMA_VERSION = 2
HUBER_TUNING = 1.345
# Six-point Huber IRLS can converge slowly when one wing is almost exactly
# linear while the other has a small coherent curvature.  A real smooth
# positive support vector required 248 updates at the pinned relative
# tolerance.  Keep the strict convergence test, allow ample deterministic
# headroom, and still fail closed if the fixed cap is reached.
MAX_ITERATIONS = 10_000
RELATIVE_COEFFICIENT_TOLERANCE = 1.0e-13


class ContinuumError(RuntimeError):
    """The six-support continuum contract failed."""


def lambda_f_lambda_to_fnu_proportional(
    wavelength_um: Sequence[float], lambda_f_lambda: Sequence[float]
) -> np.ndarray:
    """Convert a ``lambda F_lambda`` vector to an F_nu-proportional vector.

    Since ``F_nu = lambda**2 F_lambda / c``, its shape is proportional to
    ``lambda * (lambda F_lambda)``.  The common positive unit constant is
    deliberately omitted because screening fits and normalizes only shape.
    """

    wavelength = np.asarray(wavelength_um, dtype=np.float64)
    values = np.asarray(lambda_f_lambda, dtype=np.float64)
    if not (
        wavelength.ndim == values.ndim == 1
        and len(wavelength) > 0
        and wavelength.shape == values.shape
        and np.all(np.isfinite(wavelength) & (wavelength > 0.0))
        and np.all(np.isfinite(values))
    ):
        raise ContinuumError("lambda F_lambda conversion needs paired finite vectors")
    converted = wavelength * values
    if not np.all(np.isfinite(converted)):
        raise ContinuumError("lambda F_lambda conversion overflowed")
    return converted


def robust_scale_invariant_linear_continuum(
    evaluation_wavelength_um: Sequence[float],
    support_wavelength_um: Sequence[float],
    support_values: Sequence[float],
    support_side: Sequence[str],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fit a homogeneous six-point Huber line and evaluate it.

    The continuum is exactly homogeneous in intent: replacing ``support_values``
    by ``a * support_values`` for positive ``a`` replaces the returned continuum
    by ``a * continuum`` while leaving weights and dimensionless coefficients
    unchanged to floating-point precision.
    """

    evaluation = np.asarray(evaluation_wavelength_um, dtype=np.float64)
    wave = np.asarray(support_wavelength_um, dtype=np.float64)
    values = np.asarray(support_values, dtype=np.float64)
    sides = np.asarray(support_side, dtype=str)
    if not (
        evaluation.ndim == wave.ndim == values.ndim == sides.ndim == 1
        and len(evaluation) > 0
        and len(wave) == len(values) == len(sides) == 6
        and np.count_nonzero(sides == "blue") == 3
        and np.count_nonzero(sides == "red") == 3
        and np.all(np.isfinite(evaluation) & (evaluation > 0.0))
        and np.all(np.isfinite(wave) & (wave > 0.0))
        and np.all(np.isfinite(values) & (values > 0.0))
        and len(np.unique(wave)) == 6
        and float(np.max(wave[sides == "blue"]))
        < float(np.min(wave[sides == "red"]))
    ):
        raise ContinuumError(
            "Continuum needs exactly three positive, distinct samples per wing"
        )

    value_scale = float(np.median(values))
    if not math.isfinite(value_scale) or value_scale <= 0.0:
        raise ContinuumError("Continuum support normalization is invalid")
    normalized_values = values / value_scale
    if not np.all(np.isfinite(normalized_values) & (normalized_values > 0.0)):
        raise ContinuumError("Normalized continuum supports are invalid")

    wavelength_center = float(np.mean(wave))
    wavelength_scale = float(np.ptp(wave))
    if not math.isfinite(wavelength_scale) or wavelength_scale <= 0.0:
        raise ContinuumError("Continuum support has no wavelength leverage")
    normalized_x = (wave - wavelength_center) / wavelength_scale
    design = np.column_stack((np.ones(6, dtype=np.float64), normalized_x))
    coefficients, _, rank, _ = np.linalg.lstsq(
        design, normalized_values, rcond=None
    )
    if rank != 2 or not np.all(np.isfinite(coefficients)):
        raise ContinuumError("Initial six-support continuum is rank deficient")

    final_weights = np.ones(6, dtype=np.float64)
    relative_change = math.inf
    converged = False
    iterations = 0
    for iterations in range(1, MAX_ITERATIONS + 1):
        residual = normalized_values - design @ coefficients
        residual_center = float(np.median(residual))
        robust_scale = 1.4826 * float(
            np.median(np.abs(residual - residual_center))
        )
        scale_floor = np.finfo(np.float64).eps * max(
            float(np.max(np.abs(normalized_values))), 1.0
        )
        if not math.isfinite(robust_scale):
            raise ContinuumError("Huber residual scale is non-finite")
        if robust_scale <= scale_floor:
            relative_change = 0.0
            converged = True
            break

        standardized = np.abs(residual - residual_center) / (
            HUBER_TUNING * robust_scale
        )
        weights = np.ones(6, dtype=np.float64)
        downweighted = standardized > 1.0
        weights[downweighted] = 1.0 / standardized[downweighted]
        weighted_design = design * np.sqrt(weights)[:, None]
        weighted_values = normalized_values * np.sqrt(weights)
        updated, _, rank, _ = np.linalg.lstsq(
            weighted_design, weighted_values, rcond=None
        )
        if rank != 2 or not np.all(np.isfinite(updated)):
            raise ContinuumError("Huber continuum became rank deficient")
        denominator = max(
            float(np.linalg.norm(coefficients, ord=np.inf)),
            float(np.linalg.norm(updated, ord=np.inf)),
            np.finfo(np.float64).tiny,
        )
        relative_change = float(
            np.linalg.norm(updated - coefficients, ord=np.inf) / denominator
        )
        coefficients = updated
        final_weights = weights
        if relative_change <= RELATIVE_COEFFICIENT_TOLERANCE:
            converged = True
            break
    if not converged:
        raise ContinuumError(
            f"Huber continuum did not converge in {MAX_ITERATIONS} iterations"
        )

    dimensionless_continuum = coefficients[0] + coefficients[1] * (
        (evaluation - wavelength_center) / wavelength_scale
    )
    continuum = value_scale * dimensionless_continuum
    if not np.all(np.isfinite(continuum) & (continuum > 0.0)):
        raise ContinuumError("Scale-invariant continuum is not finite and positive")
    diagnostics = {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "method": "median_normalized_deterministic_huber_irls_line_on_six_supports",
        "support_count": 6,
        "support_value_normalization": "median_of_six_positive_support_values",
        "support_value_scale_in_input_units": value_scale,
        "dimensionless_intercept": float(coefficients[0]),
        "dimensionless_slope": float(coefficients[1]),
        "wavelength_center_um": wavelength_center,
        "wavelength_scale_um": wavelength_scale,
        "huber_tuning": HUBER_TUNING,
        "iterations": iterations,
        "maximum_iterations": MAX_ITERATIONS,
        "relative_coefficient_tolerance": RELATIVE_COEFFICIENT_TOLERANCE,
        "final_relative_coefficient_change": relative_change,
        "converged": converged,
        "support_weights": final_weights.tolist(),
    }
    return continuum, diagnostics
