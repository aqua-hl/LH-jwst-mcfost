#!/usr/bin/env python3
"""Numerical observation/model operator for the TMC1A H2O v3 screen.

The frozen observation is expressed as ``F_nu`` in Jy.  MCFOST Method-2
plane 0 is ``lambda F_lambda`` in W m-2, so every model LSF node is converted
to Jy *before* Gauss-Hermite aggregation.  The six-support continuum is then
fit with the scale-invariant v2 Huber kernel in the same representation for
both the observation and the model.

This module is deliberately free of filesystem I/O.  It also implements the
operational relative segment-gain sensitivity: use the zero-sum log-gain gauge
``g235=exp(-delta/2), g395=exp(+delta/2)`` on the *observation only*, then
refit the continuum and recompute all optical depths.  The canonical screen
fixes ``delta=0``; specified nonzero values are robustness scenarios, not fit
parameters.  A downstream scorer must never add 3% as an independent
per-anchor variance.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np

import h2o_sparse_continuum_kernel_v2 as continuum_kernel


KERNEL_ID = "tmc1a_h2o_observation_measurement_v3"
SCHEMA_VERSION = 3
SPEED_OF_LIGHT_M_S = 299_792_458.0
MICRON_TO_M = 1.0e-6
JY_TO_W_M2_HZ = 1.0e-26
LAMBDA_FLAMBDA_TO_FNU_JY_FACTOR = (
    MICRON_TO_M / SPEED_OF_LIGHT_M_S / JY_TO_W_M2_HZ
)
LSF_NODES_PER_TARGET = 5
TARGET_RESOLVING_POWER = 2400.0
WAVELENGTH_ABS_TOLERANCE_UM = 5.0e-12
MCFOST_FITS_WAVE_SIGNIFICANT_DIGITS = 7
EXPECTED_SEGMENTS = ("G235H", "G395H")


class MeasurementContractError(RuntimeError):
    """The v3 observation/model measurement contract was violated."""


def mcfost_fits_wave_half_unit_tolerance_um(requested_wavelength_um: float) -> float:
    """Half a unit in the seventh significant digit of a requested wavelength."""

    requested = float(requested_wavelength_um)
    if not math.isfinite(requested) or requested <= 0.0:
        raise MeasurementContractError("Requested MCFOST wavelength is invalid")
    decimal_place = math.floor(math.log10(requested)) - (
        MCFOST_FITS_WAVE_SIGNIFICANT_DIGITS - 1
    )
    return 0.5 * 10.0**decimal_place


def mcfost_fits_wave_matches_requested(
    product_fits_wave_um: float, requested_wavelength_um: float
) -> bool:
    """Whether a FITS WAVE value is compatible with seven-digit serialization."""

    product = float(product_fits_wave_um)
    requested = float(requested_wavelength_um)
    if not math.isfinite(product) or product <= 0.0:
        return False
    tolerance = mcfost_fits_wave_half_unit_tolerance_um(requested)
    floating_slack = 8.0 * np.finfo(float).eps * max(abs(requested), 1.0)
    return abs(product - requested) <= tolerance + floating_slack


def lambda_f_lambda_w_m2_to_fnu_jy(
    wavelength_um: Sequence[float], lambda_f_lambda_w_m2: Sequence[float]
) -> np.ndarray:
    """Convert MCFOST ``lambda F_lambda`` values to ``F_nu`` in Jy.

    ``F_nu = lambda * (lambda F_lambda) / c``.  Wavelength is supplied in
    microns and 1 Jy is 1e-26 W m-2 Hz-1, giving the exact factor used here.
    """

    proportional = continuum_kernel.lambda_f_lambda_to_fnu_proportional(
        wavelength_um, lambda_f_lambda_w_m2
    )
    converted = proportional * LAMBDA_FLAMBDA_TO_FNU_JY_FACTOR
    if not np.all(np.isfinite(converted)):
        raise MeasurementContractError("F_nu conversion produced non-finite values")
    return converted


def aggregate_model_lsf_nodes_to_fnu_jy(
    node_rows: Sequence[Mapping[str, Any]],
    requested_wavelength_um: Sequence[float],
    product_fits_wavelength_um: Sequence[float],
    node_lambda_f_lambda_w_m2: Sequence[float],
    *,
    expected_parent_ids: Sequence[str],
    parent_id_key: str = "parent_id",
) -> dict[str, float]:
    """Validate the frozen roster, convert node fluxes, and aggregate in Jy.

    Flux-only positional aggregation is deliberately forbidden.  The caller
    must provide (1) exact requested wavelengths read from a separately
    hash-bound lambda/task roster, (2) the lower-precision ``WAVE`` values read
    from the MCFOST FITS products in flux order, and (3) the exact ordered
    parent roster expected by the scorer.  MCFOST 4.1.13 serializes ``WAVE`` at
    about seven significant digits, so the product value is checked against
    the exact request using an explicit half-unit-at-seven-digits bound.  Jy
    conversion uses the exact request, never the rounded FITS value.
    """

    expected = list(expected_parent_ids)
    if (
        not expected
        or any(
            not isinstance(parent, str)
            or not parent.strip()
            or parent.strip().casefold() in {"none", "null", "nan"}
            for parent in expected
        )
        or len(set(expected)) != len(expected)
    ):
        raise MeasurementContractError("Expected parent roster is empty or malformed")
    expected_count = len(expected) * LSF_NODES_PER_TARGET
    if not (
        len(node_rows)
        == len(requested_wavelength_um)
        == len(product_fits_wavelength_um)
        == len(node_lambda_f_lambda_w_m2)
        == expected_count
    ):
        raise MeasurementContractError(
            "LSF rows, MCFOST wavelength/flux vectors, and parent roster differ"
        )
    exact_wavelengths = np.asarray(requested_wavelength_um, dtype=float)
    product_wavelengths = np.asarray(product_fits_wavelength_um, dtype=float)
    model_flux = np.asarray(node_lambda_f_lambda_w_m2, dtype=float)
    if not (
        exact_wavelengths.ndim == product_wavelengths.ndim == model_flux.ndim == 1
        and np.all(np.isfinite(exact_wavelengths) & (exact_wavelengths > 0.0))
        and np.all(np.diff(exact_wavelengths) > 0.0)
        and np.all(np.isfinite(product_wavelengths) & (product_wavelengths > 0.0))
        and np.all(np.diff(product_wavelengths) > 0.0)
        and np.all(np.isfinite(model_flux) & (model_flux > 0.0))
    ):
        raise MeasurementContractError("Malformed or non-monotonic MCFOST node product")

    frozen_wavelengths: list[float] = []
    grouped: dict[str, list[tuple[int, float, float, float]]] = {}
    hermite_nodes, hermite_raw_weights = np.polynomial.hermite.hermgauss(
        LSF_NODES_PER_TARGET
    )
    hermite_weights = hermite_raw_weights / math.sqrt(math.pi)
    sigma_log = 1.0 / (2.354820045 * TARGET_RESOLVING_POWER)
    for position, row in enumerate(node_rows):
        raw_parent = row.get(parent_id_key)
        if (
            not isinstance(raw_parent, str)
            or not raw_parent.strip()
            or raw_parent.strip().casefold() in {"none", "null", "nan"}
        ):
            raise MeasurementContractError(f"Malformed parent ID at node row {position}")
        parent = raw_parent
        expected_parent = expected[position // LSF_NODES_PER_TARGET]
        if parent != expected_parent:
            raise MeasurementContractError(
                f"LSF parent roster/order changed at row {position}: {parent!r}"
            )
        raw_index = row.get("lsf_node_index")
        raw_global_index = row.get("lsf_node_global_index")
        if (
            isinstance(raw_index, (bool, np.bool_))
            or not isinstance(raw_index, (int, np.integer))
            or isinstance(raw_global_index, (bool, np.bool_))
            or not isinstance(raw_global_index, (int, np.integer))
        ):
            raise MeasurementContractError(f"LSF indices must be integers at row {position}")
        index = int(raw_index)
        if index != position % LSF_NODES_PER_TARGET or int(raw_global_index) != position:
            raise MeasurementContractError(f"LSF node index/order changed at row {position}")
        try:
            frozen = float(row["lsf_node_wavelength_um"])
            center = float(row["target_wavelength_um"])
            weight = float(row["lsf_quadrature_weight"])
        except (KeyError, TypeError, ValueError) as error:
            raise MeasurementContractError(f"Malformed numeric LSF row {position}") from error
        expected_node = center * math.exp(
            math.sqrt(2.0) * sigma_log * float(hermite_nodes[index])
        )
        if not (
            math.isfinite(frozen)
            and frozen > 0.0
            and math.isfinite(center)
            and center > 0.0
            and math.isclose(
                frozen,
                expected_node,
                rel_tol=0.0,
                abs_tol=WAVELENGTH_ABS_TOLERANCE_UM,
            )
            and math.isclose(
                weight,
                float(hermite_weights[index]),
                rel_tol=0.0,
                abs_tol=1.0e-14,
            )
            and math.isclose(
                float(exact_wavelengths[position]),
                frozen,
                rel_tol=0.0,
                abs_tol=WAVELENGTH_ABS_TOLERANCE_UM,
            )
        ):
            raise MeasurementContractError(f"LSF wavelength/weight identity failed at row {position}")
        product_value = float(product_wavelengths[position])
        if not mcfost_fits_wave_matches_requested(product_value, frozen):
            raise MeasurementContractError(
                f"MCFOST FITS WAVE does not match exact requested node at row {position}"
            )
        frozen_wavelengths.append(frozen)
        grouped.setdefault(parent, []).append((index, weight, frozen, center))

    if np.any(np.diff(np.asarray(frozen_wavelengths, dtype=float)) <= 0.0):
        raise MeasurementContractError("Frozen LSF-node wavelength roster is not monotonic")
    if list(grouped) != expected or set(grouped) != set(expected):
        raise MeasurementContractError("Frozen LSF parent roster differs from expectation")

    converted = lambda_f_lambda_w_m2_to_fnu_jy(
        exact_wavelengths, model_flux
    )
    result: dict[str, float] = {}
    for parent_index, parent in enumerate(expected):
        values = grouped[parent]
        if (
            len(values) != LSF_NODES_PER_TARGET
            or [value[0] for value in values] != list(range(LSF_NODES_PER_TARGET))
            or any(values[index][2] >= values[index + 1][2] for index in range(4))
            or len({value[3] for value in values}) != 1
        ):
            raise MeasurementContractError(f"Incomplete or unordered LSF nodes for {parent}")
        if not math.isclose(
            math.fsum(value[1] for value in values),
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ):
            raise MeasurementContractError(f"Non-normalized LSF weights for {parent}")
        offset = parent_index * LSF_NODES_PER_TARGET
        result[parent] = math.fsum(
            values[index][1] * float(converted[offset + index])
            for index in range(LSF_NODES_PER_TARGET)
        )
    return result


def _validated_continuum_arrays(
    target_wavelength_um: Sequence[float],
    target_fnu_jy: Sequence[float],
    target_segment: Sequence[str],
    support_wavelength_um: Sequence[float],
    support_fnu_jy: Sequence[float],
    support_segment: Sequence[str],
    support_side: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    target_wave = np.asarray(target_wavelength_um, dtype=float)
    target_flux = np.asarray(target_fnu_jy, dtype=float)
    support_wave = np.asarray(support_wavelength_um, dtype=float)
    support_flux = np.asarray(support_fnu_jy, dtype=float)
    if any(not isinstance(item, str) for item in target_segment) or any(
        not isinstance(item, str) for item in support_segment
    ) or any(not isinstance(item, str) for item in support_side):
        raise MeasurementContractError("Segment and support-side labels must be strings")
    target_groups = np.asarray(target_segment, dtype=str)
    support_groups = np.asarray(support_segment, dtype=str)
    sides = np.asarray(support_side, dtype=str)
    if not (
        target_wave.ndim
        == target_flux.ndim
        == target_groups.ndim
        == support_wave.ndim
        == support_flux.ndim
        == support_groups.ndim
        == sides.ndim
        == 1
        and target_wave.shape == target_flux.shape == target_groups.shape
        and support_wave.shape == support_flux.shape == support_groups.shape == sides.shape
        and len(target_wave) > 0
        and len(support_wave) == 6
        and np.all(np.isfinite(target_wave) & (target_wave > 0.0))
        and np.all(np.diff(target_wave) > 0.0)
        and np.all(np.isfinite(target_flux) & (target_flux > 0.0))
        and np.all(np.isfinite(support_wave) & (support_wave > 0.0))
        and np.all(np.diff(support_wave) > 0.0)
        and np.all(np.isfinite(support_flux) & (support_flux > 0.0))
    ):
        raise MeasurementContractError("Malformed target/support arrays")
    if set(target_groups) != set(EXPECTED_SEGMENTS) or set(support_groups) != set(
        EXPECTED_SEGMENTS
    ):
        raise MeasurementContractError(
            f"Target and support rosters must each contain {list(EXPECTED_SEGMENTS)}"
        )
    blue = sides == "blue"
    red = sides == "red"
    if not (
        np.count_nonzero(blue) == np.count_nonzero(red) == 3
        and np.all(support_groups[blue] == "G235H")
        and np.all(support_groups[red] == "G395H")
    ):
        raise MeasurementContractError(
            "Supports must map exactly: three blue=G235H and three red=G395H"
        )
    return (
        target_wave,
        target_flux,
        target_groups,
        support_wave,
        support_flux,
        support_groups,
        sides,
    )


def compute_fnu_continuum_and_tau(
    target_wavelength_um: Sequence[float],
    target_fnu_jy: Sequence[float],
    target_segment: Sequence[str],
    support_wavelength_um: Sequence[float],
    support_fnu_jy: Sequence[float],
    support_segment: Sequence[str],
    support_side: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Apply the identical neutral Fnu continuum/tau operator to either side.

    This is the sole implementation used for both the observation and the
    model after model LSF nodes have been validated, converted, and aggregated
    to Jy.  It performs no gain fitting or rescaling.
    """

    (
        target_wave,
        target_flux,
        _,
        support_wave,
        support_flux,
        _,
        sides,
    ) = _validated_continuum_arrays(
        target_wavelength_um,
        target_fnu_jy,
        target_segment,
        support_wavelength_um,
        support_fnu_jy,
        support_segment,
        support_side,
    )
    continuum, diagnostics = continuum_kernel.robust_scale_invariant_linear_continuum(
        target_wave, support_wave, support_flux, sides
    )
    tau = -np.log(target_flux / continuum)
    if not np.all(np.isfinite(tau)):
        raise MeasurementContractError("Optical depths are non-finite")
    return continuum, tau, {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "operator": "neutral_identical_Fnu_continuum_and_tau_for_observation_or_model",
        "flux_representation": "F_nu_Jy",
        "continuum": diagnostics,
    }


def recompute_observation_relative_segment_gain_continuum_and_tau(
    target_wavelength_um: Sequence[float],
    target_fnu_jy: Sequence[float],
    target_segment: Sequence[str],
    support_wavelength_um: Sequence[float],
    support_fnu_jy: Sequence[float],
    support_segment: Sequence[str],
    support_side: Sequence[str],
    relative_log_gain_delta: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Apply a gauge-fixed observation relative gain, refit, and recompute tau.

    ``delta = ln(g395/g235)`` with the zero-sum gauge
    ``g235=exp(-delta/2)``, ``g395=exp(+delta/2)``.  The canonical screen fixes
    delta to zero.  Nonzero calls are predeclared observation-only robustness
    scenarios; this function is not an optimizer or a fitted likelihood
    nuisance.  Model fluxes must never be modified by these gains.
    """

    (
        target_wave,
        target_flux,
        target_groups,
        support_wave,
        support_flux,
        support_groups,
        sides,
    ) = _validated_continuum_arrays(
        target_wavelength_um,
        target_fnu_jy,
        target_segment,
        support_wavelength_um,
        support_fnu_jy,
        support_segment,
        support_side,
    )
    delta = float(relative_log_gain_delta)
    if not math.isfinite(delta) or abs(delta) > 1.0:
        raise MeasurementContractError("Relative log gain is non-finite or outside screen range")
    gains = {"G235H": math.exp(-0.5 * delta), "G395H": math.exp(0.5 * delta)}
    target_scaled = target_flux * np.asarray([gains[item] for item in target_groups])
    support_scaled = support_flux * np.asarray([gains[item] for item in support_groups])
    continuum, tau, neutral_diagnostics = compute_fnu_continuum_and_tau(
        target_wave,
        target_scaled,
        target_groups,
        support_wave,
        support_scaled,
        support_groups,
        sides,
    )
    return continuum, tau, {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "flux_representation": "F_nu_Jy",
        "application_side": "observation_only_never_model",
        "relative_log_gain_delta_ln_g395_over_g235": delta,
        "segment_gains": gains,
        "gauge": "zero_sum_log_gain_g235_exp_minus_delta_over_2_g395_exp_plus_delta_over_2",
        "canonical_delta": 0.0,
        "nonzero_delta_role": "predeclared_robustness_scenario_not_fit_or_optimized",
        "nuisance_application": (
            "multiply_all_targets_and_six_supports_by_parent_segment_gain_then_"
            "refit_scale_invariant_huber_continuum_and_recompute_tau"
        ),
        "neutral_operator": neutral_diagnostics,
    }
