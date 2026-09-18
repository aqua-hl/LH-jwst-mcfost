#!/usr/bin/env python3
"""Method-2 signed-plane0 aperture kernel with Monte-Carlo diagnostics.

Version 2 deliberately changes only the two checks that were demonstrated to
be packet-noise diagnostics in faint near-IR cells:

* signed-versus-zero-clipped aperture disagreement; and
* native-pixel component negativity.

Neither diagnostic can accept or reject a v2 realization.  The signed,
full-Mueller total-I plane 0 remains the sole estimator.  All finite-flux,
coeval closure, geometry, aperture support, edge-loss, PSF-convolution-loss,
component-finiteness, and component-aperture-fraction checks remain gates.

The numerical implementation delegates to the hash-pinned v1 kernel while a
process-local lock temporarily removes only those two diagnostic thresholds.
The v1 globals are restored even on failure, and the returned record is then
revalidated under the explicit v2 contract.  This avoids duplicating the
expensive aperture/convolution algorithm while preserving its bitwise path.
"""

from __future__ import annotations

import contextlib
import copy
import json
import math
import os
import sys
import tempfile
import threading
from pathlib import Path
from typing import Any, Iterator, Mapping

import strict_continuum_measurement_kernel_v1 as v1


KERNEL_ID = "strict_continuum_method2_signed_aperture_v2"
SCHEMA_VERSION = 2
BASE_KERNEL_SHA256 = (
    "8cfb9aa37e1a4eb6410d84ee2b66b58a00d574f43d94720ab0e513d69a4d951a"
)
BASE_KERNEL_PATH = Path(v1.__file__).resolve()

PHYSICAL_CHECK_KEYS = (
    "finite_positive_plane0_flux",
    "coeval_plane0_abs_relative_difference_le_1e-5",
    "edge_positive_fraction_lt_1e-3",
    "convolution_flux_loss_lt_1e-4",
    "finite_positive_signed_aperture",
    "aperture_plus_6sigma_support",
)
_RELAXATION_LOCK = threading.Lock()

MeasurementContractError = v1.MeasurementContractError
MeasurementGateError = v1.MeasurementGateError


def fixed_policy_v2() -> dict[str, Any]:
    """Return a fresh JSON-safe copy of the immutable v2 policy."""
    return {
        "model_distance_pc": v1.MODEL_DISTANCE_PC,
        "target_distance_pc": v1.TARGET_DISTANCE_PC,
        "distance_flux_factor_140pc_to_147pc": v1.DISTANCE_FACTOR,
        "aperture_radius_arcsec": v1.APERTURE_RADIUS_ARCSEC,
        "aperture_subpixels_per_axis": v1.APERTURE_SUBPIXELS,
        "psf_truncation_sigma": v1.PSF_TRUNCATION_SIGMA,
        "significant_component_coeval_relative_difference_max": v1.COEVAL_TOLERANCE,
        "component_significant_fraction": v1.COMPONENT_SIGNIFICANT_FRACTION,
        "component_active_fraction": v1.COMPONENT_ACTIVE_FRACTION,
        "positive_edge_fraction_max_exclusive": v1.EDGE_FRACTION_MAX,
        "psf_convolution_flux_loss_fraction_max_exclusive": v1.CONVOLUTION_LOSS_MAX,
        "component_aperture_fraction_min_inclusive": v1.COMPONENT_APERTURE_FRACTION_MIN,
        "component_aperture_fraction_max_inclusive": v1.COMPONENT_APERTURE_FRACTION_MAX,
        "canonical_plane_zero_based": 0,
        "canonical_estimator": "signed_full_mueller_total_i_plane0",
        "signed_flux_used_for_adjudication": True,
        "zero_clipped_flux_used_for_adjudication": False,
        "componentwise_flux_used_for_adjudication": False,
        "full_plane_negative_fraction_role": "diagnostic_only",
        "native_pixel_component_negative_fraction_role": "diagnostic_only",
        "signed_vs_zero_clipped_aperture_difference_role": "diagnostic_only",
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "physical_component_integrity_gates_applied": True,
        "base_numerical_kernel": {
            "path": "jwst_mcfost/strict_continuum_measurement_kernel_v1.py",
            "sha256": BASE_KERNEL_SHA256,
        },
    }


def _require_base_kernel() -> None:
    actual = v1.sha256(BASE_KERNEL_PATH)
    if actual != BASE_KERNEL_SHA256:
        raise MeasurementContractError(
            f"Base v1 kernel SHA-256 changed: {actual}; expected {BASE_KERNEL_SHA256}"
        )


@contextlib.contextmanager
def _diagnostic_only_relaxation() -> Iterator[None]:
    """Temporarily remove exactly the two v1 diagnostic thresholds."""
    with _RELAXATION_LOCK:
        _require_base_kernel()
        before = {
            "APERTURE_CLIP_MAX": v1.APERTURE_CLIP_MAX,
            "NEGATIVE_FRACTION_MAX": v1.NEGATIVE_FRACTION_MAX,
        }
        if before != {"APERTURE_CLIP_MAX": 1.0e-3, "NEGATIVE_FRACTION_MAX": 1.0e-3}:
            raise MeasurementContractError("Base v1 diagnostic constants changed")
        # The v1 serializer forbids non-finite JSON values.  The largest
        # finite binary64 value is therefore used as the implementation-only
        # ceiling; v2 applies no scientific acceptance threshold to either
        # diagnostic.
        v1.APERTURE_CLIP_MAX = sys.float_info.max
        v1.NEGATIVE_FRACTION_MAX = sys.float_info.max
        try:
            yield
        finally:
            v1.APERTURE_CLIP_MAX = before["APERTURE_CLIP_MAX"]
            v1.NEGATIVE_FRACTION_MAX = before["NEGATIVE_FRACTION_MAX"]
            if (
                v1.APERTURE_CLIP_MAX != 1.0e-3
                or v1.NEGATIVE_FRACTION_MAX != 1.0e-3
            ):
                raise RuntimeError("Failed to restore v1 measurement constants")


def _component_physical_pass(row: Mapping[str, Any]) -> bool:
    active = row.get("active_ge_1pct_component_l1") is True
    significant = row.get("significant_ge_1e-8_total_scale") is True
    finite = row.get("finite_map_contract_pass") is True
    zero = row.get("zero_map_requires_zero_sed_component_pass") is True
    aperture = row.get("aperture_fraction_contract_pass") is True
    edge = float(row.get("edge_positive_flux_fraction", math.inf))
    loss = float(row.get("convolution_flux_loss_fraction", math.inf))
    coeval = abs(float(row.get("coeval_sed_minus_map_over_map", math.inf)))
    return bool(
        finite
        and zero
        and aperture
        and ((not active) or (0.0 <= edge < v1.EDGE_FRACTION_MAX))
        and ((not active) or (0.0 <= loss < v1.CONVOLUTION_LOSS_MAX))
        and ((not significant) or coeval <= v1.COEVAL_TOLERANCE)
    )


def _convert_relaxed_v1_record(record: Mapping[str, Any]) -> dict[str, Any]:
    converted = copy.deepcopy(dict(record))
    measurement = converted.get("measurement")
    if not isinstance(measurement, dict):
        raise MeasurementContractError("Base measurement record is absent")
    checks = measurement.get("quality_checks")
    if not isinstance(checks, Mapping):
        raise MeasurementContractError("Base quality checks are absent")
    physical_checks = {name: checks.get(name) is True for name in PHYSICAL_CHECK_KEYS}
    if not all(physical_checks.values()):
        failed = tuple(sorted(name for name, passed in physical_checks.items() if not passed))
        raise MeasurementGateError(
            f"V2 physical quality gate failed: {list(failed)}",
            failed_checks=failed,
        )

    integrity = measurement.get("intrinsic_component_integrity")
    if not isinstance(integrity, dict) or not isinstance(
        integrity.get("components"), dict
    ):
        raise MeasurementContractError("Base component integrity is absent")
    component_pass: dict[str, bool] = {}
    component_negative: dict[str, float] = {}
    for name, row_value in integrity["components"].items():
        if not isinstance(row_value, dict):
            raise MeasurementContractError(f"Malformed component row: {name}")
        passed = _component_physical_pass(row_value)
        component_pass[str(name)] = passed
        component_negative[str(name)] = float(row_value["negative_flux_fraction"])
        row_value["physical_integrity_pass"] = passed
        row_value["native_pixel_negative_fraction_role"] = "diagnostic_only"
        row_value["quality_pass"] = passed
    if not all(component_pass.values()):
        failed = tuple(
            sorted(f"component_{name}" for name, passed in component_pass.items() if not passed)
        )
        raise MeasurementGateError(
            f"V2 component physical gate failed: {list(failed)}",
            failed_checks=failed,
        )

    dclip = float(measurement["aperture_clip_symmetric_fraction"])
    full_negative = float(measurement["full_plane_negative_fraction_diagnostic"])
    if not (
        math.isfinite(dclip)
        and dclip >= 0.0
        and math.isfinite(full_negative)
        and full_negative >= 0.0
        and all(math.isfinite(value) and value >= 0.0 for value in component_negative.values())
    ):
        raise MeasurementContractError("V2 diagnostic values are not finite/non-negative")

    integrity["intrinsic_component_product_integrity_gates_applied"] = False
    integrity["physical_component_integrity_gates_applied"] = True
    integrity["native_pixel_component_negativity_gate_applied"] = False
    integrity["native_pixel_component_negativity_role"] = "diagnostic_only"
    integrity["component_physical_quality_checks"] = component_pass
    integrity["quality_pass"] = True
    measurement["schema_version"] = SCHEMA_VERSION
    measurement["kernel_id"] = KERNEL_ID
    measurement["canonical_estimator"] = "signed_full_mueller_method2_total_i_plane0"
    measurement["policy"] = fixed_policy_v2()
    measurement["quality_checks"] = physical_checks
    measurement["diagnostics"] = {
        "signed_vs_zero_clipped_aperture_symmetric_fraction": dclip,
        "signed_vs_zero_clipped_role": "diagnostic_only",
        "zero_clipped_flux_admissible_for_adjudication": False,
        "full_plane_native_pixel_negative_fraction": full_negative,
        "full_plane_native_pixel_negative_fraction_role": "diagnostic_only",
        "component_native_pixel_negative_fraction": component_negative,
        "component_native_pixel_negative_fraction_role": "diagnostic_only",
    }
    converted["schema_version"] = SCHEMA_VERSION
    converted["kernel_id"] = KERNEL_ID
    validate_task_measurement_record_v2(converted)
    return converted


def measure_method2_fits_task_v2(**kwargs: Any) -> dict[str, Any]:
    """Read and measure one product with v1 numerics and the v2 gate policy."""
    with _diagnostic_only_relaxation():
        base = v1.measure_method2_fits_task_v1(**kwargs)
    return _convert_relaxed_v1_record(base)


def validate_measurement_v2(measurement: Mapping[str, Any]) -> None:
    if not (
        isinstance(measurement, Mapping)
        and measurement.get("schema_version") == SCHEMA_VERSION
        and measurement.get("kernel_id") == KERNEL_ID
        and measurement.get("valid") is True
        and measurement.get("canonical_estimator")
        == "signed_full_mueller_method2_total_i_plane0"
        and measurement.get("policy") == fixed_policy_v2()
    ):
        raise MeasurementContractError("V2 measurement identity/policy changed")
    checks = measurement.get("quality_checks")
    if not (
        isinstance(checks, Mapping)
        and set(checks) == set(PHYSICAL_CHECK_KEYS)
        and all(value is True for value in checks.values())
    ):
        raise MeasurementContractError("V2 physical quality-check table changed")
    signed140 = float(measurement.get("signed_aperture_flux_140pc_w_m2", math.nan))
    signed147 = float(measurement.get("signed_aperture_flux_147pc_w_m2", math.nan))
    clipped = float(measurement.get("zero_clipped_aperture_flux_140pc_w_m2", math.nan))
    dclip = float(measurement.get("aperture_clip_symmetric_fraction", math.nan))
    if not (
        math.isfinite(signed140)
        and signed140 > 0.0
        and math.isclose(signed147, signed140 * v1.DISTANCE_FACTOR, rel_tol=1.0e-14)
        and math.isfinite(clipped)
        and math.isclose(dclip, v1.symmetric_fraction(signed140, clipped), rel_tol=1.0e-14)
        and dclip >= 0.0
        and abs(float(measurement.get("coeval_plane0_relative_difference", math.inf)))
        <= v1.COEVAL_TOLERANCE
        and 0.0 <= float(measurement.get("plane0_edge_positive_flux_fraction", math.inf))
        < v1.EDGE_FRACTION_MAX
        and 0.0
        <= float(measurement.get("plane0_convolution_flux_loss_fraction", math.inf))
        < v1.CONVOLUTION_LOSS_MAX
    ):
        raise MeasurementContractError("V2 measurement scalar/physical gates changed")
    diagnostics = measurement.get("diagnostics")
    if not (
        isinstance(diagnostics, Mapping)
        and diagnostics.get("signed_vs_zero_clipped_aperture_symmetric_fraction") == dclip
        and diagnostics.get("signed_vs_zero_clipped_role") == "diagnostic_only"
        and diagnostics.get("zero_clipped_flux_admissible_for_adjudication") is False
        and diagnostics.get("component_native_pixel_negative_fraction_role")
        == "diagnostic_only"
    ):
        raise MeasurementContractError("V2 diagnostic role changed")
    integrity = measurement.get("intrinsic_component_integrity")
    components = integrity.get("components") if isinstance(integrity, Mapping) else None
    if not (
        isinstance(integrity, Mapping)
        and integrity.get("physical_component_integrity_gates_applied") is True
        and integrity.get("native_pixel_component_negativity_gate_applied") is False
        and integrity.get("aggregate_componentwise_flux_gate_applied") is False
        and integrity.get("quality_pass") is True
        and isinstance(components, Mapping)
        and set(components) == set(v1.COMPONENTS.values())
        and all(
            isinstance(row, Mapping)
            and row.get("physical_integrity_pass") is True
            and row.get("native_pixel_negative_fraction_role") == "diagnostic_only"
            and _component_physical_pass(row)
            for row in components.values()
        )
    ):
        raise MeasurementContractError("V2 component physical integrity changed")
    geometry = measurement.get("geometry")
    if not (
        isinstance(geometry, Mapping)
        and int(geometry.get("nx", 0)) > 0
        and int(geometry.get("ny", 0)) > 0
        and float(geometry.get("available_half_field_arcsec", -math.inf))
        >= float(geometry.get("required_half_field_arcsec", math.inf))
    ):
        raise MeasurementContractError("V2 geometry/support changed")
    v1._require_finite_json_tree(measurement)


def validate_task_measurement_record_v2(record: Mapping[str, Any]) -> None:
    if not (
        isinstance(record, Mapping)
        and record.get("schema_version") == SCHEMA_VERSION
        and record.get("kernel_id") == KERNEL_ID
        and record.get("valid") is True
        and isinstance(record.get("task_identity"), Mapping)
        and bool(record["task_identity"])
    ):
        raise MeasurementContractError("V2 task-measurement identity changed")
    inputs = record.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != {"sed_rt", "image_rt", "closure"}:
        raise MeasurementContractError("V2 task-measurement input table changed")
    for name, item in inputs.items():
        if not (
            isinstance(item, Mapping)
            and isinstance(item.get("path"), str)
            and item.get("path")
            and isinstance(item.get("sha256"), str)
            and len(item["sha256"]) == 64
        ):
            raise MeasurementContractError(f"Malformed v2 task input: {name}")
    validate_measurement_v2(record.get("measurement"))
    v1._require_finite_json_tree(record)


def validate_task_measurement_input_files_v2(record: Mapping[str, Any]) -> None:
    validate_task_measurement_record_v2(record)
    inputs = record["inputs"]
    for name in ("sed_rt", "image_rt", "closure"):
        item = inputs[name]
        actual = v1.require_regular(Path(item["path"]), f"v2 task input {name}")
        if actual != item["sha256"]:
            raise MeasurementContractError(f"V2 task input {name} SHA-256 changed")
    closure = v1._load_json_object(Path(inputs["closure"]["path"]), "v2 closure")
    v1._validate_external_closure_identity(
        closure,
        sed_path=Path(inputs["sed_rt"]["path"]),
        image_path=Path(inputs["image_rt"]["path"]),
        sed_sha256=inputs["sed_rt"]["sha256"],
        image_sha256=inputs["image_rt"]["sha256"],
    )


def atomic_write_task_measurement_json_v2(
    path: Path, record: Mapping[str, Any], *, replace_existing: bool = False
) -> str:
    validate_task_measurement_input_files_v2(record)
    path = Path(path)
    if not replace_existing and (path.exists() or path.is_symlink()):
        raise MeasurementContractError(f"Refusing pre-existing v2 measurement: {path}")
    payload = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    decoded = json.loads(payload)
    validate_task_measurement_record_v2(decoded)
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
        temporary.unlink(missing_ok=True)
        raise
    return v1.sha256(path)


def load_task_measurement_json_v2(
    path: Path, *, expected_sha256: str | None = None, verify_input_files: bool = True
) -> dict[str, Any]:
    path = Path(path)
    actual = v1.require_regular(path, "v2 task measurement JSON")
    if expected_sha256 is not None and actual != expected_sha256:
        raise MeasurementContractError(
            f"V2 task measurement SHA-256 changed: {actual}; expected {expected_sha256}"
        )
    record = v1._load_json_object(path, "v2 task measurement JSON")
    validate_task_measurement_record_v2(record)
    if verify_input_files:
        validate_task_measurement_input_files_v2(record)
    if v1.require_regular(path, "post-read v2 measurement JSON") != actual:
        raise MeasurementContractError("V2 task measurement changed while being read")
    return record
