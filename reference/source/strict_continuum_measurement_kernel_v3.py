#!/usr/bin/env python3
"""Secure filesystem adapter for the frozen Method-2 v2 measurement.

The v1 and v2 kernels predate a strict same-byte input contract: they hash a
pathname and later ask JSON/Astropy to reopen that pathname.  A pathname can
therefore be rebound between those operations.  This adapter deliberately
does not change the v2 numerical estimator or any v2 scientific gate.  It
changes only how the three input artifacts are acquired:

* lstat the requested pathname and require a non-empty regular file;
* open the final component with ``O_NOFOLLOW`` where the platform provides it;
* bind lstat/open/fstat identity and metadata;
* read the complete file twice through that same descriptor and require the
  two byte sequences to be identical;
* require unchanged fstat and pathname identity around both reads;
* hash those exact immutable bytes; and
* parse JSON/FITS only from those same bytes.

Gzip-compressed FITS files are decompressed from the captured bytes and then
opened through ``io.BytesIO(...), memmap=False``.  At no point does Astropy
receive an input pathname.  The nested measurement remains an exact v2
measurement record; the outer v3 envelope identifies and audits the secure
snapshot adapter.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import os
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np
from astropy.io import fits

import strict_continuum_measurement_kernel_v1 as v1
import strict_continuum_measurement_kernel_v2 as v2


KERNEL_ID = "strict_continuum_method2_secure_snapshot_adapter_v3"
SCHEMA_VERSION = 3
SNAPSHOT_CONTRACT_ID = "stable_regular_file_double_read_byte_snapshot_v2"
V1_KERNEL_SHA256 = (
    "8cfb9aa37e1a4eb6410d84ee2b66b58a00d574f43d94720ab0e513d69a4d951a"
)
V2_KERNEL_SHA256 = (
    "89ccb15d8712785eb8cb8c598a277ab71deec698e153c04c44d9572e829ed1af"
)
V1_KERNEL_PATH = Path(v1.__file__).resolve()
V2_KERNEL_PATH = Path(v2.__file__).resolve()
READ_BLOCK_BYTES = 1024 * 1024
JSON_MAX_BYTES = 16 * 1024 * 1024
SED_FITS_MAX_UNCOMPRESSED_BYTES = 64 * 1024 * 1024

_V2_MEASUREMENT_KEYS = frozenset(
    {
        "schema_version",
        "kernel_id",
        "valid",
        "canonical_estimator",
        "policy",
        "wavelength_um",
        "instrument",
        "signed_aperture_flux_140pc_w_m2",
        "signed_aperture_flux_147pc_w_m2",
        "zero_clipped_aperture_flux_140pc_w_m2",
        "aperture_clip_symmetric_fraction",
        "aperture_clip_evidence_kind",
        "full_plane_negative_fraction_diagnostic",
        "plane0_edge_positive_flux_fraction",
        "plane0_convolution_flux_loss_fraction",
        "coeval_plane0_relative_difference",
        "coeval_max_significant_relative_difference",
        "intrinsic_component_integrity",
        "geometry",
        "quality_checks",
        "diagnostics",
    }
)
_V2_INTEGRITY_KEYS = frozenset(
    {
        "evidence_kind",
        "intrinsic_component_product_integrity_gates_applied",
        "physical_component_integrity_gates_applied",
        "native_pixel_component_negativity_gate_applied",
        "native_pixel_component_negativity_role",
        "plane0_vs_sum_planes4to7_gate_applied",
        "aggregate_componentwise_flux_gate_applied",
        "significant_fraction_of_total_scale",
        "active_fraction_of_component_l1",
        "component_l1_sed_flux_140pc_w_m2",
        "significant_component_count",
        "active_component_count",
        "max_abs_significant_all_plane_coeval_relative_difference",
        "max_abs_significant_component_coeval_relative_difference",
        "max_active_component_negative_flux_fraction",
        "max_active_component_edge_positive_flux_fraction",
        "max_active_component_convolution_flux_loss_fraction",
        "components",
        "component_physical_quality_checks",
        "quality_pass",
    }
)
_V2_COMPONENT_ROW_KEYS = frozenset(
    {
        "plane_zero_based",
        "significant_ge_1e-8_total_scale",
        "active_ge_1pct_component_l1",
        "sed_flux_140pc_w_m2",
        "map_total_fixed_fsum_140pc_w_m2",
        "map_total_numpy_140pc_w_m2",
        "map_min_pixel_flux_w_m2",
        "map_max_pixel_flux_w_m2",
        "map_positive_flux_140pc_w_m2",
        "map_negative_magnitude_140pc_w_m2",
        "coeval_sed_minus_map_over_map",
        "negative_flux_fraction",
        "edge_positive_flux_fraction",
        "convolution_flux_loss_fraction",
        "psf_map_aperture_flux_140pc_w_m2",
        "psf_map_aperture_fraction",
        "map_pixels_finite",
        "zero_map_requires_zero_sed_component_pass",
        "finite_map_contract_pass",
        "aperture_fraction_contract_pass",
        "active_morphology_contract_pass",
        "physical_integrity_pass",
        "native_pixel_negative_fraction_role",
        "quality_pass",
    }
)
_V2_GEOMETRY_KEYS = frozenset(
    {
        "nx",
        "ny",
        "cdelt1_deg",
        "cdelt2_deg",
        "crpix1_fits",
        "crpix2_fits",
        "centre_x_zero_based",
        "centre_y_zero_based",
        "target_pixel_arcsec",
        "aperture_radius_pixels",
        "psf_fwhm_arcsec",
        "psf_sigma_pixels",
        "available_half_field_arcsec",
        "required_half_field_arcsec",
    }
)
_V2_DIAGNOSTIC_KEYS = frozenset(
    {
        "signed_vs_zero_clipped_aperture_symmetric_fraction",
        "signed_vs_zero_clipped_role",
        "zero_clipped_flux_admissible_for_adjudication",
        "full_plane_native_pixel_negative_fraction",
        "full_plane_native_pixel_negative_fraction_role",
        "component_native_pixel_negative_fraction",
        "component_native_pixel_negative_fraction_role",
    }
)

MeasurementContractError = v1.MeasurementContractError
MeasurementGateError = v1.MeasurementGateError


def _require_exact_mapping_keys(
    value: Any, expected: frozenset[str], label: str
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != set(expected):
        raise MeasurementContractError(f"Exact v2 mapping schema changed: {label}")
    return value


def _require_exact_v2_measurement_schema(measurement: Any) -> None:
    """Forbid unreviewed extension fields around the frozen v2 semantics."""
    root = _require_exact_mapping_keys(
        measurement, _V2_MEASUREMENT_KEYS, "measurement"
    )
    integrity = _require_exact_mapping_keys(
        root.get("intrinsic_component_integrity"),
        _V2_INTEGRITY_KEYS,
        "measurement.intrinsic_component_integrity",
    )
    component_names = frozenset(v1.COMPONENTS.values())
    components = _require_exact_mapping_keys(
        integrity.get("components"),
        component_names,
        "measurement.intrinsic_component_integrity.components",
    )
    for name in sorted(component_names):
        _require_exact_mapping_keys(
            components.get(name),
            _V2_COMPONENT_ROW_KEYS,
            f"measurement.intrinsic_component_integrity.components.{name}",
        )
    _require_exact_mapping_keys(
        integrity.get("component_physical_quality_checks"),
        component_names,
        "measurement.intrinsic_component_integrity.component_physical_quality_checks",
    )
    _require_exact_mapping_keys(
        root.get("geometry"), _V2_GEOMETRY_KEYS, "measurement.geometry"
    )
    _require_exact_mapping_keys(
        root.get("quality_checks"),
        frozenset(v2.PHYSICAL_CHECK_KEYS),
        "measurement.quality_checks",
    )
    diagnostics = _require_exact_mapping_keys(
        root.get("diagnostics"), _V2_DIAGNOSTIC_KEYS, "measurement.diagnostics"
    )
    _require_exact_mapping_keys(
        diagnostics.get("component_native_pixel_negative_fraction"),
        component_names,
        "measurement.diagnostics.component_native_pixel_negative_fraction",
    )


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _stat_signature(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    """Metadata that must remain identical throughout one snapshot read."""
    return (
        int(stat.S_IFMT(value.st_mode)),
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


def _require_regular_stat(
    value: os.stat_result, *, path: Path, label: str
) -> None:
    if not stat.S_ISREG(value.st_mode) or int(value.st_size) <= 0:
        raise MeasurementContractError(
            f"Expected a non-empty regular non-symlink {label}: {path}"
        )


@dataclass(frozen=True)
class StableByteSnapshot:
    """Immutable bytes and audit metadata from one securely opened inode."""

    requested_path: Path
    resolved_path: Path
    payload: bytes
    sha256: str
    st_mode_type: int
    st_dev: int
    st_ino: int
    size_bytes: int
    st_mtime_ns: int
    st_ctime_ns: int

    def __post_init__(self) -> None:
        if not (
            isinstance(self.payload, bytes)
            and self.payload
            and self.size_bytes == len(self.payload)
            and self.sha256 == _sha256_bytes(self.payload)
            and self.st_mode_type == stat.S_IFREG
            and self.requested_path.is_absolute()
            and self.resolved_path.is_absolute()
        ):
            raise MeasurementContractError("Malformed stable byte snapshot")

    @property
    def stat_signature(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.st_mode_type,
            self.st_dev,
            self.st_ino,
            self.size_bytes,
            self.st_mtime_ns,
            self.st_ctime_ns,
        )

    def input_record(self) -> dict[str, str]:
        return {"path": str(self.resolved_path), "sha256": self.sha256}

    def audit_record(self) -> dict[str, Any]:
        return {
            "snapshot_contract_id": SNAPSHOT_CONTRACT_ID,
            "requested_path": str(self.requested_path),
            "resolved_path": str(self.resolved_path),
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "st_mode_type": self.st_mode_type,
            "st_dev": self.st_dev,
            "st_ino": self.st_ino,
            "st_mtime_ns": self.st_mtime_ns,
            "st_ctime_ns": self.st_ctime_ns,
        }

    def assert_current_path_identity(self, label: str) -> None:
        """Require the pathname still denotes the inode that was captured."""
        try:
            current = os.lstat(self.requested_path)
        except OSError as error:
            raise MeasurementContractError(
                f"Cannot re-lstat {label}: {self.requested_path}: {error}"
            ) from error
        if _stat_signature(current) != self.stat_signature:
            raise MeasurementContractError(
                f"{label} pathname/inode/metadata changed after its byte snapshot: "
                f"{self.requested_path}"
            )
        try:
            resolved = self.requested_path.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise MeasurementContractError(
                f"Cannot re-resolve {label}: {self.requested_path}: {error}"
            ) from error
        if resolved != self.resolved_path:
            raise MeasurementContractError(
                f"{label} resolved pathname changed after its byte snapshot: "
                f"{self.requested_path}"
            )
        try:
            final = os.lstat(self.requested_path)
        except OSError as error:
            raise MeasurementContractError(
                f"Cannot final-lstat {label}: {self.requested_path}: {error}"
            ) from error
        if _stat_signature(final) != self.stat_signature:
            raise MeasurementContractError(
                f"{label} changed while its resolved identity was checked: "
                f"{self.requested_path}"
            )


def _read_descriptor_pass(
    descriptor: int,
    expected_size: int,
    *,
    label: str,
    requested: Path,
    pass_label: str,
) -> bytes:
    """Read exactly one complete pass and prove EOF at the frozen size."""
    blocks: list[bytes] = []
    remaining = expected_size
    while remaining:
        try:
            block = os.read(descriptor, min(READ_BLOCK_BYTES, remaining))
        except OSError as error:
            raise MeasurementContractError(
                f"Cannot read {pass_label} of {label}: {requested}: {error}"
            ) from error
        if not block:
            raise MeasurementContractError(
                f"Unexpected EOF in {pass_label} of {label}: {requested}"
            )
        blocks.append(block)
        remaining -= len(block)
    try:
        extra = os.read(descriptor, 1)
    except OSError as error:
        raise MeasurementContractError(
            f"Cannot verify EOF after {pass_label} of {label}: "
            f"{requested}: {error}"
        ) from error
    if extra:
        raise MeasurementContractError(
            f"{label} grew during {pass_label}: {requested}"
        )
    payload = b"".join(blocks)
    if len(payload) != expected_size:
        raise MeasurementContractError(
            f"{label} {pass_label} byte count changed: {requested}"
        )
    return payload


def secure_read_bytes(
    path: Path,
    label: str,
    *,
    max_bytes: Optional[int] = None,
) -> StableByteSnapshot:
    """Return a snapshot only after two matching complete descriptor reads.

    Both passes use the same securely opened regular-file descriptor.  Exact
    byte equality and SHA-256 equality are required so correctness does not
    depend on prompt ``mtime``/``ctime`` visibility from a network filesystem.
    This is the public workflow-validation primitive.  Consumers must parse
    ``snapshot.payload`` rather than reopening ``snapshot.resolved_path``.
    """
    requested = Path(os.path.abspath(os.fspath(Path(path))))
    if max_bytes is not None and (
        isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes <= 0
    ):
        raise MeasurementContractError("max_bytes must be a positive integer")
    try:
        before_lstat = os.lstat(requested)
    except OSError as error:
        raise MeasurementContractError(
            f"Cannot lstat {label}: {requested}: {error}"
        ) from error
    _require_regular_stat(before_lstat, path=requested, label=label)
    if max_bytes is not None and before_lstat.st_size > max_bytes:
        raise MeasurementContractError(
            f"{label} exceeds the frozen byte limit: {before_lstat.st_size} > "
            f"{max_bytes}: {requested}"
        )

    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    # If a regular file is exchanged for a FIFO after lstat, O_NONBLOCK keeps
    # the open from hanging; fstat below then rejects it.
    flags |= getattr(os, "O_NONBLOCK", 0)
    try:
        descriptor = os.open(requested, flags)
    except OSError as error:
        raise MeasurementContractError(
            f"Cannot securely open {label}: {requested}: {error}"
        ) from error

    try:
        opened = os.fstat(descriptor)
        _require_regular_stat(opened, path=requested, label=label)
        if _stat_signature(opened) != _stat_signature(before_lstat):
            raise MeasurementContractError(
                f"{label} changed between lstat and open/fstat: {requested}"
            )
        expected_size = int(opened.st_size)
        payload = _read_descriptor_pass(
            descriptor,
            expected_size,
            label=label,
            requested=requested,
            pass_label="first descriptor read",
        )
        after_first_fstat = os.fstat(descriptor)
        if _stat_signature(after_first_fstat) != _stat_signature(opened):
            raise MeasurementContractError(
                f"{label} inode/size/timestamps changed during its first read: "
                f"{requested}"
            )
        try:
            os.lseek(descriptor, 0, os.SEEK_SET)
        except OSError as error:
            raise MeasurementContractError(
                f"Cannot rewind {label} for its second descriptor read: "
                f"{requested}: {error}"
            ) from error
        repeated_payload = _read_descriptor_pass(
            descriptor,
            expected_size,
            label=label,
            requested=requested,
            pass_label="second descriptor read",
        )
        # Direct equality is the authoritative no-hybrid check.  Requiring
        # equal digests as well is intentionally redundant: it makes the
        # hash-repeatability contract explicit for provenance and guards
        # against a future refactor that might replace the exact comparison.
        if repeated_payload != payload or _sha256_bytes(
            repeated_payload
        ) != _sha256_bytes(payload):
            raise MeasurementContractError(
                f"{label} was not byte-for-byte repeatable across two reads of "
                f"the same descriptor: {requested}"
            )
        after_second_fstat = os.fstat(descriptor)
        if _stat_signature(after_second_fstat) != _stat_signature(opened):
            raise MeasurementContractError(
                f"{label} inode/size/timestamps changed across its two reads: "
                f"{requested}"
            )
    finally:
        os.close(descriptor)

    try:
        after_lstat = os.lstat(requested)
    except OSError as error:
        raise MeasurementContractError(
            f"Cannot post-read lstat {label}: {requested}: {error}"
        ) from error
    if _stat_signature(after_lstat) != _stat_signature(before_lstat):
        raise MeasurementContractError(
            f"{label} pathname identity changed while being read: {requested}"
        )
    try:
        resolved = requested.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise MeasurementContractError(
            f"Cannot resolve stable {label}: {requested}: {error}"
        ) from error
    try:
        final_lstat = os.lstat(requested)
    except OSError as error:
        raise MeasurementContractError(
            f"Cannot final-lstat {label}: {requested}: {error}"
        ) from error
    if _stat_signature(final_lstat) != _stat_signature(before_lstat):
        raise MeasurementContractError(
            f"{label} changed while its resolved path was captured: {requested}"
        )

    signature = _stat_signature(opened)
    return StableByteSnapshot(
        requested_path=requested,
        resolved_path=resolved,
        payload=payload,
        sha256=_sha256_bytes(payload),
        st_mode_type=signature[0],
        st_dev=signature[1],
        st_ino=signature[2],
        size_bytes=signature[3],
        st_mtime_ns=signature[4],
        st_ctime_ns=signature[5],
    )


def secure_read_json_object(
    path: Path,
    label: str,
    *,
    expected_sha256: Optional[str] = None,
    max_bytes: int = JSON_MAX_BYTES,
) -> tuple[dict[str, Any], StableByteSnapshot]:
    """Snapshot, hash, and JSON-decode one object from the same bytes."""
    snapshot = secure_read_bytes(path, label, max_bytes=max_bytes)
    if expected_sha256 is not None and snapshot.sha256 != expected_sha256:
        raise MeasurementContractError(
            f"{label} SHA-256 changed: {snapshot.sha256}; expected "
            f"{expected_sha256}"
        )
    try:
        value = json.loads(snapshot.payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementContractError(
            f"Cannot decode {label} from its stable byte snapshot: {error}"
        ) from error
    if not isinstance(value, dict):
        raise MeasurementContractError(f"{label} root is not an object")
    snapshot.assert_current_path_identity(label)
    return value, snapshot


def fixed_snapshot_contract_v3() -> dict[str, Any]:
    return {
        "snapshot_contract_id": SNAPSHOT_CONTRACT_ID,
        "regular_nonempty_final_component_required": True,
        "symlink_final_component_rejected": True,
        "nonregular_final_component_rejected": True,
        "o_nofollow_used_where_available": True,
        "o_nonblock_used_to_prevent_fifo_race_hang_where_available": True,
        "lstat_open_fstat_identity_required": True,
        "two_identical_full_reads_of_same_descriptor_required": True,
        "descriptor_read_repeatability_comparison": (
            "exact_byte_equality_and_sha256_equality_both_required"
        ),
        "unchanged_size_mtime_ns_ctime_ns_required": True,
        "post_read_path_identity_required": True,
        "hash_basis": "exact_original_artifact_bytes_from_open_descriptor",
        "json_parse_basis": "same_exact_original_artifact_bytes",
        "fits_parse_basis": (
            "io.BytesIO_of_same_snapshot_or_gzip_decompression_derived_only_"
            "from_same_snapshot"
        ),
        "fits_memmap": False,
        "numerical_and_scientific_gate_contract": {
            "id": v2.KERNEL_ID,
            "sha256": V2_KERNEL_SHA256,
        },
    }


def _require_base_kernels() -> None:
    for path, label, expected in (
        (V1_KERNEL_PATH, "frozen v1 numerical kernel", V1_KERNEL_SHA256),
        (V2_KERNEL_PATH, "frozen v2 gate kernel", V2_KERNEL_SHA256),
    ):
        snapshot = secure_read_bytes(path, label, max_bytes=4 * 1024 * 1024)
        if snapshot.sha256 != expected:
            raise MeasurementContractError(
                f"{label} SHA-256 changed: {snapshot.sha256}; expected {expected}"
            )


def _fits_bytes(
    snapshot: StableByteSnapshot, *, label: str, max_uncompressed_bytes: int
) -> bytes:
    payload = snapshot.payload
    if payload.startswith(b"\x1f\x8b"):
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(payload), mode="rb") as stream:
                decoded = stream.read(max_uncompressed_bytes + 1)
        except (OSError, EOFError) as error:
            raise MeasurementContractError(
                f"Cannot decompress {label} stable byte snapshot: {error}"
            ) from error
        if len(decoded) > max_uncompressed_bytes:
            raise MeasurementContractError(
                f"Decompressed {label} exceeds the frozen byte limit"
            )
        return decoded
    if len(payload) > max_uncompressed_bytes:
        raise MeasurementContractError(f"{label} exceeds the frozen FITS byte limit")
    return payload


def _parse_sed_snapshot(
    snapshot: StableByteSnapshot,
) -> tuple[np.ndarray, tuple[int, ...], Any, int, str, np.ndarray]:
    payload = _fits_bytes(
        snapshot,
        label="coeval Method-2 SED",
        max_uncompressed_bytes=SED_FITS_MAX_UNCOMPRESSED_BYTES,
    )
    try:
        with fits.open(io.BytesIO(payload), memmap=False) as hdus:
            sed_raw = np.array(hdus[0].data, dtype=np.float64, copy=True)
            sed_shape = tuple(sed_raw.shape)
            sed_header = hdus[0].header.copy()
            sed_bitpix = int(sed_header.get("BITPIX", 0))
            sed_unit = str(sed_header.get("BUNIT", "")).lower().replace(" ", "")
            wave = np.asarray(hdus[1].data, dtype=np.float64).ravel().copy()
    except (OSError, IndexError, TypeError, ValueError) as error:
        raise MeasurementContractError(
            f"Cannot parse Method-2 SED from its stable byte snapshot: {error}"
        ) from error
    return sed_raw, sed_shape, sed_header, sed_bitpix, sed_unit, wave


def _parse_image_snapshot(
    snapshot: StableByteSnapshot, *, expected_nx: int, expected_ny: int
) -> tuple[np.ndarray, tuple[int, ...], Any, int, str]:
    image_data_bytes = 8 * int(expected_nx) * int(expected_ny) * 4
    limit = max(64 * 1024 * 1024, 2 * image_data_bytes + 32 * 1024 * 1024)
    payload = _fits_bytes(
        snapshot,
        label="coeval Method-2 image",
        max_uncompressed_bytes=limit,
    )
    try:
        with fits.open(io.BytesIO(payload), memmap=False) as hdus:
            image_raw = np.array(hdus[0].data, copy=True)
            image_shape = tuple(image_raw.shape)
            image_bitpix = int(hdus[0].header.get("BITPIX", 0))
            image_header = hdus[0].header.copy()
            image_unit = str(image_header.get("BUNIT", "")).lower().replace(" ", "")
    except (OSError, IndexError, TypeError, ValueError) as error:
        raise MeasurementContractError(
            f"Cannot parse Method-2 image from its stable byte snapshot: {error}"
        ) from error
    return image_raw, image_shape, image_header, image_bitpix, image_unit


def _reported_path(value: Any, *, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise MeasurementContractError(f"Closure {label} path is absent")
    path = Path(value)
    if not path.is_absolute():
        raise MeasurementContractError(f"Closure {label} path is not absolute")
    # Closure producers store Path.resolve() output.  Lexical normalization is
    # sufficient here and intentionally does not reopen or resolve a pathname.
    return Path(os.path.normpath(value))


def _validate_closure_snapshot_identity(
    closure: Mapping[str, Any],
    *,
    sed_snapshot: StableByteSnapshot,
    image_snapshot: StableByteSnapshot,
) -> None:
    inputs = closure.get("inputs")
    if not isinstance(inputs, Mapping):
        raise MeasurementContractError("Closure report input table is absent")
    if not (
        _reported_path(inputs.get("custom_sed"), label="SED")
        == sed_snapshot.resolved_path
        and _reported_path(inputs.get("image"), label="image")
        == image_snapshot.resolved_path
        and inputs.get("custom_sed_sha256") == sed_snapshot.sha256
        and inputs.get("image_sha256") == image_snapshot.sha256
        and inputs.get("stock_sed") is None
        and inputs.get("stock_sed_sha256") is None
    ):
        raise MeasurementContractError(
            "Closure report product identity differs from stable input snapshots"
        )


def _fits_gate_and_measure_v2(
    *,
    sed_snapshot: StableByteSnapshot,
    image_snapshot: StableByteSnapshot,
    closure_snapshot: StableByteSnapshot,
    closure: Mapping[str, Any],
    wavelength_um: float,
    instrument: str,
    expected_nx: int,
    expected_ny: int,
    expected_cdelt_deg: float,
    expected_psf_fwhm_arcsec: float,
    task_identity: Mapping[str, Any],
) -> dict[str, Any]:
    (
        sed_raw,
        sed_shape,
        sed_header,
        sed_bitpix,
        sed_unit,
        wave,
    ) = _parse_sed_snapshot(sed_snapshot)
    (
        image_raw,
        image_shape,
        image_header,
        image_bitpix,
        image_unit,
    ) = _parse_image_snapshot(
        image_snapshot, expected_nx=expected_nx, expected_ny=expected_ny
    )

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
        and math.isclose(
            float(wave[0]), wavelength_um, rel_tol=1.0e-6, abs_tol=1.0e-8
        )
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

    base_inputs = {
        "sed_rt": sed_snapshot.input_record(),
        "image_rt": image_snapshot.input_record(),
        "closure": closure_snapshot.input_record(),
    }
    with v2._diagnostic_only_relaxation():
        measurement = v1.measure_method2_signed_aperture_v1(
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
        base_record = {
            "schema_version": v1.SCHEMA_VERSION,
            "kernel_id": v1.KERNEL_ID,
            "valid": True,
            "task_identity": dict(task_identity),
            "inputs": base_inputs,
            "measurement": measurement,
        }
        v1.validate_task_measurement_record_v1(base_record)
    converted = v2._convert_relaxed_v1_record(base_record)
    v2.validate_task_measurement_record_v2(converted)
    return converted


def measure_method2_fits_task_v3(
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
    """Securely snapshot all inputs, then run the exact v2 measurement."""
    _require_base_kernels()
    image_data_bytes = 8 * int(expected_nx) * int(expected_ny) * 4
    image_file_limit = max(
        64 * 1024 * 1024, 2 * image_data_bytes + 32 * 1024 * 1024
    )
    snapshots = {
        "sed_rt": secure_read_bytes(
            sed_path,
            "coeval Method-2 SED",
            max_bytes=SED_FITS_MAX_UNCOMPRESSED_BYTES,
        ),
        "image_rt": secure_read_bytes(
            image_path, "coeval Method-2 image", max_bytes=image_file_limit
        ),
        "closure": secure_read_bytes(
            closure_path, "coeval closure report", max_bytes=JSON_MAX_BYTES
        ),
    }
    try:
        closure = json.loads(snapshots["closure"].payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementContractError(
            f"Cannot decode closure from its stable byte snapshot: {error}"
        ) from error
    if not isinstance(closure, dict):
        raise MeasurementContractError("Closure report root is not an object")
    _validate_closure_snapshot_identity(
        closure,
        sed_snapshot=snapshots["sed_rt"],
        image_snapshot=snapshots["image_rt"],
    )
    v2_record = _fits_gate_and_measure_v2(
        sed_snapshot=snapshots["sed_rt"],
        image_snapshot=snapshots["image_rt"],
        closure_snapshot=snapshots["closure"],
        closure=closure,
        wavelength_um=wavelength_um,
        instrument=instrument,
        expected_nx=expected_nx,
        expected_ny=expected_ny,
        expected_cdelt_deg=expected_cdelt_deg,
        expected_psf_fwhm_arcsec=expected_psf_fwhm_arcsec,
        task_identity=task_identity,
    )
    for name, snapshot in snapshots.items():
        snapshot.assert_current_path_identity(f"post-measurement input {name}")
    record = {
        "schema_version": SCHEMA_VERSION,
        "kernel_id": KERNEL_ID,
        "valid": True,
        "task_identity": dict(task_identity),
        "snapshot_contract": fixed_snapshot_contract_v3(),
        "numerical_kernel": {
            "id": v2.KERNEL_ID,
            "sha256": V2_KERNEL_SHA256,
            "measurement_schema_version": v2.SCHEMA_VERSION,
        },
        "inputs": v2_record["inputs"],
        "input_snapshots": {
            name: snapshot.audit_record() for name, snapshot in snapshots.items()
        },
        "measurement": v2_record["measurement"],
    }
    validate_task_measurement_record_v3(record)
    return record


def validate_task_measurement_record_v3(record: Mapping[str, Any]) -> None:
    expected_keys = {
        "schema_version",
        "kernel_id",
        "valid",
        "task_identity",
        "snapshot_contract",
        "numerical_kernel",
        "inputs",
        "input_snapshots",
        "measurement",
    }
    if not (
        isinstance(record, Mapping)
        and set(record) == expected_keys
        and record.get("schema_version") == SCHEMA_VERSION
        and record.get("kernel_id") == KERNEL_ID
        and record.get("valid") is True
        and isinstance(record.get("task_identity"), Mapping)
        and bool(record["task_identity"])
        and record.get("snapshot_contract") == fixed_snapshot_contract_v3()
        and record.get("numerical_kernel")
        == {
            "id": v2.KERNEL_ID,
            "sha256": V2_KERNEL_SHA256,
            "measurement_schema_version": v2.SCHEMA_VERSION,
        }
    ):
        raise MeasurementContractError("V3 secure measurement identity changed")
    inputs = record.get("inputs")
    evidence = record.get("input_snapshots")
    names = {"sed_rt", "image_rt", "closure"}
    if not (
        isinstance(inputs, Mapping)
        and set(inputs) == names
        and isinstance(evidence, Mapping)
        and set(evidence) == names
    ):
        raise MeasurementContractError("V3 secure input tables changed")
    evidence_keys = {
        "snapshot_contract_id",
        "requested_path",
        "resolved_path",
        "sha256",
        "size_bytes",
        "st_mode_type",
        "st_dev",
        "st_ino",
        "st_mtime_ns",
        "st_ctime_ns",
    }
    for name in sorted(names):
        item = inputs[name]
        row = evidence[name]
        if not (
            isinstance(item, Mapping)
            and set(item) == {"path", "sha256"}
            and isinstance(item.get("path"), str)
            and Path(item["path"]).is_absolute()
            and isinstance(item.get("sha256"), str)
            and len(item["sha256"]) == 64
            and all(character in "0123456789abcdef" for character in item["sha256"])
            and isinstance(row, Mapping)
            and set(row) == evidence_keys
            and row.get("snapshot_contract_id") == SNAPSHOT_CONTRACT_ID
            and row.get("resolved_path") == item["path"]
            and row.get("sha256") == item["sha256"]
            and isinstance(row.get("requested_path"), str)
            and Path(row["requested_path"]).is_absolute()
            and isinstance(row.get("size_bytes"), int)
            and not isinstance(row.get("size_bytes"), bool)
            and row["size_bytes"] > 0
            and row.get("st_mode_type") == stat.S_IFREG
            and all(
                isinstance(row.get(key), int) and not isinstance(row.get(key), bool)
                for key in ("st_dev", "st_ino", "st_mtime_ns", "st_ctime_ns")
            )
        ):
            raise MeasurementContractError(f"Malformed v3 secure input: {name}")
    _require_exact_v2_measurement_schema(record.get("measurement"))
    v2.validate_measurement_v2(record.get("measurement"))
    v1._require_finite_json_tree(record)


def secure_snapshot_task_inputs_v3(
    record: Mapping[str, Any],
) -> dict[str, StableByteSnapshot]:
    """Re-snapshot and hash-bind all raw inputs referenced by a v3 record."""
    validate_task_measurement_record_v3(record)
    inputs = record["inputs"]
    evidence = record["input_snapshots"]
    snapshots: dict[str, StableByteSnapshot] = {}
    for name in ("sed_rt", "image_rt", "closure"):
        item = inputs[name]
        recorded_snapshot = evidence[name]
        recorded_size = int(recorded_snapshot["size_bytes"])
        limit = (
            min(recorded_size, JSON_MAX_BYTES)
            if name == "closure"
            else recorded_size
        )
        snapshot = secure_read_bytes(
            Path(recorded_snapshot["requested_path"]),
            f"v3 task input {name}",
            max_bytes=limit,
        )
        if snapshot.input_record() != dict(item):
            raise MeasurementContractError(
                f"V3 task input {name} resolved path/SHA-256 changed"
            )
        if snapshot.audit_record() != dict(recorded_snapshot):
            raise MeasurementContractError(
                f"V3 task input {name} stable snapshot evidence changed"
            )
        snapshots[name] = snapshot
    try:
        closure = json.loads(snapshots["closure"].payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementContractError(
            f"Cannot decode v3 closure from its stable byte snapshot: {error}"
        ) from error
    if not isinstance(closure, dict):
        raise MeasurementContractError("V3 closure root is not an object")
    _validate_closure_snapshot_identity(
        closure,
        sed_snapshot=snapshots["sed_rt"],
        image_snapshot=snapshots["image_rt"],
    )
    for name, snapshot in snapshots.items():
        snapshot.assert_current_path_identity(f"verified v3 input {name}")
    return snapshots


def validate_task_measurement_input_files_v3(record: Mapping[str, Any]) -> None:
    secure_snapshot_task_inputs_v3(record)


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write_task_measurement_json_v3(
    path: Path,
    record: Mapping[str, Any],
    *,
    replace_existing: bool = False,
) -> str:
    """Publish one strict record; no-clobber publication is atomic via link."""
    validate_task_measurement_record_v3(record)
    validate_task_measurement_input_files_v3(record)
    try:
        payload = (
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
        decoded = json.loads(payload.decode("utf-8"))
    except (TypeError, ValueError, UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementContractError(
            f"V3 measurement is not strict JSON: {error}"
        ) from error
    validate_task_measurement_record_v3(decoded)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary = Path(temporary_name)
    published = False
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        if replace_existing:
            os.replace(temporary, path)
            published = True
        else:
            try:
                os.link(temporary, path, follow_symlinks=False)
            except FileExistsError as error:
                raise MeasurementContractError(
                    f"Refusing pre-existing v3 measurement: {path}"
                ) from error
            except OSError as error:
                raise MeasurementContractError(
                    f"Cannot atomically no-clobber publish v3 measurement: "
                    f"{path}: {error}"
                ) from error
            published = True
            temporary.unlink()
        _fsync_directory(path.parent)
    except BaseException:
        if not published:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        raise
    snapshot = secure_read_bytes(
        path, "published v3 measurement JSON", max_bytes=JSON_MAX_BYTES
    )
    if snapshot.payload != payload:
        raise MeasurementContractError(
            f"Published v3 measurement bytes differ from serialized bytes: {path}"
        )
    loaded = json.loads(snapshot.payload.decode("utf-8"))
    validate_task_measurement_record_v3(loaded)
    snapshot.assert_current_path_identity("published v3 measurement JSON")
    return snapshot.sha256


def load_task_measurement_json_v3(
    path: Path,
    *,
    expected_sha256: Optional[str] = None,
    verify_input_files: bool = True,
) -> dict[str, Any]:
    """Hash and decode a durable v3 record from one stable byte snapshot."""
    record, snapshot = secure_read_json_object(
        path,
        "v3 task measurement JSON",
        expected_sha256=expected_sha256,
        max_bytes=JSON_MAX_BYTES,
    )
    validate_task_measurement_record_v3(record)
    if verify_input_files:
        validate_task_measurement_input_files_v3(record)
    snapshot.assert_current_path_identity("post-validation v3 measurement JSON")
    return record
