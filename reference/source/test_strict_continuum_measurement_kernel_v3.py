#!/usr/bin/env python3
"""Security and numerical-contract tests for the secure v3 adapter."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import os
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any
from unittest import mock

import numpy as np
from astropy.io import fits

import strict_continuum_measurement_kernel_v1 as v1
import strict_continuum_measurement_kernel_v2 as v2
import strict_continuum_measurement_kernel_v3 as kernel


def _closure_for_arrays(
    sed: np.ndarray, image: np.ndarray, wavelength_um: float
) -> dict[str, Any]:
    totals = {
        plane: float(np.sum(image[plane], dtype=np.float64))
        for plane in v1.CLOSURE_LABELS
    }
    total_scale = max(abs(float(sed[0])), abs(totals[0]))
    significant_cut = v1.COMPONENT_SIGNIFICANT_FRACTION * total_scale
    rows: list[dict[str, Any]] = []
    errors: list[float] = []
    for plane, label in v1.CLOSURE_LABELS.items():
        sed_value = float(sed[plane])
        image_value = totals[plane]
        significant = (
            plane == 0
            or max(abs(sed_value), abs(image_value)) > significant_cut
        )
        relative = v1.relative_difference(sed_value, image_value)
        if significant:
            errors.append(abs(relative))
        rows.append(
            {
                "plane": plane,
                "label": label,
                "significant": significant,
                "image_flux": image_value,
                "custom_sed_flux": sed_value,
                "custom_relative_difference": relative,
            }
        )
    maximum = max(errors, default=math.inf)
    return {
        "schema_version": 1,
        "wavelength_um": wavelength_um,
        "tolerance": v1.COEVAL_TOLERANCE,
        "significant_fraction": v1.COMPONENT_SIGNIFICANT_FRACTION,
        "stock_component_closure_relative": None,
        "max_significant_custom_image_relative_difference": maximum,
        "passes_requested_closure_tolerance": maximum <= v1.COEVAL_TOLERANCE,
        "planes": rows,
    }


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_fixture(root: Path, *, scale: float = 1.0) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    ny = nx = 161
    centre = nx // 2
    wavelength = 1.5
    image = np.zeros((8, ny, nx), dtype=np.float32)
    image[0, centre, centre] = scale
    image[4:8, centre, centre] = scale / 4.0
    sed = np.sum(image, axis=(-2, -1), dtype=np.float64)
    target_pixel_arcsec = 0.02
    cdelt = (
        target_pixel_arcsec
        / 3600.0
        * v1.TARGET_DISTANCE_PC
        / v1.MODEL_DISTANCE_PC
    )
    sed_path = root / "sed_rt.fits.gz"
    image_path = root / "RT.fits.gz"
    closure_path = root / "closure.json"

    sed_hdu = fits.PrimaryHDU(data=sed.reshape(8, 1, 1, 1))
    sed_hdu.header["BUNIT"] = "W.m-2"
    wave_hdu = fits.ImageHDU(data=np.asarray([wavelength], dtype=np.float64))
    fits.HDUList([sed_hdu, wave_hdu]).writeto(sed_path)

    image_hdu = fits.PrimaryHDU(data=image[:, None, None])
    image_hdu.header["BUNIT"] = "W.m-2.pixel-1"
    image_hdu.header["WAVE"] = np.float32(wavelength)
    image_hdu.header["CDELT1"] = -cdelt
    image_hdu.header["CDELT2"] = cdelt
    image_hdu.header["CRPIX1"] = 81.0
    image_hdu.header["CRPIX2"] = 81.0
    image_hdu.header["FLUX_1"] = "I = total flux"
    image_hdu.header["FLUX_5"] = "direct star light"
    image_hdu.header["FLUX_6"] = "scattered star light"
    image_hdu.header["FLUX_7"] = "direct thermal emission"
    image_hdu.header["FLUX_8"] = "scattered thermal emission"
    image_hdu.writeto(image_path)

    closure = _closure_for_arrays(sed, image, wavelength)
    closure["inputs"] = {
        "custom_sed": str(sed_path.resolve()),
        "custom_sed_sha256": _sha(sed_path),
        "image": str(image_path.resolve()),
        "image_sha256": _sha(image_path),
        "stock_sed": None,
        "stock_sed_sha256": None,
    }
    closure_path.write_text(
        json.dumps(closure, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "sed_path": sed_path,
        "image_path": image_path,
        "closure_path": closure_path,
        "wavelength_um": wavelength,
        "instrument": "NIRSpec",
        "expected_nx": nx,
        "expected_ny": ny,
        "expected_cdelt_deg": cdelt,
        "expected_psf_fwhm_arcsec": 0.033 * wavelength,
        "task_identity": {"task_index": 7, "anchor_tag": "w001"},
        "expected_flux_140pc": scale,
    }


def _measurement_kwargs(fixture: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in fixture.items()
        if key != "expected_flux_140pc"
    }


class SecureMeasurementKernelV3Tests(unittest.TestCase):
    def test_frozen_v1_v2_bytes_are_unchanged(self) -> None:
        self.assertEqual(_sha(Path(v1.__file__)), kernel.V1_KERNEL_SHA256)
        self.assertEqual(_sha(Path(v2.__file__)), kernel.V2_KERNEL_SHA256)

    def test_exact_snapshot_hash_drives_bytesio_parse_and_v2_measurement(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-") as directory:
            fixture = _write_fixture(Path(directory) / "good")
            exact_payloads = {
                name: fixture[path_key].read_bytes()
                for name, path_key in (
                    ("sed_rt", "sed_path"),
                    ("image_rt", "image_path"),
                    ("closure", "closure_path"),
                )
            }
            real_fits_open = kernel.fits.open
            observed: list[tuple[type[Any], Any]] = []

            def recording_open(file_object: Any, *args: Any, **kwargs: Any) -> Any:
                observed.append((type(file_object), kwargs.get("memmap")))
                return real_fits_open(file_object, *args, **kwargs)

            with mock.patch.object(kernel.fits, "open", side_effect=recording_open):
                record = kernel.measure_method2_fits_task_v3(
                    **_measurement_kwargs(fixture)
                )
            legacy_v2_record = v2.measure_method2_fits_task_v2(
                **_measurement_kwargs(fixture)
            )

            kernel.validate_task_measurement_record_v3(record)
            self.assertEqual(record["kernel_id"], kernel.KERNEL_ID)
            self.assertEqual(record["measurement"]["kernel_id"], v2.KERNEL_ID)
            self.assertEqual(record["measurement"]["policy"], v2.fixed_policy_v2())
            self.assertEqual(record["measurement"], legacy_v2_record["measurement"])
            self.assertAlmostEqual(
                record["measurement"]["signed_aperture_flux_140pc_w_m2"],
                fixture["expected_flux_140pc"],
                places=15,
            )
            self.assertEqual(len(observed), 2)
            self.assertTrue(all(kind is io.BytesIO for kind, _ in observed))
            self.assertTrue(all(memmap is False for _, memmap in observed))
            for name, payload in exact_payloads.items():
                digest = hashlib.sha256(payload).hexdigest()
                self.assertEqual(record["inputs"][name]["sha256"], digest)
                self.assertEqual(record["input_snapshots"][name]["sha256"], digest)
                self.assertEqual(
                    record["input_snapshots"][name]["size_bytes"], len(payload)
                )

            output = Path(directory) / "measurement.json"
            digest = kernel.atomic_write_task_measurement_json_v3(output, record)
            self.assertEqual(
                kernel.load_task_measurement_json_v3(
                    output, expected_sha256=digest
                ),
                record,
            )
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.atomic_write_task_measurement_json_v3(output, record)

    def test_transient_inode_swap_during_parse_cannot_inject_alternate_fits(self) -> None:
        """Reproduce the old parse-time swap while isolating same-byte parsing."""
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-swap-") as directory:
            root = Path(directory)
            good = _write_fixture(root / "good", scale=1.0)
            alternate = _write_fixture(root / "alternate", scale=2.0)
            original_image_sha = _sha(good["image_path"])
            saved = root / "saved-original-image.fits.gz"
            real_fits_open = kernel.fits.open
            calls = 0

            def swapping_open(file_object: Any, *args: Any, **kwargs: Any) -> Any:
                nonlocal calls
                calls += 1
                self.assertIsInstance(file_object, io.BytesIO)
                if calls != 2:
                    return real_fits_open(file_object, *args, **kwargs)
                os.replace(good["image_path"], saved)
                os.replace(alternate["image_path"], good["image_path"])
                try:
                    return real_fits_open(file_object, *args, **kwargs)
                finally:
                    os.replace(good["image_path"], alternate["image_path"])
                    os.replace(saved, good["image_path"])

            # Disable only the redundant post-parse pathname check in this
            # test.  It isolates the key property: Astropy parses the captured
            # good bytes, never the 2x file temporarily installed at the path.
            with mock.patch.object(kernel.fits, "open", side_effect=swapping_open), mock.patch.object(
                kernel.StableByteSnapshot,
                "assert_current_path_identity",
                autospec=True,
                return_value=None,
            ):
                record = kernel.measure_method2_fits_task_v3(
                    **_measurement_kwargs(good)
                )
            self.assertEqual(calls, 2)
            self.assertEqual(record["inputs"]["image_rt"]["sha256"], original_image_sha)
            self.assertAlmostEqual(
                record["measurement"]["signed_aperture_flux_140pc_w_m2"],
                1.0,
                places=15,
            )
            self.assertNotAlmostEqual(
                record["measurement"]["signed_aperture_flux_140pc_w_m2"],
                2.0,
                places=15,
            )

    def test_transient_inode_swap_during_json_parse_keeps_original_snapshot(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-json-") as directory:
            root = Path(directory)
            path = root / "input.json"
            alternate = root / "alternate.json"
            saved = root / "saved-original.json"
            path.write_text('{"value":"original"}\n', encoding="utf-8")
            alternate.write_text('{"value":"alternate"}\n', encoding="utf-8")
            expected_sha = _sha(path)
            real_loads = kernel.json.loads

            def swapping_loads(document: Any, *args: Any, **kwargs: Any) -> Any:
                os.replace(path, saved)
                os.replace(alternate, path)
                try:
                    return real_loads(document, *args, **kwargs)
                finally:
                    os.replace(path, alternate)
                    os.replace(saved, path)

            with mock.patch.object(kernel.json, "loads", side_effect=swapping_loads), mock.patch.object(
                kernel.StableByteSnapshot,
                "assert_current_path_identity",
                autospec=True,
                return_value=None,
            ):
                value, snapshot = kernel.secure_read_json_object(path, "JSON fixture")
            self.assertEqual(value, {"value": "original"})
            self.assertEqual(snapshot.sha256, expected_sha)
            self.assertEqual(hashlib.sha256(snapshot.payload).hexdigest(), expected_sha)

    def test_concurrent_no_clobber_publish_has_exactly_one_winner(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-publish-") as directory:
            root = Path(directory)
            fixture = _write_fixture(root / "fixture")
            record = kernel.measure_method2_fits_task_v3(
                **_measurement_kwargs(fixture)
            )
            output = root / "measurement.json"
            barrier = threading.Barrier(2)

            def publish() -> tuple[str, str]:
                barrier.wait()
                try:
                    return (
                        "published",
                        kernel.atomic_write_task_measurement_json_v3(output, record),
                    )
                except kernel.MeasurementContractError as error:
                    return "refused", str(error)

            with ThreadPoolExecutor(max_workers=2) as executor:
                outcomes = list(executor.map(lambda _: publish(), range(2)))
            self.assertEqual([state for state, _ in outcomes].count("published"), 1)
            self.assertEqual([state for state, _ in outcomes].count("refused"), 1)
            winner_sha = next(value for state, value in outcomes if state == "published")
            self.assertEqual(
                kernel.load_task_measurement_json_v3(
                    output, expected_sha256=winner_sha
                ),
                record,
            )

    def test_fabricated_snapshot_metadata_cannot_be_published(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-evidence-") as directory:
            root = Path(directory)
            fixture = _write_fixture(root / "fixture")
            record = kernel.measure_method2_fits_task_v3(
                **_measurement_kwargs(fixture)
            )
            forged = copy.deepcopy(record)
            row = forged["input_snapshots"]["image_rt"]
            row["size_bytes"] += 1
            row["st_ino"] += 1
            row["requested_path"] = str(root / "fabricated-image-path.fits.gz")
            # The pure record validator can only type/cross-bind serialized
            # values.  Any operation that trusts or publishes the record must
            # re-snapshot the live files and reject this fabricated evidence.
            kernel.validate_task_measurement_record_v3(forged)
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.secure_snapshot_task_inputs_v3(forged)
            output = root / "forged-measurement.json"
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.atomic_write_task_measurement_json_v3(output, forged)
            self.assertFalse(output.exists())

    def test_unreviewed_v2_measurement_extension_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-schema-") as directory:
            root = Path(directory)
            fixture = _write_fixture(root / "fixture")
            record = kernel.measure_method2_fits_task_v3(
                **_measurement_kwargs(fixture)
            )
            forged = copy.deepcopy(record)
            forged["measurement"]["unvetted_extra"] = "finite-but-not-v2-schema"
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.validate_task_measurement_record_v3(forged)
            output = root / "extended-measurement.json"
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.atomic_write_task_measurement_json_v3(output, forged)
            self.assertFalse(output.exists())

    def test_lstat_to_open_regular_inode_swap_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-open-race-") as directory:
            root = Path(directory)
            path = root / "input.bin"
            replacement = root / "replacement.bin"
            saved = root / "saved.bin"
            path.write_bytes(b"first-regular-inode")
            replacement.write_bytes(b"second-regular-inode")
            real_open = os.open
            swapped = False

            def swapping_open(open_path: Any, flags: int, *args: Any, **kwargs: Any) -> int:
                nonlocal swapped
                if Path(open_path) == path and not swapped:
                    os.replace(path, saved)
                    os.replace(replacement, path)
                    swapped = True
                return real_open(open_path, flags, *args, **kwargs)

            try:
                with mock.patch.object(kernel.os, "open", side_effect=swapping_open):
                    with self.assertRaises(kernel.MeasurementContractError):
                        kernel.secure_read_bytes(path, "lstat/open regular swap")
            finally:
                if swapped:
                    os.replace(path, replacement)
                    os.replace(saved, path)

    def test_lstat_regular_to_open_fifo_race_is_rejected_without_hanging(self) -> None:
        if not hasattr(os, "mkfifo"):
            self.skipTest("FIFO unsupported")
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-fifo-race-") as directory:
            root = Path(directory)
            path = root / "input.bin"
            fifo = root / "replacement.fifo"
            saved = root / "saved.bin"
            path.write_bytes(b"regular-before-lstat")
            os.mkfifo(fifo)
            real_open = os.open
            swapped = False

            def swapping_open(open_path: Any, flags: int, *args: Any, **kwargs: Any) -> int:
                nonlocal swapped
                if Path(open_path) == path and not swapped:
                    os.replace(path, saved)
                    os.replace(fifo, path)
                    swapped = True
                return real_open(open_path, flags, *args, **kwargs)

            try:
                with mock.patch.object(kernel.os, "open", side_effect=swapping_open):
                    with self.assertRaises(kernel.MeasurementContractError):
                        kernel.secure_read_bytes(path, "lstat regular/open FIFO race")
            finally:
                if swapped:
                    os.replace(path, fifo)
                    os.replace(saved, path)

    def test_persistent_inode_swap_during_parse_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-rebind-") as directory:
            root = Path(directory)
            good = _write_fixture(root / "good", scale=1.0)
            alternate = _write_fixture(root / "alternate", scale=2.0)
            saved = root / "saved-original-image.fits.gz"
            real_fits_open = kernel.fits.open
            calls = 0
            swapped = False

            def swapping_open(file_object: Any, *args: Any, **kwargs: Any) -> Any:
                nonlocal calls, swapped
                calls += 1
                if calls == 2:
                    os.replace(good["image_path"], saved)
                    os.replace(alternate["image_path"], good["image_path"])
                    swapped = True
                return real_fits_open(file_object, *args, **kwargs)

            try:
                with mock.patch.object(kernel.fits, "open", side_effect=swapping_open):
                    with self.assertRaises(kernel.MeasurementContractError):
                        kernel.measure_method2_fits_task_v3(
                            **_measurement_kwargs(good)
                        )
            finally:
                if swapped:
                    os.replace(good["image_path"], alternate["image_path"])
                    os.replace(saved, good["image_path"])

    def test_stable_snapshot_requires_two_complete_descriptor_reads(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-double-read-") as directory:
            path = Path(directory) / "stable.bin"
            payload = b"stable-byte-snapshot"
            path.write_bytes(payload)
            original_read = os.read
            calls = 0
            returned_bytes = 0

            def counting_read(descriptor: int, count: int) -> bytes:
                nonlocal calls, returned_bytes
                block = original_read(descriptor, count)
                calls += 1
                returned_bytes += len(block)
                return block

            with mock.patch.object(kernel.os, "read", side_effect=counting_read):
                snapshot = kernel.secure_read_bytes(path, "stable double-read fixture")
            self.assertEqual(snapshot.payload, payload)
            self.assertEqual(returned_bytes, 2 * len(payload))
            self.assertGreaterEqual(calls, 4)  # data + EOF for both passes

    def test_same_size_already_read_mutation_rejects_with_frozen_metadata(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-mutate-") as directory:
            path = Path(directory) / "large.bin"
            path.write_bytes(b"a" * (2 * kernel.READ_BLOCK_BYTES + 17))
            original_read = os.read
            original_lstat = os.lstat
            original_fstat = os.fstat
            frozen_stat = original_lstat(path)
            mutated = False

            def mutating_read(descriptor: int, count: int) -> bytes:
                nonlocal mutated
                block = original_read(descriptor, count)
                if not mutated:
                    with path.open("r+b") as stream:
                        # This byte was already returned in pass one's first
                        # block.  A one-pass reader would preserve the old byte
                        # in its payload while the live inode contains the new
                        # byte; pass two must detect that mismatch even if NFS
                        # reports stale timestamps.
                        stream.seek(3)
                        stream.write(b"b")
                        stream.flush()
                        os.fsync(stream.fileno())
                    mutated = True
                return block

            def frozen_lstat(candidate: Any, *args: Any, **kwargs: Any) -> Any:
                value = original_lstat(candidate, *args, **kwargs)
                if Path(candidate) == path:
                    return frozen_stat
                return value

            def frozen_fstat(descriptor: int) -> Any:
                value = original_fstat(descriptor)
                if (value.st_dev, value.st_ino) == (
                    frozen_stat.st_dev,
                    frozen_stat.st_ino,
                ):
                    return frozen_stat
                return value

            with mock.patch.object(
                kernel.os, "read", side_effect=mutating_read
            ), mock.patch.object(
                kernel.os, "lstat", side_effect=frozen_lstat
            ), mock.patch.object(
                kernel.os, "fstat", side_effect=frozen_fstat
            ):
                with self.assertRaisesRegex(
                    kernel.MeasurementContractError, "not byte-for-byte repeatable"
                ):
                    kernel.secure_read_bytes(
                        path, "same-size mutation with frozen metadata"
                    )
            self.assertTrue(mutated)
            self.assertEqual(path.stat().st_size, 2 * kernel.READ_BLOCK_BYTES + 17)

    def test_symlink_and_fifo_are_rejected_without_following_or_blocking(self) -> None:
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v3-types-") as directory:
            root = Path(directory)
            target = root / "target.bin"
            target.write_bytes(b"regular")
            symlink = root / "symlink.bin"
            try:
                symlink.symlink_to(target)
            except (OSError, NotImplementedError) as error:
                self.skipTest(f"symlink unsupported: {error}")
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.secure_read_bytes(symlink, "symlink fixture")

            if not hasattr(os, "mkfifo"):
                self.skipTest("FIFO unsupported")
            fifo = root / "input.fifo"
            os.mkfifo(fifo)
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.secure_read_bytes(fifo, "FIFO fixture")


if __name__ == "__main__":
    unittest.main()
