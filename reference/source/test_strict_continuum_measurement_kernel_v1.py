#!/usr/bin/env python3
"""Regression, gate-boundary, and serialization tests for the v1 kernel."""

from __future__ import annotations

import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

import analyze_strict_continuum_refinement_v1 as legacy
import parity_strict_continuum_measurement_kernel_v1 as parity
import strict_continuum_measurement_kernel_v1 as kernel


def synthetic_inputs(*, plane0_negative: float = 0.0) -> dict[str, object]:
    ny = nx = 161
    image = np.zeros((8, ny, nx), dtype=np.float32)
    centre = nx // 2
    image[0, centre, centre] = 1.0
    if plane0_negative:
        image[0, centre, centre + 20] = -float(plane0_negative)
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
    return {
        "sed_flux_140pc_w_m2": sed,
        "image_planes_140pc_w_m2_pixel": image,
        "closure": parity.closure_for_arrays(sed, image, wavelength),
        "wavelength_um": wavelength,
        "instrument": "NIRSpec",
        "cdelt1_deg": -cdelt,
        "cdelt2_deg": cdelt,
        "crpix1_fits": 81.0,
        "crpix2_fits": 81.0,
        "expected_cdelt_deg": cdelt,
        "expected_psf_fwhm_arcsec": 0.033 * wavelength,
    }


class MeasurementKernelV1Tests(unittest.TestCase):
    def test_scientific_constants_are_exact_legacy_values(self) -> None:
        names = (
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
        for name in names:
            self.assertEqual(getattr(kernel, name), getattr(legacy, name), name)
        self.assertEqual(kernel.COMPONENTS, legacy.COMPONENTS)
        self.assertEqual(kernel.CLOSURE_LABELS, legacy.CLOSURE_LABELS)

    def test_fractional_aperture_is_bitwise_legacy_equal(self) -> None:
        old = legacy.fractional_circle_weights((47, 53), 25.2, 22.7, 14.35)
        new = kernel.fractional_circle_weights((47, 53), 25.2, 22.7, 14.35)
        self.assertTrue(np.array_equal(old, new))

    def test_deterministic_old_new_structural_parity(self) -> None:
        result = parity.synthetic_case()
        self.assertEqual(result["status"], "parity_pass")
        self.assertEqual(result["differences"], [])
        self.assertEqual(result["legacy_outcome"], "pass")
        self.assertEqual(result["new_outcome"], "pass")

        negative = parity.synthetic_negative_rejection_case()
        self.assertEqual(negative["status"], "parity_pass")
        self.assertEqual(negative["differences"], [])
        self.assertEqual(negative["legacy_outcome"], "reject")
        self.assertEqual(negative["new_outcome"], "reject")

    def test_pure_kernel_passes_and_uses_signed_plane0(self) -> None:
        measurement = kernel.measure_method2_signed_aperture_v1(
            **synthetic_inputs()
        )
        kernel.validate_measurement_v1(measurement)
        self.assertTrue(measurement["valid"])
        self.assertGreater(measurement["signed_aperture_flux_140pc_w_m2"], 0.0)
        self.assertEqual(measurement["aperture_clip_symmetric_fraction"], 0.0)
        self.assertAlmostEqual(
            measurement["signed_aperture_flux_147pc_w_m2"],
            measurement["signed_aperture_flux_140pc_w_m2"]
            * (140.0 / 147.0) ** 2,
            places=15,
        )
        self.assertEqual(
            measurement["geometry"]["centre_x_zero_based"],
            measurement["geometry"]["crpix1_fits"] - 1.0,
        )
        self.assertEqual(
            measurement["geometry"]["centre_y_zero_based"],
            measurement["geometry"]["crpix2_fits"] - 1.0,
        )

    def test_negative_plane0_clip_gate_rejects(self) -> None:
        with self.assertRaises(kernel.MeasurementGateError) as context:
            kernel.measure_method2_signed_aperture_v1(
                **synthetic_inputs(plane0_negative=0.01)
            )
        self.assertIn(
            "direct_aperture_clip_symmetric_fraction_le_1e-3",
            context.exception.failed_checks,
        )

    def test_active_component_negative_gate_rejects(self) -> None:
        values = synthetic_inputs()
        image = np.array(values["image_planes_140pc_w_m2_pixel"], copy=True)
        image[4, 80, 100] = -0.001
        sed = np.sum(image, axis=(-2, -1), dtype=np.float64)
        values["image_planes_140pc_w_m2_pixel"] = image
        values["sed_flux_140pc_w_m2"] = sed
        values["closure"] = parity.closure_for_arrays(sed, image, 1.5)
        with self.assertRaises(kernel.MeasurementGateError) as context:
            kernel.measure_method2_signed_aperture_v1(**values)
        self.assertIn("component_direct_star", context.exception.failed_checks)

    def test_closure_report_tamper_rejects(self) -> None:
        values = synthetic_inputs()
        closure = copy.deepcopy(values["closure"])
        closure["planes"][0]["custom_sed_flux"] *= 1.01
        values["closure"] = closure
        with self.assertRaises(kernel.MeasurementContractError):
            kernel.measure_method2_signed_aperture_v1(**values)

    def test_serialized_gate_boundaries_are_inclusive_and_exclusive(self) -> None:
        measurement = kernel.measure_method2_signed_aperture_v1(
            **synthetic_inputs()
        )
        at_clip = copy.deepcopy(measurement)
        signed = float(at_clip["signed_aperture_flux_140pc_w_m2"])
        ratio = (2.0 + kernel.APERTURE_CLIP_MAX) / (
            2.0 - kernel.APERTURE_CLIP_MAX
        )
        at_clip["zero_clipped_aperture_flux_140pc_w_m2"] = signed * ratio
        at_clip["aperture_clip_symmetric_fraction"] = kernel.symmetric_fraction(
            signed, signed * ratio
        )
        kernel.validate_measurement_v1(at_clip)

        above_clip = copy.deepcopy(at_clip)
        above_clip["aperture_clip_symmetric_fraction"] = math.nextafter(
            kernel.APERTURE_CLIP_MAX, math.inf
        )
        with self.assertRaises(kernel.MeasurementContractError):
            kernel.validate_measurement_v1(above_clip)

        at_edge = copy.deepcopy(measurement)
        at_edge["plane0_edge_positive_flux_fraction"] = kernel.EDGE_FRACTION_MAX
        with self.assertRaises(kernel.MeasurementContractError):
            kernel.validate_measurement_v1(at_edge)

    def test_fits_adapter_and_atomic_hash_pinned_roundtrip(self) -> None:
        values = synthetic_inputs()
        sed = np.asarray(values["sed_flux_140pc_w_m2"], dtype=np.float64)
        image = np.asarray(
            values["image_planes_140pc_w_m2_pixel"], dtype=np.float32
        )
        wavelength = float(values["wavelength_um"])
        with tempfile.TemporaryDirectory(prefix="strict-kernel-v1-") as directory:
            root = Path(directory)
            sed_path = root / "sed_rt.fits.gz"
            image_path = root / "RT.fits.gz"
            closure_path = root / "closure.json"
            output_path = root / "measurement.json"

            sed_hdu = fits.PrimaryHDU(data=sed.reshape(8, 1, 1, 1))
            sed_hdu.header["BUNIT"] = "W.m-2"
            wave_hdu = fits.ImageHDU(data=np.asarray([wavelength], dtype=np.float64))
            fits.HDUList([sed_hdu, wave_hdu]).writeto(sed_path)

            map_hdu = fits.PrimaryHDU(data=image[:, None, None])
            map_hdu.header["BUNIT"] = "W.m-2.pixel-1"
            map_hdu.header["WAVE"] = np.float32(wavelength)
            map_hdu.header["CDELT1"] = values["cdelt1_deg"]
            map_hdu.header["CDELT2"] = values["cdelt2_deg"]
            map_hdu.header["CRPIX1"] = values["crpix1_fits"]
            map_hdu.header["CRPIX2"] = values["crpix2_fits"]
            map_hdu.header["FLUX_1"] = "I = total flux"
            map_hdu.header["FLUX_5"] = "direct star light"
            map_hdu.header["FLUX_6"] = "scattered star light"
            map_hdu.header["FLUX_7"] = "direct thermal emission"
            map_hdu.header["FLUX_8"] = "scattered thermal emission"
            map_hdu.writeto(image_path)

            closure = copy.deepcopy(values["closure"])
            closure["inputs"] = {
                "custom_sed": str(sed_path.resolve()),
                "custom_sed_sha256": kernel.sha256(sed_path),
                "image": str(image_path.resolve()),
                "image_sha256": kernel.sha256(image_path),
                "stock_sed": None,
                "stock_sed_sha256": None,
            }
            closure_path.write_text(
                json.dumps(closure, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            record = kernel.measure_method2_fits_task_v1(
                sed_path=sed_path,
                image_path=image_path,
                closure_path=closure_path,
                wavelength_um=wavelength,
                instrument="NIRSpec",
                expected_nx=161,
                expected_ny=161,
                expected_cdelt_deg=float(values["expected_cdelt_deg"]),
                expected_psf_fwhm_arcsec=float(
                    values["expected_psf_fwhm_arcsec"]
                ),
                task_identity={"task_index": 7, "anchor_tag": "w001"},
            )
            digest = kernel.atomic_write_task_measurement_json_v1(
                output_path, record
            )
            loaded = kernel.load_task_measurement_json_v1(
                output_path, expected_sha256=digest
            )
            self.assertEqual(loaded, record)
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.atomic_write_task_measurement_json_v1(output_path, record)

            output_path.write_text(output_path.read_text() + " ", encoding="utf-8")
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.load_task_measurement_json_v1(
                    output_path, expected_sha256=digest
                )
            digest = kernel.atomic_write_task_measurement_json_v1(
                output_path, record, replace_existing=True
            )
            with image_path.open("ab") as stream:
                stream.write(b"changed-after-measurement")
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.load_task_measurement_json_v1(
                    output_path, expected_sha256=digest
                )

    def test_invalid_record_is_not_published(self) -> None:
        measurement = kernel.measure_method2_signed_aperture_v1(
            **synthetic_inputs()
        )
        record = {
            "schema_version": 1,
            "kernel_id": kernel.KERNEL_ID,
            "valid": True,
            "task_identity": {"task_index": 0},
            "inputs": {
                name: {"path": f"/{name}", "sha256": "0" * 64}
                for name in ("sed_rt", "image_rt", "closure")
            },
            "measurement": measurement,
        }
        record["measurement"]["signed_aperture_flux_147pc_w_m2"] *= 1.1
        with tempfile.TemporaryDirectory(prefix="strict-kernel-invalid-") as directory:
            output = Path(directory) / "measurement.json"
            with self.assertRaises(kernel.MeasurementContractError):
                kernel.atomic_write_task_measurement_json_v1(output, record)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
