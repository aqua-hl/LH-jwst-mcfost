"""Small analytical/numerical tests, not a scientific convergence campaign."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.photometry import (  # noqa: E402
    FWHM_TO_SIGMA, PhotometryError, SPEED_OF_LIGHT_M_S, _adjoint_weights,
    fractional_circle_weights, gaussian_psf_fwhm_arcsec, measure_image,
)


class PhotometryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.path = self.root / "RT.fits.gz"
        self.n = 129
        self.wave = 3.0
        self.pixel = 0.05
        self.plane = np.zeros((self.n, self.n))
        self.plane[self.n // 2, self.n // 2] = 2e-14

    def tearDown(self):
        self.temp.cleanup()

    def write_image(self, plane=None, updates=None, nplanes=8):
        plane = self.plane if plane is None else plane
        data = np.zeros((nplanes, 1, 1, *plane.shape), dtype=np.float64)
        data[0, 0, 0] = plane
        if nplanes > 4:
            data[4, 0, 0] = plane
        hdu = fits.PrimaryHDU(data)
        hdu.header.update({"BUNIT": "W.m-2.pixel-1", "WAVE": self.wave,
                           "CDELT1": -self.pixel / 3600, "CDELT2": self.pixel / 3600,
                           "CRPIX1": (plane.shape[1] + 1) / 2,
                           "CRPIX2": (plane.shape[0] + 1) / 2,
                           "CUNIT1": "deg", "CUNIT2": "deg",
                           "FLUX_1": "Total flux", "FLUX_5": "Direct star flux",
                           "FLUX_6": "Scattered star flux", "FLUX_7": "Direct thermal flux",
                           "FLUX_8": "Scattered thermal flux"})
        if updates:
            hdu.header.update(updates)
        hdu.writeto(self.path, overwrite=True)

    def measure(self, **kwargs):
        options = dict(path=self.path, wavelength_um=self.wave, instrument="NIRSpec",
                       aperture_radius_arcsec=1.0, model_distance_pc=140,
                       target_distance_pc=147)
        options.update(kwargs)
        return measure_image(**options)

    def write_sed(self, total):
        data = np.zeros((8, 1, 1, 1))
        data[0, 0, 0, 0] = total
        primary = fits.PrimaryHDU(data)
        primary.header["BUNIT"] = "W.m-2"
        path = self.root / "sed_rt.fits.gz"
        fits.HDUList([primary, fits.ImageHDU(np.array([self.wave]))]).writeto(path, overwrite=True)
        return path

    def test_point_source_units_and_distance(self):
        self.write_image()
        result = self.measure(psf_fwhm_arcsec=0)
        expected = 2e-14 * (140 / 147)**2
        self.assertAlmostEqual(result["flux_w_m2"] / expected, 1.0, places=13)
        self.assertAlmostEqual(result["flux_jy"] / (expected * 3e-6 / SPEED_OF_LIGHT_M_S * 1e26),
                               1.0, places=13)
        self.assertAlmostEqual(result["diagnostics"]["target_pixel_arcsec"], self.pixel * 140 / 147)
        self.assertEqual(result["quality_status"], "passed_without_coeval_closure")
        json.dumps(result, allow_nan=False)

    def test_distance_changes_aperture_mapping_as_well_as_flux(self):
        self.plane[64, 79] = 6e-14  # 0.75 arcsec at model distance.
        self.write_image()
        near = self.measure(aperture_radius_arcsec=.5, target_distance_pc=140, psf_fwhm_arcsec=0)
        far = self.measure(aperture_radius_arcsec=.5, target_distance_pc=280, psf_fwhm_arcsec=0)
        self.assertAlmostEqual(near["flux_w_m2"] / 2e-14, 1)
        self.assertAlmostEqual(far["flux_w_m2"] / (8e-14 / 4), 1)

    def test_fast_boundary_equals_frozen_midpoint_kernel(self):
        reference = ROOT / "reference/source/strict_continuum_measurement_kernel_v1.py"
        spec = importlib.util.spec_from_file_location("frozen_photometry_reference", reference)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        shape, cx, cy, radius = (47, 53), 25.2, 22.6, 10.83
        fast = fractional_circle_weights(shape, cx, cy, radius, subpixels=64)
        expected = module.fractional_circle_weights(shape, cx, cy, radius)
        np.testing.assert_array_equal(fast, expected)

    def test_adjoint_equals_direct_gaussian_for_signed_image(self):
        y, x = np.indices(self.plane.shape)
        plane = 1e-15 * np.exp(-((x - 64)**2 + (y - 64)**2) / (2 * 5**2))
        plane[64, 66] = -1e-15  # Signed packets must be retained.
        self.write_image(plane)
        result = self.measure(aperture_subpixels=64)
        pixel = self.pixel * 140 / 147
        sigma = (0.033 * self.wave) / FWHM_TO_SIGMA / pixel
        weights = fractional_circle_weights(plane.shape, 64, 64, 1 / pixel, 64)
        blurred = gaussian_filter(plane, sigma, mode="constant", cval=0, truncate=6)
        expected = np.sum(blurred * weights) * (140 / 147)**2
        self.assertAlmostEqual(result["flux_w_m2"] / expected, 1, places=13)
        clipped = np.sum(gaussian_filter(np.maximum(plane, 0), sigma, mode="constant",
                                         cval=0, truncate=6) * weights) * (140 / 147)**2
        self.assertLess(result["flux_w_m2"], clipped)
        self.assertGreater(result["diagnostics"]["full_plane_negative_fraction_diagnostic"], 0)
        self.assertAlmostEqual(
            result["diagnostics"]["aperture_clip_symmetric_fraction_diagnostic"],
            abs(expected - clipped) / (0.5 * (abs(expected) + abs(clipped))), places=13)
        self.assertFalse(result["diagnostics"]["zero_clipping_used_for_scoring"])

    def test_adjoint_field_loss_equals_direct_convolution(self):
        rng = np.random.default_rng(7)
        plane = rng.normal(size=(37, 41))
        weights, sy, sx = _adjoint_weights(plane.shape, 20, 18, 8, 1.3, 16)
        direct = gaussian_filter(plane, 1.3, mode="constant", cval=0, truncate=6)
        self.assertAlmostEqual(float(np.einsum("ij,i,j->", plane, sy, sx)),
                               float(direct.sum()), places=12)
        self.assertFalse(weights.flags.writeable)

    def test_psf_conventions(self):
        self.assertAlmostEqual(gaussian_psf_fwhm_arcsec("NIRSpec", 2), .066)
        self.assertAlmostEqual(gaussian_psf_fwhm_arcsec("MIRI", 10), .436)

    def test_coeval_closure_passes_and_does_not_rescale(self):
        self.write_image()
        original = self.measure()
        sed = self.write_sed(2e-14 * (1 + 1e-6))
        result = self.measure(coeval_sed_path=sed)
        self.assertEqual(result["quality_status"], "passed")
        self.assertEqual(result["flux_w_m2"], original["flux_w_m2"])
        self.assertTrue(result["diagnostics"]["coeval_closure_checked"])

    def test_coeval_closure_failure_surfaces(self):
        self.write_image()
        sed = self.write_sed(2.1e-14)
        with self.assertRaises(PhotometryError) as caught:
            self.measure(coeval_sed_path=sed)
        self.assertIn("coeval_plane0_abs_relative_difference_le_1e-5", caught.exception.failed_checks)

    def test_invalid_metadata_and_shape(self):
        for updates in ({"BUNIT": "Jy/pixel"}, {"WAVE": 4}, {"CRPIX1": 64},
                        {"CDELT2": self.pixel / 1800}, {"PC1_1": 1},
                        {"FLUX_1": "Stokes Q"}, {"CUNIT1": "radian"}):
            with self.subTest(updates=updates):
                self.write_image(updates=updates)
                with self.assertRaises(PhotometryError):
                    self.measure()
        self.write_image(nplanes=4)
        with self.assertRaises(PhotometryError):
            self.measure()

    def test_nonfinite_and_nonpositive(self):
        for bad in (np.nan, -1e-14, 0):
            self.plane[64, 64] = bad
            self.write_image()
            with self.subTest(value=bad), self.assertRaises(PhotometryError):
                self.measure()

    def test_support_and_edge_checks(self):
        self.write_image()
        with self.assertRaises(PhotometryError) as caught:
            self.measure(aperture_radius_arcsec=4)
        self.assertIn("aperture_plus_6sigma_support", caught.exception.failed_checks)
        self.plane[0, 64] = 1e-14
        self.write_image()
        with self.assertRaises(PhotometryError) as caught:
            self.measure()
        self.assertIn("edge_positive_fraction_lt_1e-3", caught.exception.failed_checks)

    def test_aperture_policy_records_boundary_warnings_without_changing_aperture(self):
        self.write_image()
        strict = self.measure()
        aperture = self.measure(quality_policy="aperture_v2")
        self.assertEqual(strict["flux_jy"], aperture["flux_jy"])
        self.plane[0, 64] = 1e-14
        self.write_image()
        with self.assertRaises(PhotometryError):
            self.measure()
        boundary = self.measure(quality_policy="aperture_v2")
        self.assertEqual(strict["flux_jy"], boundary["flux_jy"])
        self.assertTrue(boundary["quality_pass"])
        self.assertIn("edge_positive_fraction_lt_1e-3", boundary["quality_warnings"])
        self.assertFalse(boundary["diagnostics"]["quality_checks"]["edge_positive_fraction_lt_1e-3"])
        self.assertTrue(all(boundary["diagnostics"]["blocking_quality_checks"].values()))

    def test_aperture_policy_keeps_support_closure_and_units_strict(self):
        self.write_image()
        with self.assertRaises(PhotometryError):
            self.measure(quality_policy="aperture_v2", aperture_radius_arcsec=4)
        sed = self.write_sed(3e-14)
        with self.assertRaises(PhotometryError) as caught:
            self.measure(quality_policy="aperture_v2", coeval_sed_path=sed)
        self.assertIn("coeval_plane0_abs_relative_difference_le_1e-5", caught.exception.failed_checks)
        self.write_image(updates={"BUNIT": "Jy"})
        with self.assertRaises(PhotometryError):
            self.measure(quality_policy="aperture_v2")
        with self.assertRaisesRegex(ValueError, "Unknown quality_policy"):
            self.measure(quality_policy="skip_checks")

    def test_subpixel_argument_and_psf_validation(self):
        self.write_image()
        self.measure(aperture_subpixels=16)  # A populated cache must not skip input validation.
        for kwargs in ({"aperture_subpixels": 0}, {"aperture_subpixels": 1.5},
                       {"aperture_subpixels": 16.0}, {"aperture_subpixels": True},
                       {"target_distance_pc": 0}, {"psf_fwhm_arcsec": -1},
                       {"psf_fwhm_arcsec": "bad"}, {"instrument": None},
                       {"instrument": "unknown"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(PhotometryError):
                self.measure(**kwargs)


if __name__ == "__main__":
    unittest.main()
