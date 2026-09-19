"""Check boundary audit against signed synthetic images and clipped support."""
import importlib.util
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid import photometry

spec = importlib.util.spec_from_file_location("boundary_audit", ROOT / "scripts/audit_aperture_boundaries.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class BoundaryAuditTests(unittest.TestCase):
    def test_bright_border_and_signed_interior_preserve_aperture(self):
        d = {"centre_x_zero_based": 70., "centre_y_zero_based": 60.,
             "aperture_radius_pixels": 10.3, "psf_sigma_pixels": 1.7,
             "psf_truncation_sigma": 6., "aperture_subpixels": 64}
        y, x = np.indices((121, 141))
        plane = 1e-14*np.exp(-((x-70)**2+(y-60)**2)/32)
        plane[60, 70] = -1e-14
        initial = audit.compare_crop(photometry, plane, d)
        plane[0] = 1e-10  # Deliberately overwhelms whole-image edge/loss tests.
        after = audit.compare_crop(photometry, plane, d)
        self.assertTrue(initial["aperture_crop_invariance_verified"])
        self.assertTrue(after["aperture_crop_invariance_verified"])
        self.assertEqual(initial["full_image_adjoint_aperture_w_m2"], after["full_image_adjoint_aperture_w_m2"])
        self.assertLess(after["relative_difference"], 1e-12)

    def test_cannot_certify_aperture_with_insufficient_support(self):
        d = {"centre_x_zero_based": 20., "centre_y_zero_based": 20.,
             "aperture_radius_pixels": 10., "psf_sigma_pixels": 4.,
             "psf_truncation_sigma": 6., "aperture_subpixels": 64}
        with self.assertRaisesRegex(ValueError, "No strictly interior crop"):
            audit.compare_crop(photometry, np.ones((41, 41)), d)

    def test_failed_fits_is_verified_read_only_and_mismatched_status_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            directory = run / "models/m0000_test"
            path = directory / "anchors/w001/attempts/attempt_001/output/seed=41001/RT.fits.gz"
            path.parent.mkdir(parents=True)
            data = np.zeros((8, 1, 1, 129, 129))
            data[0, 0, 0, 64, 64] = 2e-14
            data[0, 0, 0, 0, :] = 1e-15
            data[4] = data[0]
            hdu = fits.PrimaryHDU(data)
            hdu.header.update({"BUNIT": "W.m-2.pixel-1", "WAVE": 3.,
                               "CDELT1": -.05/3600, "CDELT2": .05/3600,
                               "CRPIX1": 65., "CRPIX2": 65.,
                               "FLUX_1": "Total flux", "FLUX_5": "Direct star flux",
                               "FLUX_6": "Scattered star flux", "FLUX_7": "Direct thermal flux",
                               "FLUX_8": "Scattered thermal flux"})
            hdu.writeto(path)
            measurement = {"aperture_radius_arcsec": 1., "target_distance_pc": 147.,
                           "aperture_subpixels": 64}
            with self.assertRaises(photometry.PhotometryError) as caught:
                photometry.measure_image(path, wavelength_um=3., instrument="NIRSpec",
                                         model_distance_pc=140., psf_fwhm_arcsec=.099, **measurement)
            status = {"state": "failed", "failed_checks": list(caught.exception.failed_checks),
                      "diagnostics": caught.exception.diagnostics}
            status_path = directory / "status.json"
            status_path.write_text(json.dumps(status))
            manifest = {"models": [{"id": directory.name, "parameters": {"distance_pc": 140.}}],
                        "anchors": [{"id": "w001", "wavelength_um": 3., "instrument": "NIRSpec",
                                     "psf_fwhm_arcsec": .099}], "measurement": measurement}
            digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            before = {p: digest(p) for p in (path, status_path)}
            result = audit.audit_model(run, manifest, 0, photometry, digest)
            self.assertTrue(result["original_diagnostics_reproduced"])
            self.assertTrue(result["crop_comparison"]["aperture_crop_invariance_verified"])
            self.assertEqual(before, {p: digest(p) for p in before})
            status["diagnostics"]["aperture_flux_at_model_distance_w_m2"] *= 2
            status_path.write_text(json.dumps(status))
            with self.assertRaisesRegex(ValueError, "latest FITS and saved status disagree"):
                audit.audit_model(run, manifest, 0, photometry, digest)


if __name__ == "__main__":
    unittest.main()
