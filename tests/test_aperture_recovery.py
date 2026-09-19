"""Exercise policy migration and checkpoint reuse with real synthetic FITS."""
from contextlib import ExitStack
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from mcfost_grid import runner, photometry
from mcfost_grid.analysis import analyze_run
from mcfost_grid.configuration import atomic_json, load_json, prepare_run, sha256


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/"scripts"/(name+".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


recovery = load_script("prepare_aperture_recovery")
audit = load_script("audit_aperture_boundaries")


class ApertureRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        config = self.root/"grid.json"
        atomic_json(config, {"schema_version": 1, "run_name": "parent", "output_dir": str(self.root),
                            "models": [{"inclination_deg": 50}, {"inclination_deg": 70}],
                            "observations": {"anchor_ids": ["w001", "w006", "w007"]},
                            "numerics": {"random_seed": 41001, "image_npix": 129}})
        self.parent = prepare_run(config)
        self.manifest = load_json(self.parent/"manifest.json")
        self.model_ids = [m["id"] for m in self.manifest["models"]]
        utils = self.root/"utils"
        for sub in ("Dust", "Lambda", "Stellar_Spectra"):
            (utils/sub).mkdir(parents=True)
            (utils/sub/"fixture").write_text("unchanged utility")
        exe = self.root/"mcfost"
        exe.write_text("#!/bin/sh\nexit 99\n")
        exe.chmod(0o755)
        self.machine = self.root/"machine.json"
        atomic_json(self.machine, {"schema_version": 1, "mcfost_executable": str(exe),
                                  "mcfost_utils": str(utils), "backend": "image_method2", "threads": 1})
        self.commands = []
        with self.patches():
            runner.run_model(self.parent, 0, self.machine)
            with self.assertRaises(photometry.PhotometryError):
                runner.run_model(self.parent, 1, self.machine)
        row = audit.audit_model(self.parent, self.manifest, 1, photometry, sha256)
        self.audit_path = self.root/"audit.json"
        atomic_json(self.audit_path, {"audit": "read_only_aperture_crop_invariance_v1", "run_id": "parent",
             "manifest_sha256": sha256(self.parent/"manifest.json"),
             "audit_script_sha256": sha256(ROOT/"scripts/audit_aperture_boundaries.py"), "models": [row],
             "original_acceptance_and_results_unchanged": True})

    def patches(self):
        stack = ExitStack()
        stack.enter_context(mock.patch.object(runner, "_invoke", side_effect=self.fake_invoke))
        stack.enter_context(mock.patch.object(runner, "check_backend", return_value=["fixture"]))
        return stack

    def fake_invoke(self, command, cwd, env, timeout, attempt):
        self.commands.append(command)
        output = cwd/command[command.index("-root_dir")+1]/"seed=41001"
        log = "Using scattering method 2\nProcessing complete\n"
        if command[1] == "temperature.para":
            product = output/"data_th/Temperature.fits.gz"
            product.parent.mkdir(parents=True)
            fits.PrimaryHDU(np.ones((3, 3))*100).writeto(product)
        else:
            product = output/"data_image/RT.fits.gz"
            product.parent.mkdir(parents=True)
            data = np.zeros((8, 1, 1, 129, 129))
            data[0, 0, 0, 64, 64] = 2e-14
            if self.model_ids[1] in str(cwd) and cwd.name == "w006":
                data[0, 0, 0, 0, :] = 1e-15
            data[4] = data[0]
            hdu = fits.PrimaryHDU(data)
            hdu.header.update({"BUNIT": "W.m-2.pixel-1", "WAVE": float(command[command.index("-img")+1]),
                               "CDELT1": -.05/3600, "CDELT2": .05/3600, "CRPIX1": 65., "CRPIX2": 65.,
                               "FLUX_1": "Total flux", "FLUX_5": "Direct star flux", "FLUX_6": "Scattered star flux",
                               "FLUX_7": "Direct thermal flux", "FLUX_8": "Scattered thermal flux"})
            hdu.writeto(product)
        atomic_json(attempt/"command.json", {"argv": command})
        atomic_json(attempt/"execution.json", {"returncode": 0, "elapsed_seconds": .01})
        (attempt/"mcfost.log").write_text(log)
        return .01, log

    def test_fork_preserves_parent_reuses_images_and_temperature_then_finishes_only_missing_anchor(self):
        before = {str(p.relative_to(self.parent)): sha256(p) for p in self.parent.rglob("*") if p.is_file()}
        receipt = recovery.prepare_recovery(self.parent, self.audit_path, "recovered")
        run = self.root/"recovered"
        manifest = load_json(run/"manifest.json")
        self.assertEqual(receipt["remaining_image_calculations"], 1)
        self.assertEqual(receipt["new_temperature_calculations"], 0)
        self.assertEqual(manifest["measurement"]["quality_policy"], "aperture_v2")
        self.assertEqual(manifest["models"], self.manifest["models"])
        self.assertEqual(manifest["anchors"], self.manifest["anchors"])
        runner.validate_inputs(run, manifest)
        before_count = len(self.commands)
        with self.patches():
            self.assertEqual(runner.run_model(run, 0, self.machine)["status"], "cached")
            self.assertEqual(runner.run_model(run, 1, self.machine)["status"], "complete")
        self.assertEqual(len(self.commands)-before_count, 1)
        self.assertEqual(self.commands[-1][1], "image.para")
        self.assertEqual(before, {str(p.relative_to(self.parent)): sha256(p) for p in self.parent.rglob("*") if p.is_file()})
        for mid in self.model_ids:
            link = run/"models"/mid/"parent_products"
            self.assertTrue(link.is_symlink())
            self.assertFalse(Path(os.readlink(link)).is_absolute())
        with mock.patch("mcfost_grid.analysis._plots", return_value=[]):
            summary = analyze_run(run)
        self.assertEqual(summary["ranked_model_count"], 2)
        self.assertEqual(summary["measurement_quality_policy"], "aperture_v2")
        self.assertEqual(len(summary["measurement_quality_warnings"]), 1)
        cached = load_json(run/"models"/self.model_ids[1]/"measurements.json")
        warned = next(m for m in cached["measurements"] if m["quality_warnings"])
        self.assertTrue(warned["quality_pass"])
        self.assertFalse(warned["diagnostics"]["quality_checks"]["edge_positive_fraction_lt_1e-3"])
        with self.assertRaises(FileExistsError):
            recovery.prepare_recovery(self.parent, self.audit_path, "recovered")
        # A parent image is never trusted merely because it was imported once.
        image = run/warned["image_path"]
        image.write_bytes(b"changed raw product")
        with self.patches(), self.assertRaisesRegex(RuntimeError, "Cached image changed"):
            runner.run_model(run, 1, self.machine)

    def test_unpassed_audit_and_active_parent_are_rejected(self):
        original = self.audit_path.read_bytes()
        value = load_json(self.audit_path)
        value["models"][0]["crop_comparison"]["relative_difference"] = .1
        atomic_json(self.audit_path, value)
        with self.assertRaisesRegex(ValueError, "audit did not pass"):
            recovery.prepare_recovery(self.parent, self.audit_path, "rejected")
        self.assertFalse((self.root/"rejected").exists())
        self.audit_path.write_bytes(original)
        with (self.parent/"models"/self.model_ids[0]/".task.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(ValueError, "still running"):
                recovery.prepare_recovery(self.parent, self.audit_path, "rejected")
        self.assertFalse((self.root/"rejected").exists())

    def test_forged_policy_cannot_bypass_mandatory_checks_or_mix_cached_policies(self):
        recovery.prepare_recovery(self.parent, self.audit_path, "recovered")
        run = self.root/"recovered"
        path = run/"models"/self.model_ids[0]/"measurements.json"
        original = path.read_bytes()
        data = load_json(path)
        data["measurements"][0]["quality_policy"] = "whole_image_v1"
        atomic_json(path, data)
        with self.patches(), self.assertRaisesRegex(RuntimeError, "quality policy differs"):
            runner.run_model(run, 0, self.machine)
        with mock.patch("mcfost_grid.analysis._plots", return_value=[]):
            self.assertEqual(analyze_run(run)["ranked_model_count"], 0)
        path.write_bytes(original)
        data = load_json(path)
        data["measurements"][0]["diagnostics"]["quality_checks"]["aperture_plus_6sigma_support"] = False
        atomic_json(path, data)
        with self.patches(), self.assertRaisesRegex(RuntimeError, "mandatory checks"):
            runner.run_model(run, 0, self.machine)


if __name__ == "__main__":
    unittest.main()
