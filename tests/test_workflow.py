"""Configuration, restart and Slurm tests; all simulator calls are mocked."""
from __future__ import annotations

from contextlib import ExitStack
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mcfost_grid import runner
from mcfost_grid.configuration import (WORKSPACE, expand_models, load_json,
                                       prepare_run, read_anchors, sha256)
from mcfost_grid.slurm import write_slurm


class CatalogueTests(unittest.TestCase):
    def test_cartesian_grid_is_explicit_and_ordered(self):
        models = expand_models({"grid": {"inclination_deg": [50, 60], "envelope_size_exponent": [2.75, 3.25]}})
        self.assertEqual(len(models), 4)
        self.assertEqual([(m["parameters"]["inclination_deg"], m["parameters"]["envelope_size_exponent"]) for m in models],
                         [(50, 2.75), (50, 3.25), (60, 2.75), (60, 3.25)])
        self.assertEqual(len({m["id"] for m in models}), 4)
        self.assertEqual(models, expand_models({"grid": {"inclination_deg": [50, 60], "envelope_size_exponent": [2.75, 3.25]}}))

    def test_duplicate_and_oversized_models_rejected(self):
        for configuration in [{"grid": {"inclination_deg": [50, 50]}},
                              {"grid": {"inclination_deg": [50, 50.0]}},
                              {"models": [{"inclination_deg": 50}, {"inclination_deg": 50}]},
                              {"grid": {"inclination_deg": [50, 60]}, "max_models": 1},
                              {"grid": {}, "models": [{}]}]:
            with self.subTest(configuration=configuration), self.assertRaises(ValueError):
                expand_models(configuration)

    def test_custom_anchors_require_matched_background_spectrum(self):
        with self.assertRaisesRegex(ValueError, "matching spectrum_file"):
            read_anchors({"observations": {"anchors_file": "custom.csv"}}, WORKSPACE)

    def test_max_models_is_positive_integer(self):
        for limit in (False, 0, -1, 1.2):
            with self.subTest(limit=limit), self.assertRaisesRegex(ValueError, "max_models"):
                expand_models({"grid": {}, "max_models": limit})

    def test_strict9_matches_frozen_wavelength_flux_and_uncertainty(self):
        anchors, _, measurement = read_anchors({}, WORKSPACE)
        old = Table.read(WORKSPACE / "reference/predictions/continuum_predictions.ecsv")
        old = old[old["model_index"] == old["model_index"][0]]
        self.assertEqual([a["id"] for a in anchors], [f"w{i:03d}" for i in range(1, 10)])
        self.assertEqual(len(anchors), 9)
        self.assertEqual(measurement["aperture_radius_arcsec"], 1.0)
        for anchor in anchors:
            rows = old[old["anchor_tag"] == anchor["id"]]
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(anchor["wavelength_um"], float(row["wavelength_um"]))
            jy_factor = anchor["wavelength_um"] * 1e-6 / 299792458.0 * 1e26
            self.assertAlmostEqual(anchor["flux_jy"] / (float(row["observed_flux_geometric_w_m2"]) * jy_factor), 1., places=12)
            self.assertAlmostEqual(anchor["uncertainty_jy"] / (float(row["observed_error_w_m2"]) * jy_factor), 1., places=12)
        self.assertGreater(min(a["wavelength_um"] for a in anchors), 1.0)

    def test_strict9_rejects_different_aperture(self):
        with self.assertRaisesRegex(ValueError, "require aperture"):
            read_anchors({"observations": {"aperture_radius_arcsec": .35}}, WORKSPACE)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        config = {"schema_version": 1, "run_name": "test_run", "output_dir": str(self.root / "runs"),
                  "smoke_test": True, "observations": {"anchor_ids": ["w001", "w008"]},
                  "numerics": {"image_npix": 31, "photons_temperature": 10,
                               "photons_image": 10, "photons_sed": 10}}
        self.config = self.root / "grid.json"
        self.config.write_text(json.dumps(config))
        self.run = prepare_run(self.config)
        self.manifest = load_json(self.run / "manifest.json")
        self.model = self.manifest["models"][0]
        self.directory = self.run / "models" / self.model["id"]
        utils = self.root / "utils"
        for name in ("Dust", "Lambda", "Stellar_Spectra"):
            (utils / name).mkdir(parents=True)
            (utils / name / "fixture.txt").write_text("frozen fixture\n")
        executable = self.root / "fake_mcfost"
        executable.write_text("#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)
        self.machine_value = {"schema_version": 1, "mcfost_executable": str(executable),
                              "mcfost_utils": str(utils), "backend": "image_method2", "threads": 1,
                              "timeout_seconds": 10, "max_memory_gb": 1,
                              "slurm": {"partition": "test", "max_parallel": 3,
                                        "setup_lines": ["module load Example/1.0"]}}
        self.machine = self.root / "machine.json"
        self.machine.write_text(json.dumps(self.machine_value))
        self.commands = []

    def tearDown(self):
        self.temporary.cleanup()

    def fake_invoke(self, command, cwd, env, timeout, attempt):
        self.commands.append((command, cwd))
        self.assertFalse(Path(command[1]).is_absolute())
        root_arg = command[command.index("-root_dir") + 1]
        self.assertFalse(Path(root_arg).is_absolute())
        output = (cwd / root_arg).resolve()
        # Upstream MCFOST nests products for -seed but still reads -Tfile
        # directly underneath root_dir. Exercise that distinction explicitly.
        products = output
        if "-seed" in command:
            products = output / ("seed=" + command[command.index("-seed") + 1])
        if command[1] == "temperature.para":
            target = products / "data_th" / "Temperature.fits.gz"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"fresh equilibrium fixture")
        else:
            self.assertEqual(command[command.index("-Tfile") + 1], "Temperature.fits.gz")
            link = output / "Temperature.fits.gz"
            self.assertTrue(link.is_symlink())
            self.assertFalse(Path(os.readlink(link)).is_absolute())
            self.assertEqual(link.read_bytes(), b"fresh equilibrium fixture")
            target = products / "data_image" / "RT.fits.gz"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"mock image " + cwd.name.encode())
            if "-rt-sed-method" in command:
                (products / "data_th").mkdir()
                (products / "data_th" / "sed_rt.fits.gz").write_bytes(b"mock coeval SED")
            self.assertIn("-rt2", command)
            self.assertNotIn("-mol", command)
        return .001, "Using scattering method 2\nProcessing complete\n"

    def patches(self, invoke=None):
        stack = ExitStack()
        stack.enter_context(mock.patch.object(runner, "_invoke", side_effect=invoke or self.fake_invoke))
        stack.enter_context(mock.patch.object(runner, "_valid_temperature"))
        stack.enter_context(mock.patch.object(runner, "check_backend", return_value=["mock MCFOST 4.1.13"]))
        stack.enter_context(mock.patch.object(runner, "cpu_capacity", return_value=4))
        stack.enter_context(mock.patch("mcfost_grid.photometry.measure_image", side_effect=lambda *args, **kwargs: {"flux_jy": 1.0, "quality_pass": True}))
        return stack

    def test_frozen_input_guard(self):
        runner.validate_inputs(self.run, self.manifest)
        with (self.directory / "temperature.para").open("a") as stream:
            stream.write("\nchanged\n")
        with self.assertRaisesRegex(RuntimeError, "Frozen run input changed"):
            runner.validate_inputs(self.run, self.manifest)

    def test_reprepare_never_overwrites_run(self):
        before = sha256(self.run / "manifest.json")
        with self.assertRaises(FileExistsError):
            prepare_run(self.config)
        self.assertEqual(before, sha256(self.run / "manifest.json"))

    def test_dotted_run_name_prepares_as_one_safe_directory(self):
        config = load_json(self.config)
        config["run_name"] = "continuum_production_v1.1"
        self.config.write_text(json.dumps(config))
        run = prepare_run(self.config)
        self.assertEqual(run, (self.root / "runs/continuum_production_v1.1").resolve())
        manifest = load_json(run / "manifest.json")
        self.assertEqual(manifest["run_id"], config["run_name"])
        runner.validate_inputs(run, manifest)

    def test_run_name_rejects_paths_empty_components_and_non_strings(self):
        config = load_json(self.config)
        before = sorted(path.name for path in (self.root / "runs").iterdir())
        for name in ("", ".", "..", ".hidden", "trailing.", "a..b", "../escape", "a/b", "a\\b", "a b", None, 1):
            with self.subTest(name=name):
                config["run_name"] = name
                self.config.write_text(json.dumps(config))
                with self.assertRaisesRegex(ValueError, "run_name"):
                    prepare_run(self.config)
        self.assertEqual(sorted(path.name for path in (self.root / "runs").iterdir()), before)

    def test_slurm_syntax_allocation_and_spool_safe_paths(self):
        result = write_slurm(self.run, self.machine)
        self.assertFalse(result["submitted"])
        for name in ("job_array.sh", "job_analysis.sh", "submit.sh"):
            parsed = subprocess.run(["bash", "-n", str(self.run / name)], capture_output=True, text=True)
            self.assertEqual(parsed.returncode, 0, parsed.stderr)
        array = (self.run / "job_array.sh").read_text()
        self.assertIn("#SBATCH --ntasks=1", array)
        self.assertIn("#SBATCH --cpus-per-task=1", array)
        self.assertIn("#SBATCH --array=0-0%3", array)
        self.assertIn('cd "$MCFOST_RUN_DIR"', array)
        self.assertIn('code/workflow.py task "$MCFOST_RUN_DIR"', array)
        self.assertNotIn("dirname", array)
        self.assertNotIn("mpirun", array)
        self.assertNotIn("srun", array)
        submit = (self.run / "submit.sh").read_text()
        self.assertIn("--export=ALL", submit)
        self.assertIn("afterany:$job", submit)
        self.assertIn("mkdir -p logs", submit)

    def test_success_and_complete_resume_do_not_rerun(self):
        with self.patches():
            first = runner.run_model(self.run, 0, self.machine)
            # Simulate interruption after measurements were committed but
            # before the complete status marker was published.
            (self.directory / "status.json").write_text(json.dumps({"state": "running"}))
            second = runner.run_model(self.run, 0, self.machine)
        self.assertEqual(first["status"], "complete")
        self.assertEqual(second["status"], "cached")
        self.assertEqual(len(self.commands), 3)
        self.assertTrue(all("-seed" not in command for command, _ in self.commands))
        self.assertTrue(load_json(self.directory / "measurements.json")["complete"])
        self.assertEqual(load_json(self.directory / "status.json")["state"], "complete")

    def test_global_runtime_binding_prevents_mixed_model_builds(self):
        runner._check_run_runtime(self.run, {"backend": "A", "binary": "sha1"}, commit=True)
        runner._check_run_runtime(self.run, {"backend": "A", "binary": "sha1"}, commit=True)
        with self.assertRaisesRegex(RuntimeError, "all models must use the same"):
            runner._check_run_runtime(self.run, {"backend": "A", "binary": "sha2"}, commit=True)

    def test_photometry_failure_preserves_named_diagnostics(self):
        from mcfost_grid.photometry import PhotometryError
        error = PhotometryError("field too small", failed_checks=("aperture_support",),
                                diagnostics={"available_half_field_arcsec": .1})
        with self.patches(), mock.patch("mcfost_grid.photometry.measure_image", side_effect=error):
            with self.assertRaises(PhotometryError):
                runner.run_model(self.run, 0, self.machine)
        state = load_json(self.directory / "status.json")
        self.assertEqual(state["state"], "failed")
        self.assertEqual(state["failed_checks"], ["aperture_support"])
        self.assertEqual(state["diagnostics"]["available_half_field_arcsec"], .1)

    def test_failed_anchor_resumes_without_repeating_temperature_or_success(self):
        failed = False
        def invoke(command, cwd, env, timeout, attempt):
            nonlocal failed
            if command[1] == "image.para" and cwd.name == "w008" and not failed:
                failed = True
                raise RuntimeError("simulated image failure")
            return self.fake_invoke(command, cwd, env, timeout, attempt)
        with self.patches(invoke):
            with self.assertRaisesRegex(RuntimeError, "simulated image failure"):
                runner.run_model(self.run, 0, self.machine)
            self.assertEqual(load_json(self.directory / "status.json")["state"], "failed")
            self.assertEqual(len(load_json(self.directory / "measurements.json")["measurements"]), 1)
            result = runner.run_model(self.run, 0, self.machine)
        self.assertEqual(result["status"], "complete")
        self.assertEqual([command[1] for command, cwd in self.commands].count("temperature.para"), 1)
        self.assertEqual(len(self.commands), 3)
        self.assertEqual(len(list((self.directory / "anchors/w008/attempts").iterdir())), 2)

    def test_completed_cache_detects_changed_temperature(self):
        with self.patches():
            runner.run_model(self.run, 0, self.machine)
            receipt = load_json(self.directory / "temperature_complete.json")
            (self.run / receipt["path"]).write_bytes(b"not the saved temperature")
            with self.assertRaisesRegex(RuntimeError, "Cached temperature changed"):
                runner.run_model(self.run, 0, self.machine)

    def test_completed_cache_detects_changed_image(self):
        with self.patches():
            runner.run_model(self.run, 0, self.machine)
            data = load_json(self.directory / "measurements.json")
            (self.run / data["measurements"][0]["image_path"]).write_bytes(b"changed image")
            with self.assertRaisesRegex(RuntimeError, "Cached image changed"):
                runner.run_model(self.run, 0, self.machine)

    def test_resume_detects_changed_external_utilities(self):
        with self.patches():
            runner.run_model(self.run, 0, self.machine)
            (self.root / "utils/Stellar_Spectra/fixture.txt").write_text("a changed stellar atmosphere\n")
            with self.assertRaisesRegex(RuntimeError, "runtime or manifest changed"):
                runner.run_model(self.run, 0, self.machine)

    def test_coeval_commands_and_saved_sed(self):
        self.machine_value["backend"] = "coeval_method2"
        self.machine.write_text(json.dumps(self.machine_value))
        with self.patches():
            runner.run_model(self.run, 0, self.machine)
            self.assertEqual(runner.run_model(self.run, 0, self.machine)["status"], "cached")
        data = load_json(self.directory / "measurements.json")
        self.assertTrue(all("coeval_sed_sha256" in m for m in data["measurements"]))
        self.assertTrue(all("-rt-sed-method" in command for command, cwd in self.commands[1:]))

    def test_seed_applies_to_temperature_and_both_image_backends(self):
        for backend in ("image_method2", "coeval_method2"):
            with self.subTest(backend=backend):
                config = load_json(self.config)
                config["run_name"] = "seeded_" + backend
                config["numerics"]["random_seed"] = 41001
                self.config.write_text(json.dumps(config))
                run = prepare_run(self.config)
                self.machine_value["backend"] = backend
                self.machine.write_text(json.dumps(self.machine_value))
                self.commands.clear()
                with self.patches():
                    self.assertEqual(runner.run_model(run, 0, self.machine)["status"], "complete")
                    self.assertEqual(runner.run_model(run, 0, self.machine)["status"], "cached")
                self.assertEqual(len(self.commands), 3)
                for command, _ in self.commands:
                    self.assertEqual(command.count("-seed"), 1)
                    self.assertEqual(command[command.index("-seed") + 1], "41001")
                manifest = load_json(run / "manifest.json")
                directory = run / "models" / manifest["models"][0]["id"]
                temperature = load_json(directory / "temperature_complete.json")
                self.assertIn("/seed=41001/", temperature["path"])
                measurements = load_json(directory / "measurements.json")["measurements"]
                self.assertTrue(all("/seed=41001/" in m["image_path"] for m in measurements))
                # A changed seed must not accept cached measurements, even if
                # the caller alters the manifest rather than the frozen JSON.
                manifest["configuration"]["numerics"]["random_seed"] = 41002
                (run / "manifest.json").write_text(json.dumps(manifest))
                with self.patches(), self.assertRaisesRegex(RuntimeError, "runtime or manifest changed"):
                    runner.run_model(run, 0, self.machine)
                self.assertEqual(len(self.commands), 3)

    def test_wrong_allocated_thread_count_fails_before_invocation(self):
        self.machine_value["threads"] = 8
        self.machine.write_text(json.dumps(self.machine_value))
        with self.patches(), self.assertRaisesRegex(ValueError, "exceeds allocated"):
            runner.run_model(self.run, 0, self.machine)
        self.assertFalse(self.commands)

    def test_fitsio_error_overrides_zero_exit_and_completion_banner(self):
        attempt = self.root / "bad_attempt"
        attempt.mkdir()
        process = mock.Mock()
        process.wait.return_value = 0
        def popen(*args, **kwargs):
            kwargs["stdout"].write(b"FITSIO Error Status = 105\nProcessing complete\n")
            return process
        with mock.patch.object(runner.subprocess, "Popen", side_effect=popen), self.assertRaisesRegex(RuntimeError, "failed or incomplete"):
            runner._invoke(["not-executed"], self.root, {"OMP_NUM_THREADS": "1"}, 10, attempt)
        self.assertEqual(load_json(attempt / "execution.json")["returncode"], 0)

    def test_runtime_limits_must_be_finite_not_booleans(self):
        for key, value in [("max_memory_gb", True), ("timeout_seconds", False), ("timeout_seconds", float("inf"))]:
            machine = dict(self.machine_value, **{key: value})
            self.machine.write_text(json.dumps(machine))
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                runner.runtime_config(self.machine)


if __name__ == "__main__":
    unittest.main()
