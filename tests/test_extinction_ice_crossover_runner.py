"""Fixed-temperature image experiments: provenance, resume and strict failures.

Fixtures are generated from tracked inputs. No downloaded simulations or actual
MCFOST executable are needed; only the external process boundary is simulated.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid import photometry, runner


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


original_builder = module("crossover_original_builder", "build_extinction_ice_numerical_package.py")


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_image(path, wave, pixel=0.05, size=129, invalid=False):
    """Sparse point-source image with independent labelled component planes."""
    data = np.zeros((8, 1, 1, size, size), dtype=np.float32)
    data[4, 0, 0, size // 2, size // 2] = 1e-14
    data[6, 0, 0, size // 2, size // 2] = 2e-15
    data[0] = data[4] + data[6]
    if invalid:
        data[0, 0, 0, size // 2, size // 2] = np.nan
    hdu = fits.PrimaryHDU(data)
    hdu.header.update({"BUNIT": "W.m-2.pixel-1", "WAVE": wave,
                       "CDELT1": -pixel / 3600, "CDELT2": pixel / 3600,
                       "CRPIX1": (size + 1) / 2, "CRPIX2": (size + 1) / 2,
                       "CUNIT1": "deg", "CUNIT2": "deg",
                       "FLUX_1": "Total flux", "FLUX_5": "Direct star flux",
                       "FLUX_6": "Scattered star flux", "FLUX_7": "Direct thermal flux",
                       "FLUX_8": "Scattered thermal flux"})
    path.parent.mkdir(parents=True, exist_ok=True)
    hdu.writeto(path, overwrite=True)


def make_source(directory, executable):
    """Construct genuine prepared inputs and synthetic successful source records."""
    original_builder.build_package(directory)
    experiment = json.loads((directory / "experiment.json").read_text())
    common = {"executable_sha256": runner.sha256(executable),
              "utilities_content_sha256": "a" * 64, "backend": "image_method2",
              "threads": 64, "python_version": "fixture-original-python",
              "packages": {"fixture": "old-runtime"}}
    write_json(directory / "runtime_binding.json",
               {"fingerprint": common, "experiment_sha256": runner.sha256(directory / "experiment.json")})
    for seed in (42004,):
        task = next(t for t in experiment["tasks"]
                    if t["seed"] == seed and t["photon_packets"] == 2048000 and t["inclination_deg"] == 50)
        run = directory / task["run_path"]
        manifest = json.loads((run / "manifest.json").read_text())
        model = manifest["models"][task["model_index"]]
        model_dir = run / "models" / model["id"]
        identity = {key: common[key] for key in ("executable_sha256", "utilities_content_sha256", "backend")}
        identity["manifest_sha256"] = runner.sha256(run / "manifest.json")
        runtime = {"mcfost_executable": str(executable), "backend": "image_method2", "threads": 64}
        write_json(model_dir / "runtime_binding.json", {"fingerprint": identity, "runtime": runtime})
        attempt = model_dir / "temperature/attempt_001"
        temperature = attempt / "output" / f"seed={seed}" / "data_th/Temperature.fits.gz"
        temperature.parent.mkdir(parents=True)
        fits.PrimaryHDU(np.full((1, 1, 70, 100), 30 + (seed - 42000) / 100)).writeto(temperature)
        write_json(attempt / "command.json",
                   {"argv": [str(executable), "temperature.para", "-seed", str(seed)],
                    "cwd": str(model_dir), "threads": "64"})
        write_json(attempt / "execution.json", {"returncode": 0, "elapsed_seconds": 1})
        (attempt / "mcfost.log").write_text("Processing complete\n")
        write_json(model_dir / "temperature_complete.json",
                   {"fingerprint": identity, "path": str(temperature.relative_to(run)),
                    "sha256": runner.sha256(temperature), "elapsed_seconds": 1})
        measurements = []
        for anchor in manifest["anchors"]:
            attempt = model_dir / "anchors" / anchor["id"] / "attempts/attempt_001"
            output = attempt / "output"
            image = output / f"seed={seed}" / f"data_{anchor['wavelength_um']}" / "RT.fits.gz"
            # The broadest MIRI PSF also fits in this small synthetic field.
            write_image(image, anchor["wavelength_um"], pixel=.1)
            shutil.copyfile(temperature, output / "Temperature.fits.gz")
            command = [str(executable), "image.para", "-img", str(anchor["wavelength_um"]),
                       "-seed", str(seed), "-rt2", "-no_T", "-Tfile", "Temperature.fits.gz"]
            write_json(attempt / "command.json", {"argv": command, "cwd": str(attempt.parents[1]), "threads": "64"})
            write_json(attempt / "execution.json", {"returncode": 0, "elapsed_seconds": 1})
            (attempt / "mcfost.log").write_text("Using ray-tracing method 2\nProcessing complete\n")
            measured = photometry.measure_image(image, wavelength_um=anchor["wavelength_um"],
                                                instrument=anchor["instrument"], model_distance_pc=140,
                                                psf_fwhm_arcsec=anchor["psf_fwhm_arcsec"], **manifest["measurement"])
            measured.update(anchor_id=anchor["id"], wavelength_um=anchor["wavelength_um"],
                            image_path=str(image.relative_to(run)), image_sha256=runner.sha256(image),
                            backend="image_method2")
            write_json(attempt / "measurement.json", measured)
            measurements.append(measured)
        write_json(model_dir / "measurements.json", {"model_id": model["id"], "model_index": 0,
                                                     "parameters": model["parameters"], "fingerprint": identity,
                                                     "complete": True, "measurements": measurements})
        write_json(model_dir / "status.json", {"state": "complete", "completed_anchors": 6, "total_anchors": 6})
    return experiment


class CrossoverRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storage = tempfile.TemporaryDirectory()
        cls.shared = Path(cls.storage.name)
        cls.executable = cls.shared / "fake_mcfost"
        cls.executable.write_text("fixed fixture executable\n")
        cls.source = cls.shared / "source"
        make_source(cls.source, cls.executable)
        cls.builder = module("crossover_builder_tests", "build_extinction_ice_crossover_v2.py")
        cls.dispatcher = module("crossover_dispatch_tests", "extinction_ice_crossover_task.py")
        cls.analyzer = module("crossover_analyzer_runner_tests", "analyze_extinction_ice_crossover_v2.py")

    @classmethod
    def tearDownClass(cls):
        cls.storage.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.bundle = self.root / "portable bundle"
        self.builder.build_package(self.source, self.bundle)
        self.experiment = self.dispatcher.validate_package(self.bundle)
        self.calls = []
        self.fail_call = None
        self.invalid_call = None
        self.runtime = {"mcfost_executable": str(self.executable), "mcfost_utils": str(self.root / "utils"),
                        "backend": "image_method2", "threads": 64, "timeout_seconds": 30,
                        "max_memory_gb": 64}
        self.patches = [mock.patch.object(self.dispatcher, "load_modules", return_value=(runner, photometry)),
                        mock.patch.object(runner, "runtime_config", return_value=self.runtime),
                        mock.patch.object(runner, "cpu_capacity", return_value=64),
                        mock.patch.object(runner, "_utilities_hash", return_value="a" * 64),
                        mock.patch.object(runner, "check_backend", return_value=["fixture method 2"]),
                        mock.patch.object(runner, "_invoke", side_effect=self.invoke)]
        for patch in self.patches:
            patch.start()

    def tearDown(self):
        for patch in reversed(self.patches):
            patch.stop()
        self.temp.cleanup()

    def invoke(self, command, cwd, env, timeout, attempt):
        """Only external execution is replaced: output photometry remains real."""
        self.assertIn("-img", command, "The crossover must never launch an equilibrium solve")
        self.assertIn("-no_T", command)
        self.assertIn("-rt2", command)
        self.assertEqual(env["OMP_NUM_THREADS"], "64")
        output = Path(cwd) / command[command.index("-root_dir") + 1]
        seed = int(command[command.index("-seed") + 1])
        temperature = output / command[command.index("-Tfile") + 1]
        wave = float(command[command.index("-img") + 1])
        parameter = (Path(cwd) / command[1]).read_text()
        map_line = next(line for line in parameter.splitlines() if "grid (nx,ny), size [AU]" in line)
        nx, ny, size_au = map(float, map_line.split()[:3])
        self.assertEqual(nx, ny)
        self.calls.append({"command": list(command), "temperature_sha256": runner.sha256(temperature),
                           "seed": seed, "size_au": size_au, "npix": int(nx), "attempt": Path(attempt)})
        failed = len(self.calls) == self.fail_call
        write_json(Path(attempt) / "command.json", {"argv": list(command), "cwd": str(cwd), "threads": "64"})
        write_json(Path(attempt) / "execution.json", {"returncode": 1 if failed else 0, "elapsed_seconds": 1})
        log = "Using ray-tracing method 2\nProcessing complete\n"
        (Path(attempt) / "mcfost.log").write_text("ERROR: fixture failure\n" if failed else log)
        if failed:
            raise RuntimeError("fixture failure")
        image = output / f"seed={seed}" / f"data_{wave}" / "RT.fits.gz"
        write_image(image, wave, pixel=size_au / 140 / nx, size=int(nx), invalid=len(self.calls) == self.invalid_call)
        return 1.0, log

    def run_task(self, index=0):
        return self.dispatcher.run_task(self.bundle, index, self.bundle / "machine.template.json")

    def task_status(self, index=0):
        task = self.experiment["tasks"][index]
        return json.loads((self.bundle / task["task_path"] / "status.json").read_text())

    def test_package_is_one_temperature_five_seeds_two_pixel_scales(self):
        experiment = self.experiment
        self.assertEqual(len(experiment["tasks"]), 5)
        self.assertEqual(experiment["schema_version"], 2)
        self.assertEqual(experiment["diagnostic_id"], "fixed_temperature_pixel_scale_v2")
        self.assertEqual(experiment["temperature_solves"], 0)
        self.assertEqual(experiment["image_requests"], 10)
        self.assertEqual({(task["temperature_seed"], task["image_seed"]) for task in experiment["tasks"]},
                         {(42004, seed) for seed in range(42001, 42006)})
        self.assertEqual(len({t["temperature_sha256"] for t in experiment["tasks"]}), 1)
        self.assertEqual(len(list((self.bundle / "inputs/temperatures").glob("*.fits.gz"))), 1)
        self.assertEqual([(g["id"], g["image_npix"], g["image_size_au"]) for g in experiment["geometries"]],
                         [("coarse", 2401, 6000), ("fine", 4801, 6000)])
        for task in experiment["tasks"]:
            self.assertEqual(runner.sha256(self.bundle / task["temperature_path"]), task["temperature_sha256"])
            work = self.bundle / task["task_path"]
            coarse = (work / "coarse/image.para").read_text().splitlines()
            fine = (work / "fine/image.para").read_text().splitlines()
            changes = [(a, b) for a, b in zip(coarse, fine) if a != b]
            self.assertEqual(len(coarse), len(fine))
            self.assertEqual(len(changes), 1)
            self.assertTrue(all("grid (nx,ny), size [AU]" in line for line in changes[0]))
            self.assertEqual(float(changes[0][0].split()[2]), float(changes[0][1].split()[2]))
        machine = json.loads((self.bundle / "machine.template.json").read_text())
        self.assertEqual(machine["max_memory_gb"], 64)
        self.assertEqual(machine["slurm"]["memory"], "96G")
        self.assertEqual(machine["slurm"]["analysis_memory"], "8G")
        with self.assertRaises(FileExistsError):
            self.builder.build_package(self.source, self.bundle)

    def test_runs_only_images_with_exact_fixed_temperature_and_caches(self):
        self.run_task()
        task = self.experiment["tasks"][0]
        self.assertEqual(len(self.calls), 2)
        self.assertEqual({call["size_au"] for call in self.calls}, {6000.})
        self.assertEqual([call["npix"] for call in self.calls], [2401, 4801])
        self.assertEqual({call["temperature_sha256"] for call in self.calls}, {task["temperature_sha256"]})
        self.assertEqual({call["seed"] for call in self.calls}, {task["image_seed"]})
        self.assertEqual(self.task_status()["state"], "complete")
        status, rows = self.dispatcher.inspect_task(self.bundle, self.experiment, task, remeasure=True)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["flux_jy"] > 0 for row in rows))
        self.run_task()
        self.assertEqual(len(self.calls), 2, "Complete results must not trigger a new external execution")
        # The separately frozen analyzer must understand real dispatcher
        # receipts and must not certify the other eight missing images.
        summary = self.analyzer.analyze_bundle(self.bundle, make_plot=False)
        self.assertEqual(summary["status"], "incomplete")
        self.assertEqual(summary["complete_tasks"], 1)
        self.assertEqual(summary["total_tasks"], 5)
        self.assertEqual(summary["validated_images"], 2)
        self.assertEqual(summary["task_state_counts"], {"complete": 1, "not_started": 4})
        self.assertFalse(summary["production_convergence_certified"])
        self.assertFalse(summary["design_complete"])
        self.assertEqual(len(summary["missing_cells"]), 8)
        self.assertEqual(summary["errors"], [])
        self.assertNotIn("temperature_changes", summary)
        self.assertFalse((self.bundle / "results/temperature_changes.csv").exists())

    def test_final_resolution_package_keeps_one_variable_and_declares_stop(self):
        bundle = self.root / "final resolution"
        result = self.builder.build_package(self.source, bundle, final_resolution=True)
        experiment = self.dispatcher.validate_package(bundle)
        self.assertEqual(result["diagnostic_id"], "fixed_temperature_final_resolution_1au_v1")
        self.assertEqual([(g["image_npix"], g["image_size_au"]) for g in experiment["geometries"]],
                         [(4801, 6000), (6001, 6000)])
        self.assertEqual(len({t["temperature_sha256"] for t in experiment["tasks"]}), 1)
        self.assertEqual(experiment["image_requests"], 10)
        policy = experiment["investigation_policy"]
        self.assertTrue(policy["close_after_completed_comparison"])
        self.assertFalse(policy["automatic_followups"])
        self.assertFalse(policy["production_dependency"])
        machine = json.loads((bundle / "machine.template.json").read_text())
        self.assertEqual(machine["max_memory_gb"], 112)
        self.assertEqual(machine["slurm"]["memory"], "160G")
        self.assertEqual(machine["slurm"]["analysis_memory"], "16G")
        for task in experiment["tasks"]:
            work = bundle / task["task_path"]
            first = (work / "coarse/image.para").read_text().splitlines()
            second = (work / "fine/image.para").read_text().splitlines()
            changes = [(a, b) for a, b in zip(first, second) if a != b]
            self.assertEqual(len(changes), 1)
            self.assertTrue(all("grid (nx,ny), size [AU]" in line for line in changes[0]))
        self.assertEqual(self.calls, [])
        # Frozen stopping policy is scientific metadata, not an editable gate.
        experiment["investigation_policy"]["production_dependency"] = True
        write_json(bundle / "experiment.json", experiment)
        (bundle / "experiment.sha256").write_text(runner.sha256(bundle / "experiment.json") + "  experiment.json\n")
        with self.assertRaisesRegex(ValueError, "stopping or reference rule"):
            self.dispatcher.validate_package(bundle)

    def test_final_resolution_rejects_wrong_common_memory_before_execution(self):
        bundle = self.root / "final resolution"
        self.builder.build_package(self.source, bundle, final_resolution=True)
        with self.assertRaisesRegex(ValueError, "112-GB"):
            self.dispatcher.run_task(bundle, 0, bundle / "machine.template.json")
        self.assertEqual(self.calls, [])

    def test_failed_second_image_resumes_without_repeating_first(self):
        self.fail_call = 2
        with self.assertRaisesRegex(RuntimeError, "fixture failure"):
            self.run_task()
        self.assertEqual(self.task_status()["state"], "failed")
        self.fail_call = None
        self.run_task()
        self.assertEqual(self.task_status()["state"], "complete")
        self.assertEqual([call["size_au"] for call in self.calls], [6000., 6000., 6000.])
        self.assertEqual([call["npix"] for call in self.calls], [2401, 4801, 4801])
        self.assertEqual(self.calls[-1]["attempt"].name, "attempt_002")

    def test_invalid_image_is_failed_and_not_accepted_as_cache(self):
        self.invalid_call = 1
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(self.task_status()["state"], "failed")
        self.assertEqual(len(self.calls), 1)

    def test_changed_frozen_temperature_is_rejected_before_execution(self):
        task = self.experiment["tasks"][0]
        path = self.bundle / task["temperature_path"]
        path.write_bytes(path.read_bytes() + b"tampered")
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(self.calls, [])

    def test_changed_image_parameters_are_rejected_before_execution(self):
        task = self.experiment["tasks"][0]
        path = self.bundle / task["task_path"] / "coarse/image.para"
        path.write_text(path.read_text() + "\nchanged after prepare\n")
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(self.calls, [])

    def test_runtime_must_match_the_original_temperature_simulator(self):
        self.runtime["threads"] = 32
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(self.calls, [])
        self.runtime["threads"] = 64
        other = self.root / "other_mcfost"
        other.write_text("different simulator\n")
        self.runtime["mcfost_executable"] = str(other)
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(self.calls, [])

    def test_modified_measurement_receipt_or_image_cache_is_rejected_without_rerunning(self):
        self.run_task()
        task = self.experiment["tasks"][0]
        measurements = self.bundle / task["task_path"] / "measurements.json"
        original = measurements.read_bytes()
        result = json.loads(original)
        result["measurements"][0]["flux_jy"] *= 1.1
        write_json(measurements, result)
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(len(self.calls), 2)
        measurements.write_bytes(original)
        command = self.calls[0]["attempt"] / "command.json"
        original = command.read_bytes()
        edited = json.loads(original)
        edited["threads"] = "32"
        write_json(command, edited)
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(len(self.calls), 2)
        command.write_bytes(original)
        image = next(self.calls[0]["attempt"].rglob("RT.fits.gz"))
        image.write_bytes(image.read_bytes() + b"tampered")
        with self.assertRaises((RuntimeError, ValueError)):
            self.run_task()
        self.assertEqual(len(self.calls), 2)

    def test_builder_accepts_downloaded_flattened_model_folder(self):
        source = self.root / "flattened source"
        shutil.copytree(self.source, source)
        subrun = source / "subruns/p2048000_s42004"
        models = list((subrun / "models").iterdir())
        model = next(path for path in models if path.name.startswith("m0000_"))
        model.rename(subrun / model.name)
        output = self.root / "from flattened download"
        self.builder.build_package(source, output)
        experiment = self.dispatcher.validate_package(output)
        self.assertEqual(len(experiment["tasks"]), 5)

    def test_old_schema_and_temperature_arm_are_rejected(self):
        path = self.bundle / "experiment.json"
        for changed in ({"schema_version": 1}, {"temperature_seeds": [42001, 42004]}):
            exp = {**self.experiment, **changed}
            write_json(path, exp)
            (self.bundle / "experiment.sha256").write_text(runner.sha256(path) + "  experiment.json\n")
            with self.assertRaises(ValueError):
                self.dispatcher.validate_package(self.bundle)


if __name__ == "__main__":
    unittest.main()
