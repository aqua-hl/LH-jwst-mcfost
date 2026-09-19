"""Numerical convergence gates and provenance tests; no simulator required."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
from astropy.io import fits
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("numerics_v2", ROOT / "scripts/analyze_extinction_ice_numerics_v2.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def fixture(bundle, flux_function=None):
    """Exercise real frozen-code/cache/raw-product validation on small fixtures."""
    probes = [{"id": f"n{i:03d}", "wavelength_um": wave} for i, wave in
              enumerate((1.2028117593024095, 2.546352514425122, 3.0764795490679795,
                         3.2961688801068254, 9.7, 18.5), 1)]
    experiment = {"schema_version": 1, "experiment_id": "extinction_ice_numerics_v2",
                  "tasks": [], "probes": probes, "resources": {"cpus_per_task": 64},
                  "acceptance": {"scatter_fraction": .01, "budget_change_fraction": .01, "confidence": .95},
                  "input_hashes": {}}
    for photons in ANALYSIS.BUDGETS:
        for seed in ANALYSIS.SEEDS:
            relative = f"subruns/p{photons:07d}_s{seed}"
            run = bundle / relative
            shutil.copytree(ROOT / "src/mcfost_grid", run / "code/src/mcfost_grid",
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            frozen = run / "inputs/frozen.txt"
            frozen.parent.mkdir()
            frozen.write_text("frozen scientific input\n")
            models = [{"index": i, "id": f"m{i:04d}_fixture", "parameters": {"inclination_deg": inclination}}
                      for i, inclination in enumerate(ANALYSIS.INCLINATIONS)]
            manifest = {"models": models, "anchors": [probe | {"score": False} for probe in probes],
                        "configuration": {"numerical_only": True, "numerics": {
                            "random_seed": seed, "photons_temperature": photons,
                            "photons_image": photons, "photons_sed": photons}},
                        "measurement": {"quality_policy": "whole_image_v1"},
                        "input_hashes": {str(path.relative_to(run)): ANALYSIS.digest(path)
                                         for path in run.rglob("*") if path.is_file()}}
            write_json(run / "manifest.json", manifest)
            manifest_hash = ANALYSIS.digest(run / "manifest.json")
            experiment["input_hashes"][f"{relative}/manifest.json"] = manifest_hash
            fingerprint = {"manifest_sha256": manifest_hash, "executable_sha256": "a" * 64,
                           "utilities_content_sha256": "b" * 64, "backend": "image_method2"}
            write_json(run / "runtime_binding.json", {"fingerprint": fingerprint})
            for model in models:
                directory = run / "models" / model["id"]
                temperature = directory / "temperature/attempt_001/output/data_th/Temperature.fits.gz"
                temperature.parent.mkdir(parents=True)
                fits.PrimaryHDU(np.array([10., 20.])).writeto(temperature)
                temp_attempt = directory / "temperature/attempt_001"
                write_json(temp_attempt / "command.json", {"argv": ["mcfost", "temperature.para", "-seed", str(seed)], "threads": "64"})
                write_json(temp_attempt / "execution.json", {"returncode": 0})
                write_json(directory / "temperature_complete.json", {
                    "fingerprint": fingerprint, "path": str(temperature.relative_to(run)), "sha256": ANALYSIS.digest(temperature)})
                write_json(directory / "runtime_binding.json", {"fingerprint": fingerprint, "runtime": {"threads": 64}})
                write_json(directory / "status.json", {"state": "complete"})
                data = {"complete": True, "model_id": model["id"], "model_index": model["index"],
                        "parameters": model["parameters"], "fingerprint": fingerprint, "measurements": []}
                for probe in probes:
                    attempt = directory / "anchors" / probe["id"] / "attempts/attempt_001"
                    image = attempt / "output/data/RT.fits.gz"
                    image.parent.mkdir(parents=True)
                    image.write_bytes(b"raw image hash fixture")
                    write_json(attempt / "command.json", {"argv": ["mcfost", "image.para", "-img", str(probe["wavelength_um"]),
                        "-seed", str(seed), "-rt2", "-no_T"], "threads": "64"})
                    write_json(attempt / "execution.json", {"returncode": 0})
                    flux = flux_function(photons, seed, model, probe) if flux_function else 1. + .0001 * (seed - 42003)
                    data["measurements"].append({"anchor_id": probe["id"], "wavelength_um": probe["wavelength_um"],
                        "flux_jy": flux, "quality_pass": True, "backend": "image_method2", "image_path": str(image.relative_to(run)),
                        "image_sha256": ANALYSIS.digest(image)})
                write_json(directory / "measurements.json", data)
                experiment["tasks"].append({"index": len(experiment["tasks"]), "run_path": relative,
                    "model_index": model["index"], "model_id": model["id"], "inclination_deg": model["parameters"]["inclination_deg"],
                    "photon_packets": photons, "seed": seed, "manifest_sha256": manifest_hash})
    write_json(bundle / "experiment.json", experiment)
    experiment_hash = ANALYSIS.digest(bundle / "experiment.json")
    (bundle / "experiment.sha256").write_text(experiment_hash + "  experiment.json\n")
    write_json(bundle / "runtime_binding.json", {"experiment_sha256": experiment_hash,
               "fingerprint": {key: fingerprint[key] for key in ANALYSIS.IDENTITY_KEYS} | {"threads": 64}})
    return experiment


class NumericalStatisticsTests(unittest.TestCase):
    def test_five_seed_interval_and_zero_scatter(self):
        result = ANALYSIS.scatter_statistics([.98, .99, 1., 1.01, 1.02])
        expected_sd = np.std([.98, .99, 1., 1.01, 1.02], ddof=1)
        self.assertEqual(result["n_seeds"], 5)
        self.assertAlmostEqual(result["sample_sd_jy"], expected_sd)
        self.assertAlmostEqual(result["sd_upper_fraction"], expected_sd * np.sqrt(4 / chi2.ppf(.025, 4)))
        self.assertAlmostEqual(result["sd_lower_fraction"], expected_sd * np.sqrt(4 / chi2.ppf(.975, 4)))
        self.assertEqual(ANALYSIS.scatter_statistics([1.] * 5)["sd_upper_fraction"], 0.)
        self.assertTrue(ANALYSIS.scatter_statistics([1.] * 5)["scatter_interval_degenerate"])
        self.assertTrue(ANALYSIS.scatter_statistics([.123456789] * 5)["scatter_interval_degenerate"])
        self.assertEqual(ANALYSIS.scatter_statistics([.123456789] * 5)["sample_sd_jy"], 0.)

    def test_pairing_is_by_seed_and_outlier_expands_interval(self):
        low = {42001: 1., 42002: 2., 42003: 3., 42004: 4., 42005: 5.}
        high = {seed: value * 1.005 for seed, value in reversed(list(low.items()))}
        result = ANALYSIS.paired_statistics(low, high)
        self.assertAlmostEqual(result["mean_fractional_change"], .005)
        self.assertLess(result["absolute_change_upper_bound"], .01)
        high[42005] = low[42005] * 1.04
        self.assertGreater(ANALYSIS.paired_statistics(low, high)["absolute_change_upper_bound"], .01)
        high.pop(42005)
        with self.assertRaisesRegex(ValueError, "same seed IDs"):
            ANALYSIS.paired_statistics(low, high)

    def test_upper_interval_gate_is_stricter_than_observed_scatter(self):
        experiment = {"probes": [{"id": "n001", "wavelength_um": 9.7}],
                      "acceptance": {"confidence": .95, "scatter_fraction": .01, "budget_change_fraction": .01}}
        rows = [{"inclination_deg": inclination, "photon_packets": budget, "seed": seed,
                 "probe_id": "n001", "flux_jy": 1 + .003 * (seed - 42003)}
                for inclination in ANALYSIS.INCLINATIONS for budget in ANALYSIS.BUDGETS for seed in ANALYSIS.SEEDS]
        scatter, changes = ANALYSIS.numerical_statistics(experiment, rows, True)
        self.assertTrue(all(r["sample_sd_fraction"] < .01 for r in scatter))
        self.assertTrue(all(r["passes_scatter"] is False for r in scatter))
        self.assertTrue(all(r["passes_budget_change"] is True for r in changes))

    def test_zero_observed_scatter_is_inconclusive(self):
        experiment = {"probes": [{"id": "n001", "wavelength_um": 9.7}],
                      "acceptance": {"confidence": .95, "scatter_fraction": .01, "budget_change_fraction": .01}}
        rows = [{"inclination_deg": inclination, "photon_packets": budget, "seed": seed,
                 "probe_id": "n001", "flux_jy": 1.}
                for inclination in ANALYSIS.INCLINATIONS for budget in ANALYSIS.BUDGETS for seed in ANALYSIS.SEEDS]
        scatter, _ = ANALYSIS.numerical_statistics(experiment, rows, True)
        self.assertTrue(all(row["passes_scatter"] is None for row in scatter))
        self.assertTrue(all(row["scatter_interval_degenerate"] for row in scatter))


class BundleIntegrityTests(unittest.TestCase):
    def test_complete_bundle_passes_and_writes_readable_reports(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            fixture(bundle)
            result = ANALYSIS.analyze_bundle(bundle)
            self.assertEqual(result["status"], "passed")
            self.assertTrue(result["numerical_gate_passed"])
            self.assertEqual(result["complete_tasks"], 20)
            self.assertEqual(len(result["scatter"]), 24)
            self.assertEqual(len(result["budget_changes"]), 12)
            self.assertFalse(result["model_ranking_performed"])
            for name in ("numerical_summary.json", "seed_fluxes.csv", "numerical_scatter.csv", "budget_changes.csv",
                         "numerical_convergence.png", "numerical_convergence.pdf"):
                self.assertGreater((bundle / "results" / name).stat().st_size, 0)

    def test_missing_task_never_passes_or_silently_drops_seed(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment = fixture(bundle)
            task = experiment["tasks"][0]
            directory = bundle / task["run_path"] / "models" / task["model_id"]
            (directory / "measurements.json").unlink()
            (directory / "status.json").unlink()
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "incomplete")
            self.assertFalse(result["numerical_gate_passed"])
            self.assertEqual(result["complete_tasks"], 19)
            self.assertEqual(len(result["tasks"]), 20)
            self.assertEqual(result["tasks"][0]["state"], "not_started")
            self.assertTrue(all(row["passes_budget_change"] is None for row in result["budget_changes"]))

    def test_mixed_runtime_across_complete_subruns_refuses_inference(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment = fixture(bundle)
            run = bundle / experiment["tasks"][0]["run_path"]
            paths = [run / "runtime_binding.json"]
            for model in ("m0000_fixture", "m0001_fixture"):
                paths.extend(run / "models" / model / name for name in
                             ("runtime_binding.json", "measurements.json", "temperature_complete.json"))
            for path in paths:
                data = ANALYSIS.read_json(path)
                data["fingerprint"]["executable_sha256"] = "c" * 64
                write_json(path, data)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["complete_tasks"], 20)
            self.assertEqual(result["status"], "integrity_error")
            self.assertFalse(result["inference_eligible"])
            self.assertTrue(any("Mixed executable" in error for error in result["errors"]))

    def test_raw_corruption_wrong_threads_and_frozen_input_are_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment = fixture(bundle)
            first, second, third = experiment["tasks"][0], experiment["tasks"][1], experiment["tasks"][2]
            directory = bundle / first["run_path"] / "models" / first["model_id"]
            measurement = ANALYSIS.read_json(directory / "measurements.json")["measurements"][0]
            (bundle / first["run_path"] / measurement["image_path"]).write_bytes(b"changed raw image")
            directory = bundle / second["run_path"] / "models" / second["model_id"]
            path = directory / "temperature/attempt_001/command.json"
            command = ANALYSIS.read_json(path)
            command["threads"] = "32"
            write_json(path, command)
            (bundle / third["run_path"] / "inputs/frozen.txt").write_text("changed input")
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertFalse(result["numerical_gate_passed"])
            self.assertIn("Cached image changed", result["tasks"][0]["reason"])
            self.assertIn("Wrong thread count", result["tasks"][1]["reason"])
            self.assertIn("Frozen input changed", result["tasks"][2]["reason"])
            self.assertEqual(result["status"], "integrity_error")

    def test_checksum_top_binding_backend_and_empty_input_contract_are_enforced(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment = fixture(bundle)
            (bundle / "runtime_binding.json").unlink()
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertFalse(result["numerical_gate_passed"])
            self.assertTrue(any("runtime_binding.json is required" in error for error in result["errors"]))
            (bundle / "experiment.sha256").write_text("0" * 64 + "  experiment.json\n")
            task = experiment["tasks"][0]
            measurement = bundle / task["run_path"] / "models" / task["model_id"] / "measurements.json"
            data = ANALYSIS.read_json(measurement)
            data["fingerprint"]["backend"] = "coeval_method2"
            write_json(measurement, data)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertTrue(any("experiment.sha256" in error for error in result["errors"]))
            self.assertIn("requires the image_method2 backend", result["tasks"][0]["reason"])
            experiment["input_hashes"] = {}
            with self.assertRaisesRegex(ValueError, "nonempty frozen input_hashes"):
                ANALYSIS.validate_experiment(experiment)

    def test_unstarted_bundle_does_not_require_runtime_binding(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment = fixture(bundle)
            (bundle / "runtime_binding.json").unlink()
            for task in experiment["tasks"]:
                directory = bundle / task["run_path"] / "models" / task["model_id"]
                (directory / "measurements.json").unlink()
                (directory / "status.json").unlink()
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["errors"], [])
            self.assertEqual(result["complete_tasks"], 0)


if __name__ == "__main__":
    unittest.main()
