"""Fixed-temperature pixel-scale accounting and frozen-code integration."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.stats import chi2, t


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("crossover_analysis", ROOT / "scripts/analyze_extinction_ice_crossover_v2.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def pixel_scale_rows():
    rows = []
    for image_index, image_seed in enumerate(ANALYSIS.IMAGE_SEEDS):
        for geometry in ANALYSIS.GEOMETRIES:
            flux = (.98 if geometry == "fine" else 1.) * (1 + .001 * (image_seed - 42003))
            rows.append({"task_index": image_index, "temperature_seed": ANALYSIS.TEMPERATURE_SEED,
                         "image_seed": image_seed, "geometry": geometry, "flux_jy": flux,
                         "components_aperture_jy": {"total_i": flux, "direct_star": .3 * flux,
                                                    "scattered_star": .4 * flux, "direct_thermal": .1 * flux,
                                                    "scattered_thermal": .2 * flux},
                         "image_npix": 2401 if geometry == "coarse" else 4801, "image_size_au": 6000,
                         "threads": 64, "host": "worker0", "quality_pass": True,
                         "temperature_sha256": "0" * 64, "image_sha256": "a" * 64,
                         "component_sum_minus_total_aperture_fraction": 0.})
    return rows


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False) + "\n")


def stub_bundle(bundle, final_resolution=False):
    """The real dispatcher is tested separately; this verifies its analyzer API.

    Pin a private dispatcher which refuses calls without remeasurement enabled,
    and delivers synthetic validated products through the real import path.
    """
    code = bundle / "code/crossover_task.py"
    code.parent.mkdir(parents=True, exist_ok=True)
    code.write_text('''import json
from pathlib import Path

def validate_package(bundle):
    return json.loads((Path(bundle) / "experiment.json").read_text())

def inspect_task(bundle, experiment, task, remeasure=True):
    if remeasure is not True:
        raise RuntimeError("Raw remeasurement must be requested")
    path = Path(bundle) / "fixture_results.json"
    results = json.loads(path.read_text())
    entry = results[str(task["index"])]
    return dict(index=task["index"], state=entry["state"]), entry["rows"]
''')
    rows = pixel_scale_rows()
    if final_resolution:
        for row in rows:
            row["image_npix"] = 4801 if row["geometry"] == "coarse" else 6001
    tasks = [{"index": image_index, "temperature_seed": ANALYSIS.TEMPERATURE_SEED,
              "image_seed": image_seed, "temperature_sha256": "0" * 64,
              "task_path": f"tasks/t{image_index:03d}"}
             for image_index, image_seed in enumerate(ANALYSIS.IMAGE_SEEDS)]
    experiment = {"schema_version": 2,
                  "diagnostic_id": ANALYSIS.FINAL_DIAGNOSTIC_ID if final_resolution else ANALYSIS.DIAGNOSTIC_ID,
                  "experiment_id": "extinction_ice_crossover_v2", "tasks": tasks,
                  "input_hashes": {"code/crossover_task.py": ANALYSIS.digest(code)},
                  "resources": {"cpus_per_task": 64}}
    write_json(bundle / "experiment.json", experiment)
    (bundle / "experiment.sha256").write_text(ANALYSIS.digest(bundle / "experiment.json") + "  experiment.json\n")
    results = {str(task["index"]): {"state": "complete", "rows": [row for row in rows if row["task_index"] == task["index"]]}
               for task in tasks}
    write_json(bundle / "fixture_results.json", results)
    return experiment, results


class CrossoverStatisticsTests(unittest.TestCase):
    def test_final_reference_rule_records_outcome_without_certifying_production(self):
        rows = pixel_scale_rows()
        for index, row in enumerate(rows):
            if row["geometry"] == "fine":
                row["flux_jy"] = rows[index - 1]["flux_jy"] * (1 + .002 + index * .00003)
        summary = ANALYSIS.summarize(rows)
        result = ANALYSIS.final_reference_result(summary)
        self.assertEqual(result["outcome"], "met")
        self.assertTrue(result["investigation_closed"])
        self.assertFalse(result["production_dependency"])
        self.assertFalse(result["automatic_followups"])
        self.assertFalse(summary["production_convergence_certified"])
        # A confidently measured 2% shift fails the reference but ends inquiry.
        result = ANALYSIS.final_reference_result(ANALYSIS.summarize(pixel_scale_rows()))
        self.assertEqual(result["outcome"], "not_met")
        self.assertTrue(result["investigation_closed"])
        incomplete = ANALYSIS.final_reference_result(ANALYSIS.summarize(rows[:-1]))
        self.assertEqual(incomplete["outcome"], "not_evaluated")
        self.assertFalse(incomplete["investigation_closed"])

    def test_degenerate_or_unmatched_pairs_do_not_establish_final_reference(self):
        rows = pixel_scale_rows()
        for row in rows:
            row["flux_jy"] = 1.
        result = ANALYSIS.final_reference_result(ANALYSIS.summarize(rows))
        self.assertEqual(result["outcome"], "not_met")
        self.assertFalse(result["criteria"]["both_total_i_scatter_upper_bounds_below_reference"])
        self.assertFalse(result["criteria"]["paired_total_i_change_interval_inside_reference"])
        rows[0]["host"] = "changed-host"
        result = ANALYSIS.final_reference_result(ANALYSIS.summarize(rows))
        self.assertFalse(result["criteria"]["all_pairs_same_host"])

    def test_sd_interval_uses_five_draws_and_zero_is_degenerate(self):
        values = [.98, .99, 1., 1.01, 1.02]
        result = ANALYSIS.scatter_statistics(values)
        sd = np.std(values, ddof=1)
        self.assertEqual(result["n_seeds"], 5)
        self.assertAlmostEqual(result["sample_sd_jy"], sd)
        self.assertAlmostEqual(result["sd_lower_fraction"], sd * np.sqrt(4 / chi2.ppf(.975, 4)))
        self.assertAlmostEqual(result["sd_upper_fraction"], sd * np.sqrt(4 / chi2.ppf(.025, 4)))
        self.assertTrue(ANALYSIS.scatter_statistics([.123456789] * 5)["scatter_interval_degenerate"])
        self.assertEqual(ANALYSIS.scatter_statistics([.123456789] * 5)["sample_sd_jy"], 0.)

    def test_pairing_keeps_seed_labels_and_outlier(self):
        reference = {42001: 1., 42002: 2., 42003: 3., 42004: 4., 42005: 5.}
        comparison = {seed: value * 1.005 for seed, value in reversed(list(reference.items()))}
        comparison[42004] = reference[42004] * .9
        result = ANALYSIS.paired_statistics(reference, comparison)
        differences = np.array([.005, .005, .005, -.1, .005])
        mean, sd = differences.mean(), differences.std(ddof=1)
        half = t.ppf(.975, 4) * sd / np.sqrt(5)
        self.assertEqual(result["n_pairs"], 5)
        self.assertAlmostEqual(result["mean_fractional_change"], mean)
        self.assertAlmostEqual(result["ci_lower"], mean - half)
        self.assertAlmostEqual(result["ci_upper"], mean + half)
        comparison.pop(42005)
        with self.assertRaisesRegex(ValueError, "same seed IDs"):
            ANALYSIS.paired_statistics(reference, comparison)

    def test_only_pixel_scale_changes_and_no_production_pass_exists(self):
        result = ANALYSIS.summarize(reversed(pixel_scale_rows()))
        self.assertEqual(result["status"], "complete_diagnostic")
        self.assertTrue(result["design_complete"])
        self.assertEqual(result["validated_images"], 10)
        self.assertEqual(len(result["scatter"]), 2)
        self.assertTrue(all(row["complete_five_seed_group"] for row in result["scatter"]))
        for row in result["spatial_changes"]:
            self.assertAlmostEqual(row["mean_fractional_change"], -.02)
            self.assertEqual(row["same_host_pairs"], 5)
            self.assertAlmostEqual(row["components"]["direct_thermal"]["mean_fractional_change"], -.02)
        self.assertEqual(len(result["spatial_changes"]), 1)
        self.assertNotIn("temperature_changes", result)
        self.assertFalse(result["production_convergence_certified"])
        self.assertEqual(result["outliers_removed"], 0)

    def test_missing_cell_remains_explicit_in_partial_statistics(self):
        rows = pixel_scale_rows()
        rows.pop()
        result = ANALYSIS.summarize(rows)
        self.assertEqual(result["status"], "incomplete")
        self.assertFalse(result["design_complete"])
        self.assertEqual(result["missing_cells"], [{"temperature_seed": 42004, "image_seed": 42005, "geometry": "fine"}])
        self.assertEqual(result["scatter"][-1]["n_seeds"], 4)
        self.assertFalse(result["scatter"][-1]["complete_five_seed_group"])
        self.assertEqual(result["spatial_changes"][-1]["n_pairs"], 4)
        self.assertFalse(result["spatial_changes"][-1]["complete_five_seed_comparison"])

    def test_duplicates_bad_factors_and_nonfinite_values_do_not_disappear(self):
        rows = pixel_scale_rows()
        rows.append(copy.deepcopy(rows[0]))
        rows[2]["temperature_seed"] = 1234
        rows[4]["flux_jy"] = float("nan")
        result = ANALYSIS.summarize(rows)
        self.assertEqual(result["status"], "integrity_error")
        self.assertEqual(len(result["rows"]), 11)
        self.assertEqual(result["validated_images"], 7)
        self.assertEqual(len(result["invalid_rows"]), 4)
        self.assertEqual(len(result["missing_cells"]), 3)
        self.assertTrue(any("Duplicate pixel-scale cell" in error for error in result["errors"]))

    def test_signed_components_and_zero_reference_are_preserved(self):
        rows = pixel_scale_rows()
        for row in rows:
            row["components_aperture_jy"]["direct_thermal"] = 0. if row["geometry"] == "coarse" else -.001
        result = ANALYSIS.summarize(rows)
        comparison = result["spatial_changes"][0]["components"]["direct_thermal"]
        self.assertEqual(comparison["mean_delta_jy"], -.001)
        self.assertFalse(comparison["fractional_comparison_defined"])
        self.assertIsNone(comparison["mean_fractional_change"])
        self.assertEqual(result["rows"][1]["components_aperture_jy"]["direct_thermal"], -.001)


class CrossoverAnalyzerIntegrationTests(unittest.TestCase):
    def test_final_resolution_report_ends_investigation_even_when_reference_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            stub_bundle(bundle, final_resolution=True)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "complete_diagnostic")
            self.assertEqual(result["diagnostic_id"], ANALYSIS.FINAL_DIAGNOSTIC_ID)
            self.assertEqual(result["final_reference_result"]["outcome"], "not_met")
            self.assertTrue(result["investigation_closed"])
            self.assertFalse(result["production_dependency"])
            self.assertFalse(result["automatic_followups"])
            self.assertFalse(result["production_convergence_certified"])
            review = (bundle / "results/REVIEW.md").read_text()
            self.assertIn("4801 to 6001", review)
            self.assertIn("0.999833 AU/pixel", review)
            self.assertNotIn("2401 and 4801", review)
            # Separate modes cannot overwrite each other's previous reports.
            stub_bundle(bundle)
            with self.assertRaisesRegex(ValueError, "another diagnostic"):
                ANALYSIS.analyze_bundle(bundle, make_plot=False)

    def test_frozen_inspector_is_used_with_remeasurement_and_reports_written(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            stub_bundle(bundle)
            result = ANALYSIS.analyze_bundle(bundle)
            self.assertEqual(result["status"], "complete_diagnostic")
            self.assertEqual(result["complete_tasks"], 5)
            self.assertTrue(result["remeasurement_performed"])
            self.assertFalse(result["production_convergence_certified"])
            self.assertEqual(result["fixed_temperature_seed"], 42004)
            self.assertEqual(result["fixed_temperature_sha256"], "0" * 64)
            self.assertEqual(result["new_temperature_solves"], 0)
            self.assertNotIn("temperature_changes", result)
            for name in ("summary.json", "seed_fluxes.csv", "scatter.csv", "spatial_changes.csv",
                         "REVIEW.md", "crossover_diagnostics.png", "crossover_diagnostics.pdf"):
                self.assertGreater((bundle / "results" / name).stat().st_size, 0)
            self.assertFalse((bundle / "results/temperature_changes.csv").exists())

    def test_incomplete_marker_and_running_task_do_not_complete_design(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            _, results = stub_bundle(bundle)
            results["0"]["state"] = "partial"
            results["1"] = {"state": "running", "rows": []}
            write_json(bundle / "fixture_results.json", results)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["complete_tasks"], 3)
            self.assertEqual(result["validated_images"], 8)
            self.assertEqual(result["tasks"][1]["state"], "running")
            self.assertEqual(result["errors"], [])

    def test_resumed_worker_change_is_reported_without_discarding_pair(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            _, results = stub_bundle(bundle)
            results["0"]["rows"][1]["host"] = "worker_after_resume"
            write_json(bundle / "fixture_results.json", results)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "complete_diagnostic")
            self.assertFalse(result["tasks"][0]["pair_host_matched"])
            self.assertEqual(result["spatial_changes"][0]["same_host_pairs"], 4)
            self.assertEqual(len(result["warnings"]), 1)
            self.assertEqual(result["validated_images"], 10)

    def test_wrong_temperature_hash_or_components_invalidate_task(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            _, results = stub_bundle(bundle)
            results["0"]["rows"][0]["temperature_sha256"] = "f" * 64
            results["1"]["rows"][0]["components_aperture_jy"].pop("direct_thermal")
            write_json(bundle / "fixture_results.json", results)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "integrity_error")
            self.assertEqual(result["complete_tasks"], 3)
            self.assertEqual(len(result["errors"]), 2)
            self.assertEqual(result["tasks"][0]["state"], "invalid")

    def test_old_schema_extra_temperature_arm_and_mixed_hashes_are_rejected(self):
        for variation, message in (("old_schema", "Expected schema 2"),
                                   ("extra_arm", "exactly five"),
                                   ("mixed_hashes", "same valid temperature")):
            with self.subTest(variation=variation), tempfile.TemporaryDirectory() as temporary:
                bundle = Path(temporary)
                experiment, _ = stub_bundle(bundle)
                if variation == "old_schema":
                    experiment["schema_version"] = 1
                elif variation == "extra_arm":
                    extra = copy.deepcopy(experiment["tasks"])
                    for task in extra:
                        task["index"] += 5
                        task["temperature_seed"] = 42001
                    experiment["tasks"].extend(extra)
                else:
                    experiment["tasks"][1]["temperature_sha256"] = "1" * 64
                write_json(bundle / "experiment.json", experiment)
                (bundle / "experiment.sha256").write_text(
                    ANALYSIS.digest(bundle / "experiment.json") + "  experiment.json\n")
                with self.assertRaisesRegex(ValueError, message):
                    ANALYSIS.analyze_bundle(bundle, make_plot=False)
                self.assertFalse((bundle / "results").exists())

    def test_changing_field_of_view_is_not_a_valid_pixel_scale_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            _, results = stub_bundle(bundle)
            results["0"]["rows"][1]["image_npix"] = 2401
            results["0"]["rows"][1]["image_size_au"] = 3000
            write_json(bundle / "fixture_results.json", results)
            result = ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(result["status"], "integrity_error")
            self.assertEqual(result["validated_images"], 8)
            self.assertIn("geometry differs", result["tasks"][0]["error"])

    def test_old_diagnostic_reports_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            stub_bundle(bundle)
            old = bundle / "results/temperature_changes.csv"
            old.parent.mkdir()
            old.write_text("old temperature-arm analysis\n")
            with self.assertRaisesRegex(ValueError, "another diagnostic"):
                ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertEqual(old.read_text(), "old temperature-arm analysis\n")
            self.assertFalse((bundle / "results/summary.json").exists())

    def test_modified_code_is_rejected_before_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            stub_bundle(bundle)
            (bundle / "code/crossover_task.py").write_text('raise RuntimeError("must never execute changed code")\n')
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                ANALYSIS.analyze_bundle(bundle, make_plot=False)
            self.assertFalse((bundle / "results").exists())

    def test_changed_manifest_and_unsafe_report_destination_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            experiment, _ = stub_bundle(bundle)
            with self.assertRaisesRegex(ValueError, "protect frozen inputs"):
                ANALYSIS.analyze_bundle(bundle, output=bundle / "code", make_plot=False)
            experiment["tasks"][0]["image_seed"] = 999
            write_json(bundle / "experiment.json", experiment)
            with self.assertRaisesRegex(ValueError, "Experiment hash differs"):
                ANALYSIS.analyze_bundle(bundle, make_plot=False)


if __name__ == "__main__":
    unittest.main()
