"""Scientific ranking-contract tests on compact fixtures, without MCFOST."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mcfost_grid.analysis import AnalysisError, analyze_run, evaluate_model, region_balanced_log_rms, validate_anchors


def anchors():
    return [
        {"id": "a", "wavelength_um": 1.2, "flux_jy": 1., "uncertainty_jy": .1, "region": "blue", "score": True},
        {"id": "b", "wavelength_um": 2., "flux_jy": 1., "uncertainty_jy": .1, "region": "blue", "score": True},
        {"id": "c", "wavelength_um": 5., "flux_jy": 1., "uncertainty_jy": .1, "region": "red", "score": True},
        {"id": "sil", "wavelength_um": 9.7, "flux_jy": 1., "uncertainty_jy": .1, "region": "silicate", "score": False},
    ]


def model(index=0):
    return {"id": f"m{index}", "index": index,
            "parameters": {"mass": float(index+1), "inclination_deg": 50.+10.*index}}


def payload(m, fluxes=(1., 1., 1., 1.)):
    return {"model_id": m["id"], "model_index": m["index"], "parameters": m["parameters"],
            "complete": True, "smoke_test": False,
            "measurements": [{"anchor_id": a["id"], "wavelength_um": a["wavelength_um"], "flux_jy": f, "quality_pass": True}
                             for a, f in zip(anchors(), fluxes)]}


def pinned_fixture(root, count=1):
    """Create compact provenance fixtures; grey spectrum parsing is mocked below."""
    (root / "inputs").mkdir()
    observed = root / "inputs/observed_spectrum.ecsv"
    observed.write_text("immutable observation fixture\n")
    manifest = {"run_id": "pinned", "anchors": anchors(), "models": [model(i) for i in range(count)],
                "input_hashes": {"inputs/observed_spectrum.ecsv": hashlib.sha256(observed.read_bytes()).hexdigest()}}
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    files = []
    for m in manifest["models"]:
        directory = root / "models" / m["id"]
        directory.mkdir(parents=True)
        data = payload(m)
        data["fingerprint"] = {"manifest_sha256": manifest_hash, "executable_sha256": "a"*64,
                               "utilities_content_sha256": "c"*64, "backend": "image_method2"}
        filename = directory / "measurements.json"
        filename.write_text(json.dumps(data))
        files.append(filename)
    return observed, files


class AnalysisContractTests(unittest.TestCase):
    def test_region_balance_not_point_count(self):
        self.assertAlmostEqual(region_balanced_log_rms({"a": [1., 1.], "b": [0.]}), math.sqrt(.5))
        result = evaluate_model(model(), anchors(), payload(model(), (10., 10., 1., 1.)))
        self.assertAlmostEqual(result["score_dex"], math.sqrt(.5))

    def test_missing_hard_anchor_excludes_even_if_complete_claimed(self):
        data = payload(model())
        data["measurements"] = data["measurements"][:2]
        result = evaluate_model(model(), anchors(), data)
        self.assertEqual(result["status"], "incomplete")
        self.assertIsNone(result["score_dex"])

    def test_unscored_missing_or_bad_anchor_does_not_change_score(self):
        data = payload(model(), (1., 1., 1., 1e9))
        self.assertEqual(evaluate_model(model(), anchors(), data)["score_dex"], 0.)
        data["measurements"].pop()
        self.assertEqual(evaluate_model(model(), anchors(), data)["score_dex"], 0.)

    def test_bad_quality_nonpositive_and_wrong_wavelength_exclude(self):
        for edit in ({"quality_pass": False}, {"flux_jy": 0.}, {"flux_jy": float("nan")}, {"wavelength_um": 1.3}):
            with self.subTest(edit=edit):
                data = payload(model())
                data["measurements"][0].update(edit)
                self.assertEqual(evaluate_model(model(), anchors(), data)["status"], "incomplete")

    def test_identity_and_duplicate_rows_rejected(self):
        data = payload(model())
        data["model_id"] = "another"
        self.assertEqual(evaluate_model(model(), anchors(), data)["status"], "invalid")
        data = payload(model())
        data["measurements"].append(data["measurements"][0])
        self.assertEqual(evaluate_model(model(), anchors(), data)["status"], "invalid")

    def test_no_chi_square_from_partial_error_subset(self):
        obs = anchors()
        obs[1]["uncertainty_jy"] = 0.
        result = evaluate_model(model(), obs, payload(model()))
        self.assertEqual(result["status"], "ranked")
        self.assertIsNone(result["chi2_diagnostic"])

    def test_common_observation_contract_invalid(self):
        obs = anchors()
        obs[0]["flux_jy"] = 0.
        with self.assertRaises(AnalysisError):
            validate_anchors(obs)

    def test_edited_pinned_observations_fail_before_report_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            observed, _ = pinned_fixture(root)
            observed.write_text("changed observation\n")
            with self.assertRaisesRegex(AnalysisError, "Observation spectrum changed"):
                analyze_run(root)
            self.assertFalse((root / "results").exists())

    @patch("mcfost_grid.analysis._plots", return_value=[])
    def test_pinned_measurement_requires_matching_fingerprint(self, _plots):
        for mutation in ("absent", "stale", "bad_backend"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                _, files = pinned_fixture(root)
                data = json.loads(files[0].read_text())
                if mutation == "absent":
                    data.pop("fingerprint")
                elif mutation == "stale":
                    data["fingerprint"]["manifest_sha256"] = "b"*64
                else:
                    data["measurements"][0]["backend"] = "coeval_method2"
                files[0].write_text(json.dumps(data))
                result = analyze_run(root)
                self.assertEqual(result["status_counts"], {"invalid": 1})
                self.assertIsNone(result["best_model"])

    @patch("mcfost_grid.analysis._plots", return_value=[])
    def test_mixed_runtime_models_are_all_excluded(self, _plots):
        for key, other in (("executable_sha256", "b"*64), ("backend", "coeval_method2"),
                           ("utilities_content_sha256", "d"*64)):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                _, files = pinned_fixture(root, count=2)
                data = json.loads(files[1].read_text())
                data["fingerprint"][key] = other
                files[1].write_text(json.dumps(data))
                result = analyze_run(root, workers=2)
                self.assertEqual(result["status_counts"], {"invalid": 2})
                self.assertIsNone(result["best_model"])
                self.assertIn("Mixed executable", result["warnings"][0])

    @patch("mcfost_grid.analysis._plots", return_value=[])
    def test_matching_pinned_fingerprints_rank_normally(self, _plots):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pinned_fixture(root, count=2)
            result = analyze_run(root)
            self.assertEqual(result["status_counts"], {"ranked": 2})

    def test_all_missing_run_is_reportable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "manifest.json").write_text(json.dumps({"run_id": "missing", "anchors": anchors(), "models": [model()]}))
            result = analyze_run(root)
            self.assertIsNone(result["best_model"])
            self.assertEqual(result["status_counts"], {"missing": 1})
            self.assertTrue(result["warnings"])

    def test_run_writes_compact_reports_and_excludes_partial_model(self):
        from astropy.table import Table
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "inputs").mkdir()
            Table({"wavelength_um": [1., 2., 4., 10.], "flux_jy": [.8, 1., 1.1, .9]}).write(
                root / "inputs/observed_spectrum.ecsv", format="ascii.ecsv")
            manifest = {"schema_version": 1, "run_id": "test", "configuration": {"smoke_test": True},
                        "models": [model(0), model(1)], "anchors": anchors()}
            (root / "manifest.json").write_text(json.dumps(manifest))
            for m in manifest["models"]:
                target = root / "models" / m["id"]
                target.mkdir(parents=True)
                data = payload(m, (1.2, 1.2, 1.2, .1))
                if m["index"] == 1:
                    data["measurements"].pop(0)
                (target / "measurements.json").write_text(json.dumps(data))
            result = analyze_run(root, workers=2)
            self.assertEqual(result["ranked_model_count"], 1)
            self.assertEqual(result["best_model"]["model_id"], "m0")
            self.assertTrue(result["smoke_test"])
            self.assertAlmostEqual(result["best_model"]["score_dex"], math.log10(1.2))
            for filename in ("summary.json", "ranking.ecsv", "predictions.ecsv", "spectrum_comparison.png",
                             "spectrum_comparison.pdf", "parameter_space.png", "parameter_space.pdf"):
                self.assertGreater((root / "results" / filename).stat().st_size, 0)
            table = Table.read(root / "results/ranking.ecsv", format="ascii.ecsv")
            self.assertEqual(list(table["rank"]), [1, 0])


if __name__ == "__main__":
    unittest.main()
