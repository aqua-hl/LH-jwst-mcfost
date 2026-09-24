"""Controlled-pair completeness, feature definitions and receipt rejection."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ANALYSIS = load("structure_analysis", ROOT / "scripts/analyze_silicate_structure_production.py")
FIXTURES = load("pilot_analysis_fixtures", ROOT / "tests/test_silicate_size_analysis.py")
OBSERVABLES = load("structure_observables", ROOT / "scripts/production_observables.py")
PRODUCTION = ANALYSIS.production


def models_fixture(catalogue, contract):
    result = []
    for record in catalogue:
        changed = record["parameters"]["envelope_amax_um"] == 1.
        vector = [b["observed_flux_jy"] * (.8 if changed else .5) for b in contract["bands"]]
        result.append({"model_index": record["index"], "model_id": f"m{record['index']}",
            "parameters": copy.deepcopy(record["parameters"]), "status": "screened",
            "band_flux_jy": vector, "baseline_score": 999. if changed else 1.,
            "quadrature_unresolved": False, "unresolved_quadrature_bands": []})
    return result


def compact_campaign(run, contract):
    manifest = FIXTURES.compact_fixture(run)
    template = PRODUCTION.read_json(run / "models/m0/measurements.json")
    catalogue = ANALYSIS.design.catalogue()
    FIXTURES.write_json(run / "production_experiment.json", {
        "catalogue": catalogue, "dust_prescription": ANALYSIS.design.DUST,
        "material_input_check": {"coverage_policy": "exact_tables_native_extrapolation",
                                  "policy_accepted": True},
        "predeclared_contrasts": ANALYSIS.predeclared_contrasts(catalogue)})
    FIXTURES.write_json(run / "inputs/production_observation_contract.json", contract)
    manifest["models"] = [{"id": f"m{r['index']}", "index": r["index"], "parameters": r["parameters"]}
                          for r in catalogue]
    manifest["anchors"] = contract["probes"]
    manifest["input_hashes"] = {relative: PRODUCTION.digest(run / relative)
                                for relative in manifest["input_hashes"]}
    FIXTURES.write_json(run / "manifest.json", manifest)
    fingerprint = {**template["fingerprint"], "manifest_sha256": PRODUCTION.digest(run / "manifest.json")}
    FIXTURES.write_json(run / "runtime_binding.json", {"fingerprint": fingerprint})
    for model in manifest["models"]:
        measurements = []
        for probe in contract["probes"]:
            receipt = copy.deepcopy(template["measurements"][0])
            receipt.update(anchor_id=probe["id"], wavelength_um=probe["wavelength_um"],
                           instrument=probe["instrument"], flux_jy=1. + .01 * model["index"])
            measurements.append(receipt)
        FIXTURES.write_json(run / "models" / model["id"] / "measurements.json", {
            "model_id": model["id"], "model_index": model["index"], "parameters": model["parameters"],
            "fingerprint": fingerprint, "complete": True, "measurements": measurements})


class MatchedCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = OBSERVABLES.build_contract(ROOT)

    def setUp(self):
        self.catalogue = ANALYSIS.design.catalogue()
        self.models = models_fixture(self.catalogue, self.contract)

    def test_all_predeclared_arms_and_pairs_ignore_scores(self):
        result = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(result["status"], "complete_campaign")
        self.assertEqual(result["pairs_complete"], 75)
        self.assertEqual(len(result["band_responses"]), 75 * 19)
        self.assertEqual(len(result["feature_responses"]), 75 * 5)
        pair = next(p for p in result["pairs"] if p["kind"] == "size")
        self.assertEqual(pair["baseline_profiled_scores"], [1., 999.])
        self.assertAlmostEqual(pair["s02"]["changed_over_reference_flux_ratio"], 1.6)
        self.assertAlmostEqual(pair["s02"]["absolute_residual_improvement_percentage_points"], 30.)
        self.assertAlmostEqual(pair["feature_shapes"][0]["absolute_shape_residual_improvement_dex"], 0.)
        for pair in (p for p in result["pairs"] if p["kind"] == "density_profile"):
            self.assertEqual(pair["reference_column_mass_msun"], pair["changed_column_mass_msun"])
            self.assertNotEqual(pair["reference_actual_mass_msun"], pair["changed_actual_mass_msun"])

    def test_power_law_continuum_and_scale_invariant_feature_diagnostic(self):
        bands = self.contract["bands"]
        flux = np.array([(.5 * (b["lower_um"] + b["upper_um"]))**2 for b in bands])
        for index, band in enumerate(bands):
            if band["id"] in {"s02", "h01", "h02", "h03", "h04"}:
                flux[index] *= .2
        shapes = ANALYSIS.feature_shapes(flux, bands)
        scaled = ANALYSIS.feature_shapes(7 * flux, bands)
        for name, values in shapes.items():
            self.assertAlmostEqual(values["core_over_continuum"], .2)
            self.assertAlmostEqual(values["core_over_continuum"], scaled[name]["core_over_continuum"])
        with self.assertRaisesRegex(ValueError, "finite and positive"):
            ANALYSIS.feature_shapes(-flux, bands)

    def test_missing_pair_and_feature_quadrature_in_shoulders_remain_explicit(self):
        self.models[0].update(status="excluded", reason="Nonpositive signed total-I")
        self.models[1].update(quadrature_unresolved=True, unresolved_quadrature_bands=["s01"])
        result = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(result["status"], "partial_campaign")
        self.assertEqual(result["pairs_complete"] + len(result["missing_pairs"]), 75)
        self.assertTrue(all(0 not in (p["reference_index"], p["changed_index"]) for p in result["pairs"]))
        affected = [row for row in result["feature_responses"]
                    if row["changed_index"] == 1 and row["feature_id"] == "silicate_9p7"]
        self.assertTrue(affected)
        self.assertTrue(all(row["changed_quadrature_unresolved"] for row in affected))
        self.catalogue[-1]["parameters"]["envelope_density_exponent"] = -1.6
        with self.assertRaisesRegex(ValueError, "catalogue differs"):
            ANALYSIS.predeclared_contrasts(self.catalogue)

    def test_compact_receipts_all_scenarios_artifacts_and_mandatory_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_campaign(run, self.contract)
            original_contrasts = PRODUCTION.controlled_contrasts

            def inherited_contrasts(*args):
                inherited = original_contrasts(*args)
                inherited["missing_contrasts"] = [{"reason": "Older ice-only comparison unavailable"}] * 17
                return inherited

            with mock.patch.object(PRODUCTION, "controlled_contrasts", side_effect=inherited_contrasts):
                result = ANALYSIS.analyze_campaign(run, make_plot=True)
            self.assertEqual(result["status"], "complete_campaign")
            self.assertFalse(result["discrete_geometry_column_validated"])
            self.assertFalse(result["posterior_intervals_reported"])
            self.assertEqual(result["material_input_check"]["coverage_policy"], "exact_tables_native_extrapolation")
            generic = PRODUCTION.read_json(run / "results/summary.json")
            self.assertEqual(generic["campaign_context"]["material_input_check"], result["material_input_check"])
            self.assertEqual(sum(len(profiles) for profiles in generic["models"][0]["profiles"].values()), 9)
            self.assertNotIn("H2O30K laboratory constants are not used", " ".join(generic["limitations"]))
            review = (run / "results/REVIEW.md").read_text()
            self.assertIn("Missing/invalid matched contrasts: 0.", review)
            self.assertNotIn("Missing/invalid matched contrasts: 17.", review)
            with (run / "results/ranking.csv").open() as stream:
                scores = json.loads(next(csv.DictReader(stream))["scenario_scores"])
            self.assertEqual(sum(len(values) for values in scores.values()), 9)
            for filename in ("campaign_summary.json", "CAMPAIGN_REVIEW.md", "matched_band_responses.csv",
                    "matched_feature_shapes.csv", "missing_pairs.csv", "matched_response_map.png",
                    "matched_response_map.pdf", "shape_tradeoffs.png", "shape_tradeoffs.pdf"):
                self.assertTrue((run / "results" / filename).is_file(), filename)
            path = run / "models/m0/measurements.json"
            payload = PRODUCTION.read_json(path)
            payload["measurements"][0]["diagnostics"]["quality_checks"]["finite_positive_plane0_flux"] = False
            FIXTURES.write_json(path, payload)
            failed = ANALYSIS.analyze_campaign(run, make_plot=False)
            self.assertEqual(failed["models_screened"], 59)
            self.assertEqual(failed["status"], "partial_campaign")
            self.assertTrue(failed["missing_pairs"])


if __name__ == "__main__":
    unittest.main()
