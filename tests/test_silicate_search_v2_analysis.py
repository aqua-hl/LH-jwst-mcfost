"""Predeclared matched-extinction contrasts and descriptive success boundaries."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ANALYSIS = load("silicate_v2_analysis", "scripts/analyze_silicate_search_v2.py")
FIXTURES = load("silicate_v2_receipt_fixtures", "tests/test_silicate_size_analysis.py")
OBSERVABLES = load("silicate_v2_observables", "scripts/production_observables.py")
PRODUCTION = ANALYSIS.production


def catalogue_fixture():
    coordinates = [("draine", carbon, ice, -1.5, inc)
                   for carbon in (0., .1, .2, .3) for ice in (.04, .08, .12) for inc in (60., 70.)]
    coordinates += [(material, .2, ice, exponent, 60.)
                    for material in ("draine", "pyroxene_mg50", "olivine_mg50")
                    for exponent in (-1.5, -1.75) for ice in (.08, .12)
                    if (material, .2, ice, exponent, 60.) not in coordinates]
    result = []
    for material, carbon, ice, exponent, inc in coordinates:
        index = len(result)
        # Synthetic Stage 0 matching; no physical opacity estimate is implied.
        mass = .000225 / (1 + 5 * carbon) * (.2393567089536 if exponent == -1.75 else 1)
        result.append({"index": index, "role": "composition" if index < 24 else "discriminator",
            "material_id": material, "carbon_mass_fraction": carbon, "ice_mass_fraction": ice,
            "density_exponent": exponent, "inclination_deg": inc,
            "parameters": {"envelope_silicate_file": material + ".dat", "envelope_density_exponent": exponent,
                "envelope_carbon_mass_fraction": carbon, "envelope_ice_mass_fraction": ice,
                "envelope_dust_mass_msun": mass, "inclination_deg": inc, "distance_pc": 140.},
            "opacity_matching": {"S": 1 / (1 + 5 * carbon), "ratio_18_9p7": .4}})
    return result


def model_fixture(catalogue, contract):
    return [{"model_index": r["index"], "model_id": f"m{r['index']}", "parameters": copy.deepcopy(r["parameters"]),
             "status": "screened", "band_flux_jy": [b["observed_flux_jy"] for b in contract["bands"]],
             "quadrature_unresolved": False, "unresolved_quadrature_bands": [],
             "baseline_score": 999. if r["carbon_mass_fraction"] > 0 else 1.} for r in catalogue]


def compact_campaign(run, catalogue, contract):
    manifest = FIXTURES.compact_fixture(run)
    template = PRODUCTION.read_json(run / "models/m0/measurements.json")
    FIXTURES.write_json(run / "production_experiment.json", {"experiment_id": "silicate_search_v2",
        "catalogue": catalogue, "predeclared_contrasts": ANALYSIS.predeclared_contrasts(catalogue),
        "stage0_opacity": {"species": {name: {"kappa_ext_2p2": float(index + 1)}
            for index, name in enumerate(("draine", "pyroxene", "olivine", "carbon", "ice"))}},
        "opacity_matching": {"method": "mass_weighted_species_DHS", "S_ref": 1.},
        "dust_prescription": {"H2O30K": "unchanged", "carbon": "separate_DHS"}})
    FIXTURES.write_json(run / "inputs/production_observation_contract.json", contract)
    manifest["models"] = [{"id": f"m{r['index']}", "index": r["index"], "parameters": r["parameters"]} for r in catalogue]
    manifest["anchors"] = contract["probes"]
    manifest["input_hashes"] = {relative: PRODUCTION.digest(run / relative) for relative in manifest["input_hashes"]}
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


class CompositionAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = OBSERVABLES.build_contract(ROOT)

    def setUp(self):
        self.catalogue = catalogue_fixture()
        self.models = model_fixture(self.catalogue, self.contract)
        self.band_index = {b["id"]: i for i, b in enumerate(self.contract["bands"])}

    def change(self, model_index, band, dex):
        self.models[model_index]["band_flux_jy"][self.band_index[band]] *= 10 ** dex

    def test_all_48_declared_pairs_are_reported_without_score_selection(self):
        result = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(result["status"], "complete_campaign")
        self.assertEqual(result["pairs_complete"], 48)
        self.assertEqual(result["pairs_complete_by_kind"], {"carbon": 18, "ice": 16, "material": 8, "profile": 6})
        self.assertEqual(len(result["band_responses"]), 48 * 19)
        self.assertEqual(len(result["feature_responses"]), 48 * 5)
        carbon = next(r for r in result["pairs"] if r["kind"] == "carbon")
        self.assertEqual(carbon["baseline_profiled_scores_diagnostic_only"], [1., 999.])
        self.assertEqual(carbon["reference_ice_mass_fraction"], carbon["changed_ice_mass_fraction"])
        self.assertNotEqual(carbon["reference_dust_mass_msun"], carbon["changed_dust_mass_msun"])
        self.assertTrue(all(r["inclination_deg"] == 60 for r in result["pairs"] if r["kind"] in {"material", "profile"}))

    def test_success_boundaries_use_raw_ice_and_do_not_require_h04(self):
        self.change(0, "s02", .10)
        for band in ("h02", "h03", "s05", "s06"):
            self.change(0, band, -.15)
        self.change(0, "h04", 2.)
        rows = ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)
        self.assertTrue(rows[0]["joint_success"])
        self.assertAlmostEqual(rows[0]["band_residuals_dex"]["h04"], 2.)
        self.assertAlmostEqual(rows[0]["feature_log_residuals_dex"]["ice_h04"], 2.)
        self.change(0, "s02", .0001)
        self.assertFalse(ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)[0]["joint_success"])

    def test_joint_success_is_individual_and_nuisance_score_cannot_rescue_failure(self):
        for model in self.models:
            model["baseline_score"] = 0.
        self.change(0, "h02", .16)
        self.change(1, "s05", .16)
        for band in [b["id"] for b in self.contract["bands"] if b["block"] == "nir_continuum"]:
            self.change(2, band, .2001)
        rows = ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)
        self.assertFalse(any(r["joint_success"] for r in rows[:3]))
        self.assertFalse(rows[0]["success_conditions"]["ice_2p95_3p11_raw"])
        self.assertFalse(rows[1]["success_conditions"]["miri_18_20"])
        self.assertFalse(rows[2]["success_conditions"]["nir_continuum"])

    def test_raw_ice_success_and_normalized_ice_failure_are_explicit(self):
        # Ice cores are exactly observed but their modeled continuum is dim.
        # Only two of nine continuum bands change, so continuum RMS still passes.
        self.change(0, "c07", -.30)
        self.change(0, "c08", -.30)
        row = ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)[0]
        self.assertTrue(row["joint_success"])
        self.assertFalse(row["normalized_ice_within_0p15_dex_diagnostic"])
        self.assertAlmostEqual(row["feature_log_residuals_dex"]["ice_h02"], .3)
        self.assertIn("RAW", ANALYSIS.SUCCESS_RULE["ice_interpretation"])

    def test_missing_receipts_and_shoulder_quadrature_are_visible_without_selection(self):
        self.models[0].update(status="excluded", reason="Nonpositive signed image flux")
        self.models[6].update(quadrature_unresolved=True, unresolved_quadrature_bands=["s01"])
        result = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(result["status"], "partial_campaign")
        self.assertEqual(result["pairs_complete"] + len(result["missing_pairs"]), 48)
        self.assertTrue(all(0 not in (r["reference_index"], r["changed_index"]) for r in result["pairs"]))
        affected = [r for r in result["feature_responses"] if r["changed_index"] == 6 and r["feature_id"] == "silicate_9p7"]
        # The carbon pair involving model 0 is missing; other pairs use 6 as a reference.
        affected += [r for r in result["feature_responses"] if r["reference_index"] == 6 and r["feature_id"] == "silicate_9p7"]
        self.assertTrue(affected)
        self.assertTrue(all(r["reference_quadrature_unresolved"] or r["changed_quadrature_unresolved"] for r in affected))
        rows = ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)
        self.assertIsNone(rows[0]["joint_success"])
        self.assertTrue(rows[6]["joint_success"])
        self.assertTrue(rows[6]["success_provisional_due_to_quadrature"])

    def test_changed_coordinates_and_receipt_parameters_rejected(self):
        bad = copy.deepcopy(self.catalogue)
        bad[-1]["density_exponent"] = -1.6
        with self.assertRaisesRegex(ValueError, "Catalogue coordinates differ"):
            ANALYSIS.predeclared_contrasts(bad)
        bad = copy.deepcopy(self.catalogue)
        bad[0]["parameters"]["envelope_ice_mass_fraction"] = .05
        with self.assertRaisesRegex(ValueError, "fraction differs"):
            ANALYSIS.predeclared_contrasts(bad)
        self.models[0]["parameters"]["distance_pc"] = 150.
        with self.assertRaisesRegex(ValueError, "Analyzed parameters differ"):
            ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)

    def test_compact_receipts_nine_scenarios_all_outputs_and_mandatory_quality(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_campaign(run, self.catalogue, self.contract)
            result = ANALYSIS.analyze_campaign(run, make_plot=True)
            self.assertEqual(result["status"], "complete_campaign")
            self.assertEqual(len(result["model_diagnostics"]), 34)
            self.assertEqual(result["pairs_complete"], 48)
            self.assertFalse(result["abundances_uniquely_identified"])
            self.assertEqual(len(result["campaign_context"]["stage0_opacity"]["species"]), 5)
            self.assertIn("S", result["catalogue_with_opacity_matching"][0]["opacity_matching"])
            generic = PRODUCTION.read_json(run / "results/summary.json")
            self.assertEqual(sum(len(profiles) for profiles in generic["models"][0]["profiles"].values()), 9)
            self.assertEqual(generic["controlled_contrasts"]["ice_responses"], [])
            with (run / "results/ranking.csv").open() as stream:
                scores = json.loads(next(csv.DictReader(stream))["scenario_scores"])
            self.assertEqual(sum(len(v) for v in scores.values()), 9)
            review = (run / "results/REVIEW.md").read_text()
            self.assertIn("RAW flux", review)
            self.assertNotIn("bare/ice ladders", review)
            for filename in ("campaign_summary.json", "CAMPAIGN_REVIEW.md", "matched_band_responses.csv",
                    "matched_feature_shapes.csv", "missing_pairs.csv", "model_success.csv", "opacity_matched_catalogue.csv",
                    "matched_response_map.png", "matched_response_map.pdf", "all_model_criteria.png", "all_model_criteria.pdf"):
                self.assertTrue((run / "results" / filename).is_file(), filename)
            path = run / "models/m0/measurements.json"
            payload = PRODUCTION.read_json(path)
            payload["measurements"][0]["diagnostics"]["quality_checks"]["finite_positive_plane0_flux"] = False
            FIXTURES.write_json(path, payload)
            failed = ANALYSIS.analyze_campaign(run, make_plot=False)
            self.assertEqual(failed["status"], "partial_campaign")
            self.assertEqual(failed["models_screened"], 33)
            self.assertIsNone(failed["model_diagnostics"][0]["joint_success"])


if __name__ == "__main__":
    unittest.main()
