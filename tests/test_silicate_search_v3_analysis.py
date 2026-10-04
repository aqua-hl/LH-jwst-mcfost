"""V3 comparisons, replica separation, bridge gating and descriptive outcomes."""
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
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ANALYSIS = load("v3_analysis", "scripts/analyze_silicate_search_v3.py")
FIXTURES = load("v3_receipts", "tests/test_silicate_size_analysis.py")
OBSERVABLES = load("v3_observables", "scripts/production_observables.py")
PRODUCTION = ANALYSIS.production


def catalogue_fixture():
    rows = []
    for config, (exponent, column) in ANALYSIS.CONFIGURATIONS.items():
        for material in ANALYSIS.MATERIALS:
            for carbon in (.15, .2, .25):
                ice = .092 if carbon == .2 else .08 if carbon == .15 else .105
                mass = .000225*column/(1+5*carbon)
                rows.append(dict(index=len(rows), role="primary", configuration_id=config, material_id=material,
                    carbon_mass_fraction=carbon, ice_mass_fraction=ice, density_exponent=exponent,
                    column_factor=column, inclination_deg=70., random_seed=43001, replica_of=None,
                    numerics={"random_seed": 43001},
                    parameters=dict(inclination_deg=70., distance_pc=140., envelope_density_exponent=exponent,
                        envelope_dust_mass_msun=mass, envelope_carbon_mass_fraction=carbon,
                        envelope_ice_mass_fraction=ice,
                        envelope_silicate_file="Draine_Si_sUV.dat" if material == "draine" else "Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat")))
    for parent in (7, 10):
        row = copy.deepcopy(rows[parent])
        row.update(index=len(rows), role="seed_replica", random_seed=43002, replica_of=parent)
        row["numerics"]["random_seed"] = 43002
        rows.append(row)
    return rows


def models_fixture(catalogue, contract):
    return [dict(model_index=r["index"], model_id=f"m{r['index']}", status="screened", parameters=copy.deepcopy(r["parameters"]),
        baseline_score=1000-r["index"], band_flux_jy=[b["observed_flux_jy"] for b in contract["bands"]],
        quadrature_unresolved=False, unresolved_quadrature_bands=[]) for r in catalogue]


def bridge_fixture(contract):
    result = []
    for index, ice in ((15, .08), (17, .12)):
        result.append(dict(model_index=index, status="screened", model_id=f"v2_m{index}",
            band_flux_jy=[b["observed_flux_jy"] for b in contract["bands"]],
            parameters=dict(inclination_deg=70., envelope_density_exponent=-1.5, envelope_carbon_mass_fraction=.2,
                            envelope_ice_mass_fraction=ice, envelope_silicate_file="Draine_Si_sUV.dat")))
    return dict(bridge_models=result)


def compact_campaign(run, catalogue, contract):
    manifest = FIXTURES.compact_fixture(run)
    template = PRODUCTION.read_json(run/"models/m0/measurements.json")
    FIXTURES.write_json(run/"inputs/production_observation_contract.json", contract)
    reference = bridge_fixture(contract)
    reference["observation_contract_sha256"] = PRODUCTION.digest(run/"inputs/production_observation_contract.json")
    FIXTURES.write_json(run/"inputs/v2_reference.json", reference)
    FIXTURES.write_json(run/"production_experiment.json", dict(experiment_id="silicate_search_v3", catalogue=catalogue,
        predeclared_contrasts=ANALYSIS.predeclared_contrasts(catalogue),
        design_provenance={"ice_crossings_sha256": "a"*64},
        ice_column_targets_msun={"draine": 9.14e-6, "pyroxene_mg50": 1.11e-5}))
    manifest["models"] = [dict(id=f"m{r['index']}", index=r["index"], parameters=r["parameters"], numerics=r["numerics"])
                          for r in catalogue]
    manifest["anchors"] = contract["probes"]
    manifest["input_hashes"]["inputs/v2_reference.json"] = PRODUCTION.digest(run/"inputs/v2_reference.json")
    manifest["input_hashes"] = {p: PRODUCTION.digest(run/p) for p in manifest["input_hashes"]}
    FIXTURES.write_json(run/"manifest.json", manifest)
    fingerprint = {**template["fingerprint"], "manifest_sha256": PRODUCTION.digest(run/"manifest.json")}
    FIXTURES.write_json(run/"runtime_binding.json", {"fingerprint": fingerprint})
    for model in manifest["models"]:
        measurements = []
        for probe in contract["probes"]:
            receipt = copy.deepcopy(template["measurements"][0])
            receipt.update(anchor_id=probe["id"], wavelength_um=probe["wavelength_um"], instrument=probe["instrument"], flux_jy=1.)
            measurements.append(receipt)
        FIXTURES.write_json(run/"models"/model["id"]/"measurements.json", dict(model_id=model["id"], model_index=model["index"],
            parameters=model["parameters"], numerics=model["numerics"], complete=True, fingerprint=fingerprint, measurements=measurements))


class V3AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = OBSERVABLES.build_contract(ROOT)

    def setUp(self):
        self.catalogue = catalogue_fixture()
        self.models = models_fixture(self.catalogue, self.contract)
        self.band_index = {b["id"]: i for i, b in enumerate(self.contract["bands"])}
        self.reference = bridge_fixture(self.contract)

    def change(self, index, band, dex):
        self.models[index]["band_flux_jy"][self.band_index[band]] *= 10**dex

    def diagnostics(self):
        result = ANALYSIS.model_diagnostics(self.catalogue, self.models, self.contract)
        scatter = ANALYSIS.seed_scatter(result)
        ANALYSIS.apply_marginality(result, scatter)
        return result, scatter

    def test_54_declarations_preserve_primary_and_replica_identity(self):
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(pairs["pairs_complete"], 54)
        self.assertEqual(pairs["pairs_complete_by_kind"], ANALYSIS.COUNTS)
        self.assertEqual(len(pairs["band_responses"]), 54*19)
        self.assertEqual(len(pairs["feature_responses"]), 54*5)
        self.assertEqual([(p["reference_index"], p["changed_index"]) for p in pairs["pairs"] if p["kind"] == "seed"], [(7, 24), (10, 25)])
        self.assertTrue(all(p["changed_index"] < 24 for p in pairs["pairs"] if p["kind"] != "seed"))
        bad = copy.deepcopy(self.catalogue)
        bad[24]["parameters"]["envelope_dust_mass_msun"] *= 2
        with self.assertRaisesRegex(ValueError, "Seed pair physical"):
            ANALYSIS.predeclared_contrasts(bad)

    def test_raw_criteria_normalized_ice_and_h04_are_separate(self):
        self.change(0, "s02", .10)
        for band in ("h02", "h03", "s05", "s06"):
            self.change(0, band, -.15)
        self.change(0, "h04", 2.)
        self.change(1, "c07", -.3)
        self.change(1, "c08", -.3)
        rows, _ = self.diagnostics()
        self.assertTrue(rows[0]["joint_success"])
        self.assertAlmostEqual(rows[0]["band_residuals_dex"]["h04"], 2.)
        self.assertTrue(rows[1]["joint_success"])
        self.assertFalse(rows[1]["normalized_ice_within_0p15_dex_diagnostic"])
        self.change(0, "s02", .0001)
        self.assertFalse(self.diagnostics()[0][0]["joint_success"])

    def test_seed_scatter_is_absolute_pair_difference_and_global_marginality(self):
        self.change(24, "s05", .02)
        self.change(25, "s05", .04)
        self.change(0, "s05", .13)
        rows, scatter = self.diagnostics()
        self.assertEqual(scatter["status"], "complete")
        self.assertAlmostEqual(scatter["pairs"][0]["criterion_absolute_paired_differences_dex"]["s05_raw_dex"], .02)
        self.assertAlmostEqual(scatter["conservative_criterion_scatter_dex"]["s05_raw_dex"], .04)
        self.assertTrue(scatter["small_responses_unresolved_warning"])
        self.assertTrue(rows[0]["joint_success"])
        self.assertTrue(rows[0]["marginal_success"])
        self.assertIn("s05_raw_dex", rows[0]["criterion_margins_within_seed_scatter"])
        self.assertIn("transfer assumption", scatter["transfer_assumption"])

    def test_bridge_range_is_widened_by_measured_seed_difference(self):
        self.change(24, "s05", .01)
        self.change(25, "s05", .02)
        self.change(1, "s05", .019)
        rows, scatter = self.diagnostics()
        bridge = ANALYSIS.bridge_check(rows, self.reference, self.contract, scatter)
        self.assertTrue(bridge["passed"])
        self.change(1, "s05", .002)
        rows, scatter = self.diagnostics()
        bridge = ANALYSIS.bridge_check(rows, self.reference, self.contract, scatter)
        self.assertFalse(bridge["passed"])
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        predictions = ANALYSIS.declared_predictions(rows)
        outcome = ANALYSIS.campaign_outcomes(rows, pairs, bridge, scatter, predictions)
        self.assertFalse(outcome["interpretive_gate_passed"])
        self.assertIsNone(outcome["success_statement_supported"])
        self.assertIsNone(outcome["stop_global_profile_search_signal"])
        self.assertTrue(rows[1]["joint_success"])

    def test_missing_seed_does_not_become_zero_scatter_or_discard_primary(self):
        self.models[25].update(status="excluded", reason="Nonpositive signed flux")
        self.models[7].update(quadrature_unresolved=True, unresolved_quadrature_bands=["s01"])
        rows, scatter = self.diagnostics()
        self.assertEqual(scatter["status"], "incomplete")
        self.assertIsNone(scatter["conservative_criterion_scatter_dex"])
        self.assertIsNone(rows[7]["marginal_success"])
        self.assertTrue(rows[7]["joint_success"])
        self.assertTrue(rows[7]["success_provisional_due_to_quadrature"])
        self.assertEqual(ANALYSIS.bridge_check(rows, self.reference, self.contract, scatter)["status"], "unavailable")
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(pairs["pairs_complete"], 53)
        self.assertEqual(pairs["missing_pairs"][0]["kind"], "seed")

    def test_ice_column_and_nir_rematch_failures_are_reported_per_configuration(self):
        self.change(6, "h02", .16)
        for band in (f"c{i:02d}" for i in range(1, 10)):
            self.change(6, band, .21)
            self.change(12, band, .21)
        rows, _ = self.diagnostics()
        prediction = ANALYSIS.declared_predictions(rows)
        self.assertFalse(prediction["ice_prediction_all_24_pass"])
        self.assertEqual(prediction["ice_column_prediction_by_configuration"][1]["raw_ice_fail_indices"], [6])
        self.assertFalse(prediction["nir_rematch_all_compositions_pass"])
        first = prediction["nir_rematch_by_composition"][0]
        self.assertEqual(first["configuration_pass"], {"shallow": False, "flat_high": False, "flat_low": True})

    def test_profile_median_thresholds_and_contingency_never_launch(self):
        # Every anchor starts 0.2 dex too bright in long MIRI. The shallow
        # change hits exactly -0.10, but breaks NIR; flat high is -0.075.
        for i, row in enumerate(self.catalogue):
            offset = {"anchor": .20, "shallow": .10, "flat_high": .125, "flat_low": .16}[row["configuration_id"]]
            for band in ("s05", "s06"):
                self.change(i, band, offset)
            if row["configuration_id"] != "anchor":
                for band in (f"c{j:02d}" for j in range(1, 10)):
                    self.change(i, band, .3)
        rows, scatter = self.diagnostics()
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        outcome = ANALYSIS.campaign_outcomes(rows, pairs, {"passed": True}, scatter, ANALYSIS.declared_predictions(rows))
        self.assertEqual(outcome["inner_flattened_contingency_discussion_configurations"], ["shallow"])
        self.assertEqual(outcome["partial_response_discussion_configurations"], ["flat_high"])
        self.assertFalse(outcome["stop_global_profile_search_signal"])
        self.assertFalse(outcome["actions_launched"])
        self.assertAlmostEqual(outcome["profile_response_by_configuration"][0]["median_long_miri_change_dex"], -.10)
        for model in self.models:
            model["band_flux_jy"] = [b["observed_flux_jy"]*10**.20 for b in self.contract["bands"]]
        rows, scatter = self.diagnostics()
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        outcome = ANALYSIS.campaign_outcomes(rows, pairs, {"passed": True}, scatter, ANALYSIS.declared_predictions(rows))
        self.assertTrue(outcome["stop_global_profile_search_signal"])

    def test_replica_only_boundary_success_is_marginal_and_prevents_null_stop(self):
        for index in range(26):
            offset = .151 if index == 7 else .149 if index == 24 else .16
            for band in ("s05", "s06"):
                self.change(index, band, offset)
        rows, scatter = self.diagnostics()
        self.assertFalse(rows[7]["joint_success"])
        self.assertTrue(rows[24]["joint_success"])
        self.assertTrue(rows[24]["marginal_success"])
        pairs = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        outcome = ANALYSIS.campaign_outcomes(rows, pairs, {"passed": True}, scatter, ANALYSIS.declared_predictions(rows))
        self.assertEqual(outcome["raw_primary_joint_success_indices"], [])
        self.assertEqual(outcome["raw_replica_joint_success_indices"], [24])
        self.assertTrue(outcome["success_statement_supported"])
        self.assertTrue(outcome["success_statement_seed_dependent"])
        self.assertFalse(outcome["stop_global_profile_search_signal"])
        self.assertEqual(outcome["inner_flattened_contingency_discussion_configurations"], [])
        replica = outcome["replica_only_success_configurations"][0]
        self.assertEqual((replica["primary_index"], replica["replica_index"]), (7, 24))
        self.assertTrue(replica["seed_dependent"])
        self.assertTrue(replica["marginal"])
        # The six primary profile changes are 0, -.009, 0, 0, 0, 0;
        # the passing replica does not enter or replace the primary response.
        shallow = outcome["profile_response_by_configuration"][0]
        self.assertAlmostEqual(shallow["median_long_miri_change_dex"], 0.)
        self.assertAlmostEqual(shallow["individual_pair_mean_changes_dex"][1], -.009)

    def test_receipts_nine_scenarios_artifacts_seed_rejection_and_bridge_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_campaign(run, self.catalogue, self.contract)
            result = ANALYSIS.analyze_campaign(run, make_plot=True)
            self.assertEqual(result["models_screened"], 26)
            self.assertEqual(result["primary_models_screened"], 24)
            self.assertEqual(result["seed_replicas_screened"], 2)
            self.assertEqual(result["pairs_complete"], 54)
            self.assertEqual(result["bridge_to_v2"]["status"], "failed")
            self.assertFalse(result["outcomes"]["interpretive_gate_passed"])
            self.assertEqual(result["campaign_context"]["design_provenance"]["ice_crossings_sha256"], "a"*64)
            with (run/"results/ranking.csv").open() as stream:
                scores = json.loads(next(csv.DictReader(stream))["scenario_scores"])
            self.assertEqual(sum(len(v) for v in scores.values()), 9)
            for filename in ("campaign_summary.json", "model_success.csv", "seed_scatter.csv", "bridge_to_v2.csv",
                    "ice_column_prediction.csv", "nir_rematch.csv", "profile_outcomes.csv", "matched_band_responses.csv",
                    "matched_feature_shapes.csv", "all_model_criteria.png", "all_model_criteria.pdf"):
                self.assertTrue((run/"results"/filename).is_file())
            path = run/"models/m24/measurements.json"
            payload = PRODUCTION.read_json(path)
            payload["numerics"]["random_seed"] = 43001
            FIXTURES.write_json(path, payload)
            result = ANALYSIS.analyze_campaign(run, make_plot=False)
            self.assertEqual(result["models_screened"], 25)
            self.assertEqual(result["primary_models_screened"], 24)
            self.assertEqual(result["bridge_to_v2"]["status"], "unavailable")
            self.assertEqual(result["pairs_complete"], 53)
            reference = PRODUCTION.read_json(run/"inputs/v2_reference.json")
            reference["observation_contract_sha256"] = "b"*64
            FIXTURES.write_json(run/"inputs/v2_reference.json", reference)
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                ANALYSIS.analyze_campaign(run, make_plot=False)


if __name__ == "__main__":
    unittest.main()
