"""Production statistics, controlled-response degeneracy and compact validation."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("production_analysis", ROOT / "scripts/analyze_extinction_ice_production.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def contract_fixture():
    probes, bands = [], []
    for index, block in enumerate(ANALYSIS.BLOCKS):
        wave = 1. + index
        instrument = "NIRSpec" if index < 2 else "MIRI"
        segment = "shared" if index < 2 else "CH2_medium"
        ids = []
        for node, wavelength in enumerate(np.linspace(wave, wave + .2, 5)):
            name = f"b{index}_p{node}"
            ids.append(name)
            probes.append({"id": name, "wavelength_um": float(wavelength), "score": False, "instrument": instrument,
                           "region": block, "image_npix": 6001, "image_size_au": 6000., "psf_fwhm_arcsec": .1})
        bands.append({"id": f"b{index}", "block": block, "instrument": instrument, "segment": segment,
                      "gain_group": instrument + "/" + segment, "lower_um": wave, "upper_um": wave + .2,
                      "observed_flux_jy": 1., "formal_sigma_jy": .05, "structure_fraction": .01,
                      "structure_sigma_jy": .01, "probe_ids": ids, "weights": [1/12, 4/12, 2/12, 4/12, 1/12],
                      "coarse_probe_ids": ids[::2], "coarse_weights": [1/6, 4/6, 1/6]})
    scenario = {"band_order": [b["id"] for b in bands], "baseline": True,
                "relative_gain_fraction_by_segment": {"NIRSpec/shared": .03, "MIRI/CH2_medium": .02}}
    base = np.eye(3) * (.05**2 + .01**2)
    design, sigmas, _ = ANALYSIS.gain_design(bands, scenario)
    gain = (design * sigmas) @ (design * sigmas).T
    scenario.update(formal_plus_structure_covariance_jy2=base.tolist(), relative_gain_covariance_jy2=gain.tolist(),
                    conditional_covariance_jy2=(base + gain).tolist())
    return {"schema_version": 1, "probes": probes, "bands": bands,
            "block_weights": ANALYSIS.WEIGHTING["baseline"], "covariance_scenarios": {"baseline": scenario},
            "common_scale_nuisance": {"kind": "gaussian_multiplicative", "mean": 1., "sigma_fraction": .03}}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False) + "\n")


def run_fixture(run):
    contract = contract_fixture()
    model = {"id": "m0000_fixture", "index": 0, "parameters": {"distance_pc": 140., "inclination_deg": 50.,
                                                               "envelope_dust_mass_msun": .00015}}
    experiment = {"catalogue": [{"index": 0, "role": "parent_control", "parameters": model["parameters"],
                                 "parent_index": 0, "ice_mass_fraction": .02}]}
    write_json(run / "production_experiment.json", experiment)
    write_json(run / "inputs/production_observation_contract.json", contract)
    quality_path = run / "code/src/mcfost_grid/quality.py"
    quality_path.parent.mkdir(parents=True)
    shutil.copy2(ROOT / "src/mcfost_grid/quality.py", quality_path)
    manifest = {"models": [model], "anchors": contract["probes"],
                "configuration": {"numerical_only": True, "numerics": {"image_npix": 6001, "image_size_au": 6000.}},
                "measurement": {"quality_policy": "aperture_v2", "aperture_radius_arcsec": 1., "target_distance_pc": 147.},
                "input_hashes": {str(path.relative_to(run)): ANALYSIS.digest(path)
                                 for path in run.rglob("*") if path.is_file()}}
    write_json(run / "manifest.json", manifest)
    fingerprint = {"manifest_sha256": ANALYSIS.digest(run / "manifest.json"), "executable_sha256": "a" * 64,
                   "utilities_content_sha256": "b" * 64, "backend": "image_method2"}
    write_json(run / "runtime_binding.json", {"fingerprint": fingerprint})
    checks = {"aperture_plus_6sigma_support": True, "finite_kernel_aperture_support": True,
              "finite_positive_plane0_flux": True, "finite_positive_signed_aperture": True,
              "edge_positive_fraction_lt_1e-3": False, "convolution_flux_loss_lt_1e-4": True}
    blocking = {k: v for k, v in checks.items() if k not in {"edge_positive_fraction_lt_1e-3", "convolution_flux_loss_lt_1e-4"}}
    pixel = 6000 / 6001 / 140
    target_pixel = pixel * 140 / 147
    measurements = []
    for probe in contract["probes"]:
        measurements.append({"anchor_id": probe["id"], "wavelength_um": probe["wavelength_um"],
            "instrument": probe["instrument"], "flux_jy": 1.01, "valid": True, "quality_pass": True,
            "quality_policy": "aperture_v2", "quality_warnings": ["edge_positive_fraction_lt_1e-3"],
            "backend": "image_method2", "image_sha256": "c" * 64,
            "canonical_estimator": "signed_direct_method2_polarized_total_i_plane0",
            "model_distance_pc": 140., "target_distance_pc": 147., "aperture_radius_arcsec": 1.,
            "diagnostics": {"quality_checks": checks, "blocking_quality_checks": blocking,
                "distance_flux_factor": (140 / 147)**2, "zero_clipping_used_for_scoring": False,
                "nx": 6001, "ny": 6001, "centre_x_zero_based": 3000., "centre_y_zero_based": 3000.,
                "aperture_radius_pixels": 1 / target_pixel, "psf_sigma_pixels": 6., "psf_truncation_sigma": 6.,
                "model_pixel_arcsec": pixel, "target_pixel_arcsec": target_pixel}})
    payload = {"model_id": model["id"], "model_index": 0, "parameters": model["parameters"],
               "fingerprint": fingerprint, "complete": True, "measurements": measurements}
    write_json(run / "models" / model["id"] / "measurements.json", payload)
    return manifest, payload


class ProductionStatisticsTests(unittest.TestCase):
    def test_block_normalization_is_independent_of_band_multiplicity(self):
        blocks = ["nir_continuum", "nir_ice", "miri"]
        residual = np.array([2., 3., 4.])
        precision = ANALYSIS.block_precision(np.eye(3), blocks, ANALYSIS.WEIGHTING["baseline"])
        score = residual @ precision @ residual
        duplicated = ["nir_continuum"] * 7 + ["nir_ice", "miri"]
        repeated = np.r_[np.repeat(2., 7), 3., 4.]
        other = ANALYSIS.block_precision(np.eye(9), duplicated, ANALYSIS.WEIGHTING["baseline"])
        self.assertAlmostEqual(score, .45 * 4 + .3 * 9 + .25 * 16)
        self.assertAlmostEqual(score, repeated @ other @ repeated)

    def test_shared_segment_gain_profile_matches_covariance_elimination(self):
        contract = contract_fixture()
        scenario = contract["covariance_scenarios"]["baseline"]
        model = np.array([1.12, 1.10, .92])
        result = ANALYSIS.profile_screen(model, contract, scenario)
        base = np.array(scenario["formal_plus_structure_covariance_jy2"])
        precision = ANALYSIS.block_precision(base, list(ANALYSIS.BLOCKS), ANALYSIS.WEIGHTING["baseline"])
        z, sigmas, _ = ANALYSIS.gain_design(contract["bands"], scenario)
        marginal_precision = np.linalg.inv(np.linalg.inv(precision) + (z * sigmas) @ (z * sigmas).T)
        residual = model - 1
        a = -(model @ marginal_precision @ residual) / (model @ marginal_precision @ model + 1/.03**2)
        adjusted = residual + a * model
        expected = adjusted @ marginal_precision @ adjusted + (a/.03)**2
        self.assertAlmostEqual(result["score"], expected)
        self.assertAlmostEqual(result["common_scale"], 1 + a)
        self.assertEqual(len(result["segment_gain_fraction"]), 2)
        self.assertGreater(result["segment_gain_fraction"]["NIRSpec/shared"], 0)
        self.assertAlmostEqual(sum(result["block_weighted_data_contribution"].values()), result["weighted_data_term"])
        for block, weight in ANALYSIS.WEIGHTING["baseline"].items():
            self.assertAlmostEqual(result["block_weighted_data_contribution"][block],
                                   weight * result["block_standardized_residual_rms"][block]**2)

    def test_nir_only_score_is_independent_of_miri_model(self):
        contract = contract_fixture()
        scenario = contract["covariance_scenarios"]["baseline"]
        first = ANALYSIS.profile_screen([1.1, .9, 1.], contract, scenario, ANALYSIS.WEIGHTING["nir_only"])
        second = ANALYSIS.profile_screen([1.1, .9, 10000.], contract, scenario, ANALYSIS.WEIGHTING["nir_only"])
        self.assertAlmostEqual(first["score"], second["score"])
        self.assertEqual(second["segment_gain_fraction"]["MIRI/CH2_medium"], 0.)

    def test_five_node_quadrature_retains_flagged_bands(self):
        contract = contract_fixture()
        values = {probe["id"]: 1. for probe in contract["probes"]}
        values[contract["bands"][0]["probe_ids"][1]] = 2.
        rows = ANALYSIS.integrate_bands(values, contract)
        self.assertEqual(len(rows), 3)
        self.assertTrue(rows[0]["quadrature_unresolved"])
        self.assertAlmostEqual(rows[0]["model_flux_jy"], 4/3)
        self.assertAlmostEqual(rows[0]["three_node_flux_jy"], 1.)

    def test_known_ice_mass_degeneracy_and_common_scale_disappear(self):
        reference = np.ones(4)
        ice = np.array([1., -1., 0., 0.])
        result = ANALYSIS.response_svd(ice, 2 * ice + .5 * reference, reference, np.eye(4))
        self.assertTrue(result["near_collinear_diagnostic"])
        self.assertLess(result["small_to_large_singular_value_ratio"], 1e-12)
        self.assertAlmostEqual(result["response_cosine"], 1.)
        self.assertIsNone(result["number_of_physical_constraints"])
        self.assertTrue(np.allclose(ANALYSIS.project_common_scale(5 * reference, reference, np.eye(4)), 0))

    def test_controlled_ladders_use_all_models_not_best_fit(self):
        contract = contract_fixture()
        catalogue, results = [], []
        direction = np.array([1., -1., 0.])
        for mass_index, mass in enumerate((.00015, .000225)):
            for fraction in (0., .02, .04, .08):
                index = len(catalogue)
                parameters = {"inclination_deg": 50., "envelope_dust_mass_msun": mass,
                              "envelope_ice_volume_fraction": fraction}
                catalogue.append({"index": index, "parameters": parameters, "role": "ice_ladder",
                                  "ice_mass_fraction": fraction, "contrast_group_id": f"group_{mass_index}"})
                flux = 100 * np.ones(3) + (fraction/.04 + mass_index) * direction
                results.append({"model_index": index, "status": "screened", "band_flux_jy": flux.tolist(),
                                "baseline_score": 1e8 if fraction else 0., "quadrature_unresolved": fraction == .08})
        output = ANALYSIS.controlled_contrasts(catalogue, results, contract)
        self.assertEqual(len(output["ice_responses"]), 6)
        self.assertEqual(len(output["ice_mass_svd"]), 3)
        self.assertTrue(all(row["near_collinear_diagnostic"] for row in output["ice_mass_svd"]))
        self.assertEqual(output["missing_contrasts"], [])
        self.assertTrue(any(row["quadrature_unresolved"] for row in output["ice_responses"]))
        self.assertIsNone(output["number_of_physical_constraints"])

    def test_covariance_order_and_double_counted_gain_are_rejected(self):
        contract = contract_fixture()
        ANALYSIS.validate_contract(contract)
        bad = copy.deepcopy(contract)
        bad["covariance_scenarios"]["baseline"]["band_order"].reverse()
        with self.assertRaisesRegex(ValueError, "band order"):
            ANALYSIS.validate_contract(bad)
        bad = copy.deepcopy(contract)
        bad["covariance_scenarios"]["baseline"]["formal_plus_structure_covariance_jy2"][0][0] *= 2
        with self.assertRaisesRegex(ValueError, "do not sum"):
            ANALYSIS.validate_contract(bad)


class CompactProductionAnalysisTests(unittest.TestCase):
    def test_compact_receipts_without_fits_write_complete_screen(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            run_fixture(run)
            result = ANALYSIS.analyze_run(run)
            self.assertEqual(result["status"], "complete_screen")
            self.assertEqual(result["models_screened"], 1)
            self.assertFalse(result["raw_fits_revalidated"])
            self.assertFalse(result["production_convergence_certified"])
            self.assertEqual(list(run.rglob("*.fits*")), [])
            for name in ("summary.json", "ranking.csv", "excluded_models.csv", "band_predictions.csv", "ice_responses.csv",
                         "ice_mass_svd.csv", "REVIEW.md", "band_comparison.png", "band_comparison.pdf", "ice_response.png", "ice_response.pdf"):
                self.assertGreater((run / "results" / name).stat().st_size, 0)

    def test_missing_duplicate_wrong_distance_and_quality_are_explicitly_excluded(self):
        for change in ("missing", "duplicate", "distance", "pixel_scale", "quality"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temporary:
                run = Path(temporary)
                manifest, payload = run_fixture(run)
                if change == "missing":
                    payload["measurements"].pop()
                elif change == "duplicate":
                    payload["measurements"].append(copy.deepcopy(payload["measurements"][0]))
                elif change == "distance":
                    payload["measurements"][0]["target_distance_pc"] = 140.
                elif change == "pixel_scale":
                    payload["measurements"][0]["diagnostics"]["model_pixel_arcsec"] *= 2
                else:
                    payload["measurements"][0]["diagnostics"]["quality_checks"]["finite_kernel_aperture_support"] = False
                path = run / "models" / manifest["models"][0]["id"] / "measurements.json"
                write_json(path, payload)
                result = ANALYSIS.analyze_run(run, make_plot=False)
                self.assertEqual(result["status"], "partial_screen")
                self.assertEqual(result["models_screened"], 0)
                self.assertEqual(result["models"][0]["status"], "excluded")
                self.assertTrue(result["models"][0]["reason"])

    def test_wrong_fingerprint_and_unstarted_runs_are_not_ranked(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            manifest, payload = run_fixture(run)
            path = run / "models" / manifest["models"][0]["id"] / "measurements.json"
            payload["fingerprint"]["executable_sha256"] = "d" * 64
            write_json(path, payload)
            result = ANALYSIS.analyze_run(run, make_plot=False)
            self.assertIn("runtime binding", result["models"][0]["reason"])
            path.unlink()
            (run / "runtime_binding.json").unlink()
            result = ANALYSIS.analyze_run(run, make_plot=False)
            self.assertEqual(result["models"][0]["status"], "missing")
            self.assertEqual(result["models_screened"], 0)

    def test_frozen_contract_tampering_aborts_before_score(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            run_fixture(run)
            path = run / "inputs/production_observation_contract.json"
            payload = ANALYSIS.read_json(path)
            payload["bands"][0]["observed_flux_jy"] *= 2
            write_json(path, payload)
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                ANALYSIS.analyze_run(run, make_plot=False)

    def test_manifest_sidecar_pins_shared_configuration(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            manifest, _ = run_fixture(run)
            (run / "manifest.sha256").write_text(ANALYSIS.digest(run / "manifest.json") + "  manifest.json\n")
            manifest["measurement"]["target_distance_pc"] = 140.
            write_json(run / "manifest.json", manifest)
            with self.assertRaisesRegex(ValueError, "Manifest checksum differs"):
                ANALYSIS.analyze_run(run, make_plot=False)


if __name__ == "__main__":
    unittest.main()
