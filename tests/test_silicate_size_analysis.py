"""Matched grain-size contrasts and mandatory compact-receipt validation."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("silicate_size_analysis", ROOT / "scripts/analyze_silicate_size_pilot.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)
PRODUCTION = ANALYSIS.production


def catalogue_fixture():
    result = []
    for inclination in (50., 70.):
        for amax in (.4, 1.):
            result.append({"index": len(result), "role": "silicate_size_pair",
                "contrast_group_id": f"inclination_{inclination:g}", "ice_mass_fraction": .08,
                "parameters": {"inclination_deg": inclination, "envelope_amax_um": amax,
                    "envelope_dust_mass_msun": .000225, "envelope_ice_volume_fraction": .2,
                    "envelope_size_exponent": 2.75, "distance_pc": 140.}})
    return result


def contract_fixture():
    bands, probes = [], []
    for prefix, block, count in (("c", "nir_continuum", 9), ("h", "nir_ice", 4), ("s", "miri", 6)):
        for j in range(count):
            lower = 1. + j * .15 if prefix == "c" else 2.7 + j * .2 if prefix == "h" else 8. + j * 1.3
            upper = lower + .1
            if prefix == "s" and j == 1:
                lower, upper = 9.3, 10.1
            band_id = f"{prefix}{j + 1:02d}"
            instrument = "MIRI" if prefix == "s" else "NIRSpec"
            ids = []
            for q, wavelength in enumerate(np.linspace(lower, upper, 5)):
                probe_id = f"{band_id}_q{q}"
                ids.append(probe_id)
                probes.append({"id": probe_id, "wavelength_um": float(wavelength), "score": False,
                    "instrument": instrument, "image_npix": 6001, "image_size_au": 6000.})
            bands.append({"id": band_id, "block": block, "instrument": instrument, "segment": "shared",
                "gain_group": instrument + "/shared", "lower_um": lower, "upper_um": upper,
                "observed_flux_jy": 1., "formal_sigma_jy": .05, "structure_sigma_jy": .01,
                "probe_ids": ids, "weights": [1/12, 4/12, 2/12, 4/12, 1/12],
                "coarse_probe_ids": ids[::2], "coarse_weights": [1/6, 4/6, 1/6]})
    for j in range(3):
        probes.append({"id": f"check{j}", "wavelength_um": 15. + j, "score": False,
                       "instrument": "MIRI", "image_npix": 6001, "image_size_au": 6000.})
    scenario = {"band_order": [b["id"] for b in bands], "baseline": True,
                "relative_gain_fraction_by_segment": {"NIRSpec/shared": .03, "MIRI/shared": .02}}
    base = np.eye(len(bands)) * .0026
    design, sigmas, _ = PRODUCTION.gain_design(bands, scenario)
    gain = (design * sigmas) @ (design * sigmas).T
    scenario.update(formal_plus_structure_covariance_jy2=base.tolist(), relative_gain_covariance_jy2=gain.tolist(),
                    conditional_covariance_jy2=(base + gain).tolist())
    return {"bands": bands, "probes": probes, "covariance_scenarios": {"baseline": scenario},
            "block_weights": PRODUCTION.WEIGHTING["baseline"],
            "common_scale_nuisance": {"kind": "gaussian_multiplicative", "mean": 1., "sigma_fraction": .03}}


def screened_fixture(catalogue, contract):
    results = []
    for record in catalogue:
        changed = record["parameters"]["envelope_amax_um"] == 1.
        flux = [.8 if b["id"] == "s02" else 1.2 for b in contract["bands"]] if changed else [.5] * 19
        results.append({"model_id": f"m{record['index']}", "model_index": record["index"], "status": "screened",
            "parameters": copy.deepcopy(record["parameters"]), "band_flux_jy": flux,
            "baseline_score": 999. if changed else 1., "quadrature_unresolved": False,
            "unresolved_quadrature_bands": []})
    return results


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False) + "\n")


def compact_fixture(run):
    catalogue = catalogue_fixture()
    contract = contract_fixture()
    write_json(run / "production_experiment.json", {"catalogue": catalogue,
        "dust_prescription": {"geometry": "DHS", "ice_material": "supplied_H2O", "ice_amax_um": .4},
        "pilot": {"size_change_scope": "silicate_only"}})
    write_json(run / "inputs/production_observation_contract.json", contract)
    shutil.copy2(ROOT / "silicate/screen_draine.json", run / "inputs/silicate_screen_draine.json")
    quality = run / "code/src/mcfost_grid/quality.py"
    quality.parent.mkdir(parents=True)
    shutil.copy2(ROOT / "src/mcfost_grid/quality.py", quality)
    manifest = {"models": [{"id": f"m{r['index']}", "index": r["index"], "parameters": r["parameters"]} for r in catalogue],
        "anchors": contract["probes"], "configuration": {"numerical_only": True,
                "numerics": {"image_npix": 6001, "image_size_au": 6000.}},
        "measurement": {"quality_policy": "aperture_v2", "target_distance_pc": 147., "aperture_radius_arcsec": 1.},
        "input_hashes": {str(p.relative_to(run)): PRODUCTION.digest(p) for p in run.rglob("*") if p.is_file()}}
    write_json(run / "manifest.json", manifest)
    fingerprint = {"manifest_sha256": PRODUCTION.digest(run / "manifest.json"), "executable_sha256": "a" * 64,
                   "utilities_content_sha256": "b" * 64, "backend": "image_method2"}
    write_json(run / "runtime_binding.json", {"fingerprint": fingerprint})
    checks = {"aperture_plus_6sigma_support": True, "finite_kernel_aperture_support": True,
              "finite_positive_plane0_flux": True, "finite_positive_signed_aperture": True,
              "edge_positive_fraction_lt_1e-3": False, "convolution_flux_loss_lt_1e-4": True}
    blocking = {k: v for k, v in checks.items() if k not in {"edge_positive_fraction_lt_1e-3", "convolution_flux_loss_lt_1e-4"}}
    pixel = 6000 / 6001 / 140
    for model in manifest["models"]:
        measurements = [{"anchor_id": p["id"], "wavelength_um": p["wavelength_um"], "instrument": p["instrument"],
            "flux_jy": 1. + .01 * model["index"], "valid": True, "quality_pass": True, "quality_policy": "aperture_v2",
            "quality_warnings": ["edge_positive_fraction_lt_1e-3"], "backend": "image_method2", "image_sha256": "c" * 64,
            "canonical_estimator": "signed_direct_method2_polarized_total_i_plane0",
            "model_distance_pc": 140., "target_distance_pc": 147., "aperture_radius_arcsec": 1.,
            "diagnostics": {"quality_checks": checks, "blocking_quality_checks": blocking,
                "distance_flux_factor": (140 / 147)**2, "zero_clipping_used_for_scoring": False,
                "nx": 6001, "ny": 6001, "centre_x_zero_based": 3000., "centre_y_zero_based": 3000.,
                "aperture_radius_pixels": 1 / (pixel * 140 / 147), "psf_sigma_pixels": 6., "psf_truncation_sigma": 6.,
                "model_pixel_arcsec": pixel, "target_pixel_arcsec": pixel * 140 / 147}} for p in contract["probes"]]
        write_json(run / "models" / model["id"] / "measurements.json", {"model_id": model["id"],
            "model_index": model["index"], "parameters": model["parameters"], "fingerprint": fingerprint,
            "complete": True, "measurements": measurements})
    return manifest


class PairedResponseTests(unittest.TestCase):
    def setUp(self):
        self.catalogue = catalogue_fixture()
        self.contract = contract_fixture()
        self.models = screened_fixture(self.catalogue, self.contract)

    def test_known_flux_response_and_absolute_residual_improvement_without_fit_selection(self):
        output = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(output["status"], "complete_pilot")
        self.assertEqual(len(output["band_responses"]), 38)
        for pair in output["pairs"]:
            s02 = pair["s02"]
            self.assertAlmostEqual(s02["changed_over_reference_flux_ratio"], 1.6)
            self.assertAlmostEqual(s02["changed_over_reference_flux_dex"], np.log10(1.6))
            self.assertAlmostEqual(s02["absolute_residual_improvement_percentage_points"], 30.)
            self.assertAlmostEqual(pair["region_residuals"]["nir_all"]["changed_fractional_residual_rms_percent"], 20.)
            self.assertEqual(len(pair["miri_bands"]), 6)
            self.assertEqual(pair["baseline_profiled_scores"], [1., 999.])

    def test_changes_in_uncontrolled_parameters_are_rejected(self):
        for key, value in (("envelope_size_exponent", 3.5), ("envelope_ice_volume_fraction", .1),
                           ("envelope_dust_mass_msun", .00015)):
            with self.subTest(parameter=key):
                bad = copy.deepcopy(self.catalogue)
                bad[3]["parameters"][key] = value
                with self.assertRaises(ValueError):
                    ANALYSIS.validate_catalogue(bad)
        bad = copy.deepcopy(self.catalogue)
        bad[3]["ice_mass_fraction"] = .04
        with self.assertRaisesRegex(ValueError, "Ice mass fraction"):
            ANALYSIS.validate_catalogue(bad)

    def test_missing_failed_pair_and_quadrature_flags_are_not_silently_dropped(self):
        self.models[0]["status"] = "excluded"
        self.models[0]["reason"] = "negative signed total-I"
        self.models[2]["quadrature_unresolved"] = True
        self.models[2]["unresolved_quadrature_bands"] = ["s02"]
        output = ANALYSIS.paired_responses(self.catalogue, self.models, self.contract)
        self.assertEqual(output["status"], "incomplete_pilot")
        self.assertEqual(output["pairs_complete"], 1)
        self.assertEqual(len(output["missing_pairs"]), 1)
        self.assertEqual(len(output["band_responses"]), 19)
        self.assertTrue(output["pairs"][0]["s02"]["reference_quadrature_unresolved"])

    def test_nonpositive_or_nonfinite_flux_is_not_a_computable_response(self):
        for value in (0., -1., float("nan"), float("inf")):
            bad = copy.deepcopy(self.models)
            bad[0]["band_flux_jy"][0] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "finite and positive"):
                ANALYSIS.paired_responses(self.catalogue, bad, self.contract)

    def test_screen_does_not_conflate_fixed_mass_and_fixed_nir_column(self):
        screen = PRODUCTION.read_json(ROOT / "silicate/screen_draine.json")
        result = ANALYSIS.opacity_screen_comparison(screen, .08)
        self.assertAlmostEqual(result["bare_grain_kappa_abs_9p7_percent_change_at_fixed_mass"], .962403, delta=.00001)
        self.assertAlmostEqual(result["bare_grain_kappa_ext_2p2_ratio_at_fixed_mass"], 3.91636, delta=.00001)
        self.assertAlmostEqual(result["bare_grain_S_9p7_per_nir_ratio_at_fixed_nir_extinction"], .2577962292)
        self.assertTrue(result["bare_screen_vs_icy_pilot_mismatch"])
        self.assertFalse(result["screen_flux_gain_dex_used_as_prediction"])


class CompactReceiptTests(unittest.TestCase):
    def test_complete_compact_run_and_required_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_fixture(run)
            result = ANALYSIS.analyze_pilot(run, make_plot=True)
            self.assertEqual(result["status"], "complete_pilot")
            self.assertFalse(result["controlling_variable_identified"])
            self.assertFalse(result["posterior_intervals_reported"])
            self.assertEqual(result["pilot"]["size_change_scope"], "silicate_only")
            self.assertEqual(result["dust_prescription"]["ice_material"], "supplied_H2O")
            inherited = PRODUCTION.read_json(run / "results/summary.json")
            self.assertNotIn("H2O30K laboratory constants are not used", " ".join(inherited["limitations"]))
            for name in ("summary.json", "pilot_summary.json", "paired_band_responses.csv", "PILOT_REVIEW.md",
                         "pilot_plot.png", "pilot_plot.pdf"):
                self.assertTrue((run / "results" / name).is_file())

    def test_mandatory_total_i_check_rejects_even_positive_aperture_flux(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_fixture(run)
            path = run / "models/m0/measurements.json"
            payload = PRODUCTION.read_json(path)
            payload["measurements"][0]["diagnostics"]["quality_checks"]["finite_positive_plane0_flux"] = False
            write_json(path, payload)
            result = ANALYSIS.analyze_pilot(run, make_plot=False)
            self.assertEqual(result["models_screened"], 3)
            self.assertEqual(result["status"], "incomplete_pilot")
            self.assertEqual(result["pairs_complete"], 1)

    def test_changed_frozen_screen_and_zero_observation_error_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            compact_fixture(run)
            screen = run / "inputs/silicate_screen_draine.json"
            screen.write_text(screen.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                ANALYSIS.analyze_pilot(run, make_plot=False)
        contract = contract_fixture()
        contract["bands"][0]["formal_sigma_jy"] = 0.
        with self.assertRaisesRegex(ValueError, "Invalid band observation"):
            PRODUCTION.validate_contract(contract)


if __name__ == "__main__":
    unittest.main()
