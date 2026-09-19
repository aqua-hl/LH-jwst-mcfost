"""Scientific observation operators and nuisance contracts for 1 AU production."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.table import Table


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("production_observables", ROOT / "scripts/production_observables.py")
observables = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(observables)


def synthetic_tables():
    table = Table({
        "instrument": ["MIRI"] * 4,
        "segment": ["CH2_medium"] * 3 + ["CH2_long"],
        "bin_index": [1, 2, 3, 1], "wavelength_um": [9., 9.6, 10.2, 9.6],
        "bin_left_um": [8.7, 9.3, 9.9, 9.3], "bin_right_um": [9.3, 9.9, 10.5, 9.9],
        "flux_jy": [2., 4., 8., 1000.], "formal_error_jy": [.2, .4, .8, 10.],
        "within_bin_scatter_jy": [.02, .04, .08, 10.],
        "n_used": [20] * 4, "n_manual_masked": [1, 2, 3, 0],
    })
    native = Table({"instrument": ["MIRI"] * 4, "segment": ["CH2_medium"] * 4,
                    "wavelength_um": [8.7, 9.3, 9.9, 10.5], "flux_jy": [2., 4., 8., 9.],
                    "error_jy": [.1] * 4, "coverage_fraction": [1.] * 4})
    return table, native


class ObservationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = observables.build_contract(ROOT)

    def test_partial_bins_are_width_means_without_stitching_or_sqrt_n_reduction(self):
        table, native = synthetic_tables()
        band = observables.select_window(table, native, ("test", "miri", "MIRI", 9., 10.2, "CH2_medium"))
        np.testing.assert_allclose([r["normalized_weight"] for r in band["selected_rows"]], [.25, .5, .25])
        self.assertAlmostEqual(band["observed_flux_jy"], 4.5)
        self.assertAlmostEqual(band["formal_sigma_jy"], .45)
        self.assertAlmostEqual(band["structure_fraction"], .01)
        self.assertEqual([r["n_manual_masked"] for r in band["selected_rows"]], [1, 2, 3])
        self.assertEqual(band["gain_group"], "MIRI/CH2_medium")

    def test_gaps_duplicates_and_native_straddles_cannot_be_silently_interpolated(self):
        table, native = synthetic_tables()
        window = ("test", "miri", "MIRI", 9., 10.2, "CH2_medium")
        with self.assertRaisesRegex(ValueError, "incomplete R100"):
            observables.select_window(table[[0, 2, 3]], native, window)
        with self.assertRaisesRegex(ValueError, "duplicate or overlapping"):
            observables.select_window(table[[0, 1, 1, 2, 3]], native, window)
        with self.assertRaisesRegex(ValueError, "exceeds native support"):
            observables.select_window(table, native, ("test", "miri", "MIRI", 9., 10.6, "CH2_medium"))
        native["coverage_fraction"][-1] = .5
        with self.assertRaisesRegex(ValueError, "exceeds native support"):
            observables.select_window(table, native, window)

    def test_empty_or_nonfinite_bin_is_rejected_instead_of_reweighted(self):
        table, native = synthetic_tables()
        window = ("test", "miri", "MIRI", 9., 10.2, "CH2_medium")
        table["n_used"][1] = 0
        with self.assertRaisesRegex(ValueError, "no retained samples"):
            observables.select_window(table, native, window)
        table["n_used"][1] = 20
        table["flux_jy"][1] = np.nan
        with self.assertRaisesRegex(ValueError, "positive finite flux"):
            observables.select_window(table, native, window)

    def test_contract_is_complete_feature_inclusive_and_keeps_legacy_checks_unweighted(self):
        contract = self.contract
        self.assertEqual(contract["probe_count"], 98)
        self.assertEqual(contract["primary_band_count"], 19)
        self.assertEqual({block: sum(b["block"] == block for b in contract["bands"])
                          for block in contract["block_weights"]},
                         {"nir_continuum": 9, "nir_ice": 4, "miri": 6})
        self.assertEqual(contract["block_weights"], {"nir_continuum": .45, "nir_ice": .30, "miri": .25})
        self.assertFalse(any(p["score"] for p in contract["probes"]))
        self.assertTrue(all(c["primary_score_weight"] == 0 for c in contract["consistency_checks"]))
        by_id = {b["id"]: b for b in contract["bands"]}
        self.assertLess(by_id["s02"]["lower_um"], 9.7)
        self.assertGreater(by_id["s02"]["upper_um"], 9.7)
        self.assertEqual((by_id["s05"]["lower_um"], by_id["s05"]["upper_um"], by_id["s05"]["segment"]),
                         (18.1, 18.8, "CH4_short"))
        self.assertTrue(all(not b["continuum_feature_mask_applied"] for b in contract["bands"]
                            if b["block"] in {"nir_ice", "miri"}))
        self.assertEqual(contract["aperture_radius_arcsec"], 1.)
        self.assertEqual((contract["model_distance_pc"], contract["target_distance_pc"]), (140., 147.))
        self.assertFalse(contract["h2o_30k_constants_required"])

    def test_gain_groups_correlate_across_nir_blocks_and_common_scale_not_double_counted(self):
        bands = self.contract["bands"]
        indices = {b["id"]: i for i, b in enumerate(bands)}
        scenario = self.contract["covariance_scenarios"]["miri_gain_0.02"]
        gain = np.asarray(scenario["relative_gain_covariance_jy2"])
        i, j = indices["c08"], indices["h02"]
        self.assertAlmostEqual(gain[i, j], .03 ** 2 * bands[i]["observed_flux_jy"] * bands[j]["observed_flux_jy"])
        self.assertEqual(gain[indices["c01"], indices["h02"]], 0.)
        i, j = indices["s05"], indices["s06"]
        self.assertAlmostEqual(gain[i, j], .02 ** 2 * bands[i]["observed_flux_jy"] * bands[j]["observed_flux_jy"])
        base = np.asarray(scenario["formal_plus_structure_covariance_jy2"])
        np.testing.assert_allclose(scenario["conditional_covariance_jy2"], base + gain)
        self.assertTrue(np.all(np.linalg.eigvalsh(base + gain) > 0))
        self.assertFalse(scenario["common_scale_in_covariance"])
        self.assertEqual(self.contract["common_scale_nuisance"]["sigma_fraction"], .03)

    def test_nested_simpson_estimates_are_means_and_the_difference_can_reveal_curvature(self):
        probes = {p["id"]: p for p in self.contract["probes"]}
        for band in self.contract["bands"]:
            lo, hi = band["lower_um"], band["upper_um"]
            for degree in range(4):
                values = [(probes[p]["wavelength_um"] - lo) ** degree for p in band["probe_ids"]]
                result = np.dot(values, band["weights"])
                exact_mean = (hi - lo) ** degree / (degree + 1)
                self.assertAlmostEqual(result / exact_mean, 1., places=11)
            self.assertEqual(band["coarse_probe_ids"], band["probe_ids"][::2])
            # x^4 is integrated differently by 3/5-node Simpson and is not its midpoint value.
            fine = np.dot(np.linspace(0., 1., 5) ** 4, band["weights"])
            coarse = np.dot(np.linspace(0., 1., 3) ** 4, band["coarse_weights"])
            self.assertLess(abs(fine - .2), abs(coarse - .2))
            self.assertGreater(abs(fine - .5 ** 4), .1)

    def test_clean_continuum_retains_five_original_nir_anchor_neighborhoods(self):
        continuum = [b for b in self.contract["bands"] if b["block"] == "nir_continuum"]
        for wave in (1.2028117593024095, 1.4987957445436022, 1.7588530996181813,
                     2.1916659199557045, 2.546352514425122):
            self.assertTrue(any(b["lower_um"] <= wave <= b["upper_um"] for b in continuum))
        for probe in self.contract["probes"]:
            self.assertEqual((probe["image_npix"], probe["image_size_au"]), (6001, 6000.))

    def test_inputs_are_immutable_and_mask_hash_is_required(self):
        before = {name: observables.digest(ROOT / name) for name in self.contract["sources"]}
        self.assertEqual(before, observables.build_contract(ROOT)["sources"])
        self.assertEqual(before, {name: observables.digest(ROOT / name) for name in before})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "reference/observations"
            directory.mkdir(parents=True)
            for name in ("continuum_sed_R100.ecsv", "tmc1a_sed_unstitched.ecsv"):
                (directory / name).symlink_to(ROOT / "reference/observations" / name)
            mask = json.loads((ROOT / "reference/observations/continuum_fit_mask_v1_manifest.json").read_text())
            mask["source"]["sha256"] = "0" * 64
            (directory / "continuum_fit_mask_v1_manifest.json").write_text(json.dumps(mask))
            with self.assertRaisesRegex(ValueError, "source hash"):
                observables.build_contract(root)


if __name__ == "__main__":
    unittest.main()
