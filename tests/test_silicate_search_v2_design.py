"""Synthetic opacities exercise catalogue arithmetic; never used to prepare a run."""
from copy import deepcopy
from itertools import product
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "src")]

from silicate_search_v2_design import (DUST, NUMERICS, REFERENCE_MASS, SPECIES,
                                      catalogue, declared_cases,
                                      matched_composition, sample,
                                      validate_opacities)


def synthetic_opacities():
    """Arbitrary positive test values, deliberately distinct in spectral shape."""
    extinctions = {
        "draine": [3, 10, 12, 30, 15, 7],
        "pyroxene_mg50": [7, 20, 19, 40, 26, 9],
        "olivine_mg50": [5, 30, 24, 55, 44, 8],
        "carbon": [13, 80, 70, 12, 9, 4],
        "ice": [1, 2, 9, 3, 4, 2],
    }
    albedos = {
        "draine": [.1, .2, .3, .4, .5, .6],
        "pyroxene_mg50": [.3, .4, .2, .15, .1, .3],
        "olivine_mg50": [.2, .1, .35, .3, .25, .4],
        "carbon": [.05, .07, .09, .11, .13, .15],
        "ice": [.6, .5, .4, .3, .2, .1],
    }
    result = {}
    for species, values in extinctions.items():
        extinction = np.asarray(values, dtype=float)
        albedo = np.asarray(albedos[species])
        result[species] = dict(wavelength_um=[1.5, 2.2, 3.11, 9.7, 18., 25.],
                               kappa_ext_cm2_g=extinction.tolist(),
                               kappa_abs_cm2_g=(extinction*(1-albedo)).tolist(),
                               kappa_sca_cm2_g=(extinction*albedo).tolist(),
                               albedo=albedo.tolist())
    return result


def key(case):
    return (case["material_id"], case["carbon_mass_fraction"],
            case["ice_mass_fraction"], case["inclination_deg"],
            case["density_exponent"])


class SilicateSearchV2DesignTests(unittest.TestCase):
    def test_exact_axes_unique_models_and_two_shared_cases(self):
        cases = declared_cases()
        stage1 = {("draine", carbon, ice, inclination, -1.5)
                  for carbon, ice, inclination in product((0., .1, .2, .3), (.04, .08, .12), (60., 70.))}
        stage2 = {(material, .2, ice, 60., exponent)
                  for material, ice, exponent in product(("draine", "pyroxene_mg50", "olivine_mg50"), (.08, .12), (-1.5, -1.75))}
        self.assertEqual(len(cases), 34)
        self.assertEqual([r["index"] for r in cases], list(range(34)))
        self.assertEqual({key(r) for r in cases}, stage1 | stage2)
        self.assertEqual({key(r) for r in cases if 1 in r["stages"]}, stage1)
        self.assertEqual({key(r) for r in cases if 2 in r["stages"]}, stage2)
        self.assertEqual({key(r) for r in cases if len(r["stages"]) == 2}, stage1 & stage2)
        self.assertEqual(len(stage1 & stage2), 2)
        self.assertEqual(sum(r["role"] == "composition" for r in cases), 24)
        self.assertEqual(sum(r["role"] == "material_structure" for r in cases), 10)
        self.assertNotIn("parameters", cases[0])  # Masses do not exist before measured DHS input.

    def test_five_synthetic_tables_validate_and_sample_without_interpolation(self):
        opacities = synthetic_opacities()
        self.assertIs(validate_opacities(opacities), opacities)
        self.assertEqual(set(opacities), set(SPECIES))
        for species, wavelength in product(SPECIES, (2.2, 9.7, 18.)):
            expected_index = {2.2: 1, 9.7: 3, 18.: 4}[wavelength]
            for field in ("kappa_ext_cm2_g", "kappa_abs_cm2_g", "kappa_sca_cm2_g", "albedo"):
                self.assertEqual(sample(opacities, species, field, wavelength), opacities[species][field][expected_index])
        with self.assertRaisesRegex(ValueError, "Missing or ambiguous"):
            sample(opacities, "draine", "kappa_ext_cm2_g", 2.25)

    def test_reference_reproduces_mass_and_all_models_match_extinction_column(self):
        opacities = synthetic_opacities()
        reference_extinction = .96*10 + .04*2
        reference_absorption = .96*30*.6 + .04*3*.7
        reference = matched_composition(opacities, "draine", 0., .04)
        self.assertEqual(reference["actual_mass_msun"], REFERENCE_MASS)
        self.assertAlmostEqual(reference["S_over_reference"], 1., places=14)
        for model in catalogue(opacities):
            match = model["opacity_matching"]
            fractions = {model["material_id"]: model["silicate_mass_fraction"],
                         "carbon": model["carbon_mass_fraction"], "ice": model["ice_mass_fraction"]}
            extinction = sum(f*opacities[s]["kappa_ext_cm2_g"][1] for s, f in fractions.items())
            abs_9p7 = sum(f*opacities[s]["kappa_abs_cm2_g"][3] for s, f in fractions.items())
            abs_18 = sum(f*opacities[s]["kappa_abs_cm2_g"][4] for s, f in fractions.items())
            profile = 1. if model["density_exponent"] == -1.5 else .2393567089536
            self.assertEqual(match["profile_mass_factor"], profile)
            self.assertAlmostEqual(match["actual_mass_msun"]*extinction/profile/
                                   (REFERENCE_MASS*reference_extinction), 1., delta=6e-13)
            self.assertAlmostEqual(match["mixture_extinction_cm2_g"], extinction, places=13)
            self.assertAlmostEqual(match["S"], abs_9p7/extinction, places=13)
            self.assertAlmostEqual(match["S_over_reference"], abs_9p7/extinction/(reference_absorption/reference_extinction), places=13)
            self.assertAlmostEqual(match["absorption_18_over_9p7"], abs_18/abs_9p7, places=13)
            self.assertEqual(match["absolute_ice_mass_msun"], match["actual_mass_msun"]*model["ice_mass_fraction"])
            self.assertFalse(match["discrete_optical_depth_match_certified"])
            self.assertEqual(match["mass_significant_digits"], 13)

    def test_materials_use_their_own_opacity_and_fixed_inputs_are_carried_over(self):
        cases = catalogue(synthetic_opacities())
        material_cases = [r for r in cases if r["carbon_mass_fraction"] == .2 and
                          r["ice_mass_fraction"] == .08 and r["inclination_deg"] == 60 and
                          r["density_exponent"] == -1.5]
        self.assertEqual(len(material_cases), 3)
        self.assertEqual(len({r["opacity_matching"]["actual_mass_msun"] for r in material_cases}), 3)
        self.assertEqual(len({r["opacity_matching"]["absorption_18_over_9p7"] for r in material_cases}), 3)
        for case in cases:
            params = case["parameters"]
            self.assertEqual(params["envelope_silicate_file"], SPECIES[case["material_id"]]["filename"])
            self.assertEqual(params["envelope_dust_mass_msun"], case["opacity_matching"]["actual_mass_msun"])
            self.assertEqual(params["envelope_carbon_mass_fraction"], case["carbon_mass_fraction"])
            self.assertEqual(params["envelope_ice_mass_fraction"], case["ice_mass_fraction"])
            self.assertEqual(params["envelope_amax_um"], .4)
            self.assertEqual(params["envelope_size_exponent"], 2.75)
            self.assertEqual(params["cavity_half_opening_deg"], 17.5)
        self.assertEqual(NUMERICS["photons_image"], 512000)
        self.assertEqual(NUMERICS["photons_temperature"], 512000)
        self.assertEqual(NUMERICS["photons_sed"], 512000)
        self.assertEqual(NUMERICS["random_seed"], 43001)
        self.assertEqual(DUST["vmax"], .1)
        self.assertEqual(DUST["species"]["ice"]["density_g_cm3"], .94)
        self.assertFalse(DUST["optical_tables_modified"])

    def test_missing_species_or_mismatched_grids_are_rejected(self):
        missing = synthetic_opacities()
        missing.pop("ice")
        with self.assertRaisesRegex(ValueError, "all five"):
            validate_opacities(missing)
        different = synthetic_opacities()
        different["ice"]["wavelength_um"][0] = 1.6
        with self.assertRaisesRegex(ValueError, "different wavelength grids"):
            validate_opacities(different)
        for value in (float("nan"), float("inf"), 0., -1.):
            bad = synthetic_opacities()
            bad["draine"]["wavelength_um"][0] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "wavelength grid"):
                validate_opacities(bad)

    def test_missing_and_ambiguous_normalization_samples_fail(self):
        for wavelength, index in ((2.2, 1), (9.7, 3), (18., 4)):
            missing = synthetic_opacities()
            for row in missing.values():
                row["wavelength_um"][index] = wavelength*1.001
            with self.subTest(missing=wavelength), self.assertRaisesRegex(ValueError, "Missing or ambiguous"):
                validate_opacities(missing)
            ambiguous = synthetic_opacities()
            for row in ambiguous.values():
                for field, values in row.items():
                    values.insert(index+1, wavelength*(1+1e-6) if field == "wavelength_um" else values[index])
            with self.subTest(ambiguous=wavelength), self.assertRaisesRegex(ValueError, "Missing or ambiguous"):
                validate_opacities(ambiguous)
            with self.assertRaisesRegex(ValueError, "Missing or ambiguous"):
                sample(ambiguous, "draine", "kappa_ext_cm2_g", wavelength)

    def test_nonfinite_negative_mismatched_arrays_and_invalid_albedo_fail(self):
        for field, value in product(("kappa_ext_cm2_g", "kappa_abs_cm2_g", "kappa_sca_cm2_g", "albedo"),
                                    (float("nan"), float("inf"), -.1)):
            bad = synthetic_opacities()
            bad["ice"][field][0] = value
            with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, "Nonfinite or negative"):
                validate_opacities(bad)
        bad = synthetic_opacities()
        bad["ice"]["albedo"][0] = 1.1
        with self.assertRaisesRegex(ValueError, "Invalid DHS albedo"):
            validate_opacities(bad)
        bad = synthetic_opacities()
        bad["ice"]["kappa_ext_cm2_g"].pop()
        with self.assertRaisesRegex(ValueError, "Nonfinite or negative"):
            validate_opacities(bad)

    def test_opacity_closure_must_agree_in_both_independent_relations(self):
        bad = synthetic_opacities()
        bad["carbon"]["kappa_ext_cm2_g"][1] *= 1.1
        with self.assertRaisesRegex(ValueError, "do not sum to extinction"):
            validate_opacities(bad)
        bad = synthetic_opacities()
        bad["carbon"]["albedo"][1] *= .5
        with self.assertRaisesRegex(ValueError, "albedo disagrees"):
            validate_opacities(bad)

    def test_zero_ice_absorption_is_allowed_but_zero_matching_extinction_is_not(self):
        opacities = synthetic_opacities()
        row = opacities["ice"]
        row["kappa_abs_cm2_g"][3] = 0.
        row["kappa_sca_cm2_g"][3] = row["kappa_ext_cm2_g"][3]
        row["albedo"][3] = 1.
        validate_opacities(opacities)
        self.assertEqual(len(catalogue(opacities)), 34)
        zero = deepcopy(opacities)
        for row in zero.values():
            for field in ("kappa_ext_cm2_g", "kappa_abs_cm2_g", "kappa_sca_cm2_g"):
                row[field][1] = 0.
        validate_opacities(zero)  # Zero is allowed locally, but cannot normalize a mass.
        with self.assertRaisesRegex(ValueError, "Invalid DHS matching opacity"):
            catalogue(zero)


if __name__ == "__main__":
    unittest.main()
