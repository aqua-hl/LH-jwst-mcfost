"""V3 mass/ice-column conservation and seed identity, with synthetic test opacities."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "src")]

from mcfost_grid.physics import render_parameter
from silicate_search_v3_design import (CONFIGURATIONS, DUST, MATERIALS, NAME,
    NUMERICS, REFERENCE_MASS, SPECIES, SUCCESS, catalogue, declared_cases, digest,
    matched_composition, model_numerics, validate_frozen_campaign, validate_opacities,
    validate_targets)

SOURCE = ROOT / "runs/silicate_search_v2_512k"
TARGETS = {"draine": 9.144280744583792e-6, "pyroxene_mg50": 1.1086300567771168e-5}


def anchors():
    return json.loads((SOURCE/"inputs/production_anchors.json").read_text())


def synthetic_opacities():
    """Arbitrary test-only spectra; these values never become production inputs."""
    wave = sorted(set([a["wavelength_um"] for a in anchors()]+[2.2, 9.7, 18.]))
    result = {}
    for species, scale, albedo in (("draine", 1500., .2), ("pyroxene_mg50", 1000., .3),
                                    ("carbon", 10000., .4), ("ice", 200., .1)):
        result[species] = dict(wavelength_um=list(wave), kappa_ext_cm2_g=[scale]*99,
            kappa_abs_cm2_g=[scale*(1-albedo)]*99, kappa_sca_cm2_g=[scale*albedo]*99,
            albedo=[albedo]*99)
    return result


def frozen_fixture(root):
    inputs = root/"inputs"
    dust = inputs/"utils/Dust"
    dust.mkdir(parents=True)
    shutil.copy2(SOURCE/"inputs/template.para", inputs/"template.para")
    for species in SPECIES.values():
        shutil.copy2(SOURCE/"inputs/utils/Dust"/species["filename"], dust)
    opacities = synthetic_opacities()
    records = catalogue(opacities, TARGETS)
    wave = opacities["draine"]["wavelength_um"]
    for filename, value in {
        "stage0_opacities.json": opacities,
        "stage0_wavelengths.json": {"wavelengths_um": wave},
        "stage0_receipt.json": {"status": "passed", "wavelengths_um": wave,
                                "species": [{"id": s} for s in SPECIES]},
        "v2_reference.json": {"ice_column_targets_msun": TARGETS},
    }.items():
        (inputs/filename).write_text(json.dumps(value))
    models = [dict(index=row["index"], id=f"m{row['index']:04d}_test",
                   parameters=row["parameters"], numerics=model_numerics(row)) for row in records]
    manifest = dict(models=models, anchors=anchors(), configuration={"numerics": NUMERICS},
                    temperature_policy="fresh equilibrium for every physical model")
    experiment = dict(experiment_id=NAME, catalogue=records, stage0_opacity=opacities,
                      ice_column_targets_msun=TARGETS, stage0_receipt_sha256=digest(inputs/"stage0_receipt.json"),
                      opacity_matching=[dict(index=r["index"], **r["opacity_matching"]) for r in records],
                      dust_prescription=DUST, success_criteria=SUCCESS, numerics=NUMERICS)
    # index=0 validates all identities and only its rendered files. This keeps
    # unit fixtures small; the package test covers full-grid materialization.
    model = models[0]
    directory = root/"models"/model["id"]
    directory.mkdir(parents=True)
    template = (inputs/"template.para").read_text()
    (directory/"temperature.para").write_text(render_parameter(template, model["parameters"], model["numerics"], "temperature"))
    for anchor in manifest["anchors"]:
        folder = directory/"anchors"/anchor["id"]
        folder.mkdir(parents=True)
        for stage in ("image", "coeval"):
            (folder/(stage+".para")).write_text(render_parameter(template, model["parameters"], model["numerics"], stage, anchor["wavelength_um"]))
    return experiment, manifest


class SilicateSearchV3DesignTests(unittest.TestCase):
    def test_ordered_24_physics_and_two_seed_replicas_are_declared_before_opacities(self):
        cases = declared_cases()
        self.assertEqual(len(cases), 26)
        self.assertEqual([r["index"] for r in cases], list(range(26)))
        expected = [(config, exponent, factor, material, carbon)
                    for config, exponent, factor in CONFIGURATIONS
                    for material in MATERIALS for carbon in (.15, .2, .25)]
        self.assertEqual([(r["configuration_id"], r["density_exponent"], r["column_factor"],
                           r["material_id"], r["carbon_mass_fraction"]) for r in cases[:24]], expected)
        self.assertTrue(all(r["inclination_deg"] == 70 and r["random_seed"] == 43001
                            and r["replica_of"] is None and r["role"] == "primary" for r in cases[:24]))
        self.assertTrue(all("parameters" not in r and "ice_mass_fraction" not in r for r in cases))
        for index, original in ((24, 7), (25, 10)):
            self.assertEqual(cases[index]["replica_of"], original)
            self.assertEqual(cases[index]["random_seed"], 43002)
            self.assertEqual(cases[index]["configuration_id"], "shallow")
            self.assertEqual(cases[index]["carbon_mass_fraction"], .2)

    def test_ice_closed_form_matches_independent_root_solve_and_exact_targets(self):
        opacities = synthetic_opacities()
        ext = {key: value["kappa_ext_cm2_g"][0] for key, value in opacities.items()}
        kref = .96*ext["draine"] + .04*ext["ice"]
        for record in catalogue(opacities, TARGETS):
            material = record["material_id"]
            target = TARGETS[material]
            carbon, c = record["carbon_mass_fraction"], record["column_factor"]

            def residual(ice):
                kmix = (1-carbon-ice)*ext[material] + carbon*ext["carbon"] + ice*ext["ice"]
                return ice*REFERENCE_MASS*c*kref/kmix-target

            solved = brentq(residual, 1e-12, 1-carbon-1e-12, xtol=1e-15)
            self.assertAlmostEqual(record["ice_mass_fraction"]/solved, 1., delta=7e-13)
            match = record["opacity_matching"]
            self.assertEqual(match["ice_column_target_msun"], target)
            self.assertAlmostEqual(match["radial_ice_column_msun"]/target, 1., delta=2e-12)
            self.assertAlmostEqual(record["silicate_mass_fraction"]+carbon+record["ice_mass_fraction"], 1., places=12)
            self.assertEqual(record["parameters"]["envelope_ice_mass_fraction"], record["ice_mass_fraction"])
            self.assertFalse(match["discrete_optical_depth_match_certified"])

    def test_dust_mass_matches_column_using_independent_profile_integrals(self):
        records = catalogue(synthetic_opacities(), TARGETS)

        def integral(power):
            return quad(lambda log_r: np.exp((power+1)*log_r), 0, np.log(3000.),
                        epsabs=1e-11, epsrel=1e-12)[0]

        base_column_per_mass = integral(-1.5)/integral(.5)
        for record in records:
            exponent = record["density_exponent"]
            match = record["opacity_matching"]
            target_column = (REFERENCE_MASS*record["column_factor"]*base_column_per_mass
                             *match["reference_extinction_cm2_g"])
            actual_column = (match["actual_mass_msun"]*integral(exponent)/integral(exponent+2)
                             *match["mixture_extinction_cm2_g"])
            self.assertAlmostEqual(actual_column/target_column, 1., delta=2e-12)
            if exponent == -1.25:
                self.assertEqual(match["profile_mass_factor"], 3.60037382695)
            if exponent == -1.:
                self.assertEqual(match["profile_mass_factor"], 10.07433956241)
            self.assertEqual(match["absolute_ice_mass_msun"], match["actual_mass_msun"]*record["ice_mass_fraction"])

    def test_seed_replicas_share_physics_and_change_only_numeric_seed(self):
        records = catalogue(synthetic_opacities(), TARGETS)
        self.assertEqual(len({json.dumps(r["parameters"], sort_keys=True) for r in records}), 24)
        self.assertEqual(len({json.dumps((r["parameters"], r["numerics"]), sort_keys=True) for r in records}), 26)
        for replica, original in ((24, 7), (25, 10)):
            self.assertEqual(records[replica]["parameters"], records[original]["parameters"])
            self.assertEqual(records[replica]["opacity_matching"], records[original]["opacity_matching"])
            self.assertEqual(records[replica]["numerics"], {**records[original]["numerics"], "random_seed": 43002})
        for record in records:
            self.assertNotIn("random_seed", record["parameters"])
            self.assertEqual(record["numerics"]["photons_temperature"], 512000)
            self.assertEqual(record["numerics"]["photons_sed"], 512000)
            self.assertEqual(record["numerics"]["photons_image"], 512000)
        for seed in (True, 43003, "43001", 43001.):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                model_numerics({"random_seed": seed})

    def test_lower_column_increases_required_ice_fraction_and_reduces_mass(self):
        records = catalogue(synthetic_opacities(), TARGETS)
        for offset in range(6):
            high, low = records[12+offset], records[18+offset]
            self.assertGreater(low["ice_mass_fraction"], high["ice_mass_fraction"])
            self.assertLess(low["parameters"]["envelope_dust_mass_msun"], high["parameters"]["envelope_dust_mass_msun"])
            self.assertAlmostEqual(low["opacity_matching"]["radial_ice_column_msun"]/
                                   high["opacity_matching"]["radial_ice_column_msun"], 1., delta=2e-12)

    def test_targets_are_positive_exact_per_material_and_impossible_mixture_rejected(self):
        for targets in ({"draine": 1e-5}, {**TARGETS, "olivine_mg50": 1e-5}):
            with self.assertRaisesRegex(ValueError, "exact Draine and pyroxene"):
                validate_targets(targets)
        for value in (0, -1., True, float("nan"), float("inf"), "0.00001"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                catalogue(synthetic_opacities(), {**TARGETS, "draine": value})
        with self.assertRaisesRegex(ValueError, "cannot be reached"):
            matched_composition(synthetic_opacities(), "draine", .2, -1.25, .86, 1.)
        for carbon in (-.1, 1., True, float("nan")):
            with self.assertRaises(ValueError):
                matched_composition(synthetic_opacities(), "draine", carbon, -1.5, 1., TARGETS["draine"])

    def test_exact_four_species_nonnegative_closure_and_99_samples_required(self):
        valid = synthetic_opacities()
        self.assertIs(validate_opacities(valid), valid)
        for mutation, message in (
            (lambda d: d.pop("ice"), "all four"),
            (lambda d: d["ice"]["wavelength_um"].pop(), "99"),
            (lambda d: d["ice"]["kappa_ext_cm2_g"].__setitem__(0, float("nan")), "Nonfinite or negative"),
            (lambda d: d["ice"]["kappa_abs_cm2_g"].__setitem__(0, -.1), "Nonfinite or negative"),
            (lambda d: d["ice"]["albedo"].__setitem__(0, 1.1), "Invalid DHS albedo"),
            (lambda d: d["ice"]["kappa_ext_cm2_g"].__setitem__(0, 300.), "do not sum"),
            (lambda d: d["ice"]["albedo"].__setitem__(0, .2), "albedo disagrees"),
        ):
            bad = deepcopy(valid)
            mutation(bad)
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                validate_opacities(bad)
        missing = deepcopy(valid)
        for row in missing.values():
            row["wavelength_um"][row["wavelength_um"].index(2.2)] = 2.201
        with self.assertRaisesRegex(ValueError, "Missing or ambiguous"):
            validate_opacities(missing)

    def test_frozen_campaign_validates_params_target_and_replica_identity(self):
        with tempfile.TemporaryDirectory(prefix="v3-design-synthetic-") as tmp:
            root = Path(tmp)
            experiment, manifest = frozen_fixture(root)
            validate_frozen_campaign(root, experiment, manifest, index=0)
            for kind in ("seed", "physics", "id", "target"):
                bad_exp, bad_man = deepcopy(experiment), deepcopy(manifest)
                if kind == "seed":
                    bad_man["models"][24]["numerics"]["random_seed"] = 43001
                elif kind == "physics":
                    bad_man["models"][24]["parameters"]["envelope_dust_mass_msun"] *= 1.01
                elif kind == "id":
                    bad_man["models"][24]["id"] = bad_man["models"][7]["id"]
                else:
                    bad_exp["ice_column_targets_msun"]["draine"] = 9.14e-6
                with self.subTest(kind=kind), self.assertRaises(ValueError):
                    validate_frozen_campaign(root, bad_exp, bad_man, index=0)
            template = root/"models"/manifest["models"][0]["id"]/"temperature.para"
            template.write_text(template.read_text()+"\n")
            with self.assertRaisesRegex(ValueError, "temperature input changed"):
                validate_frozen_campaign(root, experiment, manifest, index=0)


if __name__ == "__main__":
    unittest.main()
