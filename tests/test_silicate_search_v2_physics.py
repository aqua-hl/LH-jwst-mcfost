"""Separate DHS carbon must change envelope composition without changing disk physics."""
import hashlib
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mcfost_grid.physics import render_parameter, validate_parameters


TEMPLATE = (ROOT / "reference/parameters/ice_v02_dust.para").read_text()
MIE = (ROOT / "reference/parameters/continuum_original_dust.para").read_text()
PARAMETERS = dict(envelope_dust_mass_msun=.000225, envelope_amax_um=.4,
                  envelope_size_exponent=2.75, inclination_deg=60,
                  envelope_ice_mass_fraction=.04)
NUMERICS = dict(photons_temperature=512000, photons_sed=512000,
                photons_image=512000, image_npix=6001,
                image_size_au=6000, grains=50)


def rows(text, label):
    return [line.split(label)[0].split() for line in text.splitlines() if label in line]


def disk(text):
    block = text.split("#Grain properties\n", 1)[1].split("#Molecular RT settings", 1)[0]
    lines = block.splitlines(keepends=True)
    boundary = next(i for i, line in enumerate(lines) if "Number of species" in line and "ZONE 2" in line)
    return "".join(lines[:boundary])


BASELINE_DISK = disk(render_parameter(TEMPLATE, PARAMETERS, NUMERICS, "temperature"))


class SilicateSearchV2PhysicsTests(unittest.TestCase):
    def test_historical_render_bytes_are_unchanged(self):
        # Golden hashes captured from the unmodified renderer, before carbon
        # support: this protects old run manifests as well as numerical values.
        hashes = {
            "dhs": ("a051683a293094785356b8a072dbdd9ffe1453c9e6c36b4b486c86053748ee80",
                    "c0a69f5724a39db6493a55cc6479a2d3a6d630fe87e918ce94a0d8e0bf78f116",
                    "3ba397e1df3f1dcadbf31d0ceaca8e804af7684192e8ea5c859b3d6dbf3e1f6e"),
            "mie": ("c209cc0c2ff9d90ab78033addb7ede90da95ae4993804496d76f7ed09b9146d7",
                    "292be3b12c73a91382e04866f7df6f1851c93e53cba588c43a53a026efc05757",
                    "cb6b223b154eeabd167709e4c5685bb1ecbcddba8f5182a58f4308b33e0b1649"),
        }
        for family, template in (("dhs", TEMPLATE), ("mie", MIE)):
            parameters = dict(PARAMETERS)
            if family == "mie":
                parameters.pop("envelope_ice_mass_fraction")
            for stage, digest in zip(("temperature", "image", "coeval"), hashes[family]):
                wave = None if stage == "temperature" else 9.7
                old = render_parameter(template, parameters, NUMERICS, stage, wave)
                self.assertEqual(hashlib.sha256(old.encode()).hexdigest(), digest)
                zero = render_parameter(template, {**parameters, "envelope_carbon_mass_fraction": 0}, NUMERICS, stage, wave)
                self.assertEqual(zero, old)

    def test_all_compositions_have_separate_pure_species_and_conserve_mass(self):
        for carbon in (0, .1, .2, .3):
            for ice in (.04, .08, .12):
                with self.subTest(carbon=carbon, ice=ice):
                    parameters = {**PARAMETERS, "envelope_carbon_mass_fraction": carbon,
                                  "envelope_ice_mass_fraction": ice}
                    output = render_parameter(TEMPLATE, parameters, NUMERICS, "temperature")
                    species = rows(output, "Grain type")[2:]
                    self.assertEqual(len(species), 3 if carbon else 2)
                    self.assertEqual(rows(output, "Number of species")[-1], [str(len(species))])
                    self.assertEqual([s[:3] for s in species], [["DHS", "1", "1"]] * len(species))
                    self.assertEqual([float(s[3]) for s in species], [0] * len(species))
                    self.assertEqual([float(s[5]) for s in species], [.1] * len(species))
                    expected = [1-ice-carbon, ice] + ([carbon] if carbon else [])
                    for actual, wanted in zip(species, expected):
                        self.assertAlmostEqual(float(actual[4]), wanted, places=14)
                    self.assertAlmostEqual(sum(float(s[4]) for s in species), 1, places=14)
                    optical = rows(output, "Optical indices file")[2:]
                    self.assertEqual([r[0] for r in optical], ["Draine_Si_sUV.dat", "H2O_30K_Leiden_mcfost.dat"] + (["ac_opct.dat"] if carbon else []))
                    self.assertEqual([float(r[1]) for r in optical], [1] * len(species))
                    self.assertEqual(disk(output), BASELINE_DISK)
                    self.assertEqual(rows(output, "Grain type")[:2], rows(TEMPLATE, "Grain type")[:2])
                    self.assertEqual(rows(output, "amin, amax")[:2], rows(TEMPLATE, "amin, amax")[:2])

    def test_carbon_sizes_follow_silicate_but_ice_stays_fixed(self):
        parameters = {**PARAMETERS, "envelope_carbon_mass_fraction": .2,
                      "envelope_amax_um": 1., "envelope_size_exponent": 3.1}
        output = render_parameter(TEMPLATE, parameters, NUMERICS, "temperature")
        sizes = rows(output, "amin, amax")[2:]
        self.assertEqual(sizes, [["0.03", "1", "3.1", "50"],
                                 ["0.03", "0.4", "2.75", "50"],
                                 ["0.03", "1", "3.1", "50"]])
        changed = render_parameter(output, {**parameters, "envelope_amax_um": .4}, NUMERICS, "temperature")
        self.assertEqual(rows(changed, "amin, amax")[2], rows(changed, "amin, amax")[4])

    def test_repeated_rendering_retains_three_species_for_each_stage(self):
        parameters = {**PARAMETERS, "envelope_carbon_mass_fraction": .2,
                      "envelope_ice_mass_fraction": .08,
                      "envelope_silicate_file": "Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat"}
        for stage in ("temperature", "image", "coeval"):
            wave = None if stage == "temperature" else 9.7
            output = render_parameter(TEMPLATE, parameters, NUMERICS, stage, wave)
            self.assertEqual(render_parameter(output, parameters, NUMERICS, stage, wave), output)
            self.assertEqual(disk(output), BASELINE_DISK)
            self.assertEqual(rows(output, "Optical indices file")[2][0], parameters["envelope_silicate_file"])

    def test_zero_carbon_removes_species_instead_of_a_zero_weight_population(self):
        output = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_carbon_mass_fraction": .2}, NUMERICS, "temperature")
        zero = render_parameter(output, {**PARAMETERS, "envelope_carbon_mass_fraction": 0}, NUMERICS, "temperature")
        self.assertEqual(rows(zero, "Number of species"), [["2"], ["2"]])
        self.assertEqual(len(rows(zero, "Optical indices file")), 4)
        self.assertEqual([float(s[4]) for s in rows(zero, "Grain type")[2:]], [.96, .04])
        self.assertEqual(disk(zero), BASELINE_DISK)
        self.assertEqual(render_parameter(zero, PARAMETERS, NUMERICS, "temperature"), zero)

    def test_invalid_fractions_fail_before_rendering(self):
        for value in (-.01, 1, float("nan"), float("inf"), True, ".2"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_parameters({"envelope_carbon_mass_fraction": value})
        for carbon, ice in ((.9, .1), (.99, .04), (.3, .7)):
            with self.subTest(carbon=carbon, ice=ice), self.assertRaisesRegex(ValueError, "sum to less than one"):
                validate_parameters({"envelope_carbon_mass_fraction": carbon, "envelope_ice_mass_fraction": ice})
        with self.assertRaisesRegex(ValueError, "ice mass fractions"):
            validate_parameters({"envelope_carbon_mass_fraction": .2, "envelope_ice_volume_fraction": .1})

    def test_unsupported_templates_and_implicit_composition_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "only for the disjoint DHS"):
            render_parameter(MIE, {"envelope_carbon_mass_fraction": .2}, NUMERICS, "temperature")
        bare = render_parameter(MIE, {"envelope_ice_volume_fraction": 0}, NUMERICS, "temperature")
        with self.assertRaisesRegex(ValueError, "only for the disjoint DHS"):
            render_parameter(bare, {"envelope_carbon_mass_fraction": .2}, NUMERICS, "temperature")
        with self.assertRaisesRegex(ValueError, "positive ice mass fraction"):
            render_parameter(TEMPLATE, {**PARAMETERS, "envelope_ice_mass_fraction": 0,
                                       "envelope_carbon_mass_fraction": .2}, NUMERICS, "temperature")
        output = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_carbon_mass_fraction": .2}, NUMERICS, "temperature")
        with self.assertRaisesRegex(ValueError, "explicit envelope_carbon_mass_fraction"):
            render_parameter(output, PARAMETERS, NUMERICS, "temperature")

    def test_mutated_third_population_fails_closed(self):
        parameters = {**PARAMETERS, "envelope_carbon_mass_fraction": .2}
        output = render_parameter(TEMPLATE, parameters, NUMERICS, "temperature")
        head, carbon = output.rsplit("DHS", 1)
        mutations = [
            carbon.replace("1  1  0.0  0.2  0.1", "2  1  0.0  0.2  0.1", 1),
            carbon.replace("ac_opct.dat  1", "ice_opct.dat  1", 1),
            carbon.replace("ac_opct.dat  1", "ac_opct.dat  0.7", 1),
            carbon.replace("0.03  0.4  2.75  50", "0.03  1  2.75  50", 1),
        ]
        for changed in mutations:
            self.assertNotEqual(changed, carbon)
            with self.assertRaises(ValueError):
                render_parameter(head + "DHS" + changed, parameters, NUMERICS, "temperature")

    def test_supported_label_case_does_not_duplicate_the_silicate_population(self):
        template = TEMPLATE.replace("Grain type", "grain type").replace("Optical indices file", "optical indices file")
        parameters = {**PARAMETERS, "envelope_carbon_mass_fraction": .2}
        output = render_parameter(template, parameters, NUMERICS, "temperature")
        self.assertEqual(rows(output, "optical indices file")[-1], ["ac_opct.dat", "1"])
        self.assertEqual(float(rows(output, "grain type")[-1][4]), .2)
        self.assertEqual(render_parameter(output, parameters, NUMERICS, "temperature"), output)


if __name__ == "__main__":
    unittest.main()
