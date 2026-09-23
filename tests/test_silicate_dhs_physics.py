"""Checks for the bounded supplied-H2O silicate-size diagnostic renderer."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mcfost_grid.physics import render_parameter, validate_parameters


TEMPLATE = (ROOT / "reference/parameters/ice_v02_dust.para").read_text()
MIE_TEMPLATE = (ROOT / "reference/parameters/continuum_original_dust.para").read_text()
PARAMETERS = {"envelope_dust_mass_msun": .000225, "envelope_ice_mass_fraction": .04,
              "envelope_amax_um": .4, "envelope_size_exponent": 2.75,
              "inclination_deg": 50}


def rows(text, tag):
    return [line.split(tag)[0].split() for line in text.splitlines() if tag in line]


class SilicateDHSPhysicsTests(unittest.TestCase):
    def test_silicate_size_pair_changes_only_silicate_size(self):
        small = render_parameter(TEMPLATE, PARAMETERS, {}, "temperature")
        large = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_amax_um": 1.0}, {}, "temperature")
        differing = [(a, b) for a, b in zip(small.splitlines(), large.splitlines()) if a != b]
        self.assertEqual(len(differing), 1)
        self.assertIn("0.03  0.4  2.75  50", differing[0][0])
        self.assertIn("0.03  1  2.75  50", differing[0][1])
        self.assertEqual(rows(small, "amin, amax")[-1], ["0.03", "0.4", "2.75", "50"])
        self.assertEqual(rows(small, "Grain type")[:2], rows(TEMPLATE, "Grain type")[:2])
        self.assertEqual(rows(small, "Optical indices file"), rows(TEMPLATE, "Optical indices file"))
        self.assertEqual(rows(small, "dust mass,"), rows(TEMPLATE, "dust mass,"))
        self.assertEqual(rows(small, "compute temperature?"), [["T", "F", "T"]])

    def test_silicate_exponent_does_not_change_ice_distribution(self):
        result = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_size_exponent": 3.1}, {}, "temperature")
        sizes = rows(result, "amin, amax")
        self.assertEqual(sizes[-2], ["0.03", "0.4", "3.1", "50"])
        self.assertEqual(sizes[-1], ["0.03", "0.4", "2.75", "50"])

    def test_explicit_mass_fraction_changes_only_species_weights(self):
        result = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_ice_mass_fraction": .08}, {}, "temperature")
        species = rows(result, "Grain type")[-2:]
        self.assertEqual(species[0], ["DHS", "1", "1", "0.0", "0.92", "0.1"])
        self.assertEqual(species[1], ["DHS", "1", "1", "0.0", "0.08", "0.1"])
        self.assertEqual(rows(result, "dust mass,"), rows(TEMPLATE, "dust mass,"))
        self.assertEqual(rows(result, "Optical indices file"), rows(TEMPLATE, "Optical indices file"))

    def test_stages_render_and_can_be_rerendered(self):
        for stage, expected in (("temperature", ["T", "F", "T"]),
                                ("image", ["F", "F", "T"]),
                                ("coeval", ["F", "T", "F"])):
            wave = None if stage == "temperature" else 9.7
            with self.subTest(stage=stage):
                result = render_parameter(TEMPLATE, {**PARAMETERS, "envelope_amax_um": 1},
                                          {"photons_temperature": 2048000, "grains": 50}, stage, wave)
                self.assertEqual(rows(result, "compute temperature?"), [expected])
                self.assertEqual(result, render_parameter(result, {**PARAMETERS, "envelope_amax_um": 1},
                                                         {"photons_temperature": 2048000, "grains": 50}, stage, wave))

    def test_fraction_units_are_explicit_and_never_implicitly_converted(self):
        for fraction in (0, .04, .999):
            validate_parameters({"envelope_ice_mass_fraction": fraction})
        for fraction in (-.01, 1, float("nan"), True, "0.04"):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                validate_parameters({"envelope_ice_mass_fraction": fraction})
        with self.assertRaisesRegex(ValueError, "not both"):
            validate_parameters({"envelope_ice_mass_fraction": .04, "envelope_ice_volume_fraction": .1})
        with self.assertRaisesRegex(ValueError, "explicit envelope_ice_mass_fraction"):
            render_parameter(TEMPLATE, {"envelope_ice_volume_fraction": .04}, {}, "temperature")
        with self.assertRaisesRegex(ValueError, "not an ice mass fraction"):
            render_parameter(MIE_TEMPLATE, {"envelope_ice_mass_fraction": .04}, {}, "temperature")
        with self.assertRaisesRegex(ValueError, "positive ice mass fraction"):
            render_parameter(TEMPLATE, {"envelope_ice_mass_fraction": 0}, {}, "temperature")

    def test_ice_population_numerics_are_fixed(self):
        with self.assertRaisesRegex(ValueError, "50 grain bins"):
            render_parameter(TEMPLATE, PARAMETERS, {"grains": 100}, "temperature")

    def test_mutated_dhs_layout_is_rejected(self):
        changes = [
            ("DHS  1 1  0.0  0.96  0.1", "DHS  2 1  0.0  0.96  0.1"),
            ("DHS  1 1  0.0  0.04  0.1", "Mie  1 1  0.0  0.04  0.1"),
            ("DHS  1 1  0.0  0.04  0.1", "DHS  1 1  0.1  0.04  0.1"),
            ("DHS  1 1  0.0  0.04  0.1", "DHS  1 1  0.0  0.04  0.8"),
            ("DHS  1 1  0.0  0.04  0.1", "DHS  1 1  0.0  0.4  0.1"),
            ("DHS  1 1  0.0  0.04  0.1", "DHS  1 1  0.0  nan  0.1"),
            ("H2O_30K_Leiden_mcfost.dat  1.0", "ice_opct.dat  1.0"),
            ("H2O_30K_Leiden_mcfost.dat  1.0", "H2O_30K_Leiden_mcfost.dat  0.9"),
            ("0.03  0.4  2.75 50", "0.03  1.0  2.75 50"),
            ("0.03  0.4  2.75 50", "0.03  0.4  3.5 50"),
            ("0.03  0.4  2.75 50", "0.03  0.4  2.75 100"),
        ]
        for before, after in changes:
            with self.subTest(after=after), self.assertRaises(ValueError):
                render_parameter(TEMPLATE.replace(before, after), PARAMETERS, {}, "temperature")


if __name__ == "__main__":
    unittest.main()
