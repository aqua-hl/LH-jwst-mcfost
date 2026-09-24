"""Column conservation and isolated rendering of the three campaign axes."""
from collections import Counter
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from mcfost_grid.physics import render_parameter, validate_parameters
from silicate_structure_design import (MATERIALS, NUMERICS, catalogue, column_mass_ratio,
                                      power_integral, validate_materials)

TEMPLATE = (ROOT / "reference/parameters/ice_v02_dust.para").read_text()


def changed_lines(before, after):
    return [(a, b) for a, b in zip(before.splitlines(), after.splitlines()) if a != b]


class SilicateStructurePhysicsTests(unittest.TestCase):
    def test_analytic_column_normalization_against_independent_quadrature(self):
        # Integrate in log radius: dr=exp(t)dt. The numerical quadrature does
        # not use the analytic implementation or its power-integral helper.
        def integral(p):
            return quad(lambda t: np.exp((p + 1) * t), 0., np.log(3000.),
                        epsabs=1e-10, epsrel=1e-12)[0]
        for p in (-1.25, -1.5, -1.75):
            ratio = column_mass_ratio(p)
            for reference_mass in (1.e-4, 1.5e-4, 2.25e-4):
                baseline_column = reference_mass * integral(-1.5) / integral(.5)
                new_column = reference_mass * ratio * integral(p) / integral(p + 2)
                self.assertAlmostEqual(new_column / baseline_column, 1., places=13)
        self.assertEqual(column_mass_ratio(-1.5), 1.)
        self.assertGreater(column_mass_ratio(-1.25), 1.)
        self.assertLess(column_mass_ratio(-1.75), 1.)
        self.assertAlmostEqual(power_integral(-1), np.log(3000.), places=13)

    def test_catalogue_is_complete_unique_and_unselected(self):
        records = catalogue()
        self.assertEqual(len(records), 60)
        self.assertEqual(len({json.dumps(r["parameters"], sort_keys=True) for r in records}), 60)
        self.assertEqual(Counter(r["role"] for r in records),
                         dict(mass_size=24, material=18, density_profile=18))
        for record in records:
            normalization = record["column_normalization"]
            self.assertFalse(normalization["discrete_mcfost_column_certified"])
            self.assertEqual(record["parameters"]["envelope_ice_mass_fraction"], .04)
            if record["role"] != "mass_size":
                parent = records[record["parent_index"]]
                self.assertEqual(parent["material_id"], "draine")
                self.assertEqual(parent["reference_mass_msun"], record["reference_mass_msun"])
                self.assertEqual(parent["parameters"]["inclination_deg"], record["parameters"]["inclination_deg"])

    def test_material_change_is_confined_to_envelope_silicate_row(self):
        parameters = catalogue()[6]["parameters"]
        baseline = render_parameter(TEMPLATE, parameters, NUMERICS, "temperature")
        for material in MATERIALS.values():
            chosen = {**parameters, "envelope_silicate_file": material["filename"]}
            modified = render_parameter(TEMPLATE, chosen, NUMERICS, "temperature")
            differences = changed_lines(baseline, modified)
            if material["filename"] != "Draine_Si_sUV.dat":
                self.assertEqual(len(differences), 1)
                self.assertIn("Optical indices file, within-species", differences[0][1])
            self.assertIn("Draine_Si_sUV.dat  1.0  Optical indices file, volume fraction", modified)
            self.assertIn("H2O_30K_Leiden_mcfost.dat  1.0", modified)
            # Resume/re-render must accept a template already containing a
            # selected laboratory table without reinterpreting disk grains.
            self.assertEqual(render_parameter(modified, chosen, NUMERICS, "temperature"), modified)

    def test_profile_changes_only_envelope_mass_and_exponent(self):
        records = catalogue()
        for record in records[42:]:
            reference = records[record["parent_index"]]
            baseline = render_parameter(TEMPLATE, reference["parameters"], NUMERICS, "temperature")
            changed = render_parameter(TEMPLATE, record["parameters"], NUMERICS, "temperature")
            differences = changed_lines(baseline, changed)
            self.assertEqual(len(differences), 2)
            self.assertIn("dust mass,", differences[0][1])
            self.assertIn("surface density exponent", differences[1][1])
            self.assertIn("-1.0  0.0               surface density exponent", changed)
            self.assertEqual(render_parameter(changed, record["parameters"], NUMERICS, "temperature"), changed)

    def test_parameter_paths_and_exponents_fail_closed(self):
        for filename in ("../Draine_Si_sUV.dat", "/tmp/Draine_Si_sUV.dat", "arbitrary.dat", 1, None, True):
            with self.subTest(filename=filename), self.assertRaisesRegex(ValueError, "supported"):
                validate_parameters(dict(envelope_silicate_file=filename))
        for exponent in (-3, 0, 1, float("nan"), float("inf"), True, "-1.5"):
            with self.subTest(exponent=exponent), self.assertRaises(ValueError):
                validate_parameters(dict(envelope_density_exponent=exponent))
        mie = (ROOT / "reference/parameters/continuum_original_dust.para").read_text()
        with self.assertRaisesRegex(ValueError, "only for the disjoint DHS"):
            render_parameter(mie, {"envelope_silicate_file": "Draine_Si_sUV.dat"}, {}, "temperature")

    def test_material_receipt_preserves_density_and_unmodified_tables(self):
        receipt = validate_materials(ROOT / "laboratory_silicates_v1/Dust",
                                     ROOT / "laboratory_silicates_v1/manifest.json")
        self.assertEqual(receipt["materials"]["olivine_mg50"]["density_g_cm3"], 3.71)
        self.assertEqual(receipt["materials"]["pyroxene_mg50"]["density_g_cm3"], 3.2)
        self.assertFalse(receipt["materials"]["olivine_mg50"]["tail_physical_accuracy_certified"])
        metadata = json.loads((ROOT / "laboratory_silicates_v1/manifest.json").read_text())
        metadata["materials"]["olivine_mg50"]["density_g_cm3"] = 3.5
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            path.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, "provenance manifest changed"):
                validate_materials(ROOT / "laboratory_silicates_v1/Dust", path)


if __name__ == "__main__":
    unittest.main()
