"""Rendering tests use preserved scientific inputs but never launch MCFOST."""
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mcfost_grid.physics import numerical_args, physical_args, render_parameter, validate_parameters


TEMPLATE = (ROOT / "reference/parameters/continuum_nominal.para").read_text()


def field(text, tag):
    rows = [line for line in text.splitlines() if tag in line]
    if len(rows) != 1:
        raise AssertionError((tag, rows))
    return rows[0].split(tag)[0].split()


class PhysicsTests(unittest.TestCase):
    def test_temperature_preserves_original_physics(self):
        result = render_parameter(TEMPLATE, {}, {}, "temperature")
        self.assertEqual(field(result, "compute temperature?"), ["T", "F", "T"])
        self.assertEqual(field(result, "nbr_photons_eq_th"), ["128000"])
        self.assertIn("0.000225", result)
        self.assertIn("0.03  0.4  2.75  50", result)

    def test_stages_and_full_polarized_output(self):
        for stage, flags in [("image", ["F", "F", "T"]), ("coeval", ["F", "T", "F"])]:
            result = render_parameter(TEMPLATE, {}, {}, stage, wavelength_um=3.0)
            self.assertEqual(field(result, "compute temperature?"), flags)
            self.assertEqual(field(result, "separation of different contributions?"), ["T", "T"])
            self.assertEqual(field(result, "wavelength file ("), ["wavelength.lambda"])
            self.assertEqual(field(result, "image symmetry"), ["F"])

    def test_all_numeric_overrides(self):
        numerics = dict(photons_temperature=500, photons_image=200, photons_sed=300,
                        grid_nr=20, grid_ntheta=12, grid_n_inner=5, image_npix=31,
                        image_size_au=1000, grains=4)
        result = render_parameter(TEMPLATE, {}, numerics, "temperature")
        self.assertEqual(field(result, "nbr_photons_eq_th"), ["500"])
        self.assertEqual(field(result, "nbr_photons_lambda"), ["300"])
        self.assertEqual(field(result, "nbr_photons_image"), ["200"])
        self.assertEqual(field(result, "n_rad (log distribution)"), ["20", "12", "1", "5"])
        self.assertEqual(field(result, "grid (nx,ny), size [AU]"), ["31", "31", "1000"])
        sizes = [line.split("amin,")[0].split() for line in result.splitlines() if "amin," in line]
        self.assertEqual([x[-1] for x in sizes], ["4", "4", "4"])

    def test_physical_overrides_only_change_named_axes(self):
        p = dict(envelope_dust_mass_msun=.0001, envelope_amax_um=2, envelope_size_exponent=3.2,
                 inclination_deg=55, distance_pc=141.8, stellar_temperature_k=4250,
                 stellar_radius_rsun=3, stellar_mass_msun=.7, envelope_ice_volume_fraction=.2,
                 cavity_half_opening_deg=20, accretion_rate_msun_per_year=1e-7)
        result = render_parameter(TEMPLATE, p, {}, "image", 3.0)
        self.assertEqual(field(result, "RT: imin, imax, n_incl, centered ?"), ["55", "55", "1", "F"])
        self.assertEqual(field(result, "distance (pc)"), ["141.8"])
        self.assertEqual(field(result, "Temp, radius (solar radius),M (solar mass)"), ["4250", "3", "0.7", "0.0", "0.0", "0.0", "T"])
        densities = [line.split("dust mass,")[0].split() for line in result.splitlines() if "dust mass," in line]
        self.assertEqual(densities[0], ["1.e-4", "100."])
        self.assertEqual(densities[1], ["0.0001", "100.0"])
        self.assertIn("0.03  2  3.2  50", result)
        self.assertIn("Draine_Si_sUV.dat  0.8", result)
        self.assertIn("ice_opct.dat  0.2", result)
        # q is not the spatial density exponent. No CLI physics is inserted in .para.
        self.assertIn("-1.5  0.0", result)
        self.assertNotIn("-cavity", result)
        self.assertNotIn("-Mdot", result)

    def test_cavity_and_accretion_arguments(self):
        args = physical_args(dict(cavity_half_opening_deg=17.5, accretion_rate_msun_per_year=7.3011756855e-8))
        self.assertEqual(args[:1], ["-cavity"])
        self.assertAlmostEqual(float(args[1]), 100 / math.tan(math.radians(17.5)), places=11)
        self.assertEqual(args[2:6], ["100", "1", "-Mdot", "1"])
        self.assertEqual(float(args[6]), 7.3011756855e-8)
        self.assertNotIn("-mol", args)
        self.assertNotIn("-rt2", args)

    def test_explicit_cli_parameters_are_required(self):
        with self.assertRaisesRegex(ValueError, "Explicit CLI physics"):
            physical_args({})

    def test_zero_accretion_is_explicit(self):
        self.assertEqual(physical_args(dict(cavity_half_opening_deg=20, accretion_rate_msun_per_year=0))[-1], "0")

    def test_optional_seed_is_a_cli_numeric_not_a_parameter_row(self):
        self.assertEqual(numerical_args({}), [])
        for seed in (1, 41001, 2147483647):
            with self.subTest(seed=seed):
                self.assertEqual(numerical_args({"random_seed": seed}), ["-seed", str(seed)])
                self.assertEqual(render_parameter(TEMPLATE, {}, {"random_seed": seed}, "temperature"),
                                 render_parameter(TEMPLATE, {}, {}, "temperature"))

    def test_seed_requires_a_positive_31_bit_json_integer(self):
        for seed in (None, False, True, 0, -1, 1.5, 41001.0, "41001", 2147483648, float("inf")):
            with self.subTest(seed=seed):
                with self.assertRaisesRegex(ValueError, "random_seed"):
                    numerical_args({"random_seed": seed})
                with self.assertRaisesRegex(ValueError, "random_seed"):
                    render_parameter(TEMPLATE, {}, {"random_seed": seed}, "temperature")

    def test_nonfinite_invalid_and_unknown_values(self):
        cases = [dict(distance_pc=0), dict(envelope_amax_um=-1), dict(inclination_deg=91),
                 dict(envelope_size_exponent=-1), dict(cavity_half_opening_deg=0),
                 dict(cavity_half_opening_deg=90), dict(envelope_ice_volume_fraction=-.01),
                 dict(envelope_ice_volume_fraction=1), dict(stellar_mass_msun=float("nan")),
                 dict(distance_pc=float("inf")), dict(distance_pc=True), dict(distance_pc="147"),
                 dict(accretion_rate_msun_per_year=-1), dict(dust_family="DHS")]
        for p in cases:
            with self.subTest(parameters=p), self.assertRaises(ValueError):
                validate_parameters(p)

    def test_amax_cannot_cross_template_amin(self):
        with self.assertRaisesRegex(ValueError, "minimum grain size"):
            render_parameter(TEMPLATE, dict(envelope_amax_um=.03), {}, "temperature")

    def test_numeric_and_stage_validation(self):
        for numerics in [dict(photons_image=0), dict(grains=1), dict(grid_nr=1),
                         dict(image_npix=10.5), dict(image_size_au=-1), dict(unknown=10)]:
            with self.subTest(numerics=numerics), self.assertRaises(ValueError):
                render_parameter(TEMPLATE, {}, numerics, "temperature")
        for stage, wave in [("image", None), ("coeval", 0), ("temperature", 3.0), ("SED", None)]:
            with self.subTest(stage=stage, wave=wave), self.assertRaises(ValueError):
                render_parameter(TEMPLATE, {}, {}, stage, wave)

    def test_dhs_template_is_rejected(self):
        template = (ROOT / "reference/parameters/ice_v02_dust.para").read_text()
        with self.assertRaisesRegex(ValueError, "DHS/disjoint"):
            render_parameter(template, {}, {}, "temperature")

    def test_duplicate_and_missing_labels_fail_closed(self):
        for broken in [TEMPLATE.replace("#Scattering method", "1 distance (pc)\n#Scattering method"),
                       TEMPLATE.replace("nbr_photons_eq_th", "unlabelled_temperature_photons")]:
            with self.assertRaisesRegex(ValueError, "exactly one"):
                render_parameter(broken, dict(distance_pc=147), {}, "temperature")

    def test_blank_lines_do_not_shift_edits(self):
        padded = TEMPLATE.replace("#Density structure", "\n\n#Density structure\n\n")
        result = render_parameter(padded, dict(envelope_dust_mass_msun=.0003), {}, "temperature")
        self.assertIn("0.0003  100.0", result)

    def test_unsupported_geometry_and_phase_function(self):
        for before, after in [("2                       1 = cylindrical", "1                       1 = cylindrical"),
                              ("1                       1=exact phase function", "2                       1=exact phase function")]:
            with self.assertRaises(ValueError):
                render_parameter(TEMPLATE.replace(before, after), {}, {}, "temperature")

    def test_changing_ice_does_not_edit_disk_silicate(self):
        rendered = render_parameter(TEMPLATE, dict(envelope_ice_volume_fraction=.1), {}, "temperature")
        self.assertIn("Draine_Si_sUV.dat  1.0  Optical", rendered)
        self.assertIn("Draine_Si_sUV.dat  0.9", rendered)

    def test_zero_ice_is_a_true_bare_species_and_preserves_disk(self):
        for stage in ("temperature", "image", "coeval"):
            wave = None if stage == "temperature" else 3.0
            rendered = render_parameter(TEMPLATE, dict(envelope_ice_volume_fraction=0), {}, stage, wave)
            self.assertNotIn("ice_opct.dat", rendered)
            types = [line.split("Grain type")[0].split() for line in rendered.splitlines() if "Grain type" in line]
            self.assertEqual(types[-1][:5], ["Mie", "1", "1", "0.0", "1.0"])
            self.assertIn("ac_opct.dat", rendered)
            self.assertIn("Draine_Si_sUV.dat  1", rendered)
            rerendered = render_parameter(rendered, dict(envelope_ice_volume_fraction=0), {}, stage, wave)
            self.assertEqual(rendered, rerendered)

    def test_tiny_mantle_remains_coated_and_bare_cannot_gain_implicit_mantle(self):
        tiny = render_parameter(TEMPLATE, dict(envelope_ice_volume_fraction=1e-6), {}, "temperature")
        self.assertIn("ice_opct.dat  1e-06", tiny)
        bare = render_parameter(TEMPLATE, dict(envelope_ice_volume_fraction=0), {}, "temperature")
        with self.assertRaisesRegex(ValueError, "Bare template"):
            render_parameter(bare, dict(envelope_ice_volume_fraction=.1), {}, "temperature")


if __name__ == "__main__":
    unittest.main()
