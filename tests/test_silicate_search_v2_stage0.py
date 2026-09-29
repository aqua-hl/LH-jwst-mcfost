"""Stage 0 covers five real simulator invocations before any mass catalogue."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import silicate_search_v2_stage0 as stage0

SOURCE = ROOT / "runs/silicate_structure_production_v1_512k"


def products(root, wave, *, opacity=None, albedo=None):
    folder = Path(root) / "data_dust"
    folder.mkdir(exist_ok=True)
    arrays = {"lambda": wave, "kappa": np.ones(len(wave)) if opacity is None else opacity,
              "albedo": np.full(len(wave), .25) if albedo is None else albedo,
              "kappa_grain": np.ones((50, len(wave)))}
    for name, data in arrays.items():
        fits.PrimaryHDU(np.asarray(data)).writeto(folder / (name + ".fits.gz"), overwrite=True)


def fake_executable(path):
    """Implement the upstream lambda parser and default output-root contract.

    The child reads every numeric line. It never skips a count header, creates
    the output directory under -root_dir, and always writes under cwd, just
    like the source version that caused the earlier preflight failures.
    """
    path.write_text(f"#!{sys.executable}\n" + '''
import os
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
args = sys.argv[1:]
assert args[0] == "check.para" and "-dust_prop" in args
assert "-img" not in args and os.environ["OMP_NUM_THREADS"] == "2"
output_root = args[args.index("-root_dir") + 1] if "-root_dir" in args else "."
(Path(output_root) / "data_dust").mkdir()
print("Computation of dust properties ... Writing dust properties", flush=True)
if os.environ.get("STAGE0_FAKE_MODE") == "fitsio":
    print("FITSIO Error Status = 105: couldn't create the named file")
    print("Exiting")
    sys.exit(0)
wave = np.loadtxt("check.lambda")
data = np.ones(len(wave))
if (os.environ.get("STAGE0_FAKE_MODE") == "ice_nan"
        and "H2O_30K_Leiden_mcfost.dat" in Path("check.para").read_text()):
    data[0] = np.nan
try:
    for name, array in {"lambda": wave.astype("f4"), "kappa": data * 20,
                        "albedo": np.full(len(wave), 0.25),
                        "kappa_grain": np.ones((50, len(wave)))}.items():
        fits.PrimaryHDU(array).writeto("data_dust/" + name + ".fits.gz")
except OSError:
    print("FITSIO Error Status = 105: couldn't create the named file")
print("Exiting")
''')
    path.chmod(0o755)


def fixture(root):
    inputs = root / "inputs"
    dust = inputs / "utils/Dust"
    dust.mkdir(parents=True)
    shutil.copy2(SOURCE / "inputs/template.para", inputs / "template.para")
    definitions = stage0.species_definitions()
    for row in definitions.values():
        shutil.copy2(SOURCE / "inputs/utils/Dust" / row["filename"], dust)
    anchors = json.loads((SOURCE / "manifest.json").read_text())["anchors"]
    production = [row["wavelength_um"] for row in anchors]
    wavelengths = sorted(set(production + [2.2, 9.7, 18.]))
    (inputs / "stage0_wavelengths.json").write_text(json.dumps(
        dict(wavelengths_um=wavelengths, production_wavelengths_um=production)))
    utils = root / "system_utils"
    for folder in ("Dust", "Lambda", "Stellar_Spectra"):
        (utils / folder).mkdir(parents=True)
    (utils / "Dust/fixed.dat").write_text("test runtime utilities")
    executable = root / "fake_mcfost"
    fake_executable(executable)
    machine = root / "machine.json"
    machine.write_text(json.dumps(dict(schema_version=1, mcfost_executable=str(executable),
                                      mcfost_utils=str(utils), threads=64, backend="image_method2")))
    return machine, wavelengths


class OutputTests(unittest.TestCase):
    def test_extinction_is_split_into_absorption_and_scattering_per_gram(self):
        wave = np.array([2.2, 9.7, 18.])
        with tempfile.TemporaryDirectory() as tmp:
            products(tmp, wave, opacity=np.array([100., 200., 300.]), albedo=np.array([1., .25, 0.]))
            result = stage0.check_outputs(tmp, wave)
            self.assertEqual(result["kappa_abs_cm2_g"], [0., 150., 300.])
            self.assertEqual(result["kappa_sca_cm2_g"], [100., 50., 0.])
            self.assertEqual(result["units"], "cm^2 per g of dust")
            self.assertEqual(result["grain_absorption_shape"], [50, 3])
            # An undefined polarizability in a zero-absorbing population does
            # not enter the extinction matching or absorption screen.
            fits.PrimaryHDU(np.full((181, 3), np.nan)).writeto(Path(tmp) / "data_dust/polarizability.fits.gz")
            self.assertEqual(stage0.check_outputs(tmp, wave), result)

    def test_every_requested_wavelength_must_be_finite_and_physical(self):
        wave = np.array([2.2, 9.7, 18.])
        with tempfile.TemporaryDirectory() as tmp:
            for opacity, albedo, error in (
                ([1., np.nan, 1.], [0., 0., 0.], "Nonfinite kappa"),
                ([1., 1., -1.], [0., 0., 0.], "Negative kappa"),
                ([1., 1., 1.], [0., 1.01, 0.], "Albedo outside"),
                ([1., 1., 1.], [0., np.inf, 0.], "Nonfinite albedo"),
                ([0., 0., 0.], [0., 0., 0.], "identically zero"),
            ):
                products(tmp, wave, opacity=opacity, albedo=albedo)
                with self.assertRaisesRegex(ValueError, error):
                    stage0.check_outputs(tmp, wave)
            products(tmp, wave)
            data = np.ones((50, 3))
            data[49, 1] = np.nan
            fits.PrimaryHDU(data).writeto(Path(tmp) / "data_dust/kappa_grain.fits.gz", overwrite=True)
            with self.assertRaisesRegex(ValueError, "Nonfinite kappa_grain"):
                stage0.check_outputs(tmp, wave)

    def test_wavelength_axis_and_single_species_grain_count_are_checked(self):
        wave = np.array([2.2, 9.7, 18.])
        with tempfile.TemporaryDirectory() as tmp:
            products(tmp, wave)
            with self.assertRaisesRegex(ValueError, "wavelengths differ"):
                stage0.check_outputs(tmp, [2.2, 9.8, 18.])
            with self.assertRaisesRegex(ValueError, "Unexpected lambda shape"):
                stage0.check_outputs(tmp, [2.2, 9.7])
            fits.PrimaryHDU(np.ones((100, 3))).writeto(Path(tmp) / "data_dust/kappa_grain.fits.gz", overwrite=True)
            with self.assertRaisesRegex(ValueError, "Unexpected kappa_grain shape"):
                stage0.check_outputs(tmp, wave)


class Stage0Tests(unittest.TestCase):
    def test_isolation_keeps_one_exact_pure_DHS_species(self):
        template = (SOURCE / "inputs/template.para").read_text()
        for row in stage0.species_definitions().values():
            result = stage0.isolated_parameter(template, row["filename"])
            section = result.split("#Grain properties")[1].split("#Molecular RT settings")[0]
            rows = [line.split() for line in section.splitlines() if line.strip()]
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows[0][0], "1")
            self.assertEqual(rows[1][:6], ["DHS", "1", "1", "0.0", "1.0", "0.1"])
            self.assertEqual(rows[2][:2], [row["filename"], "1.0"])
            self.assertEqual(rows[4][:4], ["0.03", "0.4", "2.75", "50"])
            self.assertIn("8 8 1 2", result)
            self.assertIn("F T F", result)
            self.assertNotIn("ZONE 2", result)
            self.assertNotIn("Mie  1 1", result)

    def test_five_subprocesses_headerless_lambda_default_root_and_cache(self):
        with tempfile.TemporaryDirectory(prefix="stage0-") as tmp:
            root = Path(tmp)
            machine, wave = fixture(root)
            result = stage0.run_stage0(root, machine)
            self.assertEqual(result["status"], "passed")
            self.assertEqual(result["wavelengths"], 99)
            receipt = stage0.validate_stage0_receipt(root)
            self.assertEqual(len(receipt["species"]), 5)
            for check in receipt["species"]:
                attempt = root / check["attempt"]
                np.testing.assert_array_equal(np.loadtxt(attempt / "check.lambda"), wave)
                self.assertNotIn("-root_dir", check["command"])
                self.assertEqual(json.loads((attempt / "opacity.json").read_text())["kappa_abs_cm2_g"], [15.] * 99)
            self.assertEqual(stage0.run_stage0(root, machine)["status"], "cached")
            self.assertEqual(len(list((root / "stage0").glob("*/attempt_*"))), 5)
            runtime = stage0.load_runner(root).runtime_config(machine)
            stage0.validate_stage0_receipt(root, runtime)
            (root / "system_utils/Dust/fixed.dat").write_text("changed runtime utility")
            with self.assertRaisesRegex(ValueError, "runtime changed"):
                stage0.validate_stage0_receipt(root, runtime)

    def test_last_species_failure_blocks_all_and_preserves_attempts(self):
        with tempfile.TemporaryDirectory(prefix="stage0-") as tmp:
            root = Path(tmp)
            machine, _ = fixture(root)
            with patch.dict(os.environ, {"STAGE0_FAKE_MODE": "ice_nan"}):
                with self.assertRaisesRegex(ValueError, "Nonfinite kappa"):
                    stage0.run_stage0(root, machine)
            receipt = json.loads((root / stage0.RECEIPT).read_text())
            self.assertEqual(receipt["status"], "failed")
            self.assertEqual([row["status"] for row in receipt["species"]], ["passed"] * 4 + ["failed"])
            with self.assertRaisesRegex(ValueError, "not passed"):
                stage0.validate_stage0_receipt(root)
            original = (root / "stage0/ice/attempt_001/data_dust/kappa.fits.gz").read_bytes()
            self.assertEqual(stage0.run_stage0(root, machine)["status"], "passed")
            self.assertEqual((root / "stage0/ice/attempt_001/data_dust/kappa.fits.gz").read_bytes(), original)
            self.assertEqual(len(list((root / "stage0").glob("*/attempt_*"))), 10)

    def test_fitsio_with_exit_zero_does_not_release_catalogue(self):
        with tempfile.TemporaryDirectory(prefix="stage0-") as tmp:
            root = Path(tmp)
            machine, _ = fixture(root)
            with patch.dict(os.environ, {"STAGE0_FAKE_MODE": "fitsio"}):
                with self.assertRaisesRegex(ValueError, "FITSIO error despite exit 0"):
                    stage0.run_stage0(root, machine)
            receipt = json.loads((root / stage0.RECEIPT).read_text())
            self.assertEqual(receipt["status"], "failed")
            self.assertEqual(receipt["species"][0]["returncode"], 0)
            self.assertEqual(len(receipt["species"]), 1)

    def test_receipt_rejects_modified_products_even_without_cluster_runtime(self):
        with tempfile.TemporaryDirectory(prefix="stage0-") as tmp:
            root = Path(tmp)
            machine, _ = fixture(root)
            stage0.run_stage0(root, machine)
            receipt = stage0.validate_stage0_receipt(root)
            opacity = root / receipt["species"][0]["opacity_file"]
            opacity.write_text(opacity.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                stage0.validate_stage0_receipt(root)

    def test_no_stress_samples_or_modified_constants_are_allowed(self):
        with tempfile.TemporaryDirectory(prefix="stage0-") as tmp:
            root = Path(tmp)
            fixture(root)
            stage0.validate_inputs(root)
            path = root / "inputs/stage0_wavelengths.json"
            original = path.read_text()
            plan = json.loads(original)
            plan["wavelengths_um"].append(2.27227)
            path.write_text(json.dumps(plan))
            with self.assertRaisesRegex(ValueError, "exactly 99"):
                stage0.validate_inputs(root)
            path.write_text(original)
            material = root / "inputs/utils/Dust/H2O_30K_Leiden_mcfost.dat"
            material.write_text(material.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "constants changed"):
                stage0.validate_inputs(root)


if __name__ == "__main__":
    unittest.main()
