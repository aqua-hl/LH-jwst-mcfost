"""V3 pure-species outputs must reproduce the direct V2 opacity reference."""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import silicate_search_v3_stage0 as stage0

SOURCE = ROOT / "runs/silicate_search_v2_512k"


def reference_values():
    # Independently transcribed test fixtures from the direct V2 Stage-0
    # samples. Tests do not require the untracked downloaded results folder.
    # Production references are still recovered and frozen by the V3 builder.
    triples = {"draine": (1485.7525634765625, 3032.604632315978, 1114.771934215509),
               "pyroxene_mg50": (1053.412841796875, 2829.2201540675846, 1397.0458975448782),
               "carbon": (10808.2900390625, 675.7331910524235, 285.55590608214237),
               "ice": (193.22640991210938, 0., 417.5813248029947)}
    return {species: dict(zip(("ext_2p2", "abs_9p7", "abs_18"), values))
            for species, values in triples.items()}


def executable(path):
    path.write_text(f"#!{sys.executable}\n" + '''
import json, os, sys
from pathlib import Path
import numpy as np
from astropy.io import fits
args = sys.argv[1:]
assert args[0] == 'check.para' and '-dust_prop' in args
assert '-root_dir' not in args and '-img' not in args
assert os.environ['OMP_NUM_THREADS'] == '2'
root = Path.cwd().parents[2]
species = Path.cwd().parent.name
reference = json.loads((root/'inputs/v2_reference.json').read_text())['species_opacities_cm2_g'][species]
wave = np.loadtxt('check.lambda')  # MCFOST consumes every numeric line, no header.
assert len(wave) == 99
kappa, albedo = np.full(len(wave), 20.), np.full(len(wave), .25)
for key, wavelength in (('ext_2p2',2.2),('abs_9p7',9.7),('abs_18',18.)):
    found, = np.flatnonzero(np.isclose(wave, wavelength, rtol=2e-6, atol=0))
    target = reference[key]
    if key == 'ext_2p2':
        kappa[found] = target
    elif target:
        kappa[found], albedo[found] = target, 0.
    else:
        kappa[found], albedo[found] = 1., 1.
if species == 'ice' and os.environ.get('V3_FAKE_MODE') == 'drift':
    kappa[np.isclose(wave,18.,rtol=2e-6,atol=0)] *= 1.000002
if species == 'ice' and os.environ.get('V3_FAKE_MODE') == 'nonzero_ice':
    albedo[np.isclose(wave,9.7,rtol=2e-6,atol=0)] = 1-1e-12
Path('data_dust').mkdir()
print('Computation of dust properties ... Writing dust properties')
if os.environ.get('V3_FAKE_MODE') == 'fitsio':
    print('FITSIO Error Status = 105')
    print('Exiting')
    sys.exit(0)
for name,array in {'lambda':wave.astype('f4'),'kappa':kappa,'albedo':albedo,
                   'kappa_grain':np.ones((50,len(wave)))}.items():
    fits.PrimaryHDU(array).writeto('data_dust/'+name+'.fits.gz')
print('Exiting')
''')
    path.chmod(0o755)


def fixture(root):
    inputs = root / "inputs"
    dust = inputs / "utils/Dust"
    dust.mkdir(parents=True)
    shutil.copy2(SOURCE / "inputs/template.para", inputs)
    shutil.copy2(SOURCE / "inputs/stage0_wavelengths.json", inputs)
    for species in stage0.species_definitions().values():
        shutil.copy2(SOURCE / "inputs/utils/Dust" / species["filename"], dust)
    (inputs / "v2_reference.json").write_text(json.dumps({"species_opacities_cm2_g": reference_values()}))
    utils = root / "system_utils"
    for name in ("Dust", "Lambda", "Stellar_Spectra"):
        (utils / name).mkdir(parents=True)
    (utils / "Dust/test.dat").write_text("unchanged mock utility")
    fake = root / "mcfost"
    executable(fake)
    machine = root / "machine.json"
    machine.write_text(json.dumps(dict(schema_version=1, mcfost_executable=str(fake),
        mcfost_utils=str(utils), threads=64, backend="image_method2")))
    return machine


class ReproductionTests(unittest.TestCase):
    def test_relative_tolerance_and_exact_zero_do_not_floor_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "inputs").mkdir()
            reference = reference_values()
            (root / "inputs/v2_reference.json").write_text(json.dumps({"species_opacities_cm2_g": reference}))
            opacity = dict(wavelength_um=[2.2, 9.7, 18.],
                           kappa_ext_cm2_g=[reference["ice"]["ext_2p2"], 1., 1.],
                           kappa_abs_cm2_g=[1., 0., reference["ice"]["abs_18"]])
            check = stage0.reference_comparison(root, "ice", opacity)
            self.assertTrue(check["passed"])
            zero = check["comparisons"]["abs_9p7"]
            self.assertEqual(zero["policy"], "exact_zero")
            self.assertEqual(zero["absolute_tolerance_cm2_g"], 0.)
            self.assertIsNone(zero["relative_difference"])
            opacity["kappa_ext_cm2_g"][0] *= 1 + .9e-6
            self.assertTrue(stage0.reference_comparison(root, "ice", opacity)["passed"])
            opacity["kappa_ext_cm2_g"][0] = reference["ice"]["ext_2p2"] * (1 + 1.1e-6)
            self.assertFalse(stage0.reference_comparison(root, "ice", opacity)["passed"])
            opacity["kappa_ext_cm2_g"][0] = reference["ice"]["ext_2p2"]
            opacity["kappa_abs_cm2_g"][1] = 1e-30
            self.assertFalse(stage0.reference_comparison(root, "ice", opacity)["passed"])
            self.assertEqual(opacity["kappa_abs_cm2_g"][1], 1e-30)

    def test_four_actual_processes_headerless_table_default_root_and_passed_cache(self):
        with tempfile.TemporaryDirectory(prefix="v3-stage0-") as tmp:
            root = Path(tmp)
            machine = fixture(root)
            result = stage0.run_stage0(root, machine)
            self.assertEqual(result["status"], "passed")
            self.assertEqual(result["species"], 4)
            self.assertEqual(result["wavelengths"], 99)
            receipt = stage0.validate_stage0_receipt(root)
            self.assertTrue(receipt["v2_reproducibility_checked"])
            self.assertIn("inputs/v2_reference.json", receipt["identity"]["input_sha256"])
            self.assertEqual([row["id"] for row in receipt["species"]], ["draine", "pyroxene_mg50", "carbon", "ice"])
            for row in receipt["species"]:
                self.assertTrue(row["v2_reproducibility"]["passed"])
                self.assertEqual(len(row["v2_reproducibility"]["comparisons"]), 3)
                self.assertNotIn("-root_dir", row["command"])
                self.assertEqual(len(np.loadtxt(root / row["attempt"] / "check.lambda")), 99)
            runtime = stage0.load_runner(root).runtime_config(machine)
            stage0.validate_stage0_receipt(root, runtime)
            self.assertEqual(stage0.run_stage0(root, machine)["status"], "cached")
            self.assertEqual(len(list((root / "stage0").glob("*/attempt_*"))), 4)

    def test_disagreement_stops_preparation_and_preserves_failed_attempt(self):
        with tempfile.TemporaryDirectory(prefix="v3-stage0-") as tmp:
            root = Path(tmp)
            machine = fixture(root)
            with patch.dict(os.environ, {"V3_FAKE_MODE": "drift"}):
                with self.assertRaisesRegex(ValueError, "V2 reproducibility failed for ice"):
                    stage0.run_stage0(root, machine)
            receipt = json.loads((root / stage0.RECEIPT).read_text())
            self.assertEqual(receipt["status"], "failed")
            self.assertFalse(receipt["v2_reproducibility_checked"])
            failed = receipt["species"][-1]["v2_reproducibility"]["comparisons"]["abs_18"]
            self.assertFalse(failed["passed"])
            self.assertGreater(failed["relative_difference"], 1e-6)
            path = root / "stage0/ice/attempt_001/data_dust/kappa.fits.gz"
            original = path.read_bytes()
            with self.assertRaisesRegex(ValueError, "not passed"):
                stage0.validate_stage0_receipt(root)
            self.assertEqual(stage0.run_stage0(root, machine)["status"], "passed")
            self.assertEqual(path.read_bytes(), original)
            self.assertEqual(len(list((root / "stage0").glob("*/attempt_*"))), 8)

    def test_nonzero_ice_absorption_and_exit_zero_fitsio_both_block(self):
        for mode, message in (("nonzero_ice", "V2 reproducibility failed for ice"),
                              ("fitsio", "FITSIO error despite exit 0")):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(prefix="v3-stage0-") as tmp:
                root = Path(tmp)
                machine = fixture(root)
                with patch.dict(os.environ, {"V3_FAKE_MODE": mode}):
                    with self.assertRaisesRegex(ValueError, message):
                        stage0.run_stage0(root, machine)
                self.assertEqual(json.loads((root / stage0.RECEIPT).read_text())["status"], "failed")

    def test_cached_comparison_and_reference_hash_cannot_be_silently_changed(self):
        with tempfile.TemporaryDirectory(prefix="v3-stage0-") as tmp:
            root = Path(tmp)
            machine = fixture(root)
            stage0.run_stage0(root, machine)
            path = root / stage0.RECEIPT
            original = path.read_text()
            receipt = json.loads(original)
            receipt["species"][0]["v2_reproducibility"]["comparisons"]["ext_2p2"]["actual_cm2_g"] += 1
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "reproducibility is absent, failed or changed"):
                stage0.validate_stage0_receipt(root)
            path.write_text(original)
            reference = root / "inputs/v2_reference.json"
            reference.write_text(reference.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "inputs, code or runtime changed"):
                stage0.run_stage0(root, machine)


if __name__ == "__main__":
    unittest.main()
