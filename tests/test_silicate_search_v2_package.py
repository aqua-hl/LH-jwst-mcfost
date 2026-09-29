"""End-to-end portable package check using synthetic FITS only in a temp directory."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "runs/silicate_structure_production_v1_512k"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fake_runtime(root):
    """No real MCFOST calls or physical opacities; five arbitrary test scales."""
    executable = root / "synthetic_mcfost_test_only"
    executable.write_text(f"#!{sys.executable}\n" + '''
from pathlib import Path
import os
import sys
import numpy as np
from astropy.io import fits
assert sys.argv[1] == "check.para" and "-dust_prop" in sys.argv
assert "-img" not in sys.argv and "-root_dir" not in sys.argv
assert os.environ["OMP_NUM_THREADS"] == "2"
parameter = Path("check.para").read_text()
scales = {"Draine_Si_sUV.dat": 20., "Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat": 30.,
          "Olivine_MgFeSiO4_Dorschner1995_mcfost.dat": 40., "ac_opct.dat": 80.,
          "H2O_30K_Leiden_mcfost.dat": 5.}
selected = [value for name, value in scales.items() if name in parameter]
assert len(selected) == 1
wave = np.atleast_1d(np.loadtxt("check.lambda"))
assert wave.size == 99
folder = Path("data_dust")
folder.mkdir()
for name, array in {"lambda": wave.astype("f4"), "kappa": np.full(wave.size, selected[0]),
                    "albedo": np.full(wave.size, .25),
                    "kappa_grain": np.full((50, wave.size), selected[0]*.75)}.items():
    fits.PrimaryHDU(array).writeto(folder/(name+".fits.gz"))
print("Synthetic test executable: Writing dust properties")
print("Exiting")
''')
    executable.chmod(0o755)
    utils = root / "synthetic_system_utils"
    for folder in ("Dust", "Lambda", "Stellar_Spectra"):
        (utils / folder).mkdir(parents=True)
    (utils / "Dust/test_only.dat").write_text("Arbitrary utility bytes for a test runtime identity.\n")
    machine = root / "synthetic_machine.json"
    machine.write_text(json.dumps(dict(schema_version=1, mcfost_executable=str(executable),
                                     mcfost_utils=str(utils), threads=64, backend="image_method2",
                                     timeout_seconds=14400, max_memory_gb=112)))
    return machine


class SilicateSearchV2PackageTests(unittest.TestCase):
    def invoke(self, args, *, cwd, success=True):
        env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "OPENBLAS_NUM_THREADS": "1",
               "MPLCONFIGDIR": str(cwd / ".test-mplconfig")}
        result = subprocess.run([sys.executable, "-B", *map(str, args)], cwd=cwd,
                                capture_output=True, text=True, env=env, timeout=120)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_bootstrap_stage0_materialization_cache_and_tamper_gates(self):
        previous_manifest = json.loads((SOURCE / "manifest.json").read_text())
        previous_model = previous_manifest["models"][0]["id"]
        inherited_paths = [SOURCE / "manifest.json", SOURCE / "experiment.json",
                           SOURCE / "inputs/template.para", SOURCE / "models" / previous_model / "temperature.para"]
        inherited_paths.extend(sorted((SOURCE / "inputs/utils/Dust").glob("*.dat")))
        original = {path: sha256(path) for path in inherited_paths}
        with tempfile.TemporaryDirectory(prefix="silicate-v2-synthetic-package-") as tmp:
            root = Path(tmp).resolve()
            bundle = root / "bundle"
            prepared = self.invoke([ROOT / "scripts/build_silicate_search_v2.py", "--output", bundle], cwd=root)
            self.assertFalse(json.loads(prepared.stdout)["mass_estimates_used"])
            bootstrap = json.loads((bundle / "bootstrap.json").read_text())
            self.assertEqual(bootstrap["model_count"], 34)
            self.assertEqual(bootstrap["comparison_count"], 48)
            self.assertFalse(bootstrap["matched_masses_available"])
            self.assertFalse(bootstrap["production_parameters_rendered"])
            self.assertTrue(all("parameters" not in row and "opacity_matching" not in row for row in bootstrap["cases"]))
            self.assertFalse((bundle / "production").exists())
            self.assertFalse((bundle / "models").exists())
            wave_plan = json.loads((bundle / "inputs/stage0_wavelengths.json").read_text())
            self.assertEqual(len(wave_plan["wavelengths_um"]), 99)
            self.assertEqual(len(wave_plan["production_wavelengths_um"]), 98)
            self.assertEqual(wave_plan["wavelengths_um"], sorted(set(wave_plan["production_wavelengths_um"] + [2.2, 9.7, 18.])))
            self.assertFalse(wave_plan["count_header"])
            for path in inherited_paths:
                relative = path.relative_to(SOURCE)
                if relative.parts[:2] == ("inputs", "utils") or relative == Path("inputs/template.para"):
                    self.assertEqual(sha256(bundle / relative), original[path])

            builder = bundle / "code/build_silicate_search_v2.py"
            self.invoke([builder, "--validate", bundle], cwd=root)
            self.invoke([builder, "--materialize", bundle], cwd=root, success=False)
            self.assertFalse((bundle / "production").exists())

            machine = fake_runtime(root)
            result = self.invoke([bundle / "code/silicate_search_v2_stage0.py", bundle, "--machine", machine], cwd=root)
            self.assertEqual(json.loads(result.stdout)["status"], "passed")
            receipt_path = bundle / "stage0/receipt.json"
            receipt = json.loads(receipt_path.read_text())
            self.assertEqual(len(receipt["species"]), 5)
            self.assertEqual([row["status"] for row in receipt["species"]], ["passed"] * 5)
            self.assertFalse((bundle / "production").exists())

            opacity_path = bundle / receipt["species"][0]["opacity_file"]
            opacity_original = opacity_path.read_bytes()
            opacity_path.write_bytes(opacity_original + b"\n")
            failure = self.invoke([builder, "--materialize", bundle], cwd=root, success=False)
            self.assertIn("artifact changed", failure.stderr)
            self.assertFalse((bundle / "production").exists())
            opacity_path.write_bytes(opacity_original)

            result = self.invoke([builder, "--materialize", bundle], cwd=root)
            self.assertEqual(json.loads(result.stdout)["status"], "prepared")
            production = bundle / "production"
            manifest = json.loads((production / "manifest.json").read_text())
            experiment = json.loads((production / "experiment.json").read_text())
            self.assertEqual(len(manifest["models"]), 34)
            self.assertEqual(len(manifest["anchors"]), 98)
            self.assertEqual(len(experiment["predeclared_contrasts"]), 48)
            self.assertEqual(experiment["stage0_receipt_sha256"], sha256(receipt_path))
            self.assertEqual(experiment["image_requests"], 3332)
            self.assertEqual(len(list((production / "models").glob("*/temperature.para"))), 34)
            self.assertEqual(len(list((production / "models").glob("*/anchors/*/image.para"))), 3332)
            self.assertEqual(len(list((production / "models").glob("*/anchors/*/coeval.para"))), 3332)
            carbon_case = next(r for r in experiment["catalogue"] if r["carbon_mass_fraction"] == .2 and r["ice_mass_fraction"] == .08 and r["density_exponent"] == -1.5)
            expected_mass = 2.25e-4 * (.96*20+.04*5)/(.72*20+.2*80+.08*5)
            self.assertAlmostEqual(carbon_case["parameters"]["envelope_dust_mass_msun"] / expected_mass, 1., delta=6e-13)
            validate = ("import json,sys;from pathlib import Path;"
                        "sys.path.insert(0,str(Path(sys.argv[1])/'code'));"
                        "from production_task import validate_package;"
                        "print(json.dumps({'models':len(validate_package(Path(sys.argv[1]))['catalogue'])}))")
            valid = self.invoke(["-c", validate, production], cwd=root)
            self.assertEqual(json.loads(valid.stdout)["models"], 34)
            immutable_hash = sha256(production / "manifest.json")
            result = self.invoke([builder, "--materialize", bundle], cwd=root)
            self.assertEqual(json.loads(result.stdout)["status"], "cached")
            self.assertEqual(sha256(production / "manifest.json"), immutable_hash)
            self.assertEqual(len(list((bundle / "build_attempts").glob("build-*"))), 1)

            receipt_original = receipt_path.read_bytes()
            updated_receipt = json.loads(receipt_original)
            updated_receipt["test_annotation"] = "A changed receipt cannot silently reuse a frozen catalogue."
            receipt_path.write_text(json.dumps(updated_receipt))
            mismatch = self.invoke([builder, "--materialize", bundle], cwd=root, success=False)
            self.assertIn("different Stage-0", mismatch.stderr)
            receipt_path.write_bytes(receipt_original)

            temperature = production / "models" / manifest["models"][0]["id"] / "temperature.para"
            temperature_original = temperature.read_bytes()
            temperature.write_bytes(temperature_original + b"\n")
            mismatch = self.invoke([builder, "--materialize", bundle], cwd=root, success=False)
            self.assertIn("Frozen input changed", mismatch.stderr)
            temperature.write_bytes(temperature_original)
            self.invoke([builder, "--validate", bundle], cwd=root)
            self.assertEqual(sha256(production / "manifest.json"), immutable_hash)

        self.assertEqual({path: sha256(path) for path in inherited_paths}, original)


if __name__ == "__main__":
    unittest.main()
