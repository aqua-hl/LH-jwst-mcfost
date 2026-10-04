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
SOURCE = ROOT / "runs/silicate_search_v2_512k"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fake_runtime(root, reference):
    """Test-only FITS match the three archived screens; other wavelengths are synthetic."""
    executable = root / "synthetic_mcfost_test_only"
    executable.write_text(f"#!{sys.executable}\nreference_path = {str(reference)!r}\n" + '''
from pathlib import Path
import json
import os
import sys
import numpy as np
from astropy.io import fits
assert sys.argv[1] == "check.para" and "-dust_prop" in sys.argv
assert "-img" not in sys.argv and "-root_dir" not in sys.argv
assert os.environ["OMP_NUM_THREADS"] == "2"
parameter = Path("check.para").read_text()
names = {"Draine_Si_sUV.dat": "draine", "Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat": "pyroxene_mg50",
         "ac_opct.dat": "carbon", "H2O_30K_Leiden_mcfost.dat": "ice"}
selected = [value for name, value in names.items() if name in parameter]
assert len(selected) == 1
wave = np.atleast_1d(np.loadtxt("check.lambda"))
assert wave.size == 99
reference = json.loads(Path(reference_path).read_text())["species_opacities_cm2_g"][selected[0]]
extinction = np.full(wave.size, 1000.)
albedo = np.full(wave.size, .25)
extinction[np.isclose(wave, 2.2, rtol=2e-6)] = reference["ext_2p2"]
for wavelength, key in [(9.7, "abs_9p7"), (18., "abs_18")]:
    location = np.flatnonzero(np.isclose(wave, wavelength, rtol=2e-6))
    assert len(location) == 1
    absorption = reference[key]
    extinction[location[0]] = max(2*absorption, 1.)
    albedo[location[0]] = 1-absorption/extinction[location[0]]
folder = Path("data_dust")
folder.mkdir()
for name, array in {"lambda": wave.astype("f4"), "kappa": extinction,
                    "albedo": albedo,
                    "kappa_grain": np.tile(extinction*(1-albedo), (50, 1))}.items():
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


class SilicateSearchV3PackageTests(unittest.TestCase):
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
        inherited_paths = [SOURCE / "bootstrap.json", SOURCE / "inputs/production_anchors.json",
                           SOURCE / "inputs/template.para"]
        inherited_paths.extend(sorted((SOURCE / "inputs/utils/Dust").glob("*.dat")))
        original = {path: sha256(path) for path in inherited_paths}
        with tempfile.TemporaryDirectory(prefix="silicate-v3-synthetic-package-") as tmp:
            root = Path(tmp).resolve()
            bundle = root / "bundle"
            prepared = self.invoke([ROOT / "scripts/build_silicate_search_v3.py", "--output", bundle], cwd=root)
            self.assertFalse(json.loads(prepared.stdout)["mass_estimates_used"])
            bootstrap = json.loads((bundle / "bootstrap.json").read_text())
            self.assertEqual(bootstrap["model_count"], 26)
            self.assertEqual(bootstrap["comparison_count"], 54)
            self.assertFalse(bootstrap["matched_masses_available"])
            self.assertFalse(bootstrap["derived_ice_fractions_available"])
            self.assertEqual(bootstrap["physical_models"], 24)
            self.assertEqual(bootstrap["seed_replicas"], 2)
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

            builder = bundle / "code/build_silicate_search_v3.py"
            self.invoke([builder, "--validate", bundle], cwd=root)
            self.invoke([builder, "--materialize", bundle], cwd=root, success=False)
            self.assertFalse((bundle / "production").exists())

            machine = fake_runtime(root, bundle / "inputs/v2_reference.json")
            result = self.invoke([bundle / "code/silicate_search_v3_stage0.py", bundle, "--machine", machine], cwd=root)
            self.assertEqual(json.loads(result.stdout)["status"], "passed")
            receipt_path = bundle / "stage0/receipt.json"
            receipt = json.loads(receipt_path.read_text())
            self.assertEqual(len(receipt["species"]), 4)
            self.assertEqual([row["status"] for row in receipt["species"]], ["passed"] * 4)
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
            self.assertEqual(len(manifest["models"]), 26)
            self.assertEqual(len(manifest["anchors"]), 98)
            self.assertEqual(len(experiment["predeclared_contrasts"]), 54)
            self.assertEqual(experiment["stage0_receipt_sha256"], sha256(receipt_path))
            self.assertEqual(experiment["image_requests"], 2548)
            self.assertEqual(len(list((production / "models").glob("*/temperature.para"))), 26)
            self.assertEqual(len(list((production / "models").glob("*/anchors/*/image.para"))), 2548)
            self.assertEqual(len(list((production / "models").glob("*/anchors/*/coeval.para"))), 2548)
            self.assertEqual(len({json.dumps(m["parameters"], sort_keys=True) for m in manifest["models"]}), 24)
            self.assertEqual(len({m["id"] for m in manifest["models"]}), 26)
            reference = json.loads((bundle/"inputs/v2_reference.json").read_text())
            self.assertEqual(experiment["ice_column_targets_msun"], reference["ice_column_targets_msun"])
            for replica, original_index in ((24, 7), (25, 10)):
                models = manifest["models"]
                self.assertEqual(models[replica]["parameters"], models[original_index]["parameters"])
                self.assertNotEqual(models[replica]["id"], models[original_index]["id"])
                self.assertEqual(models[replica]["numerics"], {**models[original_index]["numerics"], "random_seed": 43002})
                self.assertEqual(models[original_index]["numerics"]["random_seed"], 43001)
                self.assertNotIn("random_seed", models[replica]["parameters"])
            validate = ("import json,sys;from pathlib import Path;"
                        "sys.path.insert(0,str(Path(sys.argv[1])/'code'));"
                        "from production_task import validate_package;"
                        "print(json.dumps({'models':len(validate_package(Path(sys.argv[1]))['catalogue'])}))")
            valid = self.invoke(["-c", validate, production], cwd=root)
            self.assertEqual(json.loads(valid.stdout)["models"], 26)
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

            # Exercise the real worker and frozen runner for replica index 25.
            # Simulations and photometry are replaced only in this test child;
            # one temperature and one image receipt are enough to check seeds.
            replica_probe = '''
import json
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np
from astropy.io import fits
production, machine = map(Path, sys.argv[1:3])
sys.path[:0] = [str(production/'code'), str(production/'code/src')]
from production_task import run_task
from mcfost_grid import runner, photometry
manifest = json.loads((production/'manifest.json').read_text())
model = manifest['models'][25]
assert runner.model_numerics(manifest, model)['random_seed'] == 43002
calls = []
def invoke(command, cwd, env, timeout, attempt):
    calls.append(command)
    assert command[command.index('-seed')+1] == '43002'
    if len(calls) == 3:
        raise RuntimeError('synthetic test stop after one image')
    output = Path(cwd)/command[command.index('-root_dir')+1]
    output.mkdir(parents=True, exist_ok=True)
    name = 'Temperature.fits.gz' if command[1] == 'temperature.para' else 'RT.fits.gz'
    fits.PrimaryHDU(np.ones((3, 3))).writeto(output/name)
    return .01, 'Using scattering method 2; synthetic test only'
with patch.object(runner, 'cpu_capacity', return_value=64), \
     patch.object(runner, 'check_backend', return_value=['synthetic test only']), \
     patch.object(runner, '_invoke', side_effect=invoke), \
     patch.object(photometry, 'measure_image', return_value={'quality_pass':True, 'flux_jy':1., 'quality_policy':'aperture_v2'}):
    try:
        run_task(production, 25, machine)
    except RuntimeError as error:
        assert str(error) == 'synthetic test stop after one image', str(error)
    else:
        raise AssertionError('Expected bounded synthetic stop')
    folder = production/'models'/model['id']
    for name in ['runtime_binding.json','temperature_complete.json','measurements.json']:
        assert json.loads((folder/name).read_text())['numerics']['random_seed'] == 43002
    temperature = folder/'temperature_complete.json'
    original = temperature.read_bytes()
    receipt = json.loads(original)
    receipt['numerics']['random_seed'] = 43001
    temperature.write_text(json.dumps(receipt))
    try:
        runner.run_model(production, 25, machine)
    except RuntimeError as error:
        assert 'Cached temperature has different model numerics or photon seed' in str(error), str(error)
    else:
        raise AssertionError('Mixed-seed temperature was accepted')
    temperature.write_bytes(original)
assert len(calls) == 3
print(json.dumps({'replica_index':25,'seed':43002,'receipts_checked':3,'mixed_seed_cache_rejected':True}))
'''
            probe = self.invoke(["-c", replica_probe, production, machine], cwd=root)
            self.assertEqual(json.loads(probe.stdout.strip().splitlines()[-1]),
                             dict(replica_index=25, seed=43002, receipts_checked=3,
                                  mixed_seed_cache_rejected=True))

        self.assertEqual({path: sha256(path) for path in inherited_paths}, original)


if __name__ == "__main__":
    unittest.main()
