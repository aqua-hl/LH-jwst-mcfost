"""Launch dependencies and the numerical material gate for the 60-model run."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
sys.path.insert(0, str(ROOT/"src"))
import check_silicate_structure_materials as check
from configure_extinction_ice_numerics_v2 import (SILICATE_STRUCTURE, bundle_layout,
    launch_scripts, validate_settings)


class LaunchTests(unittest.TestCase):
    def test_preflight_dependency_before_array_and_afterany_analysis(self):
        with tempfile.TemporaryDirectory(prefix="structure-cluster-") as tmp:
            root = Path(tmp)
            experiment = dict(schema_version=1, experiment_id=SILICATE_STRUCTURE,
                              tasks=[dict(index=i) for i in range(60)])
            machine = dict(schema_version=1, mcfost_executable="/unused/mcfost",
                           mcfost_utils="/unused/utils", slurm={"python": sys.executable})
            validate_settings(machine, experiment)
            self.assertEqual(machine["threads"], 64)
            self.assertEqual(machine["max_memory_gb"], 112)
            self.assertEqual(machine["slurm"]["max_parallel"], 16)
            self.assertEqual(machine["slurm"]["memory"], "160G")
            configured = root/"machine.cluster.json"
            configured.write_text(json.dumps(machine, indent=2, allow_nan=False)+"\n")
            scripts = launch_scripts(root, experiment, machine, configured)
            for name, lines in scripts.items():
                script = "\n".join(lines)+"\n"
                subprocess.run(["bash", "-n"], input=script, text=True, check=True)
                (root/name).write_text(script)
            self.assertIn("#SBATCH --array=0-59%16", scripts["job_array.sh"])
            self.assertIn("#SBATCH --cpus-per-task=2", scripts["job_materials.sh"])
            self.assertIn("code/check_materials.py", "\n".join(scripts["job_materials.sh"]))
            self.assertIn("code/analyze_silicate_structure_production.py", "\n".join(scripts["job_analysis.sh"]))
            submit = "\n".join(scripts["submit.sh"])
            self.assertIn('--dependency="afterok:$materials_job" --kill-on-invalid-dep=yes', submit)
            self.assertIn('--dependency="afterany:$job"', submit)
            self.assertLess(submit.index("job_materials.sh"), submit.index("job_array.sh"))

    def test_layout_rejects_missing_gate_or_wrong_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"code").mkdir()
            experiment = dict(schema_version=1, experiment_id=SILICATE_STRUCTURE,
                              tasks=[dict(index=i) for i in range(60)])
            (root/"experiment.json").write_text(json.dumps(experiment))
            for name in ("manifest.json", "manifest.sha256", "production_experiment.json",
                         "code/production_task.py", "code/analyze_extinction_ice_production.py",
                         "code/analyze_silicate_structure_production.py"):
                (root/name).touch()
            with self.assertRaisesRegex(ValueError, "material preflight"):
                bundle_layout(root)
            (root/"code/check_materials.py").touch()
            self.assertEqual(bundle_layout(root), experiment)
            experiment["tasks"].pop()
            (root/"experiment.json").write_text(json.dumps(experiment))
            with self.assertRaisesRegex(ValueError, "60 tasks"):
                bundle_layout(root)


def write_outputs(root, wave, *, opacity=None, albedo=None):
    root = Path(root)/"data_dust"
    root.mkdir(exist_ok=True)
    ones = np.ones(len(wave))
    arrays = dict(lambda_=wave, kappa=ones if opacity is None else opacity,
                  albedo=.5*ones if albedo is None else albedo, g=.2*ones,
                  kappa_grain=np.ones((2, len(wave))), phase_function=np.ones((181, len(wave))),
                  polarizability=np.zeros((181, len(wave))))
    for name, values in arrays.items():
        fits.PrimaryHDU(np.asarray(values)).writeto(root/(name.rstrip("_")+".fits.gz"), overwrite=True)


class MaterialTests(unittest.TestCase):
    def test_all_materials_must_pass_before_receipt_can_release_workers(self):
        import shutil
        from silicate_structure_design import catalogue
        records = catalogue()
        unique = {}
        for row in records:
            p = row["parameters"]
            unique.setdefault((p["envelope_silicate_file"], p["envelope_amax_um"]), row)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"inputs/utils/Dust").mkdir(parents=True)
            shutil.copy2(ROOT/"reference/dust/H2O_30K_Leiden_mcfost.dat", root/"inputs/utils/Dust")
            models = [dict(id=f"m{i}", parameters=row["parameters"]) for i, row in enumerate(unique.values())]
            (root/"manifest.json").write_text(json.dumps(dict(models=models, anchors=[dict(wavelength_um=9.7)])))
            for model in models:
                folder = root/"models"/model["id"]
                folder.mkdir(parents=True)
                shutil.copy2(ROOT/"reference/parameters/ice_v02_dust.para", folder/"temperature.para")
            runner = SimpleNamespace(
                runtime_config=lambda _: {"mcfost_executable": "/fake/mcfost"},
                environment=lambda _, runtime: {"OMP_NUM_THREADS": str(runtime["threads"])},
                now=lambda: "2026-09-24T00:00:00+00:00",
                atomic_json=lambda path, payload: Path(path).write_text(json.dumps(payload)))
            dispatch = SimpleNamespace(validate_package=lambda _: {"experiment_id": SILICATE_STRUCTURE},
                                       load_runner=lambda _: runner)
            def dust_command(command, **kwargs):
                self.assertIn("-dust_prop", command)
                self.assertNotIn("-img", command)
                self.assertEqual(kwargs["timeout"], 120)
                self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "2")
                wave = np.loadtxt(kwargs["cwd"]/"check.lambda", skiprows=1)
                write_outputs(kwargs["cwd"], wave)
                return SimpleNamespace(returncode=0)
            with patch.dict(sys.modules, {"production_task": dispatch}), \
                    patch.object(check, "identity", return_value={"runtime": "same"}), \
                    patch.object(check.subprocess, "run", side_effect=dust_command) as invocation:
                self.assertEqual(check.run_checks(root, "unused")["status"], "passed")
                self.assertEqual(invocation.call_count, 4)
                check.validate_receipt(root, {})
                self.assertEqual(check.run_checks(root, "unused")["status"], "cached")
                self.assertEqual(invocation.call_count, 4)
                receipt = json.loads((root/check.RECEIPT).read_text())
                receipt["status"] = "failed"
                (root/check.RECEIPT).write_text(json.dumps(receipt))
                with self.assertRaisesRegex(ValueError, "not passed"):
                    check.validate_receipt(root, {})
            with patch.dict(sys.modules, {"production_task": dispatch}), \
                    patch.object(check, "identity", return_value={"runtime": "same"}), \
                    patch.object(check.subprocess, "run", return_value=SimpleNamespace(returncode=-11)):
                with self.assertRaisesRegex(ValueError, "initialization failed"):
                    check.run_checks(root, "unused")
                receipt = json.loads((root/check.RECEIPT).read_text())
                self.assertEqual(receipt["status"], "failed")
                self.assertEqual(receipt["checks"][0]["returncode"], -11)
                self.assertEqual(len(list((root/"material_preflight"/"m0").glob("attempt_*"))), 2)

    def test_isolated_probe_preserves_exact_envelope_species(self):
        original = (ROOT/"reference/parameters/ice_v02_dust.para").read_text()
        result = check.isolated_parameter(original)
        def grain_rows(text):
            section = text.split("#Grain properties", 1)[1].split("#Molecular RT settings", 1)[0]
            return [line for line in section.splitlines() if line.strip()]
        self.assertEqual(grain_rows(result)[1:], grain_rows(original)[10:])
        self.assertNotIn("ac_opct.dat", result)
        self.assertIn("F T F  no temperature", result)
        self.assertIn("8 8 1 2", result)
        self.assertIn("3 3 6000.", result)
        self.assertEqual(result.count("DHS  1 1"), 2)

    def test_full_range_grid_contains_all_zero_k_and_probes(self):
        import shutil
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"inputs/utils/Dust").mkdir(parents=True)
            ice = ROOT/"reference/dust/H2O_30K_Leiden_mcfost.dat"
            shutil.copy2(ice, root/"inputs/utils/Dust"/ice.name)
            probes = [dict(wavelength_um=v) for v in (1.74, 2.546, 9.7, 27.51)]
            wave = check.wavelength_grid(root, {"anchors": probes})
            self.assertEqual((wave[0], wave[-1]), (.1, 3000.))
            self.assertTrue(np.all(np.diff(wave) > 0))
            rows = [line.split("#", 1)[0].strip() for line in ice.read_text().splitlines()]
            samples = np.loadtxt([row for row in rows if row][1:])
            knots = samples[samples[:, 2] == 0, 0]
            self.assertEqual(len(knots), 206)
            self.assertTrue(np.isin(knots, wave).all())
            self.assertTrue(np.isin([r["wavelength_um"] for r in probes], wave).all())

    def test_dust_outputs_reject_invalid_numbers_and_wavelengths(self):
        wave = np.array([.1, 1., 9.7, 3000.])
        with tempfile.TemporaryDirectory() as tmp:
            write_outputs(tmp, wave)
            result = check.check_outputs(tmp, wave)
            self.assertEqual(result["wavelength_count"], 4)
            self.assertEqual(result["absorption_min_cm2_g"], .5)
            for opacity, albedo, message in (
                (np.array([1., np.nan, 1., 1.]), None, "Nonfinite"),
                (np.array([1., -1., 1., 1.]), None, "Negative"),
                (None, np.array([.5, 1.01, .5, .5]), "Albedo outside"),
            ):
                write_outputs(tmp, wave, opacity=opacity, albedo=albedo)
                with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                    check.check_outputs(tmp, wave)
            write_outputs(tmp, wave)
            with self.assertRaisesRegex(ValueError, "wavelengths"):
                check.check_outputs(tmp, wave*2)

    def test_receipt_binds_runtime_prescriptions_and_output_bytes(self):
        models = [dict(id=f"m{i}", parameters=dict(envelope_silicate_file=material, envelope_amax_um=size))
                  for i, (material, size) in enumerate((("Draine", .4), ("Draine", 1.), ("Olivine", .4), ("Pyroxene", .4)))]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"manifest.json").write_text(json.dumps({"models": models}))
            (root/"material_preflight").mkdir()
            artifact = root/"material_preflight/kappa.fits"
            artifact.write_bytes(b"a checked product")
            identity = {"runtime": "original"}
            receipt = dict(status="passed", identity=identity, checks=[
                dict(model_id=m["id"], status="passed", artifact_sha256={str(artifact.relative_to(root)): check.digest(artifact)})
                for m in models])
            (root/check.RECEIPT).write_text(json.dumps(receipt))
            with patch.object(check, "identity", return_value=identity):
                self.assertEqual(check.validate_receipt(root, {}), receipt)
                artifact.write_bytes(b"changed")
                with self.assertRaisesRegex(ValueError, "artifacts changed"):
                    check.validate_receipt(root, {})
            with patch.object(check, "identity", return_value={"runtime": "different"}):
                with self.assertRaisesRegex(ValueError, "runtime or frozen inputs changed"):
                    check.validate_receipt(root, {})


if __name__ == "__main__":
    unittest.main()
