"""Numerical-only mode and immutable, portable 20-task package contracts."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.configuration import prepare_run, read_anchors
from mcfost_grid.analysis import AnalysisError, analyze_run


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


builder = module("numerical_builder", "build_extinction_ice_numerical_package.py")
dispatcher = module("numerical_dispatch", "extinction_ice_numerical_task.py")


class NumericalOnlyTests(unittest.TestCase):
    def test_unscored_probes_need_no_observed_placeholders_and_cannot_be_ranked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            Table(rows=[dict(id="n001", wavelength_um=9.7, instrument="MIRI", region="probe", score=False)]).write(root / "probes.ecsv")
            config = {"schema_version": 1, "run_name": "test", "output_dir": str(root / "runs"),
                      "numerical_only": True, "observations": {"anchors_file": "probes.ecsv",
                      "spectrum_file": str(ROOT / "reference/observations/continuum_sed_R100.ecsv")}}
            anchors, _, _ = read_anchors(config, root)
            self.assertNotIn("flux_jy", anchors[0])
            self.assertFalse(anchors[0]["score"])
            path = root / "config.json"
            path.write_text(json.dumps(config))
            run = prepare_run(path)
            with self.assertRaisesRegex(AnalysisError, "Numerical-only"):
                analyze_run(run)
            del config["numerical_only"]
            with self.assertRaisesRegex(ValueError, "requires"):
                read_anchors(config, root)
            config["numerical_only"] = "yes"
            with self.assertRaisesRegex(ValueError, "boolean"):
                read_anchors(config, root)


class PackageTests(unittest.TestCase):
    def test_twenty_distinct_replicates_relocate_and_detect_input_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = subprocess.run([sys.executable, "-B", str(ROOT / "scripts/build_extinction_ice_numerical_package.py"),
                                       "--no-archive", "--output", str(root / "prepared")],
                                      capture_output=True, text=True)
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            result = json.loads(prepared.stdout)
            self.assertEqual(result["tasks"], 20)
            self.assertIsNone(result["archive"])
            bundle = root / "portable folder"
            (root / "prepared").rename(bundle)
            experiment = dispatcher.validate_package(bundle)
            self.assertEqual(dispatcher.status(bundle)["counts"], {"not_started": 20})
            self.assertEqual(len({t["run_path"] for t in experiment["tasks"]}), 10)
            for task in experiment["tasks"]:
                manifest = json.loads((bundle / task["run_path"] / "manifest.json").read_text())
                self.assertEqual(manifest["configuration"]["numerics"]["random_seed"], task["seed"])
                self.assertTrue(all(a["score"] is False and "flux_jy" not in a for a in manifest["anchors"]))
            validated = subprocess.run([sys.executable, "-B", str(bundle / "code/numerical_task.py"), str(bundle), "--validate"], capture_output=True, text=True)
            self.assertEqual(validated.returncode, 0, validated.stderr)
            # Machine partition/account settings are deliberately editable.
            (bundle / "machine.template.json").write_text('{"editable": true}')
            dispatcher.validate_package(bundle)
            # Shared package identity must prevent mixed binaries even when
            # tasks live in different otherwise-independent frozen subruns.
            executable = root / "fake_mcfost"
            executable.write_text("fixture binary A")
            runtime = {"mcfost_executable": str(executable), "threads": 64,
                       "backend": "image_method2"}
            fake = SimpleNamespace(runtime_config=lambda path: runtime,
                                   cpu_capacity=lambda: 64,
                                   _utilities_hash=lambda value: "a" * 64,
                                   now=lambda: "test", run_model=mock.Mock(return_value={"status": "complete"}),
                                   atomic_json=lambda path, value: path.write_text(json.dumps(value)))
            with mock.patch.object(dispatcher, "load_runner", return_value=fake):
                dispatcher.run_task(bundle, 0, bundle / "machine.template.json")
                self.assertEqual(fake.run_model.call_count, 1)
                executable.write_text("fixture binary B")
                with self.assertRaisesRegex(RuntimeError, "runtime changed"):
                    dispatcher.run_task(bundle, 2, bundle / "machine.template.json")
                self.assertEqual(fake.run_model.call_count, 1)
                runtime["threads"] = 32
                with self.assertRaisesRegex(ValueError, "64-thread"):
                    dispatcher.run_task(bundle, 2, bundle / "machine.template.json")
            task = experiment["tasks"][0]
            para = bundle / task["run_path"] / "models" / task["model_id"] / "temperature.para"
            para.write_text(para.read_text() + "\nchanged")
            with self.assertRaisesRegex(ValueError, "Frozen package file"):
                dispatcher.validate_package(bundle)
            with self.assertRaises(FileExistsError):
                builder.build_package(bundle)


if __name__ == "__main__":
    unittest.main()
