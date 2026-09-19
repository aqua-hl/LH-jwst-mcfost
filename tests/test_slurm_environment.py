"""Cluster environment selection must fail before submission or simulation."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "configure_slurm_environment", ROOT / "scripts/configure_slurm_environment.py")
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.slurm import write_slurm


class SlurmEnvironmentTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.run = self.root / "run"
        (self.run / "code").mkdir(parents=True)
        (self.run / "code/workflow.py").write_text("# frozen fixture\n")
        self.manifest = {"run_id": "continuum_production_v1.1",
                         "models": [{"id": str(i)} for i in range(96)]}
        (self.run / "manifest.json").write_text(json.dumps(self.manifest))
        self.executable = self.root / "mcfost"
        self.executable.write_text("#!/bin/sh\nexit 99\n")
        self.executable.chmod(0o755)
        self.utilities = self.root / "utils"
        for name in ("Dust", "Lambda", "Stellar_Spectra"):
            (self.utilities / name).mkdir(parents=True)
        self.machine = self.root / "machine.json"
        self.value = {"schema_version": 1, "mcfost_executable": "mcfost",
                      "mcfost_utils": "utils", "threads": 64,
                      "slurm": {"python": "python", "max_parallel": 16,
                                "setup_lines": ["module load working-environment"],
                                "analysis_cpus": 1}}
        self.machine.write_text(json.dumps(self.value))
        self.original_manifest = (self.run / "manifest.json").read_bytes()
        self.calls = []
        self.real_subprocess_run = subprocess.run

    def tearDown(self):
        self.temporary.cleanup()

    def fake_process(self, command, **kwargs):
        self.calls.append((command, kwargs))
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0)
        if "slurm" in command:
            write_slurm(command[4], command[6])
            return subprocess.CompletedProcess(command, 0)
        if command[:2] == ["bash", "-n"]:
            return self.real_subprocess_run(command, **kwargs)
        self.fail(f"Unexpected external command: {command}")

    def test_pins_venv_python_and_generates_array_analysis_without_submission(self):
        venv = self.root / "venv/bin"
        venv.mkdir(parents=True)
        selected_python = venv / "python"
        selected_python.symlink_to(sys.executable)
        with mock.patch.object(helper.sys, "executable", str(selected_python)), \
             mock.patch.object(helper.shutil, "which", return_value=str(self.executable)), \
             mock.patch.object(helper.subprocess, "run", side_effect=self.fake_process):
            configured = helper.configure(self.run, self.machine)
        machine = json.loads(configured.read_text())
        self.assertEqual(configured, self.root / "machine.cluster.json")
        self.assertEqual(machine["slurm"]["python"], str(selected_python))
        self.assertNotEqual(machine["slurm"]["python"], str(selected_python.resolve()))
        self.assertEqual(machine["mcfost_executable"], str(self.executable))
        self.assertEqual(machine["mcfost_utils"], str(self.utilities))
        self.assertEqual(json.loads(self.machine.read_text()), self.value)
        self.assertEqual((self.run / "manifest.json").read_bytes(), self.original_manifest)
        self.assertIn("module load working-environment", machine["slurm"]["setup_lines"])
        for name in ("job_array.sh", "job_analysis.sh"):
            script = (self.run / name).read_text()
            self.assertIn(f"exec {selected_python} -B code/workflow.py", script)
            self.assertIn("from mcfost_grid.photometry import measure_image", script)
            self.assertIn("from mcfost_grid.analysis import analyze_run", script)
            self.assertIn("socket.gethostname()", script)
            self.assertIn("p.__version__", script)
            self.assertIn("p.__file__", script)
            self.assertIn("MCFOST_AUTO_UPDATE=0", script)
        array = (self.run / "job_array.sh").read_text()
        self.assertIn("#SBATCH --array=0-95%16", array)
        self.assertIn("#SBATCH --cpus-per-task=64", array)
        self.assertEqual(self.calls[0][1]["cwd"], self.run)
        self.assertEqual(self.calls[0][1]["env"]["MCFOST_AUTO_UPDATE"], "0")
        self.assertEqual(self.calls[0][1]["env"]["OPENBLAS_NUM_THREADS"], "1")
        self.assertTrue(all(command[0] not in ("sbatch", str(self.executable))
                            for command, _ in self.calls))

    def test_missing_dependency_fails_before_writing_any_launch_files(self):
        with mock.patch.object(helper.shutil, "which", return_value=str(self.executable)), \
             mock.patch.object(helper.subprocess, "run", side_effect=
                               subprocess.CalledProcessError(1, ["python", "-c", "import astropy"])):
            with self.assertRaises(subprocess.CalledProcessError):
                helper.configure(self.run, self.machine)
        self.assertFalse((self.root / "machine.cluster.json").exists())
        self.assertFalse((self.run / "job_array.sh").exists())

    def test_different_existing_cluster_configuration_is_preserved(self):
        existing = self.root / "machine.cluster.json"
        original = '{"notes": "in use by an existing job"}\n'
        existing.write_text(original)
        with mock.patch.object(helper.shutil, "which", return_value=str(self.executable)), \
             mock.patch.object(helper.subprocess, "run", side_effect=self.fake_process):
            with self.assertRaisesRegex(FileExistsError, "Existing jobs may use it"):
                helper.configure(self.run, self.machine)
        self.assertEqual(existing.read_text(), original)
        self.assertFalse((self.run / "job_array.sh").exists())

    def test_missing_utilities_fails_before_processes_or_launch_files(self):
        (self.utilities / "Lambda").rmdir()
        with mock.patch.object(helper.shutil, "which", return_value=str(self.executable)), \
             mock.patch.object(helper.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "missing Lambda"):
                helper.configure(self.run, self.machine)
        run.assert_not_called()
        self.assertFalse((self.root / "machine.cluster.json").exists())

    def test_machine_relative_paths_survive_sibling_output(self):
        self.value["mcfost_executable"] = "./mcfost"
        self.value["mcfost_utils"] = "$SLURM_ENV_TEST_UTILS"
        self.machine.write_text(json.dumps(self.value))
        with mock.patch.dict(os.environ, {"SLURM_ENV_TEST_UTILS": "utils"}), \
             mock.patch.object(helper.subprocess, "run", side_effect=self.fake_process):
            configured = helper.configure(self.run, self.machine)
        machine = json.loads(configured.read_text())
        self.assertEqual(machine["mcfost_executable"], str(self.executable))
        self.assertEqual(machine["mcfost_utils"], str(self.utilities))

    def test_probe_executes_frozen_imports_and_reports_package_locations(self):
        source = self.run / "code/src/mcfost_grid"
        source.mkdir(parents=True)
        (source / "__init__.py").write_text("")
        (source / "photometry.py").write_text("def measure_image(): pass\n")
        (source / "analysis.py").write_text("def analyze_run(): pass\n")
        (source / "runner.py").write_text(
            "def validate_inputs(run, manifest):\n"
            "    if manifest.get('run_id') != 'continuum_production_v1.1':\n"
            "        raise RuntimeError('Frozen input mismatch')\n")
        probe = helper.dependency_probe(str(self.executable), str(self.utilities))
        self.assertNotIn("\n", probe)
        environment = os.environ.copy()
        environment.update(MPLCONFIGDIR=str(self.root / ".mplconfig"),
                           XDG_CACHE_HOME=str(self.root / ".cache"))
        result = subprocess.run([sys.executable, "-B", "-c", probe], cwd=self.run,
                                env=environment, check=True, capture_output=True, text=True)
        receipt = json.loads(result.stdout)
        self.assertEqual(receipt["environment_check"], "passed")
        self.assertTrue(receipt["host"])
        self.assertEqual(set(receipt["packages"]), {"numpy", "scipy", "astropy", "matplotlib"})
        self.assertTrue(all(package["file"] for package in receipt["packages"].values()))
        (self.run / "manifest.json").write_text(json.dumps({"run_id": "changed"}))
        result = subprocess.run([sys.executable, "-B", "-c", probe], cwd=self.run,
                                env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Frozen input mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
