"""Submission stays explicit and uses one checked environment across all seeds."""
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
    "configure_extinction_ice_numerics_v2", ROOT / "scripts/configure_extinction_ice_numerics_v2.py")
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)


class NumericalClusterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="numerics cluster ")
        self.root = Path(self.temporary.name).resolve()
        self.bundle = self.root / "bundle space"
        code = self.bundle / "code"
        code.mkdir(parents=True)
        (code / "numerical_task.py").write_text(
            "def validate_package(bundle):\n"
            "    if not (bundle / 'package_ok').is_file():\n"
            "        raise RuntimeError('Package input mismatch')\n")
        (self.bundle / "package_ok").touch()
        (code / "analyze_extinction_ice_numerics_v2.py").write_text("# frozen analyzer fixture\n")
        tasks = []
        for pair in range(10):
            relative = f"subruns/pair{pair}"
            run = self.bundle / relative
            source = run / "code/src/mcfost_grid"
            source.mkdir(parents=True)
            (run / "code/workflow.py").write_text("# frozen workflow fixture\n")
            (source / "__init__.py").write_text("")
            (source / "photometry.py").write_text("def measure_image(): pass\n")
            (source / "runner.py").write_text(
                "def validate_inputs(run, manifest):\n"
                "    if manifest.get('valid') is not True:\n"
                "        raise RuntimeError('Frozen subrun input mismatch')\n")
            (run / "manifest.json").write_text(json.dumps({"valid": True}))
            tasks.extend({"index": pair * 2 + index, "run_path": relative, "model_index": index}
                         for index in range(2))
        (self.bundle / "experiment.json").write_text(json.dumps(
            {"experiment_id": "extinction_ice_numerics_v2", "tasks": tasks}))
        self.executable = self.root / "fake mcfost"
        self.executable.write_text("#!/bin/sh\nexit 99\n")
        self.executable.chmod(0o755)
        self.utilities = self.root / "utils"
        for name in ("Dust", "Lambda", "Stellar_Spectra"):
            (self.utilities / name).mkdir(parents=True)
        self.machine = self.root / "machine.json"
        self.value = {"schema_version": 1, "mcfost_executable": "./fake mcfost",
                      "mcfost_utils": "utils", "threads": 64,
                      "slurm": {"max_parallel": 16, "setup_lines": ["# environment already active"]}}
        self.machine.write_text(json.dumps(self.value))
        self.calls = []
        self.real_run = subprocess.run

    def tearDown(self):
        self.temporary.cleanup()

    def fake_process(self, command, **kwargs):
        self.calls.append((command, kwargs))
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0)
        if command[:2] == ["bash", "-n"]:
            return self.real_run(command, **kwargs)
        self.fail(f"Unexpected external command: {command}")

    def configure(self):
        with mock.patch.object(helper.subprocess, "run", side_effect=self.fake_process):
            return helper.configure(self.bundle, self.machine)

    def test_pins_venv_and_checks_scripts_without_submission(self):
        venv = self.root / "venv space/bin"
        venv.mkdir(parents=True)
        selected_python = venv / "python"
        selected_python.symlink_to(sys.executable)
        with mock.patch.object(helper.sys, "executable", str(selected_python)):
            configured = self.configure()
        machine = json.loads(configured.read_text())
        self.assertEqual(machine["slurm"]["python"], str(selected_python))
        self.assertNotEqual(str(selected_python), str(selected_python.resolve()))
        self.assertEqual(machine["mcfost_executable"], str(self.executable))
        self.assertEqual(machine["mcfost_utils"], str(self.utilities))
        self.assertEqual(json.loads(self.machine.read_text()), self.value)
        array = (self.bundle / "job_array.sh").read_text()
        analysis = (self.bundle / "job_analysis.sh").read_text()
        self.assertIn("#SBATCH --array=0-19%16", array)
        self.assertIn("#SBATCH --cpus-per-task=64", array)
        self.assertIn("#SBATCH --time=24:00:00", array)
        self.assertIn("#SBATCH --mem=16G", array)
        self.assertIn("code/numerical_task.py", array)
        self.assertIn("code/analyze_extinction_ice_numerics_v2.py", analysis)
        for script in (array, analysis):
            self.assertIn("validate_package(bundle)", script)
            self.assertIn("validate_inputs(run", script)
            self.assertIn("socket.gethostname()", script)
            self.assertIn("p.__file__", script)
            self.assertIn("Machine configuration changed", script)
        self.assertEqual(len(self.calls), 4)
        self.assertEqual(self.calls[0][1]["env"]["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual(self.calls[0][1]["env"]["MCFOST_AUTO_UPDATE"], "0")
        self.assertTrue(all(command[0] not in ("sbatch", str(self.executable))
                            for command, _ in self.calls))

    def test_failure_does_not_write_launch_files_and_existing_config_is_preserved(self):
        with mock.patch.object(helper.subprocess, "run", side_effect=
                               subprocess.CalledProcessError(1, ["python", "-c", "import scipy"])):
            with self.assertRaises(subprocess.CalledProcessError):
                helper.configure(self.bundle, self.machine)
        configured = self.root / "machine.cluster.json"
        self.assertFalse(configured.exists())
        self.assertFalse((self.bundle / "submit.sh").exists())
        configured.write_text('{"already": "in use"}\n')
        with self.assertRaisesRegex(FileExistsError, "Existing jobs may use it"):
            self.configure()
        self.assertEqual(configured.read_text(), '{"already": "in use"}\n')
        self.assertFalse((self.bundle / "submit.sh").exists())

    def test_resource_and_path_injection_rejected_before_probe(self):
        for key, value in (("memory", "16G\n#SBATCH --nodes=40"), ("partition", "debug;bad"),
                           ("time", "24:99:00"), ("max_parallel", 17), ("analysis_cpus", True)):
            machine = json.loads(json.dumps(self.value))
            machine["slurm"][key] = value
            self.machine.write_text(json.dumps(machine))
            with self.assertRaises(ValueError):
                self.configure()
        self.assertEqual(self.calls, [])
        experiment = json.loads((self.bundle / "experiment.json").read_text())
        experiment["tasks"][0]["run_path"] = "../../outside"
        (self.bundle / "experiment.json").write_text(json.dumps(experiment))
        with self.assertRaisesRegex(ValueError, "Invalid or incomplete task subrun"):
            self.configure()

    def test_real_probe_checks_all_subruns_and_top_package(self):
        probe = helper.dependency_probe(str(self.executable), str(self.utilities))
        environment = os.environ.copy()
        environment.update(MPLCONFIGDIR=str(self.root / ".mplconfig"))
        arguments = [sys.executable, "-B", "-c", probe]
        result = self.real_run(arguments, cwd=self.bundle, env=environment,
                               capture_output=True, text=True, check=True)
        receipt = json.loads(result.stdout)
        self.assertEqual(receipt["frozen_subruns_checked"], 10)
        self.assertEqual(set(receipt["packages"]), {"numpy", "scipy", "astropy", "matplotlib"})
        self.assertTrue(all(package["file"] for package in receipt["packages"].values()))
        last = self.bundle / "subruns/pair9/manifest.json"
        last.write_text('{"valid": false}')
        result = self.real_run(arguments, cwd=self.bundle, env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Frozen subrun input mismatch", result.stderr)
        last.write_text('{"valid": true}')
        (self.bundle / "package_ok").unlink()
        result = self.real_run(arguments, cwd=self.bundle, env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Package input mismatch", result.stderr)

    def fake_scheduler(self, fail_analysis=False):
        fake_bin = self.root / "scheduler"
        fake_bin.mkdir(exist_ok=True)
        scheduler = fake_bin / "sbatch"
        scheduler.write_text("#!/bin/bash\nset -eu\n"
                             "printf '%s\\n' \"$*\" >> \"$FAKE_SBATCH_LOG\"\n"
                             "if [[ $* == *job_analysis.sh* ]]; then\n"
                             + ("  echo 'analysis unavailable' >&2; exit 1\n" if fail_analysis
                                else "  echo '222;cluster'\n")
                             + "else\n  echo '111;cluster'\nfi\n")
        scheduler.chmod(0o755)
        return {**os.environ, "PATH": str(fake_bin) + os.pathsep + os.environ["PATH"],
                "FAKE_SBATCH_LOG": str(self.root / "scheduler.log")}

    def test_explicit_submission_records_both_jobs_with_afterany(self):
        configured = self.configure()
        result = self.real_run(["bash", str(self.bundle / "submit.sh"), str(configured)],
                               env=self.fake_scheduler(), capture_output=True, text=True, check=True)
        receipt = list((self.bundle / "submissions").iterdir())
        self.assertEqual(len(receipt), 1)
        self.assertEqual((receipt[0] / "array_job_id.txt").read_text(), "111\n")
        self.assertEqual((receipt[0] / "analysis_job_id.txt").read_text(), "222\n")
        calls = (self.root / "scheduler.log").read_text().splitlines()
        self.assertEqual(len(calls), 2)
        self.assertIn("--dependency=afterany:111", calls[1])
        self.assertIn("Submitted analysis 222", result.stdout)

    def test_analysis_submission_failure_keeps_array_receipt_and_changed_machine_blocks(self):
        configured = self.configure()
        environment = self.fake_scheduler(fail_analysis=True)
        arguments = ["bash", str(self.bundle / "submit.sh"), str(configured)]
        result = self.real_run(arguments, env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        receipt = list((self.bundle / "submissions").iterdir())[0]
        self.assertEqual((receipt / "array_job_id.txt").read_text(), "111\n")
        self.assertFalse((receipt / "analysis_job_id.txt").exists())
        self.assertIn("do not resubmit the array", result.stderr)
        configured.write_text(configured.read_text() + "\n")
        result = self.real_run(arguments, env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Machine configuration changed", result.stderr)
        self.assertEqual(len((self.root / "scheduler.log").read_text().splitlines()), 2)


if __name__ == "__main__":
    unittest.main()
