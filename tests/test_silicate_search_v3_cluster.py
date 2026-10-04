"""Cluster dependencies, immutable environments and staged dispatch contracts."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import configure_silicate_search_v3 as configure
import silicate_search_v3_task as task


def machine_settings(root):
    value = dict(schema_version=1, mcfost_executable=str(root / "mcfost"),
                 mcfost_utils=str(root / "utils"), backend="image_method2",
                 slurm=dict(python=sys.executable, partition="compute", account="science",
                            reservation="reserved", setup_lines=["export TEST_SETUP=preserved"],
                            analysis_time="02:00:00"))
    configure.validate_settings(value, {"experiment_id": "silicate_search_v3"})
    return value


def scripts_fixture(root):
    machine = machine_settings(root)
    path = root / "machine.cluster.json"
    path.write_text(json.dumps(machine, indent=2, allow_nan=False) + "\n")
    scripts = configure.launch_scripts(root, machine, path,
                                       {"executable_sha256": "fake", "utilities_content_sha256": "fake"})
    for name, lines in scripts.items():
        text = "\n".join(lines) + "\n"
        subprocess.run(["bash", "-n"], input=text, text=True, check=True)
        (root / name).write_text(text)
    return machine, path, scripts


class LaunchTests(unittest.TestCase):
    def test_resources_setup_and_all_job_environment_checks(self):
        with tempfile.TemporaryDirectory(prefix="v3 cluster ") as tmp:
            root = Path(tmp)
            machine, _, scripts = scripts_fixture(root)
            self.assertEqual(machine["max_memory_gb"], 112)
            self.assertIn("#SBATCH --array=0-25%16", scripts["job_array.sh"])
            self.assertIn("#SBATCH --cpus-per-task=64", scripts["job_array.sh"])
            self.assertIn("#SBATCH --mem=160G", scripts["job_array.sh"])
            self.assertIn("#SBATCH --time=24:00:00", scripts["job_array.sh"])
            self.assertIn("#SBATCH --cpus-per-task=2", scripts["job_stage0.sh"])
            self.assertIn("#SBATCH --mem=8G", scripts["job_stage0.sh"])
            self.assertIn("#SBATCH --time=00:30:00", scripts["job_stage0.sh"])
            self.assertIn("#SBATCH --time=02:00:00", scripts["job_analysis.sh"])
            for name in ("job_stage0.sh", "job_array.sh", "job_analysis.sh"):
                code = "\n".join(scripts[name])
                self.assertIn("#SBATCH --partition=compute", code)
                self.assertIn("#SBATCH --account=science", code)
                self.assertIn("#SBATCH --reservation=reserved", code)
                self.assertIn("export TEST_SETUP=preserved", code)
                self.assertIn("Machine configuration changed", code)
                self.assertIn("utilities_content_sha256", code)
                self.assertIn("import numpy, scipy, astropy, matplotlib", code)
                self.assertIn("OMP_NUM_THREADS", code)

    def test_fake_sbatch_records_stage0_then_dependent_array_then_analysis(self):
        with tempfile.TemporaryDirectory(prefix="v3 cluster ") as tmp:
            root = Path(tmp)
            _, machine, _ = scripts_fixture(root)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            fake = bin_dir / "sbatch"
            fake.write_text(f"#!{sys.executable}\n" + '''import json, os
from pathlib import Path
import sys
path = Path(os.environ['FAKE_SBATCH_RECEIPT'])
calls = json.loads(path.read_text()) if path.exists() else []
calls.append(sys.argv[1:])
path.write_text(json.dumps(calls))
print(str(1000 + len(calls)) + ';test-cluster')
''')
            fake.chmod(0o755)
            receipt = root / "calls.json"
            env = {**os.environ, "PATH": str(bin_dir) + os.pathsep + os.environ["PATH"],
                   "FAKE_SBATCH_RECEIPT": str(receipt)}
            subprocess.run(["bash", str(root / "submit.sh"), str(machine)],
                           cwd="/tmp", env=env, check=True, text=True, capture_output=True)
            calls = json.loads(receipt.read_text())
            self.assertEqual([row[-1] for row in calls], ["job_stage0.sh", "job_array.sh", "job_analysis.sh"])
            self.assertIn("--dependency=afterok:1001", calls[1])
            self.assertIn("--kill-on-invalid-dep=yes", calls[1])
            self.assertIn("--dependency=afterany:1002", calls[2])
            folder, = (root / "submissions").iterdir()
            self.assertEqual((folder / "stage0_job_id.txt").read_text(), "1001\n")
            self.assertEqual((folder / "array_job_id.txt").read_text(), "1002\n")
            self.assertEqual((folder / "analysis_job_id.txt").read_text(), "1003\n")
            # Editing the machine after configuration must stop before sbatch.
            machine.write_text(machine.read_text() + "\n")
            failed = subprocess.run(["bash", str(root / "submit.sh"), str(machine)],
                                    cwd="/tmp", env=env, text=True, capture_output=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("Machine configuration changed", failed.stderr)
            self.assertEqual(len(json.loads(receipt.read_text())), 3)

    def test_configure_pins_venv_symlink_without_executing_simulator(self):
        with tempfile.TemporaryDirectory(prefix="v3 config ") as tmp:
            root = Path(tmp)
            code = root / "code"
            code.mkdir()
            (code / "build_silicate_search_v3.py").write_text(
                "def validate_bootstrap(bundle):\n    return {'model_count': 26}\n")
            shutil.copytree(ROOT / "src/mcfost_grid", code / "src/mcfost_grid",
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            for folder in ("Dust", "Lambda", "Stellar_Spectra"):
                (root / "utils" / folder).mkdir(parents=True)
            (root / "utils/Dust/example.dat").write_text("immutable utility")
            executable = root / "mcfost"
            marker = root / "SIMULATOR_WAS_EXECUTED"
            executable.write_text(f"#!{sys.executable}\nfrom pathlib import Path\nPath({str(marker)!r}).touch()\n")
            executable.chmod(0o755)
            python = root / "venv/bin/python"
            python.parent.mkdir(parents=True)
            python.symlink_to(sys.executable)
            source = root / "machine.template.json"
            source.write_text(json.dumps(machine_settings(root)))
            import build_silicate_search_v3 as builder
            with patch.object(builder, "validate_bootstrap", return_value={"model_count": 26}), \
                    patch.object(configure.sys, "executable", str(python)):
                configured = configure.configure(root, source)
                settings = json.loads(configured.read_text())
                self.assertEqual(settings["slurm"]["python"], str(python))
                self.assertFalse(marker.exists())
                self.assertEqual(configure.configure(root, source), configured)
                changed = json.loads(source.read_text())
                changed["slurm"]["partition"] = "other"
                source.write_text(json.dumps(changed))
                with self.assertRaisesRegex(FileExistsError, "different cluster configuration"):
                    configure.configure(root, source)
            snapshot = json.loads((root / "cluster_environment.json").read_text())
            self.assertEqual(snapshot["python"], str(python))
            self.assertEqual(set(snapshot["packages"]), {"numpy", "scipy", "astropy", "matplotlib"})
            self.assertIn("utilities_content_sha256", snapshot["runtime_identity"])


class DispatcherTests(unittest.TestCase):
    def test_stage0_failure_cannot_materialize_and_success_preserves_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = []
            builder = SimpleNamespace(validate_bootstrap=lambda bundle: calls.append("validate") or {"model_count": 26},
                                      materialize=Mock(side_effect=lambda bundle: calls.append("materialize") or {"models": 26}))
            stage = SimpleNamespace(run_stage0=Mock(side_effect=ValueError("invalid opacity")))
            with patch.dict(sys.modules, {"build_silicate_search_v3": builder, "silicate_search_v3_stage0": stage}):
                with self.assertRaisesRegex(ValueError, "invalid opacity"):
                    task.dispatch(tmp, stage0=True, machine=Path(tmp) / "machine.json")
                builder.materialize.assert_not_called()
                calls.clear()
                stage.run_stage0.side_effect = lambda *_: calls.append("stage0") or {"status": "cached"}
                result = task.dispatch(tmp, stage0=True, machine=Path(tmp) / "machine.json")
                self.assertEqual(calls, ["validate", "stage0", "materialize"])
                self.assertEqual(result["production"]["models"], 26)

    def test_children_use_production_code_fresh_venv_and_exit_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            (root / "production").mkdir()
            builder = SimpleNamespace(validate_bootstrap=Mock(return_value={"model_count": 26}))
            with patch.dict(sys.modules, {"build_silicate_search_v3": builder}), \
                    patch.object(task.subprocess, "run") as child, \
                    patch.object(task.sys, "executable", "/shared/venv/bin/python"):
                task.dispatch(root, index=25, machine=root / "machine.json")
                args = child.call_args.args[0]
                self.assertEqual(args[0], "/shared/venv/bin/python")
                self.assertEqual(args[2], str(root / "production/code/production_task.py"))
                self.assertIn("--index", args)
                self.assertIn("25", args)
                self.assertEqual(child.call_args.kwargs, {"check": True})
                task.dispatch(root, status=True)
                self.assertIn("--status", child.call_args.args[0])
                task.dispatch(root, analyze=True)
                self.assertEqual(child.call_args.args[0][-2:], ["--output", str(root / "results")])
                child.side_effect = subprocess.CalledProcessError(7, "worker")
                with self.assertRaises(subprocess.CalledProcessError) as caught:
                    task.dispatch(root, index=0, machine=root / "machine.json")
                self.assertEqual(caught.exception.returncode, 7)
                self.assertEqual(builder.validate_bootstrap.call_count, 4)

    def test_absent_production_reports_stage0_without_creating_placeholders(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            builder = SimpleNamespace(validate_bootstrap=Mock(return_value={"model_count": 26}))
            with patch.dict(sys.modules, {"build_silicate_search_v3": builder}):
                result = task.dispatch(root, status=True)
                self.assertFalse(result["production_materialized"])
                self.assertEqual(result["stage0"]["status"], "not_started")
                result = task.dispatch(root, validate=True)
                self.assertFalse(result["mcfost_invoked"])
                with self.assertRaisesRegex(ValueError, "Production is absent"):
                    task.dispatch(root, index=0, machine=root / "machine.json")
                self.assertFalse((root / "production").exists())


if __name__ == "__main__":
    unittest.main()
