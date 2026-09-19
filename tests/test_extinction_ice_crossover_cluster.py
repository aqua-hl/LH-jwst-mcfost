"""The fixed-temperature image experiment uses the checked cluster environment."""
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
    "configure_crossover_cluster", ROOT / "scripts/configure_extinction_ice_numerics_v2.py")
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)


class CrossoverClusterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="crossover cluster ")
        self.root = Path(self.temporary.name).resolve()
        self.bundle = self.root / "bundle space"
        source = self.bundle / "code/src/mcfost_grid"
        source.mkdir(parents=True)
        (source / "__init__.py").write_text("")
        (source / "runner.py").write_text("# frozen runner fixture\n")
        (source / "photometry.py").write_text("def measure_image(): pass\n")
        (self.bundle / "code/crossover_task.py").write_text(
            "import sys\n"
            "def validate_package(bundle):\n"
            "    if (bundle / 'frozen_temperature').read_bytes() != b'checked source temperature':\n"
            "        raise RuntimeError('Frozen crossover temperature mismatch')\n"
            "def load_modules(bundle):\n"
            "    sys.path.insert(0, str(bundle / 'code/src'))\n"
            "    from mcfost_grid import runner, photometry\n"
            "    return runner, photometry\n")
        (self.bundle / "code/analyze_extinction_ice_crossover_v2.py").write_text("# analyzer fixture\n")
        (self.bundle / "frozen_temperature").write_bytes(b"checked source temperature")
        self.experiment = {"experiment_id": "extinction_ice_crossover_v2", "schema_version": 2,
                           "diagnostic_id": "fixed_temperature_pixel_scale_v2",
                           "tasks": [{"index": index} for index in range(5)]}
        (self.bundle / "experiment.json").write_text(json.dumps(self.experiment))
        self.executable = self.root / "fake mcfost"
        self.executable.write_text("#!/bin/sh\nexit 99\n")
        self.executable.chmod(0o755)
        self.utilities = self.root / "utilities"
        for name in ("Dust", "Lambda", "Stellar_Spectra"):
            (self.utilities / name).mkdir(parents=True)
        self.machine = self.root / "machine.json"
        self.settings = {"schema_version": 1, "mcfost_executable": str(self.executable),
                         "mcfost_utils": str(self.utilities), "threads": 64,
                         "slurm": {"max_parallel": 16}}
        self.machine.write_text(json.dumps(self.settings))

    def tearDown(self):
        self.temporary.cleanup()

    def configure_with_real_probe(self):
        # The executable intentionally exits 99: successful configuration also
        # establishes that neither a temperature solve nor image was launched.
        return helper.configure(self.bundle, self.machine)

    def test_real_configuration_checks_frozen_inputs_and_routes_five_pairs(self):
        configured = self.configure_with_real_probe()
        machine = json.loads(configured.read_text())
        self.assertEqual(machine["slurm"]["python"], os.path.abspath(sys.executable))
        array = (self.bundle / "job_array.sh").read_text()
        analysis = (self.bundle / "job_analysis.sh").read_text()
        self.assertIn("#SBATCH --array=0-4%5", array)
        self.assertIn("#SBATCH --cpus-per-task=64", array)
        self.assertIn("#SBATCH --mem=96G", array)
        self.assertEqual(machine["max_memory_gb"], 64)
        self.assertIn("code/crossover_task.py", array)
        self.assertIn("code/analyze_extinction_ice_crossover_v2.py", analysis)
        self.assertIn("#SBATCH --mem=8G", analysis)
        self.assertNotIn("code/numerical_task.py", array)
        for content in (array, analysis):
            self.assertIn("validate_package(bundle)", content)
            self.assertIn("load_modules(bundle)", content)
            self.assertIn("Machine configuration changed", content)
            self.assertIn("OPENBLAS_NUM_THREADS=1", content)
        self.assertIn("--dependency=\"afterany:$job\"", (self.bundle / "submit.sh").read_text())
        self.assertEqual(json.loads(self.machine.read_text()), self.settings)

    def test_corrupted_temperature_blocks_configuration_before_launch_files(self):
        (self.bundle / "frozen_temperature").write_bytes(b"changed temperature")
        with self.assertRaises(subprocess.CalledProcessError):
            self.configure_with_real_probe()
        self.assertFalse(self.machine.with_name("machine.cluster.json").exists())
        self.assertFalse((self.bundle / "submit.sh").exists())

    def test_incomplete_or_misindexed_pair_grid_rejected_before_probe(self):
        for tasks in (self.experiment["tasks"][:-1], [{"index": 1}] * 5):
            experiment = {**self.experiment, "tasks": tasks}
            (self.bundle / "experiment.json").write_text(json.dumps(experiment))
            with mock.patch.object(helper.subprocess, "run") as process:
                with self.assertRaisesRegex(ValueError, "five consecutively indexed"):
                    helper.configure(self.bundle, self.machine)
                process.assert_not_called()
        (self.bundle / "experiment.json").write_text(json.dumps(self.experiment))
        (self.bundle / "code/src/mcfost_grid/photometry.py").unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete crossover bundle"):
            helper.bundle_layout(self.bundle)

    def test_old_temperature_crossover_and_mislabeled_design_rejected_before_probe(self):
        for experiment in (
            {"schema_version": 1, "experiment_id": "extinction_ice_crossover_v2",
             "tasks": [{"index": index} for index in range(10)]},
            {**self.experiment, "schema_version": 1},
            {**self.experiment, "diagnostic_id": "temperature_crossover"},
        ):
            (self.bundle / "experiment.json").write_text(json.dumps(experiment))
            with mock.patch.object(helper.subprocess, "run") as process:
                with self.assertRaisesRegex(ValueError, "schema 2 fixed-temperature pixel-scale"):
                    helper.configure(self.bundle, self.machine)
                process.assert_not_called()
            self.assertFalse((self.bundle / "submit.sh").exists())

    def test_smaller_allocation_and_worker_probe_keep_the_same_frozen_contract(self):
        self.settings["slurm"]["max_parallel"] = 4
        self.machine.write_text(json.dumps(self.settings))
        self.configure_with_real_probe()
        self.assertIn("#SBATCH --array=0-4%4", (self.bundle / "job_array.sh").read_text())
        probe = helper.dependency_probe(str(self.executable), str(self.utilities), self.experiment)
        environment = {**os.environ, "MPLCONFIGDIR": str(self.bundle / ".mplconfig")}
        result = subprocess.run([sys.executable, "-B", "-c", probe], cwd=self.bundle,
                                env=environment, capture_output=True, text=True, check=True)
        receipt = json.loads(result.stdout)
        self.assertEqual(receipt["frozen_crossover_tasks_checked"], 5)
        self.assertEqual(set(receipt["packages"]), {"numpy", "scipy", "astropy", "matplotlib"})
        self.assertNotIn("frozen_subruns_checked", receipt)
        self.assertTrue(all(value["file"] for value in receipt["packages"].values()))
        (self.bundle / "frozen_temperature").write_bytes(b"changed after configuration")
        result = subprocess.run([sys.executable, "-B", "-c", probe], cwd=self.bundle,
                                env=environment, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Frozen crossover temperature mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
