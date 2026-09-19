"""Recovery submission tests use fake scheduler responses and never launch jobs."""
from contextlib import redirect_stderr, redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("retry_slurm_models", ROOT / "scripts/retry_slurm_models.py")
retry = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(retry)


class RetryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.run = self.root / "continuum_production_v1.1"
        self.run.mkdir()
        (self.run / "manifest.json").write_text(json.dumps({
            "run_id": self.run.name, "models": [{"index": i} for i in range(96)],
        }))
        self.machine = self.root / "machine.json"
        self.machine.write_text(json.dumps({"slurm": {"python": "python", "max_parallel": 16}}))
        for name, directive, operation in (
            ("job_array.sh", "#SBATCH --array=0-95%16\n#SBATCH --cpus-per-task=64", 'task "$MCFOST_RUN_DIR"'),
            ("job_analysis.sh", "#SBATCH --cpus-per-task=1", 'analyze "$MCFOST_RUN_DIR"'),
        ):
            (self.run / name).write_text(
                f'#!/bin/bash\n{directive}\nset -euo pipefail\nmodule load Example\n'
                f'cd "$MCFOST_RUN_DIR"\nexec python -B code/workflow.py {operation}\n'
            )
        self.python = self.root / "shared-venv/bin/python"
        self.python.parent.mkdir(parents=True)
        self.python.symlink_to("/usr/bin/python3")
        self.queue = "10514475\n10514475_0\n10514476\n10514476\n"
        self.calls = []
        self.sbatch = []
        self.analysis_fails = False

    def invoke(self, command, **kwargs):
        self.calls.append((command, kwargs))
        stdout = ""
        if command[0] == "squeue":
            stdout = self.queue
        elif command[0] == "sbatch":
            self.sbatch.append(command)
            if len(self.sbatch) == 2 and self.analysis_fails:
                raise subprocess.CalledProcessError(1, command, stderr="scheduler temporarily unavailable")
            stdout = f"{10514500 + len(self.sbatch) - 1};cluster\n"
        elif command[0] not in (str(self.python), "bash"):
            self.fail(f"Unexpected external command: {command!r}")
        return subprocess.CompletedProcess(command, 0, stdout=stdout, stderr="")

    def execute(self, *options):
        with mock.patch.object(retry.sys, "executable", str(self.python)), \
                mock.patch.object(retry.subprocess, "run", side_effect=self.invoke), \
                redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return retry.main([str(self.run), "--machine", str(self.machine), *options])

    def receipt(self):
        receipts = list((self.run / "recovery").glob("*/recovery.json"))
        self.assertEqual(len(receipts), 1)
        return receipts[0], json.loads(receipts[0].read_text())

    def test_83_retries_wait_for_original_array_and_analysis_waits_for_all(self):
        self.execute("--indices", "8-10,16-95", "--after-job", "10514475",
                     "--after-job", "10514475", "--submit")
        self.assertEqual(len(self.sbatch), 2)
        array, analysis = self.sbatch
        self.assertIn("--dependency=afterany:10514475", array)
        dependency = next(arg for arg in analysis if arg.startswith("--dependency="))
        jobs = dependency.removeprefix("--dependency=afterany:").split(":")
        self.assertEqual(jobs, ["10514475", "10514475_0", "10514476", "10514500"])
        path, receipt = self.receipt()
        expected = [8, 9, 10, *range(16, 96)]
        self.assertEqual(receipt["indices"], expected)
        self.assertEqual(len(receipt["indices"]), 83)
        self.assertEqual(receipt["after_jobs"], ["10514475"])
        self.assertEqual(receipt["analysis_dependencies"], jobs)
        self.assertTrue(receipt["submitted"])
        self.assertEqual(receipt["analysis_job_id"], "10514501")
        generated = (path.parent / "job_array.sh").read_text()
        self.assertIn("#SBATCH --array=" + ",".join(map(str, expected)) + "%16", generated)
        self.assertIn("#SBATCH --cpus-per-task=64", generated)

    def test_dry_run_pins_venv_symlink_preserves_setup_and_probes_both_jobs(self):
        originals = {p: p.read_bytes() for p in (self.machine, self.run / "job_array.sh", self.run / "job_analysis.sh")}
        self.execute("--indices", "10,8-9,8", "--after-job", "10514475", "--after-job", "10514476")
        self.assertFalse(self.sbatch)
        self.assertFalse(any(command[0] == "squeue" for command, _ in self.calls))
        path, receipt = self.receipt()
        self.assertEqual(receipt["indices"], [8, 9, 10])
        self.assertEqual(receipt["after_jobs"], ["10514475", "10514476"])
        self.assertEqual(receipt["python"], str(self.python))
        self.assertFalse(receipt["submitted"])
        self.assertEqual(self.calls[0][0][0], str(self.python))
        for name in ("job_array.sh", "job_analysis.sh"):
            text = (path.parent / name).read_text()
            self.assertIn("module load Example", text)
            self.assertIn(f"exec {self.python} -B code/workflow.py", text)
            self.assertIn("numpy, scipy, astropy, matplotlib", text)
            self.assertIn("from mcfost_grid.photometry import measure_image", text)
            self.assertLess(text.index("module load Example"), text.index("from mcfost_grid.photometry"))
        self.assertIn("#SBATCH --array=8,9,10%3", (path.parent / "job_array.sh").read_text())
        for path, before in originals.items():
            self.assertEqual(path.read_bytes(), before)

    def test_invalid_selection_rejected_before_external_commands_or_recovery_creation(self):
        for indices in ("", ",", "8,", ",8", "8,,9", "8-", "-1", "8--10", "10-8", "8 9", "96", "8-96", "1-999999999999999999"):
            with self.subTest(indices=indices), self.assertRaises(SystemExit):
                self.execute("--indices", indices, "--submit")
            self.assertEqual(self.calls, [])
            self.assertFalse((self.run / "recovery").exists())

    def test_invalid_parent_job_rejected_before_external_commands(self):
        for job in ("", "0", "-1", "1.5", "10514475_0", "10514475,10514476", "abc"):
            with self.subTest(job=job), self.assertRaises(SystemExit):
                self.execute("--indices", "8", "--after-job", job, "--submit")
            self.assertEqual(self.calls, [])
            self.assertFalse((self.run / "recovery").exists())

    def test_analysis_submission_failure_keeps_retry_receipt_and_never_resubmits_array(self):
        self.analysis_fails = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.execute("--indices", "8-10,16-95", "--after-job", "10514475", "--submit")
        self.assertEqual(len(self.sbatch), 2)
        self.assertEqual(sum(command[-1].endswith("job_array.sh") for command in self.sbatch), 1)
        _, receipt = self.receipt()
        self.assertTrue(receipt["submitted"])
        self.assertEqual(receipt["retry_array_job_id"], "10514500")
        self.assertNotIn("analysis_job_id", receipt)
        self.assertEqual(receipt["analysis_submission_error"], "scheduler temporarily unavailable")
        self.assertIn("10514475", receipt["analysis_dependencies"])


if __name__ == "__main__":
    unittest.main()
