"""Per-task validation retains shared/selected integrity without all-grid reads."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid import runner


class SelectedInputValidationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.run = Path(self.temporary.name).resolve()
        self.shared = ["code/src/mcfost_grid/runner.py", "inputs/utils/Dust/ice_opct.dat",
                       "inputs/configuration.json", "inputs/production_observation_contract.json",
                       "experiment.json"]
        self.paths = self.shared + ["models/m0000_a/temperature.para",
                                   "models/m0000_a/anchors/c01_q0/image.para",
                                   "models/m0001_b/temperature.para"]
        for relative in self.paths:
            path = self.run / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(relative)
        self.manifest = {"models": [{"id": "m0000_a"}, {"id": "m0001_b"}],
                         "input_hashes": {p: runner.sha256(self.run / p) for p in self.paths}}

    def test_selected_validation_hashes_all_shared_and_only_its_model_files(self):
        with patch.object(runner, "sha256", wraps=runner.sha256) as hashed:
            runner.validate_inputs(self.run, self.manifest, model_id="m0000_a")
        paths = {str(call.args[0].relative_to(self.run)) for call in hashed.call_args_list}
        self.assertEqual(paths, set(self.paths) - {"models/m0001_b/temperature.para"})
        with patch.object(runner, "sha256", wraps=runner.sha256) as hashed:
            runner.validate_inputs(self.run, self.manifest)
        self.assertEqual(hashed.call_count, len(self.paths))

    def test_every_shared_input_and_selected_parameter_mutation_still_fails(self):
        for relative in self.shared + ["models/m0000_a/anchors/c01_q0/image.para"]:
            with self.subTest(relative=relative):
                path = self.run / relative
                original = path.read_text()
                path.write_text("changed")
                try:
                    with self.assertRaisesRegex(RuntimeError, "Frozen run input changed"):
                        runner.validate_inputs(self.run, self.manifest, model_id="m0000_a")
                finally:
                    path.write_text(original)

    def test_unrelated_mutation_is_detected_by_its_own_task_and_full_audit(self):
        (self.run / "models/m0001_b/temperature.para").write_text("changed")
        runner.validate_inputs(self.run, self.manifest, model_id="m0000_a")
        with self.assertRaisesRegex(RuntimeError, "Frozen run input changed"):
            runner.validate_inputs(self.run, self.manifest, model_id="m0001_b")
        with self.assertRaisesRegex(RuntimeError, "Frozen run input changed"):
            runner.validate_inputs(self.run, self.manifest)

    def test_unknown_or_unsafe_selected_ids_are_rejected(self):
        for model_id in ("m9999_unknown", "../m0000_a", "/m0000_a", "", True):
            with self.subTest(model_id=model_id), self.assertRaises(ValueError):
                runner.validate_inputs(self.run, self.manifest, model_id=model_id)

    def test_unsafe_or_unknown_other_model_paths_are_not_skipped(self):
        for relative in ("models/m0001_b/../../escape", "/absolute/escape", ".",
                         "models/m9999_unknown/temperature.para"):
            with self.subTest(relative=relative):
                manifest = {**self.manifest, "input_hashes": {**self.manifest["input_hashes"], relative: "0" * 64}}
                with self.assertRaises(RuntimeError):
                    runner.validate_inputs(self.run, manifest, model_id="m0000_a")

    def test_shared_symlink_cannot_substitute_even_identical_input_bytes(self):
        path = self.run / self.shared[0]
        target = self.run / "identical.py"
        target.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(target)
        with self.assertRaisesRegex(RuntimeError, "Frozen run input changed"):
            runner.validate_inputs(self.run, self.manifest, model_id="m0000_a")


if __name__ == "__main__":
    unittest.main()
