"""Per-model seeds reach every simulator command and cannot share caches."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import shutil
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "src"))
import test_workflow
from mcfost_grid import runner


class ModelNumericsTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_workflow.WorkflowTests(methodName="runTest")
        self.fixture.setUp()

    def tearDown(self):
        self.fixture.tearDown()

    def seeded_manifest(self, replica=True):
        fixture = self.fixture
        manifest = deepcopy(fixture.manifest)
        manifest["configuration"]["numerics"]["random_seed"] = 43001
        manifest["models"][0]["numerics"] = {**manifest["configuration"]["numerics"], "random_seed": 43002}
        if replica:
            base = deepcopy(manifest["models"][0])
            base.update(id="m0001_seed43001", index=1)
            base["numerics"]["random_seed"] = 43001
            manifest["models"].append(base)
            target = fixture.run / "models" / base["id"]
            shutil.copytree(fixture.directory, target)
            for path in target.rglob("*"):
                if path.is_file():
                    manifest["input_hashes"][str(path.relative_to(fixture.run))] = runner.sha256(path)
        (fixture.run / "manifest.json").write_text(json.dumps(manifest))
        return manifest

    def test_replica_seed_reaches_fresh_temperature_and_both_image_backends(self):
        fixture = self.fixture
        manifest = self.seeded_manifest()
        for backend in ("image_method2", "coeval_method2"):
            with self.subTest(backend=backend):
                # Independent physical runs share the same global seed but
                # preserve their own actual CLI seed and product directory.
                if backend == "coeval_method2":
                    fixture.tearDown()
                    fixture.setUp()
                    manifest = self.seeded_manifest()
                fixture.machine_value["backend"] = backend
                fixture.machine.write_text(json.dumps(fixture.machine_value))
                with fixture.patches():
                    for index, seed in ((0, 43002), (1, 43001)):
                        before = len(fixture.commands)
                        self.assertEqual(runner.run_model(fixture.run, index, fixture.machine)["status"], "complete")
                        commands = fixture.commands[before:]
                        self.assertEqual(len(commands), 3)
                        self.assertEqual(commands[0][0][1], "temperature.para")
                        for command, _ in commands:
                            self.assertEqual(command.count("-seed"), 1)
                            self.assertEqual(command[command.index("-seed") + 1], str(seed))
                        self.assertEqual(runner.run_model(fixture.run, index, fixture.machine)["status"], "cached")
                        self.assertEqual(len(fixture.commands), before + 3)
                        directory = fixture.run / "models" / manifest["models"][index]["id"]
                        for name in ("runtime_binding", "temperature_complete", "measurements"):
                            receipt = json.loads((directory / f"{name}.json").read_text())
                            self.assertEqual(receipt["numerics"], manifest["models"][index]["numerics"])
                        temperature = json.loads((directory / "temperature_complete.json").read_text())
                        self.assertIn(f"/seed={seed}/", temperature["path"])
                # Seed differences are allowed within one grid runtime.
                shared = json.loads((fixture.run / "runtime_binding.json").read_text())["fingerprint"]
                for model in manifest["models"]:
                    receipt = json.loads((fixture.run / "models" / model["id"] / "measurements.json").read_text())
                    self.assertEqual(receipt["fingerprint"], shared)

    def test_incomplete_or_nonseed_overrides_rejected_before_simulator(self):
        fixture = self.fixture
        original = self.seeded_manifest(replica=False)
        for override in ({"random_seed": 43002}, None,
                         {**original["models"][0]["numerics"], "photons_image": 999},
                         {**original["models"][0]["numerics"], "random_seed": True},
                         {**original["models"][0]["numerics"], "random_seed": 0}):
            with self.subTest(override=override):
                manifest = deepcopy(original)
                manifest["models"][0]["numerics"] = override
                (fixture.run / "manifest.json").write_text(json.dumps(manifest))
                with fixture.patches(), self.assertRaises(ValueError):
                    runner.run_model(fixture.run, 0, fixture.machine)
                self.assertFalse(fixture.commands)

    def test_cross_seed_temperature_measurement_and_runtime_receipts_are_rejected(self):
        fixture = self.fixture
        manifest = self.seeded_manifest()
        with fixture.patches():
            for index in (0, 1):
                runner.run_model(fixture.run, index, fixture.machine)
            source = fixture.run / "models" / manifest["models"][1]["id"]
            target = fixture.directory
            for name in ("temperature_complete.json", "measurements.json", "runtime_binding.json"):
                original = (target / name).read_bytes()
                (target / name).write_bytes((source / name).read_bytes())
                with self.subTest(receipt=name), self.assertRaisesRegex(RuntimeError, "numerical|numerics|seed"):
                    runner.run_model(fixture.run, 0, fixture.machine)
                (target / name).write_bytes(original)
            self.assertEqual(len(fixture.commands), 6)

    def test_changing_model_seed_changes_manifest_fingerprint_and_rejects_cache(self):
        fixture = self.fixture
        manifest = self.seeded_manifest(replica=False)
        with fixture.patches():
            runner.run_model(fixture.run, 0, fixture.machine)
            manifest["models"][0]["numerics"]["random_seed"] = 43003
            (fixture.run / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "runtime or manifest changed"):
                runner.run_model(fixture.run, 0, fixture.machine)
            self.assertEqual(len(fixture.commands), 3)

    def test_legacy_global_numerics_and_receipt_schema_are_unchanged(self):
        fixture = self.fixture
        with fixture.patches():
            runner.run_model(fixture.run, 0, fixture.machine)
            self.assertEqual(runner.run_model(fixture.run, 0, fixture.machine)["status"], "cached")
        self.assertTrue(all("-seed" not in command for command, _ in fixture.commands))
        for name in ("runtime_binding", "temperature_complete", "measurements"):
            self.assertNotIn("numerics", json.loads((fixture.directory / f"{name}.json").read_text()))
        self.assertEqual(runner.model_numerics(fixture.manifest, fixture.model),
                         fixture.manifest["configuration"]["numerics"])


if __name__ == "__main__":
    unittest.main()
