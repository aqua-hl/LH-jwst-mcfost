"""Freeze per-model material identity without implicit or missing dust inputs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from mcfost_grid.configuration import expand_models, prepare_run
from silicate_structure_design import catalogue, ICE_NAME


class MaterialConfigurationTests(unittest.TestCase):
    def test_material_is_preserved_in_identity_and_changes_model_hash(self):
        base = catalogue()[0]["parameters"]
        other = {**base, "envelope_silicate_file": "Olivine_MgFeSiO4_Dorschner1995_mcfost.dat"}
        models = expand_models({"models": [base, other]})
        self.assertEqual(models[0]["parameters"], base)
        self.assertEqual(models[1]["parameters"], other)
        self.assertNotEqual(models[0]["id"].split("_")[1], models[1]["id"].split("_")[1])
        with self.assertRaisesRegex(ValueError, "Duplicate physical model"):
            expand_models({"models": [other, other]})

    def test_all_sixty_records_survive_freezing_without_coercing_material_names(self):
        records = catalogue()
        models = expand_models({"models": [r["parameters"] for r in records]})
        self.assertEqual([m["parameters"] for m in models], [r["parameters"] for r in records])
        self.assertEqual(len({m["id"] for m in models}), 60)

    def test_missing_selected_material_fails_before_any_run_is_created(self):
        lab_model = next(r["parameters"] for r in catalogue() if r["role"] == "material")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = dict(schema_version=1, run_name="missing_material", output_dir=str(root),
                          template=str(ROOT / "reference/parameters/ice_v02_dust.para"),
                          models=[lab_model], dust_files={ICE_NAME: str(ROOT / "reference/dust" / ICE_NAME)})
            path = root / "grid.json"
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "Models require supplied dust_files"):
                prepare_run(path)
            self.assertFalse((root / "missing_material").exists())


if __name__ == "__main__":
    unittest.main()
