"""Production catalogue, frozen-input isolation and Slurm handoff contracts."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
from build_extinction_ice_production import catalogue, ice_volume_fraction
from extinction_ice_production_task import digest, validate_package, NAME
from configure_extinction_ice_numerics_v2 import launch_scripts, validate_settings
from mcfost_grid.configuration import atomic_json, expand_models


class ProductionCatalogueTests(unittest.TestCase):
    def test_previous_96_are_preserved_and_ice_axis_is_controlled(self):
        parent = json.loads((ROOT / "config/grid.continuum-production-v1.1.json").read_text())
        records = catalogue(parent)
        self.assertEqual(len(records), 144)
        self.assertEqual([r["parameters"] for r in records[:96]],
                         [m["parameters"] for m in expand_models(parent)])
        groups = {}
        for row in records[96:]:
            groups.setdefault(row["contrast_group_id"], []).append(row)
        self.assertEqual(len(groups), 12)
        for group in groups.values():
            self.assertEqual([r["ice_mass_fraction"] for r in group], [0., .02, .04, .08])
            self.assertEqual(group[0]["dust_family"], "bare_silicate_Mie")
            base = {k: v for k, v in group[0]["parameters"].items() if k != "envelope_ice_volume_fraction"}
            for row in group:
                self.assertEqual(base, {k: v for k, v in row["parameters"].items() if k != "envelope_ice_volume_fraction"})
                f = row["parameters"]["envelope_ice_volume_fraction"]
                self.assertAlmostEqual(f*1.2/(f*1.2+(1-f)*3.5), row["ice_mass_fraction"])

    def test_mass_fraction_is_not_silently_used_as_volume(self):
        self.assertEqual(ice_volume_fraction(0), 0)
        self.assertAlmostEqual(ice_volume_fraction(.04), .10835913312693499)
        self.assertNotEqual(ice_volume_fraction(.04), .04)

    def test_array_caps_and_custom_afterany_analysis(self):
        experiment = {"schema_version": 1, "experiment_id": NAME,
                      "tasks": [{"index": i} for i in range(144)]}
        machine = dict(schema_version=1, mcfost_executable="/tmp/fake mcfost",
                       mcfost_utils="/tmp/utilities", slurm={"python": "/tmp/shared venv/bin/python"})
        validate_settings(machine, experiment)
        self.assertEqual(machine["max_memory_gb"], 112)
        self.assertEqual(machine["slurm"]["memory"], "160G")
        scripts = launch_scripts(Path("/tmp/run space"), experiment, machine, Path("/tmp/machine space.json"))
        array, analysis, submit = ("\n".join(scripts[name]) for name in ("job_array.sh", "job_analysis.sh", "submit.sh"))
        self.assertIn("--array=0-143%16", array)
        self.assertIn("--cpus-per-task=64", array)
        self.assertIn("code/production_task.py", array)
        self.assertIn("code/analyze_extinction_ice_production.py", analysis)
        self.assertIn('--dependency="afterany:$job"', submit)
        self.assertNotIn("final_resolution", submit)
        for lines in scripts.values():
            subprocess.run(["bash", "-n"], input="\n".join(lines)+"\n", text=True, check=True)

    def test_wrong_memory_is_rejected_before_worker_allocation(self):
        with self.assertRaisesRegex(ValueError, "112"):
            validate_settings(dict(schema_version=1, max_memory_gb=12), {"experiment_id": NAME})


class ProductionIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = Path(self.tmp.name)
        models = [dict(index=i, id=f"m{i}", parameters={"distance_pc": 140.}) for i in range(2)]
        probes = [{"id": "n1", "wavelength_um": 3., "score": False}]
        experiment = dict(schema_version=1, experiment_id=NAME,
                          tasks=[dict(index=m["index"], model_id=m["id"]) for m in models],
                          catalogue=[dict(index=m["index"], parameters=m["parameters"]) for m in models])
        for name in ("experiment.json", "production_experiment.json"):
            atomic_json(self.run / name, experiment)
        (self.run / "experiment.sha256").write_text(digest(self.run / "experiment.json")+"  experiment.json\n")
        atomic_json(self.run / "inputs/production_observation_contract.json", {"probes": probes})
        for model in models:
            path = self.run / "models" / model["id"] / "temperature.para"
            path.parent.mkdir(parents=True)
            path.write_text("immutable model parameter")
        manifest = dict(models=models, anchors=probes,
                        configuration=dict(numerical_only=True, numerics=dict(
                            photons_temperature=2048000, photons_image=2048000, photons_sed=2048000,
                            image_npix=6001, image_size_au=6000.)),
                        measurement=dict(aperture_radius_arcsec=1., target_distance_pc=147.,
                                         aperture_subpixels=64, quality_policy="aperture_v2"),
                        input_hashes={str(p.relative_to(self.run)): digest(p) for p in self.run.rglob("*") if p.is_file()})
        self.write_manifest(manifest)

    def tearDown(self):
        self.tmp.cleanup()

    def write_manifest(self, manifest):
        atomic_json(self.run / "manifest.json", manifest)
        (self.run / "manifest.sha256").write_text(digest(self.run / "manifest.json")+"  manifest.json\n")

    def test_unrelated_model_is_deferred_but_shared_inputs_are_always_checked(self):
        validate_package(self.run)
        (self.run / "models/m1/temperature.para").write_text("changed")
        validate_package(self.run, 0)
        for index in (None, 1):
            with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                validate_package(self.run, index)
        (self.run / "inputs/production_observation_contract.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "Frozen input changed"):
            validate_package(self.run, 0)

    def test_bad_index_manifest_and_path_never_reach_simulator(self):
        for index in (-1, 2, True, "0"):
            with self.assertRaisesRegex(ValueError, "index"):
                validate_package(self.run, index)
        manifest = json.loads((self.run / "manifest.json").read_text())
        manifest["input_hashes"]["models/m1/../../../outside"] = "0"*64
        self.write_manifest(manifest)
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            validate_package(self.run, 0)
        (self.run / "manifest.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "checksum"):
            validate_package(self.run)


if __name__ == "__main__":
    unittest.main()
