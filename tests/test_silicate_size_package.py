"""Four-model freezing, input gates and cluster scripts; never run MCFOST.

Default packages preserve the actual archived H2O table, with explicit native
extrapolation acceptance. The temporary H2O extension below is deliberately
synthetic. It exercises the optional full-range coverage/provenance contract
only and is not a physical material or usable scientific extension.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from build_silicate_size_pilot import build_package
from configure_extinction_ice_numerics_v2 import bundle_layout, configure, launch_scripts, validate_settings
from extinction_ice_production_task import validate_package
from mcfost_grid.configuration import expand_models
from silicate_pilot_design import (ICE_NAME, LAB_SHA256, NAME, NUMERICS, catalogue,
                                   digest, optical_table, plan, validate_ice)

LAB = ROOT / "reference/dust" / ICE_NAME
HISTORICAL_PROVENANCE = ROOT / "silicate/constants/H2O_30K_historical.provenance.json"


def fixture(directory):
    """Add two test-only endpoint rows while retaining all actual lab samples."""
    lines = LAB.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith("#"))
    table = directory / "synthetic_test_only.dat"
    table.write_text("\n".join([
        "# SYNTHETIC TEST FIXTURE: NOT PHYSICAL OPTICAL CONSTANTS",
        *lines[:header + 1], "0.1 1.3 1e-8", *lines[header + 1:], "3000 1.3 1e-8", ""]))
    provenance = directory / "synthetic_test_only.provenance.json"
    metadata = dict(schema_version=1, table_sha256=digest(table), laboratory_sha256=LAB_SHA256,
                    laboratory_temperature_k=30., measured_region_treatment="preserved_original_samples",
                    sources=["Synthetic unit-test fixture, not a scientific material"],
                    extension_method="Two artificial endpoint rows for input-gate testing only",
                    join_description="Original laboratory samples retained verbatim",
                    tail_temperature_assumptions="No physical temperature assumption: test fixture only")
    provenance.write_text(json.dumps(metadata))
    return table, provenance


def fields(text, tag):
    return [line.split(tag)[0].split() for line in text.splitlines() if tag in line]


@contextmanager
def replaced(path, content):
    previous = path.read_bytes()
    path.write_bytes(content if isinstance(content, bytes) else content.encode())
    try:
        yield
    finally:
        path.write_bytes(previous)


class SilicateInputGateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="test-only-silicate-inputs-")
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.table, self.provenance = fixture(self.directory)

    def test_plan_accepts_historical_policy_but_still_requires_preparation_and_configuration(self):
        design = plan()
        self.assertTrue(design["input_policy_ready"])
        self.assertFalse(design["launch_ready"])
        self.assertEqual(design["temperature_solves"], 4)
        self.assertEqual(design["image_requests"], 392)
        rows = catalogue()
        self.assertEqual([(r["parameters"]["inclination_deg"], r["parameters"]["envelope_amax_um"])
                          for r in rows], [(50, .4), (50, 1), (70, .4), (70, 1)])
        for row in rows:
            self.assertEqual(row["parameters"]["envelope_dust_mass_msun"], .000225)
            self.assertEqual(row["parameters"]["envelope_ice_mass_fraction"], .04)
            self.assertNotIn("envelope_ice_volume_fraction", row["parameters"])

    def test_strict_full_range_policy_rejects_short_or_missing_table_before_preparation(self):
        destination = self.directory / "must_not_exist"
        for table in (LAB, self.directory / "absent.dat"):
            metadata = json.loads(self.provenance.read_text())
            if table.exists():
                metadata["table_sha256"] = digest(table)
            self.provenance.write_text(json.dumps(metadata))
            with self.subTest(table=table), self.assertRaises(ValueError):
                build_package(destination, table, self.provenance)
            self.assertFalse(destination.exists())

    def test_full_range_supplied_option_remains_available_and_reports_no_extrapolation(self):
        receipt = validate_ice(self.table, self.provenance, LAB)
        self.assertTrue(receipt["coverage_passed"])
        self.assertTrue(receipt["policy_accepted"])
        self.assertFalse(receipt["native_extrapolation_required"])
        bundle = self.directory / "strict_full_range_test_only"
        build_package(bundle, self.table, self.provenance)
        self.assertEqual((bundle / "inputs/utils/Dust" / ICE_NAME).read_bytes(), self.table.read_bytes())
        validate_package(bundle)

    def test_historical_policy_requires_explicit_acceptance_and_the_exact_archived_table(self):
        original = json.loads(HISTORICAL_PROVENANCE.read_text())
        for field, value in (("coverage_policy", "full_range_supplied"),
                             ("native_extrapolation_accepted_for_pilot", False),
                             ("zero_k_policy", "floor_zero_samples"),
                             ("table_sha256", "0" * 64)):
            metadata = self.directory / "historical.provenance.json"
            metadata.write_text(json.dumps({**original, field: value}))
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_ice(LAB, metadata, LAB)
        metadata.write_text(json.dumps({k: v for k, v in original.items() if k != "coverage_policy"}))
        with self.assertRaises(ValueError):
            validate_ice(LAB, metadata, LAB)
        changed = self.directory / "changed_short_table.dat"
        changed.write_bytes(LAB.read_bytes() + b"\n# Different archived input\n")
        metadata.write_text(json.dumps({**original, "table_sha256": digest(changed)}))
        with self.assertRaises(ValueError):
            validate_ice(changed, metadata, LAB)

    def test_provenance_is_required_and_binds_exact_material_and_reference(self):
        with self.assertRaisesRegex(ValueError, "Missing H2O .*provenance"):
            validate_ice(self.table, self.directory / "missing.json", LAB)
        original = json.loads(self.provenance.read_text())
        for field, value, message in (
            ("table_sha256", "0" * 64, "table_sha256"),
            ("laboratory_sha256", "0" * 64, "reference hash"),
            ("laboratory_temperature_k", 100, "30 K"),
            ("measured_region_treatment", "smoothed", "preserved original"),
            ("extension_method", "TODO", "placeholder"),
            ("sources", [], "placeholder"),
        ):
            with self.subTest(field=field), replaced(self.provenance, json.dumps({**original, field: value})):
                with self.assertRaisesRegex(ValueError, message):
                    validate_ice(self.table, self.provenance, LAB)

    def test_invalid_numeric_tables_and_changed_measured_region_are_rejected(self):
        text = self.table.read_text()
        for changed, message in (
            (text.replace("0.1 1.3 1e-8", "0.1 1.3 -1"), "Invalid wavelength"),
            (text.replace("0.1 1.3 1e-8", "0.1 nan 1e-8"), "Nonfinite"),
            (text.replace("0.1 1.3 1e-8", "0.1 1.3 1e-8\n0.1 1.3 1e-8"), "strictly increasing"),
            (text.replace("0.1 1.3 1e-8", "0.1 1.3"), "three columns"),
        ):
            with self.subTest(message=message), replaced(self.table, changed):
                with self.assertRaisesRegex(ValueError, message):
                    optical_table(self.table)
        header, samples = optical_table(self.table)
        samples[10][1] += .01
        changed = " ".join(map(str, header)) + "\n" + "\n".join(" ".join(map(str, row)) for row in samples) + "\n"
        with replaced(self.table, changed):
            metadata = json.loads(self.provenance.read_text())
            metadata["table_sha256"] = digest(self.table)
            with replaced(self.provenance, json.dumps(metadata)):
                with self.assertRaisesRegex(ValueError, "changes the archived measured region"):
                    validate_ice(self.table, self.provenance, LAB)

    def test_mass_fraction_omits_only_the_implicit_volume_default(self):
        configs = [dict(models=[{"envelope_ice_mass_fraction": .04}]),
                   dict(fixed={"envelope_ice_mass_fraction": .04}),
                   dict(grid={"envelope_ice_mass_fraction": [.04]})]
        for config in configs:
            with self.subTest(config=config):
                params = expand_models(config)[0]["parameters"]
                self.assertEqual(params["envelope_ice_mass_fraction"], .04)
                self.assertNotIn("envelope_ice_volume_fraction", params)
                with self.assertRaisesRegex(ValueError, "not both"):
                    expand_models({**config, "fixed": {**config.get("fixed", {}), "envelope_ice_volume_fraction": .1}})
        self.assertIn("envelope_ice_volume_fraction", expand_models({})[0]["parameters"])


class SilicateFrozenPackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="test-only-silicate-package-")
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.directory = Path(cls.tmp.name)
        cls.table, cls.provenance = LAB, HISTORICAL_PROVENANCE
        cls.bundle = cls.directory / NAME
        with mock.patch("subprocess.run", side_effect=AssertionError("Preparation must not execute a simulator")):
            cls.result = build_package(cls.bundle)
        cls.manifest = json.loads((cls.bundle / "manifest.json").read_text())
        cls.experiment = json.loads((cls.bundle / "experiment.json").read_text())

    def test_four_independent_temperature_models_and_98_probes_are_frozen(self):
        self.assertEqual(self.result["models"], 4)
        self.assertEqual(self.result["temperature_solves"], 4)
        self.assertEqual(self.result["probes_per_model"], 98)
        self.assertFalse(self.result["mcfost_invoked"])
        self.assertFalse(self.result["jobs_submitted"])
        self.assertEqual(len(self.manifest["models"]), 4)
        self.assertEqual(len(self.manifest["anchors"]), 98)
        self.assertEqual(self.manifest["temperature_policy"], "fresh equilibrium for every physical model")
        self.assertEqual(self.manifest["configuration"]["numerics"], NUMERICS)
        self.assertEqual(self.manifest["measurement"]["quality_policy"], "aperture_v2")
        temperatures = []
        for model in self.manifest["models"]:
            folder = self.bundle / "models" / model["id"]
            temperature = (folder / "temperature.para").read_text()
            temperatures.append(temperature)
            self.assertEqual(fields(temperature, "compute temperature?"), [["T", "F", "T"]])
            self.assertEqual(fields(temperature, "amin, amax")[-1], ["0.03", "0.4", "2.75", "50"])
            self.assertEqual(len(list((folder / "anchors").glob("*/image.para"))), 98)
            self.assertEqual(len(list((folder / "anchors").glob("*/coeval.para"))), 98)
            self.assertFalse(list(folder.rglob("Temperature.fits*")))
        for a, b in ((0, 1), (2, 3)):
            different = [(x, y) for x, y in zip(temperatures[a].splitlines(), temperatures[b].splitlines()) if x != y]
            self.assertEqual(len(different), 1)
            self.assertIn("amin, amax", different[0][0])
            self.assertIn("0.03  0.4  2.75  50", different[0][0])
            self.assertIn("0.03  1  2.75  50", different[0][1])
        validate_package(self.bundle)

    def test_default_archived_material_is_byte_preserved_with_an_honest_coverage_receipt(self):
        self.assertEqual(digest(self.table), LAB_SHA256)
        frozen = self.bundle / "inputs/utils/Dust" / ICE_NAME
        self.assertEqual(frozen.read_bytes(), LAB.read_bytes())
        header, samples = optical_table(frozen)
        self.assertEqual(header, [.94, 150.])
        self.assertEqual(sum(row[2] == 0 for row in samples), 206)
        receipt = self.experiment["ice_input_check"]
        self.assertEqual(receipt["coverage_policy"], "historical_native_extrapolation")
        self.assertFalse(receipt["coverage_passed"])
        self.assertTrue(receipt["policy_accepted"])
        self.assertTrue(receipt["native_extrapolation_required"])
        self.assertFalse(receipt["tail_physical_accuracy_certified"])
        self.assertEqual(receipt["zero_k_rows"], 206)

    def test_all_frozen_parameter_files_use_512k_photons(self):
        keys = ("photons_temperature", "photons_sed", "photons_image")
        for numerics in (NUMERICS, self.manifest["configuration"]["numerics"],
                         self.experiment["numerics"]):
            self.assertEqual({key: numerics[key] for key in keys},
                             {key: 512000 for key in keys})
        tags = ("nbr_photons_eq_th", "nbr_photons_lambda", "nbr_photons_image")
        for model in self.manifest["models"]:
            folder = self.bundle / "models" / model["id"]
            paths = [folder / "temperature.para", *sorted(folder.glob("anchors/*/image.para")),
                     *sorted(folder.glob("anchors/*/coeval.para"))]
            self.assertEqual(len(paths), 197)
            for path in paths:
                text = path.read_text()
                for tag in tags:
                    with self.subTest(model=model["id"], file=path.relative_to(folder), tag=tag):
                        values = fields(text, tag)
                        self.assertEqual(len(values), 1)
                        self.assertEqual(float(values[0][0]), 512000)

    def test_resealed_manifest_cannot_restore_previous_2m_pilot_photons(self):
        manifest_path = self.bundle / "manifest.json"
        checksum_path = self.bundle / "manifest.sha256"
        for key in ("photons_temperature", "photons_sed", "photons_image"):
            manifest = json.loads(manifest_path.read_text())
            manifest["configuration"]["numerics"][key] = 2048000
            with self.subTest(key=key), replaced(manifest_path, json.dumps(manifest)):
                with replaced(checksum_path, digest(manifest_path) + "  manifest.json\n"):
                    with self.assertRaisesRegex(ValueError, "Photon count differs"):
                        validate_package(self.bundle, 0)

    def test_supplied_material_and_provenance_are_frozen_and_tamper_checked(self):
        for relative, source in (("inputs/utils/Dust/" + ICE_NAME, self.table),
                                 ("inputs/ice_input_provenance.json", self.provenance),
                                 ("inputs/H2O30K_laboratory_reference.dat", LAB)):
            with self.subTest(relative=relative):
                frozen = self.bundle / relative
                self.assertEqual(digest(frozen), digest(source))
                self.assertEqual(self.manifest["input_hashes"][relative], digest(source))
                with replaced(frozen, frozen.read_bytes() + b"\n# changed\n"):
                    with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                        validate_package(self.bundle, 0)

    def test_worker_checks_shared_inputs_and_its_own_model_only(self):
        other = self.bundle / "models" / self.manifest["models"][1]["id"] / "temperature.para"
        with replaced(other, "changed model input\n"):
            validate_package(self.bundle, 0)
            for index in (None, 1):
                with self.subTest(index=index), self.assertRaisesRegex(ValueError, "Frozen input changed"):
                    validate_package(self.bundle, index)
            shared = self.bundle / "inputs/production_observation_contract.json"
            with replaced(shared, "{}"):
                with self.assertRaisesRegex(ValueError, "Frozen input changed"):
                    validate_package(self.bundle, 0)

    def test_cluster_array_and_afterany_paired_analysis(self):
        self.assertEqual(bundle_layout(self.bundle)["experiment_id"], NAME)
        machine = dict(schema_version=1, mcfost_executable="/tmp/fake mcfost",
                       mcfost_utils="/tmp/fake utilities", slurm={"python": "/tmp/shared venv/bin/python"})
        validate_settings(machine, self.experiment)
        self.assertEqual(machine["threads"], 64)
        self.assertEqual(machine["max_memory_gb"], 112)
        self.assertEqual(machine["slurm"]["memory"], "160G")
        self.assertEqual(machine["slurm"]["max_parallel"], 4)
        scripts = launch_scripts(self.bundle, self.experiment, machine, self.directory / "machine.cluster.json")
        text = {name: "\n".join(lines) + "\n" for name, lines in scripts.items()}
        self.assertIn("--array=0-3%4", text["job_array.sh"])
        self.assertIn("--cpus-per-task=64", text["job_array.sh"])
        self.assertIn("code/production_task.py", text["job_array.sh"])
        self.assertIn("code/analyze_silicate_size_pilot.py", text["job_analysis.sh"])
        self.assertIn('--dependency="afterany:$job"', text["submit.sh"])
        self.assertIn("SLURM_ARRAY_TASK_ID", text["job_array.sh"])
        for script in text.values():
            subprocess.run(["bash", "-n"], input=script, text=True, check=True)
        machine["slurm"]["max_parallel"] = 5
        with self.assertRaisesRegex(ValueError, "four concurrent"):
            validate_settings(machine, self.experiment)

    def test_frozen_cli_and_configuration_use_the_copied_code_without_simulation(self):
        result = subprocess.run([sys.executable, "-B", str(self.bundle / "code/production_task.py"),
                                 str(self.bundle), "--validate"],
                                capture_output=True, text=True, check=True, cwd=self.directory)
        self.assertEqual(json.loads(result.stdout)["models"], 4)
        # The configuration probe may stat this file, but must never execute it.
        executable = self.directory / "must-not-run-mcfost"
        marker = self.directory / "simulator-was-invoked"
        executable.write_text("#!/bin/sh\ntouch '" + str(marker) + "'\nexit 97\n")
        executable.chmod(0o755)
        utils = self.directory / "utils"
        for name in ("Dust", "Lambda", "Stellar_Spectra"):
            (utils / name).mkdir(parents=True)
        machine = json.loads((self.bundle / "machine.template.json").read_text())
        machine.update(mcfost_executable=str(executable), mcfost_utils=str(utils))
        source = self.directory / "machine.json"
        source.write_text(json.dumps(machine))
        configured = configure(self.bundle, source)
        self.assertTrue(configured.is_file())
        self.assertFalse(marker.exists())
        self.assertIn("code/analyze_silicate_size_pilot.py", (self.bundle / "job_analysis.sh").read_text())


if __name__ == "__main__":
    unittest.main()
