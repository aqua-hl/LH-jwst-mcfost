"""Launch dependencies and the numerical material gate for the 60-model run."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
sys.path.insert(0, str(ROOT/"src"))
import check_silicate_structure_materials as check
from configure_extinction_ice_numerics_v2 import (SILICATE_STRUCTURE, bundle_layout,
    launch_scripts, validate_settings)


class LaunchTests(unittest.TestCase):
    def test_preflight_dependency_before_array_and_afterany_analysis(self):
        with tempfile.TemporaryDirectory(prefix="structure-cluster-") as tmp:
            root = Path(tmp)
            experiment = dict(schema_version=1, experiment_id=SILICATE_STRUCTURE,
                              tasks=[dict(index=i) for i in range(60)])
            machine = dict(schema_version=1, mcfost_executable="/unused/mcfost",
                           mcfost_utils="/unused/utils", slurm={"python": sys.executable})
            validate_settings(machine, experiment)
            self.assertEqual(machine["threads"], 64)
            self.assertEqual(machine["max_memory_gb"], 112)
            self.assertEqual(machine["slurm"]["max_parallel"], 16)
            self.assertEqual(machine["slurm"]["memory"], "160G")
            configured = root/"machine.cluster.json"
            configured.write_text(json.dumps(machine, indent=2, allow_nan=False)+"\n")
            scripts = launch_scripts(root, experiment, machine, configured)
            for name, lines in scripts.items():
                script = "\n".join(lines)+"\n"
                subprocess.run(["bash", "-n"], input=script, text=True, check=True)
                (root/name).write_text(script)
            self.assertIn("#SBATCH --array=0-59%16", scripts["job_array.sh"])
            self.assertIn("#SBATCH --cpus-per-task=2", scripts["job_materials.sh"])
            self.assertIn("code/check_materials.py", "\n".join(scripts["job_materials.sh"]))
            self.assertIn("code/analyze_silicate_structure_production.py", "\n".join(scripts["job_analysis.sh"]))
            submit = "\n".join(scripts["submit.sh"])
            self.assertIn('--dependency="afterok:$materials_job" --kill-on-invalid-dep=yes', submit)
            self.assertIn('--dependency="afterany:$job"', submit)
            self.assertLess(submit.index("job_materials.sh"), submit.index("job_array.sh"))

    def test_layout_rejects_missing_gate_or_wrong_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"code").mkdir()
            experiment = dict(schema_version=1, experiment_id=SILICATE_STRUCTURE,
                              tasks=[dict(index=i) for i in range(60)])
            (root/"experiment.json").write_text(json.dumps(experiment))
            for name in ("manifest.json", "manifest.sha256", "production_experiment.json",
                         "code/production_task.py", "code/analyze_extinction_ice_production.py",
                         "code/analyze_silicate_structure_production.py"):
                (root/name).touch()
            with self.assertRaisesRegex(ValueError, "material preflight"):
                bundle_layout(root)
            (root/"code/check_materials.py").touch()
            self.assertEqual(bundle_layout(root), experiment)
            experiment["tasks"].pop()
            (root/"experiment.json").write_text(json.dumps(experiment))
            with self.assertRaisesRegex(ValueError, "60 tasks"):
                bundle_layout(root)


def write_outputs(root, wave, *, opacity=None, albedo=None, include_g=False):
    root = Path(root)/"data_dust"
    root.mkdir(exist_ok=True)
    ones = np.ones(len(wave))
    arrays = dict(lambda_=wave, kappa=ones if opacity is None else opacity,
                  albedo=.5*ones if albedo is None else albedo,
                  kappa_grain=np.ones((2, len(wave))), phase_function=np.ones((181, len(wave))),
                  polarizability=np.zeros((181, len(wave))))
    if include_g:
        arrays["g"] = .2*ones
    for name, values in arrays.items():
        fits.PrimaryHDU(np.asarray(values)).writeto(root/(name.rstrip("_")+".fits.gz"), overwrite=True)


def write_temperature_inputs(root, models):
    """Use real temperature-stage files, including their default-grid flags."""
    from mcfost_grid.physics import render_parameter
    from silicate_structure_design import NUMERICS, catalogue
    template = (ROOT/"reference/parameters/ice_v02_dust.para").read_text()
    defaults = catalogue()[0]["parameters"]
    for model in models:
        folder = Path(root)/"models"/model["id"]
        folder.mkdir(parents=True)
        text = render_parameter(template, model.get("parameters", defaults), NUMERICS, "temperature")
        (folder/"temperature.para").write_text(text)


class MaterialTests(unittest.TestCase):
    def test_all_materials_must_pass_before_receipt_can_release_workers(self):
        import shutil
        from silicate_structure_design import catalogue
        records = catalogue()
        unique = {}
        for row in records:
            p = row["parameters"]
            unique.setdefault((p["envelope_silicate_file"], p["envelope_amax_um"]), row)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"inputs/utils/Dust").mkdir(parents=True)
            shutil.copy2(ROOT/"reference/dust/H2O_30K_Leiden_mcfost.dat", root/"inputs/utils/Dust")
            models = [dict(id=f"m{i}", parameters=row["parameters"]) for i, row in enumerate(unique.values())]
            (root/"manifest.json").write_text(json.dumps(dict(models=models, anchors=[dict(wavelength_um=9.7)])))
            write_temperature_inputs(root, models)
            runner = SimpleNamespace(
                runtime_config=lambda _: {"mcfost_executable": "/fake/mcfost"},
                environment=lambda _, runtime: {"OMP_NUM_THREADS": str(runtime["threads"])},
                now=lambda: "2026-09-24T00:00:00+00:00",
                atomic_json=lambda path, payload: Path(path).write_text(json.dumps(payload)))
            dispatch = SimpleNamespace(validate_package=lambda _: {"experiment_id": SILICATE_STRUCTURE},
                                       load_runner=lambda _: runner)
            required_failure = {"enabled": False}
            def dust_command(command, **kwargs):
                self.assertIn("-dust_prop", command)
                self.assertNotIn("-img", command)
                self.assertNotIn("-root_dir", command)
                self.assertEqual(kwargs["timeout"], 120)
                self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "2")
                # MCFOST lect_lambda consumes every numeric row; there is no
                # count header to skip. Compare the file itself with the plan.
                wave = np.loadtxt(kwargs["cwd"]/"check.lambda")
                expected = check.wavelength_grid(root, json.loads((root/"manifest.json").read_text()))
                np.testing.assert_array_equal(wave, expected)
                write_outputs(kwargs["cwd"], wave)
                required = check.required_wavelength_grid(root, json.loads((root/"manifest.json").read_text()))
                matched = np.any(np.isclose(wave[:, None], required, rtol=2e-6, atol=0), axis=1)
                extra_index = np.flatnonzero(~matched)[0]
                opacity = np.ones(len(wave))
                opacity[extra_index] = np.nan
                if required_failure["enabled"] and kwargs["cwd"].parent.name == "m3":
                    opacity[np.flatnonzero(matched)[0]] = np.nan
                fits.PrimaryHDU(opacity).writeto(kwargs["cwd"]/"data_dust/kappa.fits.gz", overwrite=True)
                grain = np.ones((50, len(wave)))
                grain[:, extra_index] = np.nan
                fits.PrimaryHDU(grain).writeto(kwargs["cwd"]/"data_dust/kappa_grain.fits.gz", overwrite=True)
                kwargs["stdout"].write("Computing dust properties ... Writing dust properties\n Exiting\n")
                return SimpleNamespace(returncode=0)
            with patch.dict(sys.modules, {"production_task": dispatch}), \
                    patch.object(check, "identity", return_value={"runtime": "same"}), \
                    patch.object(check.subprocess, "run", side_effect=dust_command) as invocation:
                self.assertEqual(check.run_checks(root, "unused")["status"], "passed")
                self.assertEqual(invocation.call_count, 4)
                check.validate_receipt(root, {})
                self.assertEqual(check.run_checks(root, "unused")["status"], "cached")
                self.assertEqual(invocation.call_count, 4)
                receipt = json.loads((root/check.RECEIPT).read_text())
                for item in receipt["checks"]:
                    extra = item["diagnostics"]["extra_wavelength_diagnostics"]
                    self.assertEqual(extra["kappa"]["nonfinite_count"], 1)
                    self.assertEqual(extra["kappa_grain"]["nonfinite_count"], 50)
                receipt["status"] = "failed"
                (root/check.RECEIPT).write_text(json.dumps(receipt))
                with self.assertRaisesRegex(ValueError, "not passed"):
                    check.validate_receipt(root, {})
                required_failure["enabled"] = True
                with self.assertRaisesRegex(ValueError, "Nonfinite"):
                    check.run_checks(root, "unused")
                failed = json.loads((root/check.RECEIPT).read_text())
                self.assertEqual(failed["status"], "failed")
                self.assertEqual([item["status"] for item in failed["checks"][:3]], ["passed"]*3)
                with self.assertRaisesRegex(ValueError, "not passed"):
                    check.validate_receipt(root, {})
            with patch.dict(sys.modules, {"production_task": dispatch}), \
                    patch.object(check, "identity", return_value={"runtime": "same"}), \
                    patch.object(check.subprocess, "run", return_value=SimpleNamespace(returncode=-11)):
                with self.assertRaisesRegex(ValueError, "initialization failed"):
                    check.run_checks(root, "unused")
                receipt = json.loads((root/check.RECEIPT).read_text())
                self.assertEqual(receipt["status"], "failed")
                self.assertEqual(receipt["checks"][0]["returncode"], -11)
                self.assertEqual(len(list((root/"material_preflight"/"m0").glob("attempt_*"))), 3)

    def test_isolated_probe_preserves_exact_envelope_species(self):
        original = (ROOT/"reference/parameters/ice_v02_dust.para").read_text()
        result = check.isolated_parameter(original)
        def grain_rows(text):
            section = text.split("#Grain properties", 1)[1].split("#Molecular RT settings", 1)[0]
            return [line for line in section.splitlines() if line.strip()]
        self.assertEqual(grain_rows(result)[1:], grain_rows(original)[10:])
        self.assertNotIn("ac_opct.dat", result)
        self.assertIn("F T F  no temperature", result)
        self.assertIn("8 8 1 2", result)
        self.assertIn("3 3 6000.", result)
        self.assertEqual(result.count("DHS  1 1"), 2)

    def test_full_range_grid_contains_all_zero_k_and_probes(self):
        import shutil
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"inputs/utils/Dust").mkdir(parents=True)
            ice = ROOT/"reference/dust/H2O_30K_Leiden_mcfost.dat"
            shutil.copy2(ice, root/"inputs/utils/Dust"/ice.name)
            probes = [dict(wavelength_um=v) for v in (1.74, 2.546, 9.7, 27.51)]
            models = [dict(id="m0")]
            write_temperature_inputs(root, models)
            manifest = {"anchors": probes, "models": models}
            wave = check.wavelength_grid(root, manifest)
            self.assertEqual((wave[0], wave[-1]), (.1, 3000.))
            self.assertTrue(np.all(np.diff(wave) > 0))
            rows = [line.split("#", 1)[0].strip() for line in ice.read_text().splitlines()]
            samples = np.loadtxt([row for row in rows if row][1:])
            knots = samples[samples[:, 2] == 0, 0]
            self.assertEqual(len(knots), 206)
            self.assertTrue(np.isin(knots, wave).all())
            self.assertTrue(np.isin([r["wavelength_um"] for r in probes], wave).all())
            required = check.required_wavelength_grid(root, manifest)
            self.assertTrue(np.isin(required, wave).all())
            self.assertTrue(np.isin(np.geomspace(.1, 3000., 100), required).all())
            self.assertTrue(np.isin([r["wavelength_um"] for r in probes], required).all())
            self.assertLess(len(required), len(wave))

    def test_required_grid_reads_thermal_centres_from_every_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            models = [dict(id="m0"), dict(id="m1")]
            write_temperature_inputs(root, models)
            other = root/"models/m1/temperature.para"
            rows = other.read_text().splitlines()
            start = rows.index("#Wavelength")
            rows[start+1] = "20 0.2 2000.0  n_lambda, lambda_min, lambda_max"
            other.write_text("\n".join(rows)+"\n")
            manifest = dict(models=models, anchors=[dict(wavelength_um=1.74), dict(wavelength_um=27.51)])
            required = check.required_wavelength_grid(root, manifest)
            self.assertTrue(np.all(np.isfinite(required)))
            self.assertTrue(np.all(np.diff(required) > 0))
            self.assertEqual((required[0], required[-1]), (.1, 3000.))
            for count, lower, upper in ((50, .1, 3000.), (20, .2, 2000.)):
                edges = np.geomspace(lower, upper, count+1)
                centers = np.sqrt(edges[:-1]*edges[1:])
                for value in centers:
                    self.assertTrue(np.any(np.isclose(required, value, rtol=2e-12, atol=0)),
                                    f"Missing {count}-bin model thermal wavelength {value}")

    def test_dust_outputs_reject_invalid_numbers_and_wavelengths(self):
        wave = np.array([.1, 1., 9.7, 3000.])
        with tempfile.TemporaryDirectory() as tmp:
            write_outputs(tmp, wave)
            result = check.check_outputs(tmp, wave)
            self.assertEqual(result["wavelength_count"], 4)
            self.assertEqual(result["absorption_min_cm2_g"], .5)
            self.assertFalse(result["g_available"])
            self.assertFalse(result["g_applicable"])
            with self.assertRaisesRegex(ValueError, "one g FITS"):
                check.check_outputs(tmp, wave, g_applicable=True)
            for opacity, albedo, message in (
                (np.array([1., np.nan, 1., 1.]), None, "Nonfinite"),
                (np.array([1., -1., 1., 1.]), None, "Negative"),
                (None, np.array([.5, 1.01, .5, .5]), "Albedo outside"),
            ):
                write_outputs(tmp, wave, opacity=opacity, albedo=albedo)
                with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                    check.check_outputs(tmp, wave)
            write_outputs(tmp, wave)
            with self.assertRaisesRegex(ValueError, "wavelengths"):
                check.check_outputs(tmp, wave*2)

    def test_optional_asymmetry_file_is_validated_when_present(self):
        wave = np.array([.1, 1., 9.7, 3000.])
        with tempfile.TemporaryDirectory() as tmp:
            write_outputs(tmp, wave, include_g=True)
            self.assertTrue(check.check_outputs(tmp, wave)["g_available"])
            target = Path(tmp)/"data_dust/g.fits.gz"
            for value, message in ((float("nan"), "Nonfinite"), (1.01, "Asymmetry outside")):
                fits.PrimaryHDU(np.full(len(wave), value)).writeto(target, overwrite=True)
                with self.subTest(value=value), self.assertRaisesRegex(ValueError, message):
                    check.check_outputs(tmp, wave)

    def test_unused_stress_wavelength_failures_are_recorded_without_modifying_fits(self):
        wave = np.array([.1, 1.74, 9.00867, 9.1, 10., 3000.])
        required = wave[[0, 1, 4, 5]]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_outputs(root, wave, opacity=np.array([1., 1., np.nan, -1., 1., 1.]),
                          albedo=np.array([.5, .5, np.nan, 1.2, .5, .5]))
            grain = np.ones((50, len(wave)))
            grain[:, 2] = np.nan
            grain[:3, 3] = -1.
            fits.PrimaryHDU(grain).writeto(root/"data_dust/kappa_grain.fits.gz", overwrite=True)
            before = {p.name: check.digest(p) for p in (root/"data_dust").iterdir()}
            with self.assertRaises(ValueError):
                check.check_outputs(root, wave)  # Omitting the subset retains strict behavior.
            result = check.check_outputs(root, wave, required_wavelengths=required)
            self.assertEqual(result["required_wavelength_count"], 4)
            self.assertEqual(result["diagnostic_wavelength_count"], 2)
            diagnostics = result["extra_wavelength_diagnostics"]
            self.assertEqual(diagnostics["kappa"]["nonfinite_count"], 1)
            self.assertEqual(diagnostics["kappa"]["negative_count"], 1)
            self.assertEqual(diagnostics["kappa_grain"]["nonfinite_count"], 50)
            self.assertEqual(diagnostics["kappa_grain"]["negative_count"], 3)
            self.assertEqual(diagnostics["albedo"]["nonfinite_count"], 1)
            self.assertEqual(diagnostics["albedo"]["out_of_range_count"], 1)
            np.testing.assert_allclose(diagnostics["kappa"]["nonfinite_wavelengths_um"], [wave[2]])
            np.testing.assert_allclose(diagnostics["kappa_grain"]["nonfinite_wavelengths_um"], [wave[2]])
            self.assertTrue(np.isfinite(result["extinction_min_cm2_g"]))
            self.assertTrue(np.isfinite(result["absorption_min_cm2_g"]))
            self.assertEqual(before, {p.name: check.digest(p) for p in (root/"data_dust").iterdir()})
            self.assertEqual(np.count_nonzero(~np.isfinite(fits.getdata(root/"data_dust/kappa_grain.fits.gz"))), 50)

    def test_invalid_required_thermal_or_image_wavelength_always_rejects(self):
        wave = np.array([.110859066, 1.74, 9.00867])
        required = wave[:2]  # Representative default thermal centre and production image.
        failures = (("kappa", np.nan), ("kappa", -1.), ("kappa_grain", np.nan),
                    ("kappa_grain", -1.), ("albedo", 1.01), ("phase_function", np.nan),
                    ("phase_function", -1.), ("polarizability", np.nan), ("g", 1.01))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for index in (0, 1):
                for name, value in failures:
                    write_outputs(root, wave, include_g=True)
                    path = root/"data_dust"/(name+".fits.gz")
                    values = fits.getdata(path)
                    values[..., index] = value
                    fits.PrimaryHDU(values).writeto(path, overwrite=True)
                    with self.subTest(index=index, name=name, value=value), self.assertRaises(ValueError):
                        check.check_outputs(root, wave, required_wavelengths=required)

    def test_required_wavelength_match_takes_priority_over_stress_label(self):
        knot = 9.00867
        wave = np.array([.1, knot, knot*(1.+1e-6), 3000.])
        with tempfile.TemporaryDirectory() as tmp:
            write_outputs(tmp, wave, opacity=np.array([1., np.nan, 1., 1.]))
            # The exact zero-k knot lies within the spectral matching tolerance
            # of a required image, so it cannot be silently exempted as stress.
            with self.assertRaises(ValueError):
                check.check_outputs(tmp, wave, required_wavelengths=wave[[0, 2, 3]])

    def test_required_set_must_be_nonempty_finite_and_present_in_output(self):
        wave = np.array([.1, 1.74, 9.00867, 3000.])
        with tempfile.TemporaryDirectory() as tmp:
            write_outputs(tmp, wave)
            for required in ([], [np.nan], [1.8], [[1.74]]):
                with self.subTest(required=required), self.assertRaises(ValueError):
                    check.check_outputs(tmp, wave, required_wavelengths=np.asarray(required))
            result = check.check_outputs(tmp, wave, required_wavelengths=wave*(1.+1e-7))
            self.assertEqual(result["required_wavelength_count"], len(wave))

    def test_spectral_axis_and_nonfinite_wavelengths_never_become_diagnostics(self):
        wave = np.array([.1, 1.74, 9.00867, 3000.])
        malformed = (("lambda", np.ones((2, len(wave)))),
                     ("kappa", np.ones((2, len(wave)))),
                     ("albedo", np.ones((2, len(wave)))),
                     ("kappa_grain", np.ones(len(wave))),
                     ("kappa_grain", np.ones((len(wave), 50))),
                     ("phase_function", np.ones((len(wave), 181))),
                     ("polarizability", np.zeros((181, len(wave)+1))))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, data in malformed:
                write_outputs(root, wave)
                fits.PrimaryHDU(data).writeto(root/"data_dust"/(name+".fits.gz"), overwrite=True)
                with self.subTest(name=name, shape=data.shape), self.assertRaises(ValueError):
                    check.check_outputs(root, wave, required_wavelengths=wave[:2])
            write_outputs(root, wave)
            broken_wave = wave.copy()
            broken_wave[2] = np.nan  # Nonfinite extra wavelength is still an invalid axis.
            fits.PrimaryHDU(broken_wave).writeto(root/"data_dust/lambda.fits.gz", overwrite=True)
            with self.assertRaises(ValueError):
                check.check_outputs(root, wave, required_wavelengths=wave[:2])

    def test_real_subprocess_default_root_and_exit_zero_failures(self):
        """Exercise the filesystem behavior omitted by the original process mock."""
        import shutil
        from silicate_structure_design import catalogue
        with tempfile.TemporaryDirectory(prefix="material-root-regression-") as tmp:
            root = Path(tmp)
            executable = root/"fake-mcfost"
            executable.write_text(f"#!{sys.executable}\n" + '''\
import os
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
print("Input file read successfully", flush=True)
mode = os.environ.get("FAKE_PREFLIGHT_FAILURE", "")
if mode == "incomplete":
    raise SystemExit(0)
root = Path(sys.argv[sys.argv.index("-root_dir")+1]) if "-root_dir" in sys.argv else Path(".")
(root/"data_dust").mkdir(parents=True)
print("Computing dust properties ... Writing dust properties", flush=True)
if mode != "missing":
    # Independently mimic lect_lambda: each numeric line is one wavelength,
    # including any integer that a caller incorrectly treats as a row count.
    samples = []
    for line in Path("check.lambda").read_text().splitlines():
        numeric = line.split("#", 1)[0].strip()
        if numeric:
            samples.append(float(numeric.split()[0]))
    # tab_lambda2 and the lambda FITS product are default REAL in 4.1.14.
    wave = np.asarray(samples, dtype=np.float32)
    ones = np.ones(len(wave))
    data = {"lambda": wave, "kappa": ones, "albedo": .5*ones,
            "kappa_grain": np.ones((2,len(wave))), "phase_function": np.ones((181,len(wave))),
            "polarizability": np.zeros((181,len(wave)))}
    try:
        # Reproduce 4.1.14: creation respects root_dir; FITS writer does not.
        for name, value in data.items():
            fits.PrimaryHDU(value).writeto(Path("data_dust")/(name+".fits.gz"))
    except OSError:
        print("FITSIO ErrorStatus = 105: failed to create data_dust/lambda.fits.gz", flush=True)
    if mode == "fitsio":
        print("FiTsIo ErRoRsTaTuS = 105", flush=True)
print(" Exiting", flush=True)
''')
            executable.chmod(0o755)
            # The previous command reproduces the user's exit-zero FITSIO error.
            legacy = root/"old_layout"
            legacy.mkdir()
            (legacy/"check.lambda").write_text(".1\n1\n9.7\n3000\n")
            outcome = subprocess.run([str(executable), "check.para", "-dust_prop", "-root_dir", "output"],
                                     cwd=legacy, capture_output=True, text=True, timeout=30)
            self.assertEqual(outcome.returncode, 0)
            self.assertIn("FITSIO ErrorStatus = 105", outcome.stdout)
            self.assertTrue((legacy/"output/data_dust").is_dir())
            self.assertFalse(list(legacy.rglob("*.fits.gz")))
            # A count-prefixed file is a different bug: the documented parser
            # emits five wavelengths, with an unintended 4-micron first row.
            bad_header = root/"old_count_header"
            bad_header.mkdir()
            (bad_header/"check.lambda").write_text("# Custom wavelengths in microns\n4\n.1\n1\n9.7\n3000\n")
            outcome = subprocess.run([str(executable), "check.para", "-dust_prop"],
                                     cwd=bad_header, capture_output=True, text=True, timeout=30)
            self.assertEqual(outcome.returncode, 0)
            actual = fits.getdata(bad_header/"data_dust/lambda.fits.gz")
            intended = np.array([.1, 1., 9.7, 3000.])
            np.testing.assert_array_equal(actual, np.r_[len(intended), intended].astype(np.float32))
            self.assertEqual(actual.shape, (len(intended)+1,))
            with self.assertRaisesRegex(ValueError, r"lambda wavelength shape.*expected \(4,\).*got \(5,\)"):
                check.check_outputs(bad_header, intended)
            (root/"inputs/utils/Dust").mkdir(parents=True)
            shutil.copy2(ROOT/"reference/dust/H2O_30K_Leiden_mcfost.dat", root/"inputs/utils/Dust")
            unique = {}
            for row in catalogue():
                p = row["parameters"]
                unique.setdefault((p["envelope_silicate_file"], p["envelope_amax_um"]), row)
            models = [dict(id=f"m{i}", parameters=row["parameters"]) for i, row in enumerate(unique.values())]
            (root/"manifest.json").write_text(json.dumps(dict(models=models, anchors=[dict(wavelength_um=9.7)])))
            write_temperature_inputs(root, models)
            fake_mode = {"value": ""}
            runner = SimpleNamespace(
                runtime_config=lambda _: {"mcfost_executable": str(executable)},
                environment=lambda _, runtime: {**os.environ, "OMP_NUM_THREADS": str(runtime["threads"]),
                                                "FAKE_PREFLIGHT_FAILURE": fake_mode["value"]},
                now=lambda: "2026-09-24T00:00:00+00:00",
                atomic_json=lambda path, payload: Path(path).write_text(json.dumps(payload)))
            dispatch = SimpleNamespace(validate_package=lambda _: {"experiment_id": SILICATE_STRUCTURE},
                                       load_runner=lambda _: runner)
            with patch.dict(sys.modules, {"production_task": dispatch}), \
                    patch.object(check, "identity", return_value={"runtime": "same"}):
                for mode, message in (("fitsio", "FITSIO error reported despite exit 0"),
                                      ("incomplete", "completion markers absent"),
                                      ("missing", "Expected one lambda FITS product")):
                    fake_mode["value"] = mode
                    with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, message) as caught:
                        check.run_checks(root, "unused")
                    self.assertIn("mcfost.log", str(caught.exception))
                    self.assertIn("Last simulator output", str(caught.exception))
                    receipt = json.loads((root/check.RECEIPT).read_text())
                    self.assertEqual(receipt["status"], "failed")
                    self.assertEqual(receipt["checks"][0]["returncode"], 0)
                    with self.assertRaisesRegex(ValueError, "not passed"):
                        check.validate_receipt(root, {})
                fake_mode["value"] = ""
                self.assertEqual(check.run_checks(root, "unused")["status"], "passed")
                receipt = check.validate_receipt(root, {})
                manifest = json.loads((root/"manifest.json").read_text())
                expected = check.wavelength_grid(root, manifest)
                for item in receipt["checks"]:
                    self.assertNotIn("-root_dir", item["command"])
                    self.assertFalse(item["diagnostics"]["g_available"])
                    self.assertFalse(item["diagnostics"]["g_applicable"])
                    self.assertEqual(len(item["diagnostics"]["product_sha256"]), 6)
                    self.assertEqual(item["diagnostics"]["wavelength_count"], len(expected))
                    model_attempt = sorted((root/"material_preflight"/item["model_id"]).glob("attempt_*"))[-1]
                    written = np.loadtxt(model_attempt/"check.lambda")
                    emitted = fits.getdata(model_attempt/"data_dust/lambda.fits.gz")
                    np.testing.assert_array_equal(written, expected)
                    np.testing.assert_array_equal(emitted, expected.astype(np.float32))
                # Failed outputs are retained; a new attempt fixes the root.
                self.assertTrue((root/"material_preflight/m0/attempt_001/mcfost.log").exists())
                self.assertTrue((root/"material_preflight/m0/attempt_004/data_dust/lambda.fits.gz").exists())

    def test_receipt_binds_runtime_prescriptions_and_output_bytes(self):
        models = [dict(id=f"m{i}", parameters=dict(envelope_silicate_file=material, envelope_amax_um=size))
                  for i, (material, size) in enumerate((("Draine", .4), ("Draine", 1.), ("Olivine", .4), ("Pyroxene", .4)))]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"manifest.json").write_text(json.dumps({"models": models}))
            (root/"material_preflight").mkdir()
            artifact = root/"material_preflight/kappa.fits"
            artifact.write_bytes(b"a checked product")
            identity = {"runtime": "original"}
            receipt = dict(status="passed", identity=identity, checks=[
                dict(model_id=m["id"], status="passed", artifact_sha256={str(artifact.relative_to(root)): check.digest(artifact)})
                for m in models])
            (root/check.RECEIPT).write_text(json.dumps(receipt))
            with patch.object(check, "identity", return_value=identity):
                self.assertEqual(check.validate_receipt(root, {}), receipt)
                artifact.write_bytes(b"changed")
                with self.assertRaisesRegex(ValueError, "artifacts changed"):
                    check.validate_receipt(root, {})
            with patch.object(check, "identity", return_value={"runtime": "different"}):
                with self.assertRaisesRegex(ValueError, "runtime or frozen inputs changed"):
                    check.validate_receipt(root, {})


if __name__ == "__main__":
    unittest.main()
