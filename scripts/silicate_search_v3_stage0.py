#!/usr/bin/env python3
"""Measure four pure DHS species with MCFOST before matching envelope masses.

The mandatory wavelength set is exactly the production image set plus 2.2,
9.7 and 18 microns. There are no extra zero-k stress samples. This is a dust
property calculation, not a temperature solution or a fit to observations.
The three matching opacities must reproduce v2 within one part per million;
an exactly zero v2 absorption must remain exactly zero. No values are floored.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

RECEIPT = "stage0/receipt.json"
POLICY = "four_pure_DHS_species_v2_reproducibility_v1"
PRODUCTS = ("lambda", "kappa", "albedo", "kappa_grain")
TIMEOUT_SECONDS = 240
REPRODUCIBILITY_RTOL = 1e-6
REFERENCE_FIELDS = {"ext_2p2": ("kappa_ext_cm2_g", 2.2),
                    "abs_9p7": ("kappa_abs_cm2_g", 9.7),
                    "abs_18": ("kappa_abs_cm2_g", 18.)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def species_definitions():
    from silicate_search_v3_design import SPECIES
    require(len(SPECIES) == 4, "Stage 0 requires four pure species")
    return SPECIES


def isolated_parameter(text, filename):
    """Keep the v02 envelope geometry and substitute one pure DHS species.

    The source template must retain its two disk and two envelope species.
    Each output is one separate population with a mass fraction of one, the
    same fixed sizes/shape as all three production envelope populations.
    """
    require(re.fullmatch(r"[A-Za-z0-9_.-]+", filename) is not None,
            "Unsafe optical-constants filename")
    lines = text.splitlines()

    def section(start, end):
        first = next(i for i, row in enumerate(lines) if row.startswith(start))
        last = next(i for i in range(first + 1, len(lines)) if lines[i].startswith(end))
        rows = [i for i in range(first + 1, last)
                if lines[i].strip() and not lines[i].lstrip().startswith("#")]
        return first, last, rows

    first, last, rows = section("#Grain properties", "#Molecular RT settings")
    require(len(rows) == 18 and lines[rows[0]].split()[0] == "2"
            and lines[rows[9]].split()[0] == "2",
            "Expected original v02 two-disk/two-envelope-species template")
    for offset in (10, 14):
        values = lines[rows[offset]].split()
        require(values[:4] == ["DHS", "1", "1", "0.0"]
                and float(values[5]) == .1
                and [float(v) for v in lines[rows[offset + 3]].split()[:4]] == [.03, .4, 2.75, 50.],
                "Stage 0 requires fixed non-porous DHS vmax=0.1, 0.03-0.4 um, q=2.75, 50 bins")
    lines[first + 1:last] = [
        "1  Number of species  ZONE 1: pure material opacity measurement",
        "DHS 1 1 0.0 1.0 0.1  type, components, mixing, porosity, mass fraction, vmax",
        f"{filename} 1.0  optical indices file, within-species volume fraction",
        "1  Heating method: RE + LTE",
        "0.03 0.4 2.75 50  amin, amax [um], exponent, size bins", ""]
    first, last, rows = section("#Density structure", "#Grain properties")
    require(len(rows) == 12 and lines[rows[6]].split()[0] == "3",
            "Expected original disk and envelope density blocks")
    lines[first + 1:last] = [lines[i].replace("ZONE 2", "ZONE 1") for i in rows[6:]] + [""]
    for start, end, replacements in (
        ("#Number of zones", "#Density structure", {0: "1  number of zones"}),
        ("#Number of photon packages", "#Wavelength", {0: "1", 1: "1", 2: "1"}),
        ("#Wavelength", "#Grid geometry", {1: "F T F  no temperature; custom dust wavelengths", 2: "check.lambda"}),
        ("#Grid geometry", "#Maps", {1: "8 8 1 2"}),
        ("#Maps", "#Scattering method", {0: "3 3 6000."}),
    ):
        _, _, rows = section(start, end)
        for offset, replacement in replacements.items():
            lines[rows[offset]] = replacement
    return "\n".join(lines) + "\n"


def validate_inputs(bundle):
    import numpy as np
    bundle = Path(bundle)
    plan = json.loads((bundle / "inputs/stage0_wavelengths.json").read_text())
    wave = np.asarray(plan["wavelengths_um"], dtype=float)
    production = np.asarray(plan["production_wavelengths_um"], dtype=float)
    require(production.shape == (98,) and bool(np.isfinite(production).all())
            and bool((production > 0).all()) and np.unique(production).size == 98,
            "Stage 0 requires 98 distinct production wavelengths")
    expected = np.unique(np.concatenate((production, [2.2, 9.7, 18.])))
    require(wave.shape == (99,) and np.array_equal(wave, expected),
            "Stage 0 wavelengths must be exactly 99 production/screen wavelengths")
    template = (bundle / "inputs/template.para").read_text()
    for species in species_definitions().values():
        path = bundle / "inputs/utils/Dust" / species["filename"]
        require(path.is_file() and not path.is_symlink() and digest(path) == species["sha256"],
                f"Supplied optical constants changed or absent: {species['filename']}")
        isolated_parameter(template, species["filename"])
    reference = json.loads((bundle / "inputs/v2_reference.json").read_text())["species_opacities_cm2_g"]
    require(set(reference) == set(species_definitions()), "V2 reference must cover the four V3 species")
    for species, values in reference.items():
        require(set(values) == set(REFERENCE_FIELDS)
                and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                        and np.isfinite(v) and v >= 0 for v in values.values())
                and values["ext_2p2"] > 0,
                f"Invalid direct V2 reference opacities: {species}")
    return wave


def check_outputs(directory, expected):
    """Read mass extinction/scattering/absorption without clipping or filling.

    MCFOST dust_prop.f90 write_dust_prop reports kappa per gram of dust and
    kappa_grain as the individual grain absorption cross section / grain mass.
    Phase/polarizability/g files do not enter mass matching and are not gated.
    """
    import numpy as np
    from astropy.io import fits
    directory = Path(directory)
    expected = np.asarray(expected, dtype=float)
    require(expected.ndim == 1 and expected.size > 0 and bool(np.isfinite(expected).all())
            and bool((expected > 0).all()) and bool((np.diff(expected) > 0).all()),
            "Invalid expected wavelengths")
    values, hashes = {}, {}
    for name in PRODUCTS:
        paths = list(directory.rglob(name + ".fits.gz")) + list(directory.rglob(name + ".fits"))
        require(len(paths) == 1, f"Expected one {name} FITS product; found {len(paths)}")
        array = np.asarray(fits.getdata(paths[0]), dtype=float).squeeze()
        shape = (50, expected.size) if name == "kappa_grain" else expected.shape
        require(array.shape == shape, f"Unexpected {name} shape: expected {shape}, got {array.shape}")
        bad = ~np.isfinite(array)
        bad_wave = expected[bad.any(axis=0) if bad.ndim == 2 else bad]
        require(not bool(bad.any()), f"Nonfinite {name} at required wavelengths (um): {bad_wave[:20].tolist()}")
        if name in ("kappa", "kappa_grain"):
            require(bool((array >= 0).all()), f"Negative {name} opacity")
        if name == "albedo":
            require(bool(((array >= 0) & (array <= 1)).all()), "Albedo outside [0, 1]")
        if name == "lambda":
            require(bool(np.allclose(array, expected, rtol=2e-6, atol=0)),
                    "Dust-property wavelengths differ from the declared wavelengths")
        values[name] = array
        hashes[str(paths[0].relative_to(directory))] = digest(paths[0])
    extinction, albedo = values["kappa"], values["albedo"]
    require(bool((extinction > 0).any()), "Dust extinction is identically zero")
    return dict(schema_version=1, wavelength_um=expected.tolist(),
                fits_wavelength_um=values["lambda"].tolist(),
                kappa_ext_cm2_g=extinction.tolist(), albedo=albedo.tolist(),
                kappa_abs_cm2_g=(extinction * (1 - albedo)).tolist(),
                kappa_sca_cm2_g=(extinction * albedo).tolist(),
                units="cm^2 per g of dust", normalization="pure species dust mass",
                grain_absorption_shape=list(values["kappa_grain"].shape), product_sha256=hashes)


def reference_comparison(bundle, species, opacity):
    """Compare direct recorded v2 values; relative tolerance implies exact zero.

    The reference's recovered least-squares coefficients are provenance only.
    Numerical negatives from that inversion are never treated as opacities.
    This function neither rounds nor replaces an actual or reference opacity.
    """
    import numpy as np
    reference = json.loads((Path(bundle) / "inputs/v2_reference.json").read_text())["species_opacities_cm2_g"][species]
    wave = np.asarray(opacity["wavelength_um"], float)
    comparisons = {}
    for key, (field, wavelength) in REFERENCE_FIELDS.items():
        matches = np.flatnonzero(np.isclose(wave, wavelength, rtol=2e-6, atol=0))
        require(len(matches) == 1, f"Missing or ambiguous reproduction wavelength {wavelength}")
        actual, expected = float(opacity[field][matches[0]]), float(reference[key])
        require(np.isfinite(actual) and actual >= 0 and np.isfinite(expected) and expected >= 0,
                f"Invalid opacity in V2 reproducibility check: {species} {key}")
        difference = abs(actual - expected)
        tolerance = REPRODUCIBILITY_RTOL * abs(expected)
        comparisons[key] = dict(wavelength_um=wavelength, actual_cm2_g=actual,
                                v2_reference_cm2_g=expected, absolute_difference_cm2_g=difference,
                                relative_difference=difference / expected if expected else None,
                                absolute_tolerance_cm2_g=tolerance,
                                policy="relative_1e-6" if expected else "exact_zero",
                                passed=bool(difference <= tolerance))
    return dict(passed=all(row["passed"] for row in comparisons.values()),
                relative_tolerance=REPRODUCIBILITY_RTOL, absolute_tolerance_floor_cm2_g=0.,
                optical_constants_or_opacities_modified=False, comparisons=comparisons)


def log_context(path):
    return f"Log: {path}\nLast simulator output:\n" + "\n".join(
        Path(path).read_text(errors="replace").splitlines()[-40:])


def validate_completion(path, returncode):
    text = Path(path).read_text(errors="replace").casefold()
    error = None
    if returncode != 0:
        error = f"exit {returncode}"
    elif re.search(r"fitsio\s*error\s*status", text):
        error = "FITSIO error despite exit 0"
    elif re.search(r"^\s*(?:error\b|fatal\b|fortran runtime error\b)", text, re.MULTILINE) \
            or "program received signal" in text or "segmentation fault" in text:
        error = "simulator error despite exit 0"
    elif "writing dust properties" not in text or not re.search(r"^\s*exiting\s*$", text, re.MULTILINE):
        error = "dust-property completion markers absent"
    if error:
        raise ValueError(f"Stage 0 initialization failed ({error}).\n{log_context(path)}")


def load_runner(bundle):
    source = Path(bundle) / "code/src"
    sys.path.insert(0, str(source))
    from mcfost_grid import runner
    return runner


def identity(bundle, runtime=None):
    bundle = Path(bundle)
    species = species_definitions()
    paths = ["inputs/template.para", "inputs/stage0_wavelengths.json", "inputs/v2_reference.json"]
    paths.extend("inputs/utils/Dust/" + row["filename"] for row in species.values())
    result = dict(policy=POLICY, code_sha256=digest(__file__),
                  input_sha256={p: digest(bundle / p) for p in paths}, species=species)
    if runtime is not None:
        runner = load_runner(bundle)
        result["runtime"] = dict(executable_sha256=digest(runtime["mcfost_executable"]),
                                 utilities_content_sha256=runner._utilities_hash(runtime),
                                 backend=runtime["backend"])
    return result


def safe_artifact(bundle, relative):
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts and bool(path.parts),
            "Unsafe Stage 0 artifact path")
    candidate = Path(bundle) / path
    require(candidate.resolve().is_relative_to(Path(bundle).resolve()) and not candidate.is_symlink()
            and candidate.is_file(), f"Stage 0 artifact missing or unsafe: {relative}")
    return candidate


def validate_stage0_receipt(bundle, runtime=None):
    """Verify products and frozen inputs; optionally bind the exact runtime."""
    bundle = Path(bundle).resolve()
    wave = validate_inputs(bundle)
    receipt = json.loads((bundle / RECEIPT).read_text())
    require(receipt.get("status") == "passed" and receipt.get("policy") == POLICY,
            "Stage 0 has not passed; all four pure species are required")
    recorded = dict(receipt.get("identity", {}))
    require(isinstance(recorded.get("runtime"), dict) and recorded["runtime"].get("executable_sha256")
            and recorded["runtime"].get("utilities_content_sha256"), "Stage 0 runtime identity absent")
    if runtime is None:
        recorded.pop("runtime")
    require(recorded == identity(bundle, runtime), "Stage 0 inputs, code or runtime changed")
    require(receipt.get("wavelengths_um") == wave.tolist(), "Stage 0 wavelength receipt changed")
    require(receipt.get("v2_reproducibility_checked") is True, "Stage 0 has not reproduced V2 opacities")
    definitions = species_definitions()
    checks = receipt.get("species", [])
    require([row.get("id") for row in checks] == list(definitions), "Stage 0 species coverage changed")
    template = (bundle / "inputs/template.para").read_text()
    for check in checks:
        species = definitions[check["id"]]
        require(check.get("status") == "passed" and check.get("filename") == species["filename"],
                "Incomplete or changed pure-species check")
        artifacts = check.get("artifact_sha256", {})
        require(bool(artifacts), "Stage 0 artifact hashes absent")
        for relative, checksum in artifacts.items():
            require(digest(safe_artifact(bundle, relative)) == checksum, "Stage 0 artifact changed")
        attempt = bundle / check["attempt"]
        require(attempt.resolve().is_relative_to((bundle / "stage0" / check["id"]).resolve()),
                "Stage 0 attempt outside its species directory")
        for filename in ("check.para", "check.lambda", "mcfost.log", "opacity.json"):
            require(str((attempt / filename).relative_to(bundle)) in artifacts,
                    f"Stage 0 required artifact not bound: {filename}")
        require((attempt / "check.para").read_text() == isolated_parameter(template, species["filename"]),
                "Stage 0 pure-species parameter file changed")
        expected_table = "\n".join(f"{v:.17g}" for v in wave) + "\n"
        require((attempt / "check.lambda").read_text() == expected_table, "Stage 0 lambda file changed")
        validate_completion(attempt / "mcfost.log", check.get("returncode"))
        opacity_path = safe_artifact(bundle, check["opacity_file"])
        require(opacity_path == attempt / "opacity.json" and digest(opacity_path) == check["opacity_sha256"],
                "Stage 0 opacity summary changed")
        derived = check_outputs(attempt, wave)
        require(json.loads(opacity_path.read_text()) == derived, "Opacity summary differs from FITS products")
        reproduced = reference_comparison(bundle, check["id"], derived)
        require(reproduced["passed"] and check.get("v2_reproducibility") == reproduced,
                "Stage 0 V2 opacity reproducibility is absent, failed or changed")
        for relative, checksum in derived["product_sha256"].items():
            require(artifacts.get(str((attempt / relative).relative_to(bundle))) == checksum,
                    "Stage 0 FITS product is not bound to the receipt")
    return receipt


def run_stage0(bundle, machine):
    bundle = Path(bundle).resolve()
    wave = validate_inputs(bundle)
    runner = load_runner(bundle)
    runtime = runner.runtime_config(machine)
    root = bundle / "stage0"
    root.mkdir(exist_ok=True)
    with (root / ".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (bundle / RECEIPT).is_file():
            previous = json.loads((bundle / RECEIPT).read_text())
            if previous.get("status") == "passed":
                validate_stage0_receipt(bundle, runtime)
                return dict(status="cached", receipt=str(bundle / RECEIPT))
        state = dict(schema_version=1, policy=POLICY, status="running", identity=identity(bundle, runtime),
                     wavelengths_um=wave.tolist(), species=[], threads=2,
                     timeout_seconds_per_species=TIMEOUT_SECONDS, started_utc=runner.now(),
                     fresh_temperatures_computed=False, images_computed=False,
                     optical_constants_modified=False, native_out_of_range_handling=True,
                     v2_reproducibility_checked=False, reproducibility_relative_tolerance=REPRODUCIBILITY_RTOL,
                     zero_reference_policy="Exactly zero absorption must reproduce as exactly zero; no floor or absolute tolerance",
                     mass_match_scope="Pure-species DHS mass opacities; no observations or fitted quantities",
                     thermal_grid_validated=False)
        runner.atomic_json(bundle / RECEIPT, state)
        try:
            for name, species in species_definitions().items():
                attempts = root / name
                attempts.mkdir(exist_ok=True)
                number = 1
                while (attempts / f"attempt_{number:03d}").exists():
                    number += 1
                attempt = attempts / f"attempt_{number:03d}"
                attempt.mkdir()
                template = (bundle / "inputs/template.para").read_text()
                (attempt / "check.para").write_text(isolated_parameter(template, species["filename"]))
                # lect_lambda reads every numeric line as a wavelength. A
                # numeric count header would become an extra sample.
                (attempt / "check.lambda").write_text("\n".join(f"{v:.17g}" for v in wave) + "\n")
                # MCFOST 4.1.14 writes data_dust relative to cwd even when
                # -root_dir points elsewhere. Every attempt uses default '.'.
                command = [runtime["mcfost_executable"], "check.para", "-dust_prop", "-max_mem", "2", "-no_backup"]
                check = dict(id=name, filename=species["filename"], status="running",
                             attempt=str(attempt.relative_to(bundle)), command=command)
                state["species"].append(check)
                runner.atomic_json(bundle / RECEIPT, state)
                started = time.monotonic()
                with (attempt / "mcfost.log").open("w") as log:
                    outcome = subprocess.run(command, cwd=attempt,
                        env=runner.environment(bundle, {**runtime, "threads": 2}),
                        stdout=log, stderr=subprocess.STDOUT, timeout=TIMEOUT_SECONDS, check=False)
                check.update(returncode=outcome.returncode, elapsed_seconds=time.monotonic() - started)
                validate_completion(attempt / "mcfost.log", outcome.returncode)
                try:
                    opacity = check_outputs(attempt, wave)
                except (ValueError, OSError) as exc:
                    raise ValueError(f"Stage 0 opacity products failed checks: {exc}.\n"
                                     f"{log_context(attempt / 'mcfost.log')}") from exc
                output = attempt / "opacity.json"
                runner.atomic_json(output, opacity)
                check["v2_reproducibility"] = reference_comparison(bundle, name, opacity)
                require(check["v2_reproducibility"]["passed"],
                        f"Stage 0 V2 reproducibility failed for {name}: "
                        + json.dumps({key: row for key, row in check["v2_reproducibility"]["comparisons"].items()
                                      if not row["passed"]}, allow_nan=False))
                check.update(status="passed", opacity_file=str(output.relative_to(bundle)),
                             opacity_sha256=digest(output), artifact_sha256={
                                 str(path.relative_to(bundle)): digest(path)
                                 for path in sorted(attempt.rglob("*")) if path.is_file()})
                runner.atomic_json(bundle / RECEIPT, state)
            require(state["identity"] == identity(bundle, runtime), "Stage 0 runtime inputs changed during calculation")
            state.update(status="passed", v2_reproducibility_checked=True, finished_utc=runner.now())
        except Exception as exc:
            if state["species"] and state["species"][-1]["status"] == "running":
                state["species"][-1].update(status="failed", error=f"{type(exc).__name__}: {exc}")
            state.update(status="failed", error=f"{type(exc).__name__}: {exc}", finished_utc=runner.now())
            runner.atomic_json(bundle / RECEIPT, state)
            raise
        runner.atomic_json(bundle / RECEIPT, state)
    return dict(status="passed", species=4, wavelengths=len(wave), receipt=str(bundle / RECEIPT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--machine", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_stage0(args.bundle, args.machine), indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
