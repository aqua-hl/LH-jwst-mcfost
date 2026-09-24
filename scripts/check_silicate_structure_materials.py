#!/usr/bin/env python3
"""Bounded isolated-DHS initialization gate; never compute temperatures/images.

MCFOST -dust_prop output contract:
https://mcfost.readthedocs.io/en/latest/outputs.html#dust-property-files
The test preserves both envelope species but removes the disk and uses an 8x8
spatial grid. It tests numerical initialization, not extrapolated-tail physics.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

NAME = "silicate_structure_production_v1"
RECEIPT = "material_preflight/receipt.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prescriptions(manifest):
    selected = {}
    for model in manifest["models"]:
        parameters = model["parameters"]
        key = (parameters["envelope_silicate_file"], parameters["envelope_amax_um"])
        selected.setdefault(key, model)
    require(len(selected) == 4, "Expected four distinct envelope material/size prescriptions")
    return list(selected.values())


def wavelength_grid(bundle, manifest):
    import numpy as np
    table = (Path(bundle) / "inputs/utils/Dust/H2O_30K_Leiden_mcfost.dat").read_text().splitlines()
    data = [row.split("#", 1)[0].strip() for row in table if row.split("#", 1)[0].strip()]
    ice = np.loadtxt(data[1:])  # First non-comment row is density/Tsub.
    knots = ice[ice[:, 2] == 0, 0]
    # Both boundaries, the actual thermal-bin centres, every production probe,
    # and the exact zero-k samples with neighbours are deliberately included.
    edges = np.geomspace(.1, 3000., 51)
    return np.unique(np.concatenate((np.geomspace(.1, 3000., 100), np.sqrt(edges[:-1]*edges[1:]),
        [a["wavelength_um"] for a in manifest["anchors"]], knots*(1.-1e-5), knots, knots*(1.+1e-5))))


def isolated_parameter(text):
    """Keep the exact envelope grain rows while removing disk-only work."""
    lines = text.splitlines()

    def section(start, end):
        first = next(i for i, row in enumerate(lines) if row.startswith(start))
        last = next(i for i in range(first+1, len(lines)) if lines[i].startswith(end))
        return first, last

    def active(first, last):
        return [i for i in range(first+1, last) if lines[i].strip() and not lines[i].lstrip().startswith("#")]

    first, last = section("#Grain properties", "#Molecular RT settings")
    rows = active(first, last)
    require(len(rows) == 18 and lines[rows[9]].split()[0] == "2"
            and lines[rows[10]].split()[0] == lines[rows[14]].split()[0] == "DHS",
            "Expected two disk species and two disjoint envelope DHS species")
    # Each single-component species occupies four rows plus one zone count.
    lines[first+1:last] = [lines[i].replace("ZONE 2", "ZONE 1") for i in rows[9:]] + [""]
    first, last = section("#Density structure", "#Grain properties")
    rows = active(first, last)
    require(len(rows) == 12 and lines[rows[6]].split()[0] == "3", "Expected disk and envelope density blocks")
    lines[first+1:last] = [lines[i].replace("ZONE 2", "ZONE 1") for i in rows[6:]] + [""]
    for start, end, updates in (
        ("#Number of zones", "#Density structure", {0: "1  number of zones"}),
        ("#Number of photon packages", "#Wavelength", {0: "1", 1: "1", 2: "1"}),
        ("#Wavelength", "#Grid geometry", {1: "F T F  no temperature; custom dust wavelengths", 2: "check.lambda"}),
        ("#Grid geometry", "#Maps", {1: "8 8 1 2"}),
        ("#Maps", "#Scattering method", {0: "3 3 6000."}),
    ):
        first, last = section(start, end)
        rows = active(first, last)
        for offset, replacement in updates.items():
            lines[rows[offset]] = replacement
    return "\n".join(lines) + "\n"


def check_outputs(directory, expected, *, g_applicable=False):
    import numpy as np
    from astropy.io import fits
    directory = Path(directory)
    values, products = {}, {}
    for name in ("lambda", "kappa", "albedo", "g", "kappa_grain", "phase_function", "polarizability"):
        paths = list(directory.rglob(name + ".fits.gz")) + list(directory.rglob(name + ".fits"))
        # Pinned MCFOST 4.1.14 writes g.fits only for HG (aniso_method=2).
        # This frozen campaign uses the exact phase function (method 1).
        # If an extra g file is present, still audit its values and identity.
        if name == "g" and not paths and not g_applicable:
            continue
        require(len(paths) == 1, f"Expected one {name} FITS product; found {len(paths)}")
        data = np.asarray(fits.getdata(paths[0]), dtype=float).squeeze()
        require(data.size > 0 and bool(np.isfinite(data).all()), f"Nonfinite or empty {name} output")
        if name in {"kappa", "kappa_grain", "phase_function"}:
            require(bool((data >= 0).all()), f"Negative {name} output")
        if name == "albedo":
            require(bool(((data >= 0) & (data <= 1)).all()), "Albedo outside [0, 1]")
        if name == "g":
            require(bool((np.abs(data) <= 1).all()), "Asymmetry outside [-1, 1]")
        values[name] = data
        products[str(paths[0].relative_to(directory))] = digest(paths[0])
    wave = values["lambda"]
    require(wave.shape == expected.shape and bool(np.allclose(wave, expected, rtol=2e-6, atol=0)),
            "Dust-property wavelengths do not match the requested full-range grid")
    opacity, albedo = values["kappa"], values["albedo"]
    require(opacity.shape == albedo.shape == wave.shape, "Unexpected integrated dust-property shape")
    require(bool((opacity > 0).any()), "Dust extinction is identically zero")
    return dict(wavelength_count=len(wave), wavelength_min_um=float(wave[0]),
                wavelength_max_um=float(wave[-1]), extinction_min_cm2_g=float(opacity.min()),
                absorption_min_cm2_g=float((opacity*(1-albedo)).min()),
                scattering_min_cm2_g=float((opacity*albedo).min()),
                g_available="g" in values, g_applicable=g_applicable, product_sha256=products)


def log_context(path):
    """Attach the simulator's explanation to failures of the outer wrapper."""
    path = Path(path)
    tail = "\n".join(path.read_text(errors="replace").splitlines()[-40:])
    return f"Log: {path}\nLast simulator output:\n{tail}"


def validate_completion(path, returncode):
    """Fortran STOP and FITSIO failures may return zero; require real completion."""
    text = Path(path).read_text(errors="replace")
    folded = text.casefold()
    reason = None
    if returncode != 0:
        reason = f"exit {returncode}"
    elif re.search(r"fitsio\s*error\s*status", folded):
        reason = "FITSIO error reported despite exit 0"
    elif re.search(r"^\s*(?:error\b|fatal\b|fortran runtime error\b)", folded, re.MULTILINE) \
            or "program received signal" in folded or "segmentation fault" in folded:
        reason = "simulator error reported despite exit 0"
    elif "writing dust properties" not in folded or not re.search(r"^\s*exiting\s*$", folded, re.MULTILINE):
        reason = "dust-property completion markers absent despite exit 0"
    if reason is not None:
        raise ValueError(f"Material initialization failed ({reason}).\n{log_context(path)}")


def identity(bundle, runtime):
    from mcfost_grid.runner import _utilities_hash
    return dict(manifest_sha256=digest(Path(bundle)/"manifest.json"),
                preflight_code_sha256=digest(__file__), executable_sha256=digest(runtime["mcfost_executable"]),
                utilities_content_sha256=_utilities_hash(runtime))


def validate_receipt(bundle, runtime):
    bundle = Path(bundle)
    path = bundle / RECEIPT
    require(path.is_file(), "Material preflight receipt absent; launch through submit.sh")
    receipt = json.loads(path.read_text())
    require(receipt.get("status") == "passed" and len(receipt.get("checks", [])) == 4,
            "Material initialization has not passed for all four prescriptions")
    require(receipt.get("identity") == identity(bundle, runtime), "Material preflight runtime or frozen inputs changed")
    manifest = json.loads((bundle/"manifest.json").read_text())
    require([row.get("model_id") for row in receipt["checks"]] == [m["id"] for m in prescriptions(manifest)],
            "Material receipt does not cover the required prescriptions")
    for check in receipt["checks"]:
        require(check.get("status") == "passed", "Incomplete material check")
        for relative, checksum in check["artifact_sha256"].items():
            path = bundle / relative
            require(not Path(relative).is_absolute() and ".." not in Path(relative).parts
                    and path.resolve().is_relative_to(bundle.resolve()) and not path.is_symlink()
                    and path.is_file() and digest(path) == checksum, "Material check artifacts changed")
    return receipt


def run_checks(bundle, machine):
    bundle = Path(bundle).resolve()
    sys.path.insert(0, str(bundle/"code"))
    from production_task import validate_package, load_runner
    experiment = validate_package(bundle)
    require(experiment["experiment_id"] == NAME, "Wrong material-check experiment")
    runner = load_runner(bundle)
    from mcfost_grid.physics import physical_args
    runtime = runner.runtime_config(machine)
    manifest = json.loads((bundle/"manifest.json").read_text())
    root = bundle / "material_preflight"
    root.mkdir(exist_ok=True)
    with (root/".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (bundle/RECEIPT).is_file():
            previous = json.loads((bundle/RECEIPT).read_text())
            if previous.get("status") == "passed":
                validate_receipt(bundle, runtime)
                return dict(status="cached", receipt=str(bundle/RECEIPT))
        wave = wavelength_grid(bundle, manifest)
        state = dict(schema_version=1, status="running", identity=identity(bundle, runtime),
                     started_utc=runner.now(), checks=[], threads=2, timeout_seconds_per_case=120,
                     scope="Isolated envelope dust initialization only; disk omitted, exact DHS species retained",
                     fresh_temperatures_computed=False, images_computed=False,
                     optical_constants_modified=False, tail_physical_accuracy_certified=False)
        runner.atomic_json(bundle/RECEIPT, state)
        try:
            for model in prescriptions(manifest):
                attempts = root/model["id"]
                attempts.mkdir(exist_ok=True)
                number = 1
                while (attempts/f"attempt_{number:03d}").exists():
                    number += 1
                attempt = attempts/f"attempt_{number:03d}"
                attempt.mkdir()
                parameter = isolated_parameter((bundle/"models"/model["id"]/"temperature.para").read_text())
                (attempt/"check.para").write_text(parameter)
                (attempt/"check.lambda").write_text(str(len(wave))+"\n"+"\n".join(f"{v:.17g}" for v in wave)+"\n")
                # In MCFOST 4.1.14, -dust_prop creates root_dir/data_dust but
                # writes FITS to cwd/data_dust. The unique attempt already
                # isolates outputs, so retain the simulator's default root '.'.
                command = [runtime["mcfost_executable"], "check.para", "-dust_prop",
                           "-max_mem", "2", "-no_backup",
                           *physical_args(model["parameters"])]
                started = time.monotonic()
                check = dict(model_id=model["id"], material=model["parameters"]["envelope_silicate_file"],
                             amax_um=model["parameters"]["envelope_amax_um"], status="running", command=command)
                state["checks"].append(check)
                runner.atomic_json(bundle/RECEIPT, state)
                with (attempt/"mcfost.log").open("w") as log:
                    outcome = subprocess.run(command, cwd=attempt,
                        env=runner.environment(bundle, {**runtime, "threads": 2}),
                        stdout=log, stderr=subprocess.STDOUT, timeout=120, check=False)
                check.update(returncode=outcome.returncode, elapsed_seconds=time.monotonic()-started)
                validate_completion(attempt/"mcfost.log", outcome.returncode)
                try:
                    check["diagnostics"] = check_outputs(attempt, wave)
                except (ValueError, OSError) as exc:
                    raise ValueError(f"Material initialization products failed checks: {exc}.\n"
                                     f"{log_context(attempt/'mcfost.log')}") from exc
                check.update(status="passed", artifact_sha256={str(p.relative_to(bundle)): digest(p)
                    for p in sorted(attempt.rglob("*")) if p.is_file()})
                runner.atomic_json(bundle/RECEIPT, state)
            require(state["identity"] == identity(bundle, runtime), "Runtime inputs changed during preflight")
            state.update(status="passed", finished_utc=runner.now())
        except Exception as exc:
            state.update(status="failed", error=f"{type(exc).__name__}: {exc}", finished_utc=runner.now())
            runner.atomic_json(bundle/RECEIPT, state)
            raise
        runner.atomic_json(bundle/RECEIPT, state)
    return dict(status="passed", prescriptions=4, receipt=str(bundle/RECEIPT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--machine", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_checks(args.bundle, args.machine), indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
