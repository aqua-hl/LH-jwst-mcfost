#!/usr/bin/env python3
"""Pin a working cluster environment and generate Slurm scripts; never submit jobs.

Run this with the working Python on the cluster, after loading its modules or
activating its environment. The interpreter, packages, MCFOST executable and
utilities must be visible from compute nodes as well as the login node.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


BLAS_VARIABLES = (
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def simulator_paths(machine, machine_path):
    """Resolve machine-relative paths before writing the sibling configuration."""
    executable = os.path.expanduser(os.path.expandvars(
        machine.get("mcfost_executable", os.environ.get("MCFOST_EXE", "mcfost"))
    ))
    if "/" in executable:
        candidate = Path(executable)
        if not candidate.is_absolute():
            candidate = machine_path.parent / candidate
        executable = str(candidate.resolve())
    else:
        executable = shutil.which(executable)
    if not executable or not Path(executable).is_file() or not os.access(executable, os.X_OK):
        raise ValueError("MCFOST executable is unavailable; load the working MCFOST environment first")
    executable = str(Path(executable).resolve())
    value = machine.get("mcfost_utils", os.environ.get("MCFOST_UTILS"))
    if not value:
        raise ValueError("Set mcfost_utils in the machine JSON or export MCFOST_UTILS")
    utilities = Path(os.path.expanduser(os.path.expandvars(value)))
    if not utilities.is_absolute():
        utilities = machine_path.parent / utilities
    utilities = utilities.resolve()
    missing = [name for name in ("Dust", "Lambda", "Stellar_Spectra")
               if not (utilities / name).is_dir()]
    if missing:
        raise ValueError(f"Incomplete MCFOST utilities: {utilities}; missing {', '.join(missing)}")
    return executable, str(utilities)


def dependency_probe(executable, utilities):
    """One probe for both the login node and every model/analysis job."""
    # No MCFOST process is invoked. Imports run in a fresh interpreter so the
    # modules checked here really belong to the selected Python environment.
    code = f"""import sys, os, socket, json
from pathlib import Path
import numpy, scipy, astropy, matplotlib
from astropy.io import fits
from scipy.ndimage import gaussian_filter
sys.path.insert(0, str(Path('code/src').resolve()))
from mcfost_grid.photometry import measure_image
from mcfost_grid.analysis import analyze_run
from mcfost_grid.runner import validate_inputs
validate_inputs(Path('.'), json.loads(Path('manifest.json').read_text()))
exe = Path({executable!r})
utilities = Path({utilities!r})
if not exe.is_file() or not os.access(exe, os.X_OK):
    raise RuntimeError(f'MCFOST unavailable on this node: {{exe}}')
if not all((utilities / name).is_dir() for name in ('Dust', 'Lambda', 'Stellar_Spectra')):
    raise RuntimeError(f'Incomplete MCFOST utilities on this node: {{utilities}}')
print(json.dumps({{'environment_check': 'passed', 'host': socket.gethostname(),
                  'python': sys.executable, 'mcfost_executable': str(exe), 'mcfost_utils': str(utilities),
                  'packages': {{p.__name__: {{'version': p.__version__, 'file': p.__file__}}
                               for p in (numpy, scipy, astropy, matplotlib)}}}}), flush=True)
"""
    # setup_lines requires a single shell line; repr preserves Python's newlines
    # without embedding literal newlines in the generated shell command.
    return f"exec({code!r})"


def configure(run, machine_path):
    run, machine_path = Path(run).resolve(), Path(machine_path).resolve()
    machine = json.loads(machine_path.read_text())
    frozen_workflow = run / "code/workflow.py"
    if not (run / "manifest.json").is_file() or not frozen_workflow.is_file():
        raise ValueError(f"Prepared run with frozen code is required: {run}")
    executable, utilities = simulator_paths(machine, machine_path)
    # Do not resolve this symlink: venv/bin/python can point to the base Python,
    # whose site-packages would lose the activated virtual environment.
    python = os.path.abspath(sys.executable)
    probe = dependency_probe(executable, utilities)
    environment = os.environ.copy()
    environment.update({name: "1" for name in BLAS_VARIABLES})
    environment.update(MCFOST_AUTO_UPDATE="0", MCFOST_UTILS=utilities,
                       MPLCONFIGDIR=str(run / ".mplconfig"))
    subprocess.run([python, "-B", "-c", probe], cwd=run, env=environment,
                   check=True, timeout=180)

    configured = machine_path.with_name(machine_path.stem + ".cluster.json")
    settings = machine.setdefault("slurm", {})
    setup = settings.get("setup_lines", [])
    if not isinstance(setup, list) or not all(isinstance(line, str) for line in setup):
        raise ValueError("slurm.setup_lines must be a list of shell command lines")
    exports = "export MCFOST_AUTO_UPDATE=0 MCFOST_UTILS=" + shlex.quote(utilities)
    probe_exports = " ".join(name + "=1" for name in BLAS_VARIABLES)
    probe_line = (
        exports + ' MPLCONFIGDIR="$MCFOST_RUN_DIR/.mplconfig"; '
        + probe_exports + " " + shlex.quote(python) + " -B -c " + shlex.quote(probe)
    )
    # Export before existing setup guards, and again immediately before the
    # check in case a module setup modifies these values. Pin Python after any
    # module/environment setup, for both the model array and dependent analysis.
    settings.update(python=python, setup_lines=[exports, *setup, probe_line])
    machine.update(mcfost_executable=executable, mcfost_utils=utilities)
    serialized = json.dumps(machine, indent=2, allow_nan=False) + "\n"
    try:
        with configured.open("x") as stream:
            stream.write(serialized)
    except FileExistsError:
        if configured.read_text() != serialized:
            raise FileExistsError(
                f"Refusing to replace a different cluster configuration: {configured}. "
                "Existing jobs may use it. Copy the template to a new filename and configure that instead."
            ) from None
    subprocess.run([python, "-B", str(frozen_workflow), "slurm", str(run),
                    "--machine", str(configured)], env=environment, check=True)
    for name in ("job_array.sh", "job_analysis.sh", "submit.sh"):
        subprocess.run(["bash", "-n", str(run / name)], check=True)
    return configured


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--machine", type=Path, required=True,
                        help="Template machine JSON; writes a sibling *.cluster.json")
    args = parser.parse_args(argv)
    configured = configure(args.run, args.machine)
    print(f"Cluster machine configuration: {configured}")
    print("No jobs submitted. The selected Python and its packages must be shared with compute nodes.")
    print("Submit from the same working cluster environment:")
    print(f"bash {shlex.quote(str(args.run.resolve() / 'submit.sh'))} {shlex.quote(str(configured))}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
