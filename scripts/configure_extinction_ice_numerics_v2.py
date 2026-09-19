#!/usr/bin/env python3
"""Pin a shared cluster environment and prepare a V2 diagnostic array.

Run with the successful shared Python environment after uploading the complete
bundle. This helper checks dependencies and immutable inputs, writes a sibling
*.cluster.json, and generates launch scripts. It never submits or runs MCFOST.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


BLAS_VARIABLES = ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS")
CROSSOVER_EXPERIMENT = "extinction_ice_crossover_v2"
PIXEL_SCALE_DIAGNOSTIC = "fixed_temperature_pixel_scale_v2"


def is_crossover(experiment):
    return experiment.get("experiment_id") == CROSSOVER_EXPERIMENT


def simulator_paths(machine, machine_path):
    executable = os.path.expanduser(os.path.expandvars(
        machine.get("mcfost_executable", os.environ.get("MCFOST_EXE", "mcfost"))))
    if "/" in executable:
        candidate = Path(executable)
        executable = str((candidate if candidate.is_absolute()
                          else machine_path.parent / candidate).resolve())
    else:
        executable = shutil.which(executable)
    if not executable or not Path(executable).is_file() or not os.access(executable, os.X_OK):
        raise ValueError("MCFOST executable is unavailable; load the working MCFOST environment first")
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
    return str(Path(executable).resolve()), str(utilities)


def bundle_layout(bundle):
    experiment = json.loads((bundle / "experiment.json").read_text())
    tasks = experiment.get("tasks", [])
    if is_crossover(experiment):
        if experiment.get("schema_version") != 2 \
                or experiment.get("diagnostic_id") != PIXEL_SCALE_DIAGNOSTIC:
            raise ValueError("Expected schema 2 fixed-temperature pixel-scale diagnostic; "
                             "prepare a fresh bundle instead of launching the old temperature crossover")
        if len(tasks) != 5 or [task.get("index") for task in tasks] != list(range(5)):
            raise ValueError("Expected exactly five consecutively indexed pixel-scale tasks")
        for relative in ("code/crossover_task.py", "code/analyze_extinction_ice_crossover_v2.py",
                         "code/src/mcfost_grid/runner.py", "code/src/mcfost_grid/photometry.py"):
            if not (bundle / relative).is_file():
                raise ValueError(f"Incomplete crossover bundle: {relative} is absent")
        # The crossover dispatcher checks every frozen temperature, image input,
        # task design and package hash during the dependency probe below.
        return experiment
    if len(tasks) != 20 or [task.get("index") for task in tasks] != list(range(20)):
        raise ValueError("Expected exactly 20 consecutively indexed numerical tasks")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", experiment.get("experiment_id", "")):
        raise ValueError("Invalid experiment_id")
    for relative in ("code/numerical_task.py", "code/analyze_extinction_ice_numerics_v2.py"):
        if not (bundle / relative).is_file():
            raise ValueError(f"Incomplete numerical bundle: {relative} is absent")
    runs = []
    for task in tasks:
        relative = task.get("run_path")
        if not isinstance(relative, str) or Path(relative).is_absolute():
            raise ValueError("Task run_path must be relative to the bundle")
        run = (bundle / relative).resolve()
        if not run.is_relative_to(bundle) or not (run / "manifest.json").is_file() \
                or not (run / "code/workflow.py").is_file():
            raise ValueError(f"Invalid or incomplete task subrun: {relative}")
        if run not in runs:
            runs.append(run)
    if len(runs) != 10:
        raise ValueError("Expected ten independent seed/photon subruns")
    return experiment


def validate_settings(machine, experiment=None):
    allowed = {"schema_version", "mcfost_executable", "mcfost_utils", "backend", "threads",
               "timeout_seconds", "max_memory_gb", "slurm", "notes"}
    if set(machine) - allowed or machine.get("schema_version") != 1:
        raise ValueError("Unsupported machine schema or keys")
    if machine.get("backend", "image_method2") != "image_method2":
        raise ValueError("This numerical package requires the image_method2 backend")
    settings = machine.setdefault("slurm", {})
    if not isinstance(settings, dict):
        raise ValueError("slurm must be an object")
    if set(settings) - {"partition", "account", "reservation", "time", "memory", "max_parallel",
                        "python", "setup_lines", "analysis_cpus", "analysis_memory", "analysis_time"}:
        raise ValueError("Unsupported Slurm setting")
    machine.setdefault("threads", 64)
    machine.setdefault("timeout_seconds", 14400)
    pixel_scale = is_crossover(experiment or {})
    machine.setdefault("max_memory_gb", 64 if pixel_scale else 12)
    settings.setdefault("max_parallel", 16)
    settings.setdefault("analysis_cpus", 2)
    settings.setdefault("time", "24:00:00")
    settings.setdefault("memory", "96G" if pixel_scale else "16G")
    settings.setdefault("analysis_time", "01:00:00")
    settings.setdefault("analysis_memory", "8G" if pixel_scale else "4G")
    for name, value in (("threads", machine["threads"]), ("max_parallel", settings["max_parallel"]),
                        ("analysis_cpus", settings["analysis_cpus"])):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if machine["threads"] != 64 or settings["max_parallel"] > 16:
        raise ValueError("The experiment uses 64 CPUs per task and at most 16 concurrent tasks")
    for name in ("timeout_seconds", "max_memory_gb"):
        value = machine[name]
        if isinstance(value, bool) or not isinstance(value, (float, int)) \
                or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive")
    for name in ("partition", "account", "reservation"):
        if name in settings and settings[name] and (not isinstance(settings[name], str)
                or not re.fullmatch(r"[A-Za-z0-9_.:,/-]+", settings[name])):
            raise ValueError(f"Invalid Slurm {name}")
    for name in ("time", "analysis_time"):
        if not isinstance(settings[name], str) or not re.fullmatch(r"(?:\d+-)?\d+:\d{2}:\d{2}", settings[name]):
            raise ValueError(f"Invalid Slurm {name}; use [days-]hours:minutes:seconds")
        parts = settings[name].split(":")
        if int(parts[-1]) > 59 or int(parts[-2]) > 59 \
                or not any(int(value) for value in re.split(r"[-:]", settings[name])):
            raise ValueError(f"Invalid Slurm {name}")
    for name in ("memory", "analysis_memory"):
        if not isinstance(settings[name], str) or not re.fullmatch(r"[1-9]\d*[KMGT]?", settings[name]):
            raise ValueError(f"Invalid Slurm {name}; use a positive integer with optional K/M/G/T")
    setup = settings.setdefault("setup_lines", [])
    if not isinstance(setup, list) or not all(isinstance(line, str) and "\n" not in line
                                            and "\r" not in line for line in setup):
        raise ValueError("setup_lines must be a list of single-line shell commands")


def dependency_probe(executable, utilities, experiment=None):
    """Read-only login/compute probe: package hashes, all subruns and imports."""
    if is_crossover(experiment or {}):
        frozen_check = """from crossover_task import validate_package, load_modules
validate_package(bundle)
runner, photometry = load_modules(bundle)
if not callable(photometry.measure_image):
    raise RuntimeError('Frozen crossover photometry is unavailable')
experiment = json.loads((bundle / 'experiment.json').read_text())
frozen_receipt = {'frozen_crossover_tasks_checked': len(experiment['tasks'])}
"""
    else:
        frozen_check = """from numerical_task import validate_package
validate_package(bundle)
experiment = json.loads((bundle / 'experiment.json').read_text())
runs = list(dict.fromkeys(task['run_path'] for task in experiment['tasks']))
sys.path.insert(0, str(bundle / runs[0] / 'code/src'))
from mcfost_grid.photometry import measure_image
from mcfost_grid.runner import validate_inputs
for relative in runs:
    run = (bundle / relative).resolve()
    if not run.is_relative_to(bundle):
        raise RuntimeError('Subrun escapes bundle')
    validate_inputs(run, json.loads((run / 'manifest.json').read_text()))
frozen_receipt = {'frozen_subruns_checked': len(runs)}
"""
    code = f"""import sys, os, socket, json
from pathlib import Path
import numpy, scipy, astropy, matplotlib
from astropy.io import fits
from scipy.ndimage import gaussian_filter
bundle = Path('.').resolve()
sys.path.insert(0, str(bundle / 'code'))
{frozen_check}
exe = Path({executable!r})
utilities = Path({utilities!r})
if not exe.is_file() or not os.access(exe, os.X_OK):
    raise RuntimeError(f'MCFOST unavailable on this node: {{exe}}')
if not all((utilities / name).is_dir() for name in ('Dust', 'Lambda', 'Stellar_Spectra')):
    raise RuntimeError(f'Incomplete MCFOST utilities on this node: {{utilities}}')
print(json.dumps({{'environment_check': 'passed', 'host': socket.gethostname(),
                  'python': sys.executable, 'mcfost_executable': str(exe), 'mcfost_utils': str(utilities),
                  **frozen_receipt,
                  'packages': {{p.__name__: {{'version': p.__version__, 'file': p.__file__}}
                               for p in (numpy, scipy, astropy, matplotlib)}}}}), flush=True)
"""
    return f"exec({code!r})"


def launch_scripts(bundle, experiment, machine, configured):
    settings = machine["slurm"]
    python = shlex.quote(settings["python"])
    utilities = shlex.quote(machine["mcfost_utils"])
    probe = shlex.quote(dependency_probe(machine["mcfost_executable"], machine["mcfost_utils"], experiment))
    crossover = is_crossover(experiment)
    dispatcher = "crossover_task.py" if crossover else "numerical_task.py"
    analyzer = ("analyze_extinction_ice_crossover_v2.py" if crossover
                else "analyze_extinction_ice_numerics_v2.py")
    last_index = len(experiment["tasks"]) - 1
    concurrency = min(settings["max_parallel"], len(experiment["tasks"]))
    machine_hash = hashlib.sha256((json.dumps(machine, indent=2, allow_nan=False) + "\n").encode()).hexdigest()
    machine_check = ("import hashlib, pathlib, sys; expected = " + repr(machine_hash)
                     + "; actual = hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest(); "
                     + "sys.exit(0 if actual == expected else "
                     + repr("Machine configuration changed: configure a new template instead") + ")")
    check_machine = f'{python} -B -c {shlex.quote(machine_check)} "$MCFOST_MACHINE"'
    common = ["#!/bin/bash", f"#SBATCH --job-name={experiment['experiment_id'][:80]}",
              "#SBATCH --nodes=1", "#SBATCH --ntasks=1"]
    for name in ("partition", "account", "reservation"):
        if settings.get(name):
            common.append(f"#SBATCH --{name}={settings[name]}")
    body = ["", "set -euo pipefail",
            ': "${MCFOST_BUNDLE_DIR:?Use submit.sh to set the absolute bundle path}"',
            ': "${MCFOST_MACHINE:?Use submit.sh to set the pinned machine JSON}"',
            f'[[ "$MCFOST_BUNDLE_DIR" == {shlex.quote(str(bundle))} ]] || {{ echo "Bundle moved: configure it again" >&2; exit 2; }}',
            f'[[ "$MCFOST_MACHINE" == {shlex.quote(str(configured))} ]] || {{ echo "Wrong machine configuration: configure it again" >&2; exit 2; }}',
            'cd "$MCFOST_BUNDLE_DIR"', f"export MCFOST_AUTO_UPDATE=0 MCFOST_UTILS={utilities}",
            *settings["setup_lines"],
            'export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1',
            'export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE',
            f"export MCFOST_AUTO_UPDATE=0 MCFOST_UTILS={utilities}",
            'export MPLCONFIGDIR="$MCFOST_BUNDLE_DIR/.mplconfig"', 'mkdir -p "$MPLCONFIGDIR"',
            check_machine, f"{python} -B -c {probe}"]
    array = common + [f"#SBATCH --cpus-per-task={machine['threads']}",
                      f"#SBATCH --array=0-{last_index}%{concurrency}",
                      f"#SBATCH --time={settings['time']}", f"#SBATCH --mem={settings['memory']}",
                      "#SBATCH --output=logs/array-%A_%a.out", "#SBATCH --error=logs/array-%A_%a.err"] + body + [
        f'exec {python} -B code/{dispatcher} "$MCFOST_BUNDLE_DIR" --index "$SLURM_ARRAY_TASK_ID" --machine "$MCFOST_MACHINE"']
    analysis = common + [f"#SBATCH --cpus-per-task={settings['analysis_cpus']}",
                         f"#SBATCH --time={settings['analysis_time']}",
                         f"#SBATCH --mem={settings['analysis_memory']}",
                         "#SBATCH --output=logs/analysis-%j.out", "#SBATCH --error=logs/analysis-%j.err"] + body + [
        f'exec {python} -B code/{analyzer} "$MCFOST_BUNDLE_DIR"']
    submit = ["#!/bin/bash", "set -euo pipefail",
              'if [[ $# -ne 1 ]]; then echo "Usage: bash submit.sh /absolute/path/to/machine.cluster.json" >&2; exit 2; fi',
              'export MCFOST_BUNDLE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"',
              f'[[ "$MCFOST_BUNDLE_DIR" == {shlex.quote(str(bundle))} ]] || {{ echo "Bundle moved: configure it again" >&2; exit 2; }}',
              f'[[ "$1" == {shlex.quote(str(configured))} ]] || {{ echo "Use the pinned cluster machine JSON shown by configure" >&2; exit 2; }}',
              'export MCFOST_MACHINE="$1"', 'test -f "$MCFOST_MACHINE"',
              'cd "$MCFOST_BUNDLE_DIR"', 'mkdir -p logs submissions',
              check_machine,
              'submission_dir="submissions/$(date -u +%Y%m%dT%H%M%SZ)_$$"',
              'mkdir "$submission_dir"', 'printf "%s\\n" "$MCFOST_MACHINE" > "$submission_dir/machine_path.txt"',
              'sbatch --parsable --export=ALL job_array.sh > "$submission_dir/array_sbatch.out"',
              'job=$(<"$submission_dir/array_sbatch.out")', 'job=${job%%;*}',
              '[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected array sbatch output; inspect $submission_dir before submitting again" >&2; exit 1; }',
              'printf "%s\\n" "$job" > "$submission_dir/array_job_id.txt"',
              'echo "Submitted numerical array $job; receipt: $MCFOST_BUNDLE_DIR/$submission_dir"',
              'if ! sbatch --parsable --export=ALL --dependency="afterany:$job" job_analysis.sh > "$submission_dir/analysis_sbatch.out" 2> "$submission_dir/analysis_sbatch.err"; then',
              '  echo "Array $job was submitted, but analysis submission failed. Keep the receipt and submit analysis separately; do not resubmit the array for this error." >&2',
              '  cat "$submission_dir/analysis_sbatch.err" >&2', '  exit 1', 'fi',
              'analysis_job=$(<"$submission_dir/analysis_sbatch.out")', 'analysis_job=${analysis_job%%;*}',
              '[[ "$analysis_job" =~ ^[0-9]+$ ]] || { echo "Unexpected analysis sbatch output; inspect $submission_dir" >&2; exit 1; }',
              'printf "%s\\n" "$analysis_job" > "$submission_dir/analysis_job_id.txt"',
              'echo "Submitted analysis $analysis_job (afterany:$job); incomplete experiments remain inconclusive."']
    return {"job_array.sh": array, "job_analysis.sh": analysis, "submit.sh": submit}


def configure(bundle, machine_path):
    bundle, machine_path = Path(bundle).resolve(), Path(machine_path).resolve()
    experiment = bundle_layout(bundle)
    machine = json.loads(machine_path.read_text())
    validate_settings(machine, experiment)
    executable, utilities = simulator_paths(machine, machine_path)
    # Preserve the venv executable path; resolving its symlink loses that venv.
    python = os.path.abspath(sys.executable)
    machine.update(mcfost_executable=executable, mcfost_utils=utilities)
    machine["slurm"]["python"] = python
    environment = os.environ.copy()
    environment.update({name: "1" for name in BLAS_VARIABLES})
    environment.update(MCFOST_AUTO_UPDATE="0", MCFOST_UTILS=utilities,
                       MPLCONFIGDIR=str(bundle / ".mplconfig"))
    subprocess.run([python, "-B", "-c", dependency_probe(executable, utilities, experiment)],
                   cwd=bundle, env=environment, check=True, timeout=180)
    configured = machine_path.with_name(machine_path.stem + ".cluster.json")
    serialized = json.dumps(machine, indent=2, allow_nan=False) + "\n"
    if configured.exists() and configured.read_text() != serialized:
        raise FileExistsError(
            f"Refusing to replace a different cluster configuration: {configured}. "
            "Existing jobs may use it. Copy the template to a new filename and configure that instead.")
    scripts = launch_scripts(bundle, experiment, machine, configured)
    # Check every generated script before publishing any launch file.
    for lines in scripts.values():
        subprocess.run(["bash", "-n"], input="\n".join(lines) + "\n", text=True, check=True)
    if not configured.exists():
        with configured.open("x") as stream:
            stream.write(serialized)
    (bundle / "logs").mkdir(exist_ok=True)
    for name, lines in scripts.items():
        (bundle / name).write_text("\n".join(lines) + "\n")
    return configured


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--machine", type=Path, required=True,
                        help="Template JSON; writes a separate sibling *.cluster.json")
    args = parser.parse_args(argv)
    configured = configure(args.bundle, args.machine)
    print(f"Cluster machine configuration: {configured}")
    print("Validated shared Python/imports and all frozen inputs; no MCFOST or Slurm job was launched.")
    print("Submit from the same working cluster environment:")
    print(f"bash {shlex.quote(str(args.bundle.resolve() / 'submit.sh'))} {shlex.quote(str(configured))}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
