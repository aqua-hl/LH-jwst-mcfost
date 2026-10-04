#!/usr/bin/env python3
"""Pin the cluster venv and prepare Stage 0 -> 26-run array -> analysis.

Configuration validates imports and frozen inputs. It does not execute MCFOST
or submit jobs. The production directory is created only after measured DHS
opacities pass on the cluster.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from configure_extinction_ice_numerics_v2 import BLAS_VARIABLES, simulator_paths, validate_settings


def dependency_probe(bundle, executable, utilities, expected_runtime=None):
    code = f"""import hashlib, json, os, socket, sys
from pathlib import Path
import numpy, scipy, astropy, matplotlib
from astropy.io import fits
from scipy.ndimage import gaussian_filter
bundle = Path({str(bundle)!r})
sys.path.insert(0, str(bundle / 'code'))
from build_silicate_search_v3 import validate_bootstrap
plan = validate_bootstrap(bundle)
sys.path.insert(0, str(bundle / 'code/src'))
from mcfost_grid.runner import _utilities_hash
from mcfost_grid.photometry import measure_image
exe, utilities = Path({executable!r}), Path({utilities!r})
if not exe.is_file() or not os.access(exe, os.X_OK):
    raise RuntimeError(f'MCFOST unavailable on this node: {{exe}}')
if not all((utilities / name).is_dir() for name in ('Dust', 'Lambda', 'Stellar_Spectra')):
    raise RuntimeError(f'Incomplete MCFOST utilities: {{utilities}}')
runtime = dict(executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
               utilities_content_sha256=_utilities_hash(dict(mcfost_utils=str(utilities))))
expected = {expected_runtime!r}
if expected is not None and runtime != expected:
    raise RuntimeError('MCFOST binary or utilities changed since cluster configuration')
print(json.dumps(dict(environment_check='passed', host=socket.gethostname(),
    python=os.path.abspath(sys.executable), mcfost_executable=str(exe), mcfost_utils=str(utilities),
    runtime_identity=runtime, frozen_declared_models=plan['model_count'],
    packages={{p.__name__: dict(version=p.__version__, file=p.__file__)
              for p in (numpy, scipy, astropy, matplotlib)}})), flush=True)
"""
    return code


def launch_scripts(bundle, machine, configured, runtime_identity):
    bundle, configured = Path(bundle), Path(configured)
    settings = machine["slurm"]
    python = shlex.quote(settings["python"])
    machine_text = json.dumps(machine, indent=2, allow_nan=False) + "\n"
    checksum = hashlib.sha256(machine_text.encode()).hexdigest()
    check_code = ("import hashlib,pathlib,sys; expected=" + repr(checksum)
                  + "; actual=hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest();"
                  + "sys.exit(0 if actual==expected else 'Machine configuration changed; configure a new template')")
    check_machine = f'{python} -B -c {shlex.quote(check_code)} "$MCFOST_MACHINE"'
    probe = shlex.quote(dependency_probe(bundle, machine["mcfost_executable"],
                                        machine["mcfost_utils"], runtime_identity))
    common = ["#!/bin/bash", "#SBATCH --job-name=silicate_search_v3", "#SBATCH --nodes=1", "#SBATCH --ntasks=1"]
    for name in ("partition", "account", "reservation"):
        if settings.get(name):
            common.append(f"#SBATCH --{name}={settings[name]}")
    body = ["", "set -euo pipefail",
        ': "${MCFOST_BUNDLE_DIR:?Use submit.sh to set the absolute bundle path}"',
        ': "${MCFOST_MACHINE:?Use submit.sh to set the pinned machine JSON}"',
        f'[[ "$MCFOST_BUNDLE_DIR" == {shlex.quote(str(bundle))} ]] || {{ echo "Bundle moved: configure again" >&2; exit 2; }}',
        f'[[ "$MCFOST_MACHINE" == {shlex.quote(str(configured))} ]] || {{ echo "Wrong machine configuration" >&2; exit 2; }}',
        'cd "$MCFOST_BUNDLE_DIR"', *settings["setup_lines"],
        f'export MCFOST_AUTO_UPDATE=0 MCFOST_UTILS={shlex.quote(machine["mcfost_utils"])}',
        'export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE',
        'export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1',
        'export MPLCONFIGDIR="$MCFOST_BUNDLE_DIR/.mplconfig"', 'mkdir -p "$MPLCONFIGDIR"',
        check_machine, f"{python} -B -c {probe}"]
    dispatcher = f'exec {python} -B code/silicate_search_v3_task.py "$MCFOST_BUNDLE_DIR"'
    scripts = {
        "job_stage0.sh": common + ["#SBATCH --cpus-per-task=2", "#SBATCH --mem=8G", "#SBATCH --time=00:30:00",
            "#SBATCH --output=logs/stage0-%j.out", "#SBATCH --error=logs/stage0-%j.err"] + body + [
            dispatcher + ' --stage0 --machine "$MCFOST_MACHINE"'],
        "job_array.sh": common + ["#SBATCH --cpus-per-task=64",
            f"#SBATCH --array=0-25%{settings['max_parallel']}", f"#SBATCH --mem={settings['memory']}",
            f"#SBATCH --time={settings['time']}", "#SBATCH --output=logs/array-%A_%a.out",
            "#SBATCH --error=logs/array-%A_%a.err"] + body + [
            dispatcher + ' --index "$SLURM_ARRAY_TASK_ID" --machine "$MCFOST_MACHINE"'],
        "job_analysis.sh": common + [f"#SBATCH --cpus-per-task={settings['analysis_cpus']}",
            f"#SBATCH --mem={settings['analysis_memory']}", f"#SBATCH --time={settings['analysis_time']}",
            "#SBATCH --output=logs/analysis-%j.out", "#SBATCH --error=logs/analysis-%j.err"] + body + [
            dispatcher + " --analyze"],
    }
    scripts["submit.sh"] = ["#!/bin/bash", "set -euo pipefail",
        'if [[ $# -ne 1 ]]; then echo "Usage: bash submit.sh /absolute/path/to/machine.cluster.json" >&2; exit 2; fi',
        'export MCFOST_BUNDLE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"',
        f'[[ "$MCFOST_BUNDLE_DIR" == {shlex.quote(str(bundle))} ]] || {{ echo "Bundle moved: configure again" >&2; exit 2; }}',
        f'[[ "$1" == {shlex.quote(str(configured))} ]] || {{ echo "Use the pinned cluster machine JSON" >&2; exit 2; }}',
        'export MCFOST_MACHINE="$1"', 'cd "$MCFOST_BUNDLE_DIR"', 'mkdir -p logs submissions', check_machine,
        'submission_dir="submissions/$(date -u +%Y%m%dT%H%M%SZ)_$$"', 'mkdir "$submission_dir"',
        'printf "%s\\n" "$MCFOST_MACHINE" > "$submission_dir/machine_path.txt"',
        'sbatch --parsable --export=ALL job_stage0.sh > "$submission_dir/stage0_sbatch.out"',
        'stage0_job=$(<"$submission_dir/stage0_sbatch.out"); stage0_job=${stage0_job%%;*}',
        '[[ "$stage0_job" =~ ^[0-9]+$ ]] || { echo "Unexpected Stage 0 sbatch output; inspect $submission_dir" >&2; exit 1; }',
        'printf "%s\\n" "$stage0_job" > "$submission_dir/stage0_job_id.txt"',
        'echo "Submitted Stage 0 $stage0_job; four DHS species must pass before the catalogue is built."',
        'sbatch --parsable --export=ALL --dependency="afterok:$stage0_job" --kill-on-invalid-dep=yes job_array.sh > "$submission_dir/array_sbatch.out"',
        'array_job=$(<"$submission_dir/array_sbatch.out"); array_job=${array_job%%;*}',
        '[[ "$array_job" =~ ^[0-9]+$ ]] || { echo "Unexpected array sbatch output; inspect $submission_dir" >&2; exit 1; }',
        'printf "%s\\n" "$array_job" > "$submission_dir/array_job_id.txt"',
        'echo "Submitted production array $array_job (afterok:$stage0_job); receipt: $MCFOST_BUNDLE_DIR/$submission_dir"',
        'if ! sbatch --parsable --export=ALL --dependency="afterany:$array_job" job_analysis.sh > "$submission_dir/analysis_sbatch.out" 2> "$submission_dir/analysis_sbatch.err"; then',
        '  echo "Array $array_job was submitted; analysis submission failed. Keep this receipt and submit analysis separately." >&2',
        '  cat "$submission_dir/analysis_sbatch.err" >&2; exit 1', 'fi',
        'analysis_job=$(<"$submission_dir/analysis_sbatch.out"); analysis_job=${analysis_job%%;*}',
        '[[ "$analysis_job" =~ ^[0-9]+$ ]] || { echo "Unexpected analysis sbatch output; inspect $submission_dir" >&2; exit 1; }',
        'printf "%s\\n" "$analysis_job" > "$submission_dir/analysis_job_id.txt"',
        'echo "Submitted analysis $analysis_job (afterany:$array_job); partial results remain marked incomplete."']
    return scripts


def configure(bundle, machine_path):
    from build_silicate_search_v3 import validate_bootstrap
    bundle, machine_path = Path(bundle).resolve(), Path(machine_path).resolve()
    validate_bootstrap(bundle)
    machine = json.loads(machine_path.read_text())
    machine.setdefault("slurm", {}).setdefault("analysis_time", "02:00:00")
    validate_settings(machine, {"experiment_id": "silicate_search_v3"})
    executable, utilities = simulator_paths(machine, machine_path)
    python = os.path.abspath(sys.executable)  # Preserve venv identity through symlinks.
    machine.update(mcfost_executable=executable, mcfost_utils=utilities)
    machine["slurm"]["python"] = python
    configured = machine_path.with_name(machine_path.stem + ".cluster.json")
    serialized = json.dumps(machine, indent=2, allow_nan=False) + "\n"
    if configured.exists() and configured.read_text() != serialized:
        raise FileExistsError(f"Refusing to replace different cluster configuration: {configured}. "
                              "Existing jobs may use it; copy the template to a new filename first.")
    environment = os.environ.copy()
    environment.pop("SLURM_ARRAY_TASK_ID", None)
    environment.update({name: "1" for name in BLAS_VARIABLES})
    environment.update(MCFOST_AUTO_UPDATE="0", MCFOST_UTILS=utilities,
                       MPLCONFIGDIR=str(bundle / ".mplconfig"))
    try:
        result = subprocess.run([python, "-B", "-c", dependency_probe(bundle, executable, utilities)],
                                cwd=bundle, env=environment, text=True, capture_output=True, timeout=180, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("Cluster dependency check failed:\n" + (exc.stderr or exc.stdout or str(exc))[-4000:]) from exc
    receipt = json.loads(result.stdout.strip().splitlines()[-1])
    if receipt.get("environment_check") != "passed" or not receipt.get("runtime_identity"):
        raise ValueError("Cluster dependency probe did not produce a valid runtime receipt")
    print(result.stdout, end="")
    scripts = launch_scripts(bundle, machine, configured, receipt["runtime_identity"])
    for lines in scripts.values():
        subprocess.run(["bash", "-n"], input="\n".join(lines) + "\n", text=True, check=True)
    if not configured.exists():
        with configured.open("x") as stream:
            stream.write(serialized)
    (bundle / "logs").mkdir(exist_ok=True)
    for name, lines in scripts.items():
        (bundle / name).write_text("\n".join(lines) + "\n")
    (bundle / "cluster_environment.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return configured


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--machine", type=Path, required=True)
    args = parser.parse_args(argv)
    configured = configure(args.bundle, args.machine)
    print("Validated frozen bootstrap and shared Python; no MCFOST calculation or job was launched.")
    print(f"bash {shlex.quote(str(args.bundle.resolve() / 'submit.sh'))} {shlex.quote(str(configured))}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
