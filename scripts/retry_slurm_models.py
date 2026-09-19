#!/usr/bin/env python3
"""Prepare retries using this Python environment; --submit explicitly submits them."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile


def positive_job_id(value):
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise argparse.ArgumentTypeError("job ID must be a positive integer (the parent array ID)")
    return value


def selected_indices(value, model_count):
    """Expand bounded, comma-separated model indices and inclusive ranges."""
    if not re.fullmatch(r"[0-9]+(?:-[0-9]+)?(?:,[0-9]+(?:-[0-9]+)?)*", value):
        raise ValueError("--indices must contain non-negative indices or inclusive ranges, e.g. 8-10,16-95")
    indices = set()
    for item in value.split(","):
        ends = item.split("-")
        start, end = int(ends[0]), int(ends[-1])
        if start > end:
            raise ValueError(f"Reversed model index range: {item}")
        if end >= model_count:
            raise ValueError(f"An index is outside this model catalogue: {item}")
        indices.update(range(start, end + 1))
    return sorted(indices)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--machine", type=Path, required=True)
    parser.add_argument("--indices", required=True, help="Comma-separated model indices or ranges, e.g. 8-10,16-95")
    parser.add_argument("--after-job", action="append", default=[], type=positive_job_id,
                        help="Wait for this job/parent array to end before retries start; repeat for multiple jobs")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    run, machine_path = args.run.resolve(), args.machine.resolve()
    manifest = json.loads((run / "manifest.json").read_text())
    machine = json.loads(machine_path.read_text())
    try:
        indices = selected_indices(args.indices, len(manifest["models"]))
    except ValueError as exc:
        parser.error(str(exc))
    after_jobs = list(dict.fromkeys(args.after_job))
    # Keep a virtual environment's bin/python symlink path: resolve() can select
    # its base interpreter and thereby lose the environment's site-packages.
    python = os.path.abspath(sys.executable)
    probe = (
        "import sys, socket, json, numpy, scipy, astropy, matplotlib; "
        "from astropy.io import fits; from scipy.ndimage import gaussian_filter; "
        "sys.path.insert(0, 'code/src'); "
        "from mcfost_grid.photometry import measure_image; "
        "from mcfost_grid.analysis import analyze_run; "
        "from mcfost_grid.runner import validate_inputs; "
        "from pathlib import Path; "
        "validate_inputs(Path('.'), json.loads(Path('manifest.json').read_text())); "
        "print(json.dumps({'host': socket.gethostname(), 'python': sys.executable, "
        "'packages': {p.__name__: {'version': p.__version__, 'file': p.__file__} "
        "for p in (numpy, scipy, astropy, matplotlib)}}), flush=True)"
    )
    env = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"
    subprocess.run([python, "-B", "-c", probe], cwd=run, env=env, check=True)
    if args.submit:
        # Snapshot every currently queued/running job with this user's run name,
        # including the ORIGINAL analysis. Waiting on arrays alone can allow that
        # older analysis to overwrite the final report after a recovery finishes.
        queue = subprocess.run([
            "squeue", "--noheader", "--array", "--states=all", f"--user={os.getuid()}",
            f"--name={manifest['run_id'][:80]}", "--format=%i",
        ], capture_output=True, text=True, check=True).stdout
        prior_jobs = sorted(set(queue.split()))
        if any(not re.fullmatch(r"\d+(_\d+)?", job) for job in prior_jobs):
            raise ValueError(f"Unexpected scheduler job IDs: {prior_jobs!r}")
    else:
        prior_jobs = []
    # Recovery files are outside the hashed scientific snapshot. Existing scripts,
    # machine configuration, model inputs and running jobs are left intact.
    recovery_root = run / "recovery"
    recovery_root.mkdir(exist_ok=True)
    recovery = Path(tempfile.mkdtemp(prefix="python-retry-", dir=recovery_root))
    env.update(MCFOST_RUN_DIR=str(run), MCFOST_MACHINE=str(machine_path))
    prefix = f"exec {shlex.quote(machine.get('slurm', {}).get('python', 'python'))} -B code/workflow.py "
    paths = {}
    for kind, source in (("array", "job_array.sh"), ("analysis", "job_analysis.sh")):
        lines = (run / source).read_text().splitlines()
        matches = [i for i, line in enumerate(lines) if line.startswith(prefix)]
        if len(matches) != 1:
            raise ValueError(f"{source}: interpreter differs from machine configuration; regenerate Slurm scripts first")
        i = matches[0]
        lines[i:i+1] = [
            f"{shlex.quote(python)} -B -c {shlex.quote(probe)}",
            f"exec {shlex.quote(python)} -B code/workflow.py " + lines[i][len(prefix):],
        ]
        if kind == "array":
            throttle = min(len(indices), int(machine.get("slurm", {}).get("max_parallel", 1)))
            selection = ",".join(map(str, indices)) + f"%{throttle}"
            array_rows = [i for i, line in enumerate(lines) if line.startswith("#SBATCH --array=")]
            if len(array_rows) != 1:
                raise ValueError("Expected one array directive")
            lines[array_rows[0]] = "#SBATCH --array=" + selection
        path = recovery / source
        path.write_text("\n".join(lines) + "\n")
        subprocess.run(["bash", "-n", str(path)], check=True)
        paths[kind] = path
    receipt_path = recovery / "recovery.json"
    receipt = dict(run=str(run), machine=str(machine_path), python=python,
                   indices=indices, prior_jobs=prior_jobs, after_jobs=after_jobs, submitted=False)
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Prepared {recovery}", flush=True)
    if not args.submit:
        print("No jobs submitted. Rerun this command with --submit on the cluster to submit.")
        return
    (run / "logs").mkdir(exist_ok=True)

    def submit(kind, dependency=None):
        command = ["sbatch", "--parsable", "--export=ALL"]
        if dependency:
            command.append("--dependency=afterany:" + ":".join(dependency))
        command.append(str(paths[kind]))
        result = subprocess.run(command, cwd=run, env=env, capture_output=True, text=True, check=True)
        job = result.stdout.strip().split(";")[0]
        if not re.fullmatch(r"\d+", job):
            raise RuntimeError(f"Unexpected sbatch response (check the queue before retrying): {result.stdout!r}")
        return job

    retry = submit("array", after_jobs)
    analysis_dependencies = list(dict.fromkeys([*prior_jobs, *after_jobs, retry]))
    receipt.update(submitted=True, retry_array_job_id=retry, analysis_dependencies=analysis_dependencies)
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Submitted retry array {retry}, indices {args.indices}", flush=True)
    try:
        analysis = submit("analysis", analysis_dependencies)
    except (subprocess.CalledProcessError, RuntimeError, OSError) as exc:
        detail = getattr(exc, "stderr", None) or str(exc)
        receipt["analysis_submission_error"] = detail
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        print(f"Retry {retry} IS submitted, but final analysis submission failed: {detail}", file=sys.stderr)
        print("Do not resubmit the retry array. After all run jobs end, run analysis using the pinned Python above.", file=sys.stderr)
        raise
    receipt["analysis_job_id"] = analysis
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Submitted final analysis {analysis}; waits for original run jobs and retry {retry}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        if isinstance(exc, subprocess.CalledProcessError) and exc.stderr:
            print(exc.stderr, file=sys.stderr)
        raise SystemExit(1)
