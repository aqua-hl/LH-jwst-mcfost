"""Generate (never submit) relocatable Slurm model arrays and analysis scripts."""
from pathlib import Path
import re
import shlex

from .configuration import load_json, _strict_keys


def write_slurm(run_dir, machine_path):
    run = Path(run_dir).resolve()
    manifest = load_json(run / "manifest.json")
    machine = load_json(machine_path)
    settings = machine.get("slurm", {})
    _strict_keys(settings, {"partition", "account", "reservation", "time", "memory", "max_parallel",
                           "python", "setup_lines", "analysis_cpus", "analysis_memory", "analysis_time"}, "slurm")
    threads = machine.get("threads", 16)
    parallel = settings.get("max_parallel", 4)
    analysis_cpus = settings.get("analysis_cpus", 2)
    for key, value in (("threads", threads), ("max_parallel", parallel), ("analysis_cpus", analysis_cpus)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{key} must be a positive integer")
    python = settings.get("python", "python")
    for value in (python, *[settings[k] for k in ("partition", "account", "reservation", "time", "memory", "analysis_time", "analysis_memory") if k in settings]):
        if not isinstance(value, str) or "\n" in value or "\r" in value:
            raise ValueError("Slurm settings must be single-line strings")
    n = len(manifest["models"])
    common = ["#!/bin/bash", f"#SBATCH --job-name={manifest['run_id'][:80]}", "#SBATCH --nodes=1", "#SBATCH --ntasks=1"]
    for key in ("partition", "account", "reservation"):
        if settings.get(key):
            if not re.fullmatch(r"[A-Za-z0-9_.:,/-]+", settings[key]):
                raise ValueError(f"Invalid Slurm {key}")
            common.append(f"#SBATCH --{key}={settings[key]}")
    setup = settings.get("setup_lines", [])
    if not isinstance(setup, list) or not all(isinstance(s, str) and "\n" not in s for s in setup):
        raise ValueError("setup_lines must be a list of shell command lines")
    # The script lives in Slurm's spool when executing. Explicit exported paths
    # from submit.sh avoid the old dirname($0) /var/lib/slurm failure.
    body = ["", "set -euo pipefail", ': "${MCFOST_RUN_DIR:?Use submit.sh, or export absolute MCFOST_RUN_DIR}"',
            ': "${MCFOST_MACHINE:?Set MCFOST_MACHINE to the cluster machine JSON}"',
            'cd "$MCFOST_RUN_DIR"', *setup,
            'export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1',
            'export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE MCFOST_AUTO_UPDATE=0',
            'export MPLCONFIGDIR="$MCFOST_RUN_DIR/.mplconfig"', 'mkdir -p "$MPLCONFIGDIR"']
    array = common + [f"#SBATCH --cpus-per-task={threads}", f"#SBATCH --array=0-{n-1}%{parallel}",
                      f"#SBATCH --time={settings.get('time', '12:00:00')}", f"#SBATCH --mem={settings.get('memory', '16G')}",
                      "#SBATCH --output=logs/array-%A_%a.out", "#SBATCH --error=logs/array-%A_%a.err"] + body + [
        f'exec {shlex.quote(python)} -B code/workflow.py task "$MCFOST_RUN_DIR" --index "$SLURM_ARRAY_TASK_ID" --machine "$MCFOST_MACHINE"']
    analysis = common + [f"#SBATCH --cpus-per-task={analysis_cpus}", f"#SBATCH --time={settings.get('analysis_time', '00:30:00')}",
                         f"#SBATCH --mem={settings.get('analysis_memory', '4G')}", "#SBATCH --output=logs/analysis-%j.out",
                         "#SBATCH --error=logs/analysis-%j.err"] + body + [
        f'exec {shlex.quote(python)} -B code/workflow.py analyze "$MCFOST_RUN_DIR" --workers "$SLURM_CPUS_PER_TASK"']
    submit = ['#!/bin/bash', 'set -euo pipefail',
              'if [[ $# -ne 1 ]]; then echo "Usage: bash submit.sh /absolute/path/to/machine.cluster.json" >&2; exit 2; fi',
              'export MCFOST_RUN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"',
              'case "$1" in /*) export MCFOST_MACHINE="$1";; *) echo "Use an absolute machine JSON path" >&2; exit 2;; esac',
              'test -f "$MCFOST_MACHINE"', 'cd "$MCFOST_RUN_DIR"', 'mkdir -p logs',
              'job=$(sbatch --parsable --export=ALL job_array.sh)', 'job=${job%%;*}',
              'echo "Submitted model array $job"',
              'sbatch --parsable --export=ALL --dependency="afterany:$job" job_analysis.sh',
              '# afterany makes the summary report failures too; incomplete models are never ranked.']
    (run / "logs").mkdir(exist_ok=True)
    for name, lines in (("job_array.sh", array), ("job_analysis.sh", analysis), ("submit.sh", submit)):
        (run / name).write_text("\n".join(lines) + "\n")
    return {"run": str(run), "models": n, "threads_per_model": threads, "max_concurrent_models": parallel,
            "max_cpus": parallel * threads, "submitted": False,
            "note": "Regenerate scripts after editing machine resources; setup_lines run as your own shell code."}
