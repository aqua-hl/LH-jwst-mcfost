#!/usr/bin/env python3
"""Entry point: python workflow.py --help. No command submits cluster jobs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import subprocess

# Set before NumPy/SciPy imports, so Python workers do not multiply BLAS pools.
for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="Freeze a JSON grid and enumerate independent physical models")
    prepare.add_argument("config", type=Path)
    for name in ("run-local", "task", "slurm"):
        p = sub.add_parser(name)
        p.add_argument("run", type=Path)
        p.add_argument("--machine", type=Path, required=True)
        if name == "run-local":
            p.add_argument("--workers", type=int, default=1)
        if name == "task":
            p.add_argument("--index", type=int, required=True)
    p = sub.add_parser("status", help="Read small status records; no scheduler access")
    p.add_argument("run", type=Path)
    p = sub.add_parser("analyze", help="Rank complete compact results and plot against grey JWST spectrum")
    p.add_argument("run", type=Path)
    p.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.command in {"run-local", "task", "analyze"}:
        frozen = args.run.resolve() / "code/workflow.py"
        if frozen.is_file() and frozen != Path(__file__).resolve():
            # Reuse the exact code frozen at preparation, even after active code
            # evolves. Paths are resolved before entering the portable snapshot.
            command = [sys.executable, "-B", str(frozen), args.command, str(args.run.resolve())]
            if hasattr(args, "machine"):
                command += ["--machine", str(args.machine.resolve())]
            if hasattr(args, "workers"):
                command += ["--workers", str(args.workers)]
            if hasattr(args, "index"):
                command += ["--index", str(args.index)]
            return subprocess.call(command)
    from mcfost_grid.configuration import prepare_run
    if args.command == "prepare":
        run = prepare_run(args.config)
        result = {"prepared_run": str(run), "next": f"python workflow.py run-local {run} --machine config/machine.local.json"}
    elif args.command == "slurm":
        from mcfost_grid.slurm import write_slurm
        result = write_slurm(args.run, args.machine)
    elif args.command == "analyze":
        from mcfost_grid.analysis import analyze_run
        result = analyze_run(args.run, workers=args.workers)
    else:
        from mcfost_grid.runner import run_local, run_model, status
        if args.command == "run-local":
            result = run_local(args.run.resolve(), args.machine.resolve(), args.workers)
        elif args.command == "task":
            result = run_model(args.run.resolve(), args.index, args.machine.resolve())
        else:
            result = status(args.run)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 1 if isinstance(result, dict) and result.get("failed", 0) else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
