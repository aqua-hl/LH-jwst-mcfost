#!/usr/bin/env python3
"""Dispatch the immutable Stage-0 bootstrap and its measured-opacity grid."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def dispatch(bundle, *, stage0=False, index=None, analyze=False, status=False,
             validate=False, machine=None):
    from build_silicate_search_v3 import validate_bootstrap
    bundle = Path(bundle).resolve()
    plan = validate_bootstrap(bundle)
    production = bundle / "production"
    if validate:
        return dict(bootstrap_valid=True, models=plan["model_count"],
                    production_materialized=production.is_dir(), mcfost_invoked=False)
    if stage0:
        if machine is None:
            raise ValueError("--stage0 requires --machine")
        from silicate_search_v3_stage0 import run_stage0
        from build_silicate_search_v3 import materialize
        receipt = run_stage0(bundle, Path(machine).resolve())
        prepared = materialize(bundle)  # Only reached after all four species pass.
        return dict(stage0=receipt, production=prepared)
    if status and not production.is_dir():
        path = bundle / "stage0/receipt.json"
        receipt = json.loads(path.read_text()) if path.is_file() else {"status": "not_started"}
        return dict(production_materialized=False, stage0=receipt,
                    note="File-based status; no live scheduler query")
    if not production.is_dir():
        raise ValueError("Production is absent; Stage 0 must pass and freeze the catalogue first")
    # Bootstrap validation imports its frozen workflow. Production's strict
    # loader intentionally rejects that module location, so dispatch children
    # in fresh processes using this same venv executable (do not resolve it).
    python = os.path.abspath(sys.executable)
    if index is not None:
        if type(index) is not int or not 0 <= index < 26:
            raise ValueError("Model index must be between 0 and 25")
        if machine is None:
            raise ValueError("--index requires --machine")
        command = [python, "-B", str(production / "code/production_task.py"),
                   str(production), "--index", str(index), "--machine", str(Path(machine).resolve())]
    elif status:
        command = [python, "-B", str(production / "code/production_task.py"), str(production), "--status"]
    elif analyze:
        # The analyzer validates all frozen inputs and compact result receipts.
        command = [python, "-B", str(production / "code/analyze_silicate_search_v3.py"),
                   str(production), "--output", str(bundle / "results")]
    else:
        raise ValueError("Choose one action")
    subprocess.run(command, check=True)  # Preserve the worker's stdout/stderr and failure code.
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--stage0", action="store_true")
    modes.add_argument("--index", type=int)
    modes.add_argument("--analyze", action="store_true")
    modes.add_argument("--status", action="store_true")
    modes.add_argument("--validate", action="store_true")
    parser.add_argument("--machine", type=Path)
    args = parser.parse_args(argv)
    result = dispatch(**vars(args))
    if result is not None:
        print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode if exc.returncode > 0 else 1)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
