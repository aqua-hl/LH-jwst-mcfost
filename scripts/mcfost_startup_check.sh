#!/bin/bash
#SBATCH --job-name=mcfost-startup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:02:00
#SBATCH --output=mcfost-startup-%j.out

# Compare the same executable on the login node and a compute node:
#   bash scripts/mcfost_startup_check.sh > mcfost-login.out 2>&1
#   sbatch scripts/mcfost_startup_check.sh
# Use the same partition/account/module environment as the failing smoke job.
# Optional: export MCFOST_EXE=/absolute/path/to/a/specific/mcfost
set -euo pipefail
ulimit -c 0
export LC_ALL=C MCFOST_AUTO_UPDATE=0 OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

printf 'Host: %s\nJob: %s\nArchitecture: ' "$(hostname)" "${SLURM_JOB_ID:-login-node}"
uname -m
printf 'OpenMP threads: %s\n' "$OMP_NUM_THREADS"
if command -v lscpu >/dev/null 2>&1; then
    lscpu
fi

binary=$(command -v "${MCFOST_EXE:-mcfost}") || {
    printf 'ERROR: MCFOST executable is not on PATH; activate the same environment as the smoke test.\n' >&2
    exit 127
}
printf '\nExecutable: %s\n' "$binary"
if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$binary"
elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$binary"
fi
if command -v file >/dev/null 2>&1; then
    file -L "$binary"
fi
if command -v ldd >/dev/null 2>&1; then
    ldd "$binary" || true
fi

printf '\nMCFOST_UTILS: %s\n' "${MCFOST_UTILS:-<unset>}"
printf 'MY_MCFOST_UTILS: %s\n' "${MY_MCFOST_UTILS:-<unset>}"
for directory in Dust Lambda Stellar_Spectra; do
    if [[ -n "${MCFOST_UTILS:-}" && -d "$MCFOST_UTILS/$directory" ]]; then
        printf 'Utilities directory accessible: %s\n' "$directory"
    else
        printf 'Utilities directory missing/inaccessible: %s\n' "$directory"
    fi
done

# Bound the startup probe on the login node as well as under Slurm.
timeout_command=$(command -v timeout || command -v gtimeout) || {
    printf 'ERROR: This diagnostic requires the coreutils timeout command.\n' >&2
    exit 127
}
help_log=$(mktemp)
trap 'rm -f -- "$help_log"' EXIT
printf '\nRunning mcfost -help with a 30-second limit and automatic updates disabled.\n'
startup_exit=0
"$timeout_command" --kill-after=5s 30s "$binary" -help > "$help_log" 2>&1 || startup_exit=$?
cat "$help_log"
printf '\nStartup exit code: %s\n' "$startup_exit"
if [[ "$startup_exit" -eq 132 ]]; then
    printf 'SIGILL reproduced. Compare CPU flags and executable checksums between hosts.\n'
elif [[ "$startup_exit" -eq 124 ]]; then
    printf 'Startup exceeded the 30-second limit.\n'
fi
if [[ "$startup_exit" -ne 0 ]]; then
    exit "$startup_exit"
fi
if ! grep -qi 'usage : mcfost' "$help_log"; then
    printf 'ERROR: The help output lacks the usage marker required by this workflow.\n' >&2
    exit 1
fi
if grep -q 'ERROR:' "$help_log"; then
    printf 'ERROR: MCFOST reported a startup error; inspect its output above.\n' >&2
    exit 1
fi
printf 'SUCCESS: MCFOST startup passed on this host; model execution still needs its smoke test.\n'
