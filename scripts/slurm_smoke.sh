#!/bin/bash
#SBATCH --job-name=tmc1a-array-test
#SBATCH --array=0-3%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128M
#SBATCH --time=00:02:00
#SBATCH --output=slurm-smoke-%A_%a.out

# Submit from the directory where you want the four output logs:
#   sbatch scripts/slurm_smoke.sh
# Add --partition=NAME and --account=NAME before the script path if required.
set -euo pipefail

: "${SLURM_ARRAY_JOB_ID:?Submit this script with sbatch}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"
: "${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}"
cd "$SLURM_SUBMIT_DIR"

printf 'Array job: %s | task: %s | host: %s\n' \
    "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$(hostname)"
printf 'Working directory: %s\n' "$PWD"
printf 'Started: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

# Tiny, distinct calculation for each task; no Python or MCFOST required.
value=$((SLURM_ARRAY_TASK_ID + 1))
printf 'Calculation: %s squared = %s\n' "$value" "$((value * value))"
printf 'SUCCESS: task %s completed\n' "$SLURM_ARRAY_TASK_ID"
