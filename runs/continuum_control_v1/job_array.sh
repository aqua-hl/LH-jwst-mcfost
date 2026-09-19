#!/bin/bash
#SBATCH --job-name=continuum_control_v1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=0-0%1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/array-%A_%a.out
#SBATCH --error=logs/array-%A_%a.err

set -euo pipefail
: "${MCFOST_RUN_DIR:?Use submit.sh, or export absolute MCFOST_RUN_DIR}"
: "${MCFOST_MACHINE:?Set MCFOST_MACHINE to the cluster machine JSON}"
cd "$MCFOST_RUN_DIR"
: "${MCFOST_UTILS:?Export MCFOST_UTILS before submitting}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE MCFOST_AUTO_UPDATE=0
export MPLCONFIGDIR="$MCFOST_RUN_DIR/.mplconfig"
mkdir -p "$MPLCONFIGDIR"
exec python -B code/workflow.py task "$MCFOST_RUN_DIR" --index "$SLURM_ARRAY_TASK_ID" --machine "$MCFOST_MACHINE"
