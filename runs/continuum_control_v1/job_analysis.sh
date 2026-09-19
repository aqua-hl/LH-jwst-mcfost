#!/bin/bash
#SBATCH --job-name=continuum_control_v1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=logs/analysis-%j.out
#SBATCH --error=logs/analysis-%j.err

set -euo pipefail
: "${MCFOST_RUN_DIR:?Use submit.sh, or export absolute MCFOST_RUN_DIR}"
: "${MCFOST_MACHINE:?Set MCFOST_MACHINE to the cluster machine JSON}"
cd "$MCFOST_RUN_DIR"
: "${MCFOST_UTILS:?Export MCFOST_UTILS before submitting}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_DYNAMIC=FALSE MCFOST_AUTO_UPDATE=0
export MPLCONFIGDIR="$MCFOST_RUN_DIR/.mplconfig"
mkdir -p "$MPLCONFIGDIR"
exec python -B code/workflow.py analyze "$MCFOST_RUN_DIR" --workers "$SLURM_CPUS_PER_TASK"
