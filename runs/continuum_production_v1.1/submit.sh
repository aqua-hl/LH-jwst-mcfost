#!/bin/bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: bash submit.sh /absolute/path/to/machine.cluster.json" >&2; exit 2; fi
export MCFOST_RUN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "$1" in /*) export MCFOST_MACHINE="$1";; *) echo "Use an absolute machine JSON path" >&2; exit 2;; esac
test -f "$MCFOST_MACHINE"
cd "$MCFOST_RUN_DIR"
mkdir -p logs
job=$(sbatch --parsable --export=ALL job_array.sh)
job=${job%%;*}
echo "Submitted model array $job"
sbatch --parsable --export=ALL --dependency="afterany:$job" job_analysis.sh
# afterany makes the summary report failures too; incomplete models are never ranked.
