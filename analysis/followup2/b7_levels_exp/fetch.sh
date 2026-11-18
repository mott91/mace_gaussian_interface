#!/bin/bash
# Follow-up 2, B7: copy outputs and results back from rune03 scratch (through tci5), without the
# large Gaussian/ORCA scratch files. Then rerun the analysis locally to check it matches:
#   bash analysis/followup2/b7_levels_exp/fetch.sh
#   cd analysis/followup2/b7_levels_exp && python3 analyse.py && diff results/b7_tables.md results_cluster/b7_tables.md
set -euo pipefail
HOST=mot@tci5
HERE=$(cd "$(dirname "$0")" && pwd)
REMOTE=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7
rsync -a --exclude 'results/' --include '*/' --include '*.log' --include '*.fchk' --include '*.out' --include '*.hess' \
  --include '*.err' --include '*_epoch.txt' --exclude '*' "$HOST:$REMOTE/" "$HERE/"
rsync -a "$HOST:$REMOTE/results/" "$HERE/results_cluster/"
