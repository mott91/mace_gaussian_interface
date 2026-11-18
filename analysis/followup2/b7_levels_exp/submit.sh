#!/bin/bash
# Follow-up 2, B7: copy the inputs to rune03 scratch (through tci5) and submit the 21 calculations
# plus the analysis job that waits for all of them. Appends "<level>/<mol> <jobid>" to submitted.txt.
#   bash analysis/followup2/b7_levels_exp/submit.sh
# Fetch results later:  bash analysis/followup2/b7_levels_exp/fetch.sh
set -euo pipefail
HOST=mot@tci5
HERE=$(cd "$(dirname "$0")" && pwd)
REMOTE=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7
SSH="ssh -n -o BatchMode=yes -o ConnectTimeout=30"   # -n: keep ssh off the while-loop's stdin

$SSH $HOST "mkdir -p $REMOTE/results"
rsync -a --exclude results --exclude __pycache__ --exclude submitted.txt "$HERE/" "$HOST:$REMOTE/"
echo "submitted $(date -Iseconds)" >> "$HERE/submitted.txt"
ids=()
while read -r j; do
  [ -z "$j" ] && continue
  id=$($SSH $HOST "cd $REMOTE/$j && sbatch --parsable job.sh")
  echo "$j $id" | tee -a "$HERE/submitted.txt"
  ids+=("$id")
done < "$HERE/jobs.txt"
dep=$(IFS=:; echo "${ids[*]}")
aid=$($SSH $HOST "cd $REMOTE && sbatch --parsable --dependency=afterok:$dep analysis_job.sh")
echo "analysis $aid (afterok:$dep)" | tee -a "$HERE/submitted.txt"
