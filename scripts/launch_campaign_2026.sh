#!/bin/zsh
# 2026 campaign launcher. Start it detached so it survives a closed terminal or lost SSH:
#
#     setsid nohup scripts/launch_campaign_2026.sh > /dev/null 2>&1 < /dev/null &
#
# Step 1 submits the B3LYP job of every molecule to the cluster at once (after the shared
# pre-optimization), so SLURM can run them whenever cores are free. Step 2 runs the ML side,
# one batch after another: they share the one GPU, and the recorded run times are thesis
# data, so they must not compete. Each batch then collects its own B3LYP results and writes
# reports and figures. A finished molecule is skipped and a submitted job is not submitted
# again, so the script can simply be started again after an interruption.
#
# Progress: tail -f campaigns/2026/logs/ml_<list>.log ; one line per step in summary.txt.

cd /home/mot/mace_gaussian
ENV=/home/mot/micromamba/envs/mace4ir_v3_accel
LOGS=campaigns/2026/logs
mkdir -p $LOGS

export PATH=$ENV/bin:$PATH CONDA_PREFIX=$ENV
# The login shell points CC/CXX at the NVIDIA HPC SDK compilers, which torch cannot use.
export CC=/usr/bin/gcc CXX=/usr/bin/g++
# cuEquivariance kernels for the energy models; compiled kernels are kept on disk.
export MACE_ENABLE_CUEQ=1 CUEQUIVARIANCE_OPS_NVRTC_CACHE_DIR=$HOME/.cache/cuequivariance_nvrtc

run() {  # run <label> <list> [extra batch options]
    local label=$1 list=$2; shift 2
    local t=$(date +%s)
    mace-gaussian batch molecules/panel_2026_$list.txt --campaign 2026 \
        --dft-on-cluster mot@tci5 "$@" > $LOGS/${label}_$list.log 2>&1
    echo "$(date '+%F %T') $label $list exit $? wall $(( $(date +%s) - t ))s" >> $LOGS/summary.txt
}

BIG=(--slurm-template templates/slurm_dft_big.sh)

# Step 1: every B3LYP job onto the cluster.
run submit big $BIG --submit-dft-only
run submit acids --submit-dft-only
run submit alcohols --submit-dft-only
run submit inorganic --submit-dft-only

# Step 2: ML runs, analysis, reports. Big last, so its long B3LYP jobs have time to finish.
run ml acids
run ml alcohols
run ml inorganic
run ml big $BIG
echo "$(date '+%F %T') DONE" >> $LOGS/summary.txt
