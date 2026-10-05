#!/bin/zsh
# 2026 campaign launcher. Start it detached so it survives a closed terminal or lost SSH:
#
#     setsid nohup scripts/launch_campaign_2026.sh > /dev/null 2>&1 < /dev/null &
#
# The batches run one after another: they share the one GPU, and the recorded run times
# are thesis data, so they must not compete. The big molecules go first because their
# B3LYP jobs take longest and are submitted when their batch starts. A finished molecule
# is skipped on restart, so the script can simply be started again after an interruption.
#
# Progress: tail -f campaigns/2026/logs/<list>.log ; one line per batch in summary.txt.

cd /home/mot/mace_gaussian
ENV=/home/mot/micromamba/envs/mace4ir_v3_accel
LOGS=campaigns/2026/logs
mkdir -p $LOGS

export PATH=$ENV/bin:$PATH CONDA_PREFIX=$ENV
# The login shell points CC/CXX at the NVIDIA HPC SDK compilers, which torch cannot use.
export CC=/usr/bin/gcc CXX=/usr/bin/g++
# cuEquivariance kernels for the energy models; compiled kernels are kept on disk.
export MACE_ENABLE_CUEQ=1 CUEQUIVARIANCE_OPS_NVRTC_CACHE_DIR=$HOME/.cache/cuequivariance_nvrtc

run() {  # run <list> [extra batch options]
    local list=$1; shift
    local t=$(date +%s)
    mace-gaussian batch molecules/panel_2026_$list.txt --campaign 2026 \
        --dft-on-cluster mot@tci5 "$@" > $LOGS/$list.log 2>&1
    echo "$(date '+%F %T') $list exit $? wall $(( $(date +%s) - t ))s" >> $LOGS/summary.txt
}

run big --slurm-template templates/slurm_dft_big.sh
run acids
run alcohols
run inorganic
echo "$(date '+%F %T') DONE" >> $LOGS/summary.txt
