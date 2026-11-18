#!/bin/bash
#SBATCH --job-name=f2b7_o_wb97mv_fix_ethanol
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/ethanol/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/ethanol/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/ethanol
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) ethanol.inp > ethanol.out
rm -f ethanol.gbw ethanol*.tmp*
date +%s > end_epoch.txt
