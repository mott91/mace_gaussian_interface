#!/bin/bash
#SBATCH --job-name=f2b7_o_wb97md3bj_fix_1-butanol
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97md3bj_fix/1-butanol/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97md3bj_fix/1-butanol/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97md3bj_fix/1-butanol
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) 1-butanol.inp > 1-butanol.out
rm -f 1-butanol.gbw 1-butanol*.tmp*
date +%s > end_epoch.txt
