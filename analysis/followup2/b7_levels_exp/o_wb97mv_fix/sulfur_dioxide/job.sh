#!/bin/bash
#SBATCH --job-name=f2b7_o_wb97mv_fix_sulfur_dioxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/sulfur_dioxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/sulfur_dioxide/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/o_wb97mv_fix/sulfur_dioxide
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) sulfur_dioxide.inp > sulfur_dioxide.out
rm -f sulfur_dioxide.gbw sulfur_dioxide*.tmp*
date +%s > end_epoch.txt
