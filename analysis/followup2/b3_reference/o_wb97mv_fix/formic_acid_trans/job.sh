#!/bin/bash
#SBATCH --job-name=f2b3_o_wb97mv_fix_formic_acid_trans
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/formic_acid_trans/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/formic_acid_trans/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/formic_acid_trans
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) formic_acid_trans.inp > formic_acid_trans.out
rm -f formic_acid_trans.gbw formic_acid_trans*.tmp*
date +%s > end_epoch.txt
