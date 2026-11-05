#!/bin/bash
#SBATCH --job-name=f2b3_o_wb97md3bj_fix_carbon_dioxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97md3bj_fix/carbon_dioxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97md3bj_fix/carbon_dioxide/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97md3bj_fix/carbon_dioxide
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) carbon_dioxide.inp > carbon_dioxide.out
rm -f carbon_dioxide.gbw carbon_dioxide*.tmp*
date +%s > end_epoch.txt
