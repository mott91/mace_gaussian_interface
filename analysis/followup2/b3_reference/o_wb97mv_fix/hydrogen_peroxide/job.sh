#!/bin/bash
#SBATCH --job-name=f2b3_o_wb97mv_fix_hydrogen_peroxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/hydrogen_peroxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/hydrogen_peroxide/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/hydrogen_peroxide
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) hydrogen_peroxide.inp > hydrogen_peroxide.out
rm -f hydrogen_peroxide.gbw hydrogen_peroxide*.tmp*
date +%s > end_epoch.txt
