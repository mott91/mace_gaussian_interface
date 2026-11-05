#!/bin/bash
#SBATCH --job-name=f2b3_g_pbe0_tzvp_fix_carbon_dioxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/carbon_dioxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/carbon_dioxide/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/carbon_dioxide
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 carbon_dioxide.gjf
formchk carbon_dioxide.chk carbon_dioxide.fchk
rm -f carbon_dioxide.chk
date +%s > end_epoch.txt
