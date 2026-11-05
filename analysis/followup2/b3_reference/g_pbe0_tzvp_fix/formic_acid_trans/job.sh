#!/bin/bash
#SBATCH --job-name=f2b3_g_pbe0_tzvp_fix_formic_acid_trans
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/formic_acid_trans/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/formic_acid_trans/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/formic_acid_trans
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 formic_acid_trans.gjf
formchk formic_acid_trans.chk formic_acid_trans.fchk
rm -f formic_acid_trans.chk
date +%s > end_epoch.txt
