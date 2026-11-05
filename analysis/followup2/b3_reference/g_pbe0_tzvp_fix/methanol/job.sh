#!/bin/bash
#SBATCH --job-name=f2b3_g_pbe0_tzvp_fix_methanol
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/methanol/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/methanol/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_pbe0_tzvp_fix/methanol
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 methanol.gjf
formchk methanol.chk methanol.fchk
rm -f methanol.chk
date +%s > end_epoch.txt
