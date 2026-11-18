#!/bin/bash
#SBATCH --job-name=f2b7_g_b3lyp_631gdp_harm_1-butanol
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_631gdp_harm/1-butanol/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_631gdp_harm/1-butanol/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_631gdp_harm/1-butanol
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 1-butanol.gjf
formchk 1-butanol.chk 1-butanol.fchk
rm -f 1-butanol.chk
date +%s > end_epoch.txt
