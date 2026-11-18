#!/bin/bash
#SBATCH --job-name=f2b7_g_b3lyp_avtz_fix_sulfur_dioxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_avtz_fix/sulfur_dioxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_avtz_fix/sulfur_dioxide/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/g_b3lyp_avtz_fix/sulfur_dioxide
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 sulfur_dioxide.gjf
formchk sulfur_dioxide.chk sulfur_dioxide.fchk
rm -f sulfur_dioxide.chk
date +%s > end_epoch.txt
