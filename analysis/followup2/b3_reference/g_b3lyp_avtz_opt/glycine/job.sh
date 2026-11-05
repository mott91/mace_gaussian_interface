#!/bin/bash
#SBATCH --job-name=f2b3_g_b3lyp_avtz_opt_glycine
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/glycine/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/glycine/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/glycine
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 glycine.gjf
formchk glycine.chk glycine.fchk
rm -f glycine.chk
date +%s > end_epoch.txt
