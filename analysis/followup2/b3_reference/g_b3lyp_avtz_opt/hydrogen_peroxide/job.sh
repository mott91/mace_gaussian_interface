#!/bin/bash
#SBATCH --job-name=f2b3_g_b3lyp_avtz_opt_hydrogen_peroxide
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/hydrogen_peroxide/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/hydrogen_peroxide/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/hydrogen_peroxide
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 hydrogen_peroxide.gjf
formchk hydrogen_peroxide.chk hydrogen_peroxide.fchk
rm -f hydrogen_peroxide.chk
date +%s > end_epoch.txt
