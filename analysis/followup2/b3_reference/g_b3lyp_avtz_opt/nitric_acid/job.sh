#!/bin/bash
#SBATCH --job-name=f2b3_g_b3lyp_avtz_opt_nitric_acid
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/nitric_acid/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/nitric_acid/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/g_b3lyp_avtz_opt/nitric_acid
date +%s > start_epoch.txt
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=8
g16 nitric_acid.gjf
formchk nitric_acid.chk nitric_acid.fchk
rm -f nitric_acid.chk
date +%s > end_epoch.txt
