#!/bin/bash
#SBATCH --job-name=f2b3_o_wb97mv_fix_carbonic_acid
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/carbonic_acid/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/carbonic_acid/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3/o_wb97mv_fix/carbonic_acid
date +%s > start_epoch.txt
module load chaos/orca/6.1.1
$(which orca) carbonic_acid.inp > carbonic_acid.out
rm -f carbonic_acid.gbw carbonic_acid*.tmp*
date +%s > end_epoch.txt
