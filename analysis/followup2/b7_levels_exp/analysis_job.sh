#!/bin/bash
#SBATCH --job-name=f2b7_analysis
#SBATCH --output=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/results/slurm_%j.out
#SBATCH --error=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7/results/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
# Follow-up 2, B7: tables and figure once all 21 calculations have finished (submit.sh adds
# --dependency=afterok:<ids>). Each calculation's job.sh ends with `date`, so a failed g16/ORCA run
# does not block this job; analyse.py lists missing outputs and abnormal terminations instead.
R=/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7
cd $R
module load rune/miniconda/tci-collection-v1
python analyse.py $R /scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3 > results/analyse_stdout.txt 2>&1
echo "exit $?" >> results/analyse_stdout.txt
