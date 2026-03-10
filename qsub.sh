#!/bin/bash

#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
#SBATCH --job-name=BSLNE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G


module purge
#module load rune/pq/0.6.0
#module load rune/gcc/13.2.0
#module load rune/dftbplus/24.1
module load g16

export OMP_NUM_THREADS=8

 g16 dft_mace_omol_20260220_113516.gjf
