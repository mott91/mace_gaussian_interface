#!/bin/bash
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune04
#SBATCH -J BF_tBu
#SBATCH -N 1                
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1

source /usr/local/_tci_software_environment_.sh
module purge

module load rune/gaussian/g16

g16 BFPtPZ_tBu_B3LYP.com
g16 BFPtPZ_tBu_MP2.com



