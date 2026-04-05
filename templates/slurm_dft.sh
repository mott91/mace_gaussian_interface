#!/bin/bash
# SLURM job template for DFT frequency calculations via Gaussian 16.
#
# Placeholders (filled by mace_gaussian/slurm.py):
#   {molecule}       - molecule name (e.g. water, aspirin)
#   {remote_dir}     - remote working directory (~/mace_gaussian_dft/<molecule>)
#   {gjf_filename}   - Gaussian input file (e.g. water_freq_anharm.gjf)
#   {chk_filename}   - checkpoint file (e.g. water_freq_anharm.chk)
#   {fchk_filename}  - formatted checkpoint file (e.g. water_freq_anharm.fchk)
#
# Edit resource requests and module loads below for your cluster.

#SBATCH --job-name=dft_{molecule}
#SBATCH --output={remote_dir}/slurm_%j.out
#SBATCH --error={remote_dir}/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03

# Load Gaussian
module load rune/gaussian/g16-c.02
export OMP_NUM_THREADS=4

cd {remote_dir}

# Run Gaussian DFT frequency calculation
g16 {gjf_filename}

# Convert checkpoint to formatted checkpoint (HPC-02)
formchk {chk_filename} {fchk_filename}
