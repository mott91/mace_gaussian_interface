# MACE-Gaussian Interface

Interface between MACE machine learning potentials and Gaussian quantum chemistry software for enhanced molecular calculations with ML-accelerated dipole moments and IR intensities.

## Features

- Multiple dipole calculation methods (MACE ML, Espaloma, xTB, geometry-based)
- Automatic geometry optimization
- ZMQ-based communication with Gaussian
- Anharmonic frequency calculations
- Diagnostic tools for environment checking

## Installation

### Prerequisites

- Python 3.8+
- Gaussian 16 (for quantum chemistry calculations)
- CUDA-capable GPU (recommended for MACE calculations)

### Option 1: Conda (Recommended)

```bash
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
conda env create -f environment.yml
conda activate mace-gaussian-interface
