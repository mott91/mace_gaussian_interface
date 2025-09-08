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
- Python 3.9+
- Gaussian 16 (for quantum chemistry calculations)
- CUDA-capable GPU (recommended for MACE calculations)

### Quick Install (Recommended)

1. **Install uv** (if not already installed):
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
# or: pip install uv
