# MACE-Gaussian Interface

Interface between MACE machine learning potentials and Gaussian quantum chemistry software for enhanced molecular calculations with ML-accelerated dipole moments and IR intensities.

## Features

- Multiple dipole calculation methods (MACE ML, Espaloma, xTB, geometry-based)
- Automatic geometry optimization using MACE potentials
- Real-time communication with Gaussian via ZMQ
- Anharmonic frequency calculations with ML-enhanced dipole derivatives

## Installation

### Prerequisites

- Python 3.9+
- Gaussian 16
- CUDA-capable GPU (recommended)

### Quick Start

1. **Install uv**:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. **Clone and install**:
```bash
git clone https://github.com/yourusername/mace-gaussian-interface.git
cd mace-gaussian-interface
uv sync
```

3. **Test installation**:
```bash
uv run python gm_main.py --diagnose
```

## Usage

### Run calculations with uv (recommended):
```bash
uv run python gm_main.py molecule.xyz
```

### Or activate environment:
```bash
source .venv/bin/activate
python gm_main.py molecule.xyz
```

## Setup Gaussian Integration

1. **Make helper script executable**:
```bash
chmod +x gm_helper.py
```

2. **Copy to your bin directory**:
```bash
cp gm_helper.py /home/bin/
```

3. **Use in Gaussian input files**:
```
# freq (anharm)
# external="/home/bin/gm_helper.py"
```

## Workflow

1. **Check environment**:
```bash
uv run python gm_main.py --diagnose
```

2. **Run calculation**:
```bash
uv run python gm_main.py your_molecule.xyz
```

This will:
- Optimize geometry using MACE
- Generate Gaussian input file (`molecule_freq_anharm.gjf`)
- Launch Gaussian with ML-enhanced dipole calculations

## Configuration

Edit `gm_main.py` to change settings:

```python
DIPOLE_METHOD = 'auto'  # Options: 'auto', 'mace_ml', 'espaloma', 'xtb', 'geometry'
CALCULATE_DIPOLE_DERIVATIVES = True
```

## Example

```bash
# Create a simple water molecule file (water.xyz)
echo "3
Water molecule
O     0.000000     0.000000     0.117000
H     0.000000     0.757000    -0.467000  
H     0.000000    -0.757000    -0.467000" > water.xyz

# Run calculation
uv run python gm_main.py water.xyz
```

## Troubleshooting

**Missing dependencies**:
```bash
uv sync  # Reinstall dependencies
```

**CUDA not available**:
```bash
nvidia-smi  # Check GPU status
```

**Gaussian not found**:
```bash
which g16  # Check if in PATH
module load gaussian  # Load if using HPC
```

**ZMQ communication errors**:
```bash
# Ensure helper script is executable and in PATH
ls -la /home/bin/gm_helper.py
```

## Dependencies

Core packages installed automatically:
- numpy, ase, torch, mace-torch, pyzmq
- espaloma_charge, rdkit (for dipole calculations)

## Quick Reference

```bash
# Install
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone <repo-url> && cd mace-gaussian-interface
uv sync

# Use
uv run python gm_main.py --diagnose  # Test
uv run python gm_main.py mol.xyz     # Calculate
```