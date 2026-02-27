# MACE-Gaussian Interface

Interface between MACE machine learning potentials and Gaussian 16 for enhanced IR spectroscopy calculations. Combines ML-accelerated dipole calculations with quantum chemistry for fast, accurate anharmonic frequency predictions.

## Features

- **ML-Enhanced Calculations**: Multiple dipole calculators (MACE, Espaloma, xTB, geometric)
- **Anharmonic Frequencies**: Full anharmonic treatment including overtones and combination bands
- **DFT Baselines**: Pure DFT calculations for rigorous comparison
- **Comprehensive Analysis**: Statistical comparison, regression plots, publication-ready HTML reports
- **Modern CLI**: Interactive command-line interface for workflow management
- **Real-time Integration**: ZMQ communication between Python and Gaussian

## Installation

### Prerequisites

- **Python 3.10** (required — matches the `mace4ir_v2` environment)
- **Gaussian 16** (must be in PATH as `g16`)
- **CUDA-capable GPU** (required for ML calculations; CPU fallback is slow)
- **micromamba or conda** (for environment management)

### Step-by-Step Installation

**1. Clone the repository**

```bash
git clone https://github.com/your-username/mace_gaussian.git
cd mace_gaussian
```

**2. Create the conda environment**

The repository includes `environment.yml` with the exact package versions used for development.

```bash
# Install micromamba if needed
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Create and activate the environment
micromamba env create -f environment.yml
micromamba activate mace4ir_v2
```

**3. Install the custom MACE packages**

This project depends on two local MACE forks included in the repository.
They must be installed in editable mode so their paths resolve correctly:

```bash
# Install standard MACE (energy/force calculator)
pip install -e mace_ML_pkg

# Install MACE with dipole support (dipole derivative calculator)
pip install -e mace_dipole_pkg
```

These packages are not on PyPI and must be installed from the local subdirectories.

**4. Install the main package**

```bash
pip install -e .
```

This installs the `mace_gaussian` package in editable mode so CLI commands and imports work.

**5. Verify the installation**

```bash
mace-gaussian diagnose
```

This checks for Gaussian 16 (`g16`), `formchk`, CUDA availability, and all dipole calculators.
A successful install shows green checkmarks for the calculators; Gaussian/CUDA checks may
show warnings if those are not available on the current machine.

### Environment Files

The repository includes:
- `environment.yml` — Complete conda/micromamba environment with exact package versions
- `requirements_mace4ir_v2.txt` — Pip package list for reference
- `install_mace_packages.sh` — Convenience script that runs steps 3 above automatically

## Quickstart

After installation, see [docs/quickstart.md](docs/quickstart.md) for a step-by-step guide
to running your first calculation and viewing the HTML analysis report.

## Usage

### Quick Example

```bash
# Run calculations on a test molecule
mace-gaussian run water.xyz

# Generate analysis report
python run_analysis.py water
```

This will:
1. Optimize geometry with MACE
2. Run ML calculations (multiple energy/dipole combos)
3. Run DFT baseline for comparison
4. Generate statistical analysis and HTML report

### Available Test Molecules

The repo includes test molecules:
- `water.xyz` - Simple (3 atoms, ~5-10 min)
- `co.xyz` - Diatomic linear molecule (2 atoms, ~2-5 min)
- `ammonia.xyz` - Small molecule (4 atoms, ~5-15 min)
- `formaldehyde.xyz` - Carbonyl test (4 atoms, ~10-20 min)
- `acoh.xyz` - Acetic acid (8 atoms, ~75 min DFT)

## Workflow Overview

```
Input XYZ → Geometry Opt (MACE) → ML Frequency Calcs → DFT Baseline → Analysis & Report
```

### Step-by-Step

#### 1. Run Calculations
```bash
# Basic run (includes DFT baseline)
mace-gaussian run molecule.xyz

# Skip DFT baseline (faster, for testing)
mace-gaussian run molecule.xyz --skip-dft-baseline

# Customize calculators
mace-gaussian run molecule.xyz \
  --energy-calculators mace_mp,mace_omol \
  --dipole-calculators espaloma,mace_ml
```

#### 2. Check Results
```bash
# List all results
mace-gaussian list

# List specific molecule
mace-gaussian list water
```

#### 3. Generate Analysis Report
```bash
# Create comprehensive analysis with plots and HTML report
python run_analysis.py water
```

The analysis generates:
- `analysis_results/water/plots/` - Spectrum and regression plots
- `analysis_results/water/data/` - CSV comparison tables
- `analysis_results/water/report.html` - **Interactive HTML report**

## Output Structure

```
comparison_results/
  └── water/
      ├── geometry_opt/
      │   ├── optimized.xyz
      │   └── results.json
      ├── wb97xd_def2tzvp/          # DFT baseline
      │   ├── gaussian_dft.log
      │   └── results.json
      └── mace_mp_espaloma/          # ML calculation
          ├── gaussian_freq.log
          └── results.json

analysis_results/
  └── water/
      ├── plots/
      │   ├── spectrum_combined.png
      │   ├── regression_combined.png
      │   └── ...
      ├── data/
      │   └── comparison_*.csv
      └── report.html               # Open this in browser!
```

## Configuration

Most settings have sensible defaults. Advanced users can edit:

**Energy Calculators** (in `cli.py` or command line):
- `mace_omol` - MACE-OMOL
- `mace_mp` - MACE-MP
- `mace_off` - MACE-OFF (not included yet)

**Dipole Calculators**:
- `espaloma` - ML charge-based
- `mace_ml` - Custom MACE dipole model
- `xtb` - Semi-empirical (not working for now)
- `geometry` - Geometric fallback

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

## For Developers

See `CLAUDE.md` for detailed architecture documentation including:
- Code structure and key design patterns
- How to add new calculators
- Parser internals
- Known issues and workarounds

## Common Issues

**"Gaussian not found"**: Make sure `g16` is in your PATH
```bash
which g16  # Should return path to Gaussian
```

**"No CUDA devices"**: ML will fallback to CPU (slower but works)

**Environment issues**: If you get import errors:
```bash
# Make sure you're in the correct environment
micromamba activate mace4ir_v2

# Verify MACE packages are installed
python -c "import mace; print('mace-torch OK')"
python -c "import mace_dipole_core; print('mace-dipole OK')"

# If MACE imports fail, reinstall:
./install_mace_packages.sh
```

**"Anharmonic intensities missing"**: This happens when:
- Molecule has imaginary frequencies (not at minimum)
- Linear molecules (Gaussian limitation)
- Parser sees `R²=N/A` in plots for limited data

**DFT baseline takes forever**: This is normal! DFT anharmonic calculations are expensive:
- Water: ~15-30 min
- Acetic acid: ~75 min  
- Consider `--skip-dft-baseline` for quick tests

