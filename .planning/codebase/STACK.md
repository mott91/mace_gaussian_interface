# Technology Stack

**Analysis Date:** 2026-02-16

## Languages

**Primary:**
- Python 3.9+ - Core application language for ML calculations, quantum chemistry interface, and analysis workflows

**Secondary:**
- Shell (bash/zsh) - CLI entry points and environment setup scripts

## Runtime

**Environment:**
- Python 3.8-3.12 supported (3.10 preferred based on environment.yml)
- Virtual environment via `uv` package manager

**Package Manager:**
- `uv` - Primary dependency manager
- `pip` - Secondary, used in environment.yml and within conda
- Lockfile: `uv.lock` present

## Frameworks

**Core Scientific Computing:**
- PyTorch 2.4.0 - Deep learning framework, used for MACE models and dipole calculations
- ASE (Atomic Simulation Environment) 3.26.0 - Molecular structure manipulation and calculator interface
- NumPy 1.26.4+ - Numerical operations, array handling
- SciPy 1.15.3 - Scientific computations (used indirectly via torch/ase)

**ML Models:**
- MACE (Moment Tensor Potential) - Custom local packages
  - `mace_torch` (custom, in `mace_ML_pkg/`) - Energy calculator
  - `mace_dipole_core` (custom, in `mace_dipole_pkg/`) - Dipole moment and polarizability calculator
- RDKit 2022.03.0+ - Molecular structure representation, SMILES parsing
- e3nn 0.4.4 - Equivariant neural networks (used by MACE)
- CUE Equivariance 0.6.0 - Equivariance operations

**Charge Calculations:**
- espaloma_charge 0.0.8+ - Charge prediction for dipole moment calculations

**Visualization & Analysis:**
- Matplotlib 3.10.6 - Static plotting, spectrum visualization
- Seaborn - Color palette management (colorblind palette)
- Pandas 2.3.2+ - Data frame operations for analysis results

**Inter-Process Communication:**
- PyZMQ (ZeroMQ) 27.1.0 - Socket-based communication between Gaussian and Python processes via IPC transport

**Testing:**
- pytest 7.0.0+ - Test framework
- (ty 0.0.1a1 - Type checking, referenced in CLAUDE.md)

**Build/Dev:**
- Hatchling - Build backend
- Ruff 0.9.0+ - Linting and code formatting (line length 100)
- GitPython 3.1.45 - Git operations tracking

**CLI & Configuration:**
- Click - Command-line interface framework (via cli.py)
- ConfigArgParse 1.7.1 - Configuration file parsing
- PyYAML 6.0.2+ - YAML configuration support

**Utilities:**
- tqdm 4.67.1 - Progress bars
- prettytable 3.16.0 - Formatted table output
- h5py 3.14.0+ - HDF5 file I/O (optional for large datasets)
- matscipy 1.1.1 - Materials science utilities
- python_hostlist 2.3.0 - HPC host parsing
- opt_einsum 3.4.0, opt-einsum-fx 0.1.4 - Einsum optimization

**GPU Acceleration:**
- CUDA Toolkit 11.8 (cu118) - GPU compute capability
  - nvidia-cuda-runtime-cu11 11.8.89
  - nvidia-cudnn-cu11 9.1.0.70
  - nvidia-cublas-cu11 11.11.3.6
  - nvidia-cusolver-cu11 11.4.1.48
- Triton 3.0.0 - Kernel compiler for GPU operations
- TorchAudio 2.4.0, TorchVision 0.19.0 - Additional PyTorch domains
- TorchMetrics 1.8.1 - Metric computation on GPU
- Torch-EMA 0.3 - Exponential moving average for models

## Key Dependencies

**Critical:**
- `torch>=1.12.0` - Neural network backbone, must support CUDA for acceptable performance
- `mace-torch>=0.3.14` - MACE-OMOL-0 support requirement, energy calculations on all molecules
- `pyzmq>=22.0.0` - Real-time communication with Gaussian external hook, essential for frequency calculations
- `ase>=3.22.0` - Atomic structure interface, geometry optimization
- `espaloma_charge>=0.0.8` - Dipole calculation baseline method

**Infrastructure:**
- `rdkit>=2022.03.0` - Molecular descriptor generation, structure parsing
- `numpy>=1.20.0` - All numerical operations throughout stack
- `matplotlib>=3.5.0` - Spectrum and regression plotting
- `pandas>=2.3.2` - Results aggregation and CSV export

## Configuration

**Environment:**
- Conda environment file: `environment.yml` (mace4ir_v2)
- Python version management via `.python-version` or pyenv
- GPU device configuration via `CUDA_VISIBLE_DEVICES` environment variable
- Model paths via environment variables:
  - `MACE_DIPOLE_MODEL_PATH` - Path to dipole model file (default: `~/mace_gaussian/dipole_model/model_1.model`)
  - `MACE_HELPER_SCRIPT_PATH` - Path to ZMQ helper script (default: `./gm_helper.py`)

**Build:**
- `pyproject.toml` - Project metadata and dependency declarations
- `[tool.ruff]` - Linting configuration (line-length: 100)
- `[tool.ruff.lint.isort]` - Import organization (known-first-party modules defined)

## Platform Requirements

**Development:**
- Linux (x86_64 verified)
- CUDA 11.8 capable GPU (NVIDIA)
- 8GB+ GPU VRAM recommended for MACE calculations
- ~4GB RAM minimum for DFT baseline

**Production:**
- Gaussian 16 must be installed and in PATH as `g16` command
- formchk utility (Gaussian checkpoint converter) required in PATH
- CUDA-capable GPU for acceptable performance (CPU fallback available but ~10x slower)
- Linux HPC environment (module system support via ConfigArgParse)

**Python-only dependencies:**
- All dependencies in `pyproject.toml` installable via `uv sync`
- MACE packages installed separately via custom installation script (see environment.yml notes)

---

*Stack analysis: 2026-02-16*
