# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MACE-Gaussian Interface bridges Machine Learning potentials (MACE) with Gaussian 16 quantum chemistry software for enhanced molecular IR spectroscopy calculations. The system uses ZMQ to inject ML-calculated dipole derivatives into Gaussian's anharmonic frequency calculations in real-time, enabling fast, accurate spectral predictions.

## Development Commands

### Environment Setup
```bash
# Install dependencies (uses uv package manager)
uv sync

# Activate virtual environment
source .venv/bin/activate

# Check installation and available calculators
uv run python gm_main.py --diagnose
python cli.py diagnose
```

### Running Calculations
```bash
# Single molecule calculation (legacy interface)
uv run python gm_main.py molecule.xyz

# CLI workflow (recommended)
python cli.py run water.xyz
python cli.py run water.xyz --skip-dft-baseline
python cli.py run water.xyz --energy-calculators mace_mp --dipole-calculators espaloma

# List results
python cli.py list
python cli.py list water

# Run comprehensive analysis workflow
python run_analysis.py water  # or acoh
```

### Code Quality
```bash
# Format code (line length 100)
black .
isort .

# Type checking
mypy gm_main.py

# Linting
flake8 .
```

## Architecture

### Core Workflow Pipeline

```
Input (xyz)
  ↓
[Geometry Optimization] → MACE (mace_omol/mace_mp/mace_off)
  ↓
[Frequency Calculation] → Gaussian + ML dipoles (via ZMQ)
  ↓
[Anharmonic Calculation] → Gaussian (overtones + combination bands)
  ↓
[Parse Results] → JSON (gaussian_parser.py)
  ↓
[Statistical Analysis] → Regression, KDE broadening (analyze_spectra.py)
  ↓
[HTML Report] → Publication-ready visualizations (html_report_generator.py)
```

### Module System Architecture

#### 1. **Modular Dipole Calculator System** (`gm_main.py`)

Abstract base class `DipoleCalculatorBase` with auto-fallback:
- `MACEMLDipoleCalculator` - Custom MACE fork with dipole prediction (primary)
- `EspalomaCalculator` - ML charge-based dipole (reliable fallback)
- `XTBCalculator` - Semi-empirical xTB method
- `GeometryBasedCalculator` - Geometric estimation (last resort)

**Critical**: The system uses module monkey-patching to switch between standard MACE and dipole-enabled MACE. See `mace_calculators.py` for `fake_module_from_real()` mechanism. Always call `cleanup_mace_modules()` after dipole calculations.

#### 2. **Dual MACE Package Architecture**

Two separate MACE installations coexist:
- `mace_ML_pkg/mace/` - Standard MACE-torch (energy/forces)
- `mace_dipole_pkg/mace_dipole_core/` - Custom fork with dipole moments

Import mechanism in `mace_calculators.py`:
- `load_standard_mace_calculator()` - Swaps in standard MACE
- `load_dipole_mace_calculator()` - Swaps in dipole-enabled MACE
- Must invalidate caches between switches

#### 3. **ZMQ Bridge for Gaussian** (`gm_helper.py`)

IPC file-based communication:
1. Gaussian launches external program (gm_helper.py) via `External="gm_helper.py"`
2. Helper script connects to ZMQ socket bound to `.ipc_file`
3. Main script receives request: `"infile|outfile"`
4. Python calculates dipoles, writes to outfile
5. Helper receives "done" message, Gaussian continues

**Path Requirements**: Helper script must be executable and in PATH or specified via `MACE_HELPER_SCRIPT_PATH` environment variable.

#### 4. **Results Management** (`results_manager.py`)

Hierarchical directory structure:
```
comparison_results/
  └── {molecule_name}/
      ├── geometry_opt/
      │   ├── results.json
      │   ├── {molecule}_initial.xyz
      │   └── {molecule}_optimized.xyz
      └── {energy_calc}_{dipole_calc}/
          ├── results.json
          ├── {molecule}_freq_anharm.gjf
          └── {molecule}_freq_anharm.log
```

**Directory Naming Convention**:
- DFT: `wb97xd` (not `wb97xd_wb97xd`)
- ML: `mace_mp_espaloma`, `mace_omol_mace_ml`

#### 5. **Comparison & Analysis Framework**

**Entry Point**: `comparison_workflow.py` orchestrates full analysis
- `find_dft_baseline()` - Auto-detects wb97xd DFT reference
- `load_calculation_results()` - Parses all JSON results
- `generate_comprehensive_report()` - Creates HTML with embedded plots

**Spectral Analysis** (`analyze_spectra.py`):
- Gaussian KDE broadening (FWHM configurable, default 8.0 cm⁻¹)
- Statistical metrics: MAE, RMSE, R², slope/intercept
- Supports fundamentals, overtones, combination bands
- Regression plots: ML frequencies vs DFT baseline

**Key Methods**:
- `SpectrumAnalyzer.compare_spectra()` - Returns `ComparisonMetrics`
- `SpectrumAnalyzer.plot_combined_regression()` - Multi-method comparison plots
- `HTMLReportGenerator.generate()` - Combines plots + data into report

### Configuration System

**Main Configuration** (`gm_main.py` lines 62-78):
```python
DEFAULT_MACE_DIPOLE_MODEL  # Set via MACE_DIPOLE_MODEL_PATH env var
DEFAULT_HELPER_SCRIPT      # Set via MACE_HELPER_SCRIPT_PATH env var
DIPOLE_METHOD              # 'auto', 'mace_ml', 'espaloma', 'xtb', 'geometry'
CALCULATE_DIPOLE_DERIVATIVES = True
```

**CLI Configuration** (`cli.py`):
- `--optimization-calculator`: mace_omol, mace_off, mace_mp (default: mace_omol)
- `--energy-calculators`: Comma-separated list (default: mace_mp,mace_omol)
- `--dipole-calculators`: Comma-separated list (default: espaloma,mace_ml)
- `--skip-dft-baseline`: Skip DFT anharmonic baseline
- `--force-optimization`: Re-optimize even if exists

### Key Data Structures

**Gaussian Parser Output** (`gaussian_parser.py`):
```python
{
  "calculator_type": "dft" | "ml",
  "energy_calculator": str,
  "dipole_calculator": str,
  "energy_eV": float,
  "frequencies": {
    "harmonic": [...],
    "anharmonic": [...],
    "overtones": [...],
    "combination_bands": [...]
  },
  "ir_intensities": {
    "harmonic": [...],
    "anharmonic": [...]
  },
  "runtime_s": float
}
```

**SpectrumData** (`analyze_spectra.py`):
```python
@dataclass
class SpectrumData:
    frequencies: np.ndarray
    intensities: np.ndarray
    labels: List[str]  # ['fundamental', 'overtone', 'combination']
```

## Critical Implementation Details

### Frequency Parsing Issues (Recent Fixes)

When adding parsers for overtones/combination bands in `gaussian_parser.py`:
1. Overtones: Match `"2(X)"` pattern where X is mode number
2. Combination bands: Match `"X(1) + Y(1)"` pattern
3. Must handle both harmonic and anharmonic sections separately
4. wb97xd sometimes missing anharmonic data - handle gracefully

### MACE Model Loading

Standard MACE models (MACE-OFF, MACE-MP, MACE-OMOL-0):
```python
calc = mace_off(model="small", device="cuda", default_dtype="float64")
```

Custom dipole model:
```python
calc = load_dipole_mace_calculator(
    model_paths="/path/to/model_1.model",
    device="cuda",
    model_type="DipolePolarizabilityMACE",
    default_dtype="float64"
)
```

### Gaussian Input Generation

Template in `gm_main.py` ~line 800:
- Route: `# wb97xd/6-31G** Freq=(Anharmonic) External="gm_helper.py"`
- Must include `%Chk=molecule.chk`
- Multiplicity: Auto-detect from atoms (singlet for even electrons, doublet for odd)
- Charge: Defaults to 0 (modify if needed)

### Plot Generation Best Practices

- Use PNG format for HTML embedding (not PDF)
- DPI=300 for publication quality
- Color palette: seaborn "colorblind" palette
- Font: Arial/Helvetica, size 10pt body, 11pt labels
- Always include R² and RMSE in regression plots
- Grid frequency range: 400-4000 cm⁻¹ (typical IR range)

## Known Issues & Workarounds

1. **Acetic Acid Parsing Bug** (commit a4384c4): Frequency matching for regression plots failing with acoh DFT calculations. Check overtone/combination band parsing logic.

2. **CUDA Device Placement**: Some calculators default to CPU. Always explicitly pass `device="cuda"` when GPU available.

3. **ZMQ Socket Cleanup**: If IPC file exists from previous crash, must manually delete `.ipc_file` before new run.

4. **Path Sensitivity**: Gaussian 16 requires absolute paths for external scripts. Never use relative paths in `.gjf` files.

5. **Module Contamination**: After dipole calculations, standard MACE imports may fail. Solution: Always call `cleanup_mace_modules()` in finally blocks.

## Development Workflow

### Adding New Dipole Calculator

1. Subclass `DipoleCalculatorBase` in `gm_main.py`
2. Implement `_check_availability()` and `calculate_dipole()`
3. Add to `get_dipole_calculator()` function
4. Return format: `(dipole_vector_3d, partial_charges_or_None)`

### Adding New Analysis Metric

1. Add field to `ComparisonMetrics` dataclass (`analyze_spectra.py`)
2. Calculate in `SpectrumAnalyzer.compare_spectra()`
3. Update `HTMLReportGenerator` to display in report
4. Add to CSV export in `comparison_workflow.py`

### Testing Calculations

Use provided test molecules:
- `water.xyz` - Simple 3-atom molecule (fast testing)
- `acoh.xyz` - Acetic acid (complex overtones/combinations)

Always run diagnostics first:
```bash
python gm_main.py --diagnose
```

## External Dependencies

- **Gaussian 16**: Must be in PATH as `g16` or loaded via HPC module system
- **CUDA**: Required for GPU acceleration (CPU fallback available but slow)
- **Pre-trained Models**: Dipole model at `~/mace_gaussian/dipole_model/model_1.model` (17 MB)
- **MACE-OMOL-0**: Requires mace-torch >= 0.3.14

## File Naming Conventions

- Molecules: `{name}.xyz` (lowercase)
- Gaussian input: `{name}_freq_anharm.gjf`
- Gaussian output: `{name}_freq_anharm.log`
- Optimized geometry: `{name}_optimized.xyz`
- Analysis plots: `{molecule}_{plot_type}.png`
- Reports: `{molecule}_spectral_analysis.html`
