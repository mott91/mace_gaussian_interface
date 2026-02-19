# Architecture

## Modular Dipole Calculator System (`gm_main.py`)

Abstract base class `DipoleCalculatorBase` with three implementations:
- `MACEMLDipoleCalculator` - Custom MACE fork with dipole prediction (primary)
- `EspalomaCalculator` - ML charge-based dipole (reliable fallback)
- `XTBCalculator` - Semi-empirical xTB method

`DipoleCalculatorFactory` manages selection with preferred order: `["mace_ml", "espaloma", "xtb"]`. Auto-selection returns first available.

## Dual MACE Package Architecture

Two separate MACE installations coexist:
- `mace_ML_pkg/mace/` - Standard MACE-torch (energy/forces)
- `mace_dipole_pkg/mace_dipole_core/` - Custom fork with dipole moments

`calculators/mace_loader.py` handles loading the dipole model via `torch.load(pickle_module=...)` class remapping, which remaps `mace.modules.models` class paths to `mace_dipole_core.modules.models` during deserialization. No `sys.modules` mutation or cache invalidation is needed.

## ZMQ Bridge for Gaussian (`gm_helper.py`)

IPC file-based communication:
1. Gaussian launches `gm_helper.py` via `External="gm_helper.py"`
2. Helper connects to ZMQ socket bound to `.ipc_file`
3. Main script receives request: `"infile|outfile"`
4. Python calculates dipoles, writes to outfile
5. Helper receives "done" message, Gaussian continues

## Mode Matching System (`mode_matching.py` + `fchk_parser.py`)

Matches vibrational modes between DFT and ML using eigenvector dot products (scalar overlap) instead of frequency ordering. Handles mode reordering between methods.

`fchk_parser.py` parses Gaussian formatted checkpoint files and converts `.chk` → `.fchk` via `formchk`.

`mode_matching.py` computes overlap matrices and generates pastel-styled heatmap PNGs.

Integration: `ComparisonWorkflow.extract_mode_mapping()` → `SpectrumAnalyzer.match_by_mode()` for rigorous frequency pairing.

## Comparison & Analysis Framework

**Entry point**: `comparison_workflow.py` orchestrates full analysis via `ComparisonWorkflow`.

`analyze_molecule(name, use_harmonic=False)` is the one-call entry point. `analyze_molecule_harmonic(name)` is the harmonic-only shortcut.

**Spectral analysis** (`analyze_spectra.py`): Gaussian KDE broadening (FWHM default 8.0 cm⁻¹), statistical metrics (MAE, RMSE, R², slope/intercept), supports fundamentals, overtones, combination bands.

**HTML reports** (`html_report_generator.py`): `HTMLReportGenerator.generate_report()` combines plots + data.

## Results Directory Structure

```
comparison_results/
  └── {molecule_name}/
      ├── geometry_opt/          # results.json, initial.xyz, optimized.xyz
      └── {energy}_{dipole}/     # results.json, .gjf, .log, .chk, .fchk

analysis_results/                # Anharmonic output
analysis_results_harmonic/       # Harmonic-only output
  └── {molecule_name}/
      ├── plots/
      ├── data/
      └── report.html
```

## Configuration

- `--optimization-calculator`: mace_omol (default), mace_off, mace_mp
- `--energy-calculators`: Comma-separated (default: mace_mp,mace_omol)
- `--dipole-calculators`: Comma-separated (default: espaloma,mace_ml)
- `--skip-dft-baseline`: Skip DFT anharmonic baseline
- `--force-optimization`: Re-optimize even if cached

Env vars: `MACE_DIPOLE_MODEL_PATH`, `MACE_HELPER_SCRIPT_PATH`

## Gaussian Input Generation

**ML runs** (`gm_main.py` `ase_to_gjf()`): Route is `# freq (anharm)` + `# external="<path>"`. No DFT functional (energy/forces from external ML). Must include `%Chk=molecule.chk`. Auto-detects multiplicity.

**DFT baseline** (`dft_baseline.py`): Route is `# opt freq(anharm) b3lyp/6-31G(d,p)`. Pure Gaussian, no external interface.
