# Architecture

## Dipole Calculator System (`mace_gaussian/calculators/`)

Abstract base class `DipoleCalculatorBase` (`base.py`) with three implementations:
- `MACEMLDipoleCalculator` (`mace_ml.py`) - Custom MACE fork with dipole prediction (primary)
- `EspalomaDipoleCalculator` (`espaloma.py`) - ML charge-based dipole (reliable fallback)
- `XTBDipoleCalculator` (`xtb.py`) - Semi-empirical xTB method

`DipoleCalculatorFactory` (`factory.py`) manages selection with preferred order: `["mace_ml", "espaloma", "xtb"]`.

## Dual MACE Package Architecture

Two separate MACE installations coexist:
- `mace_ML_pkg/mace/` - Standard MACE-torch (energy/forces)
- `mace_dipole_pkg/mace_dipole_core/` - Custom fork with dipole moments

`calculators/mace_loader.py` handles loading the dipole model via `torch.load(pickle_module=...)` class remapping (`mace.modules.models` → `mace_dipole_core.modules.models`). No `sys.modules` mutation needed.

## Gaussian Interface (`mace_gaussian/gaussian/`)

- `io.py` - Input/output: `ase_to_gjf()` generates `.gjf` files, parser functions
- `runner.py` - Gaussian 16 execution
- `parser.py` - Log file parsing (frequencies, intensities, overtones, combinations)
- `fchk.py` - Formatted checkpoint parsing, `.chk` → `.fchk` conversion via `formchk`
- `zmq_server.py` - ZMQ IPC server for real-time dipole injection

### ZMQ Bridge

1. Gaussian launches `gm_helper.py` via `External="<path>"`
2. Helper connects to ZMQ socket bound to `.ipc_file`
3. Main script receives request, calculates dipoles, writes output
4. Helper receives "done", Gaussian continues

## Analysis Framework (`mace_gaussian/analysis/`)

- `analysis_workflow.py` - `ComparisonWorkflow` orchestrator; `analyze_molecule()` / `analyze_molecule_harmonic()` entry points
- `analyze_spectra.py` - `SpectrumAnalyzer`: KDE broadening, statistical metrics (MAE, RMSE, R², slope/intercept)
- `mode_matching.py` - Eigenvector dot-product mode matching, overlap heatmaps
- `html_report_generator.py` - HTML report combining plots + data
- `batch_report.py` - Multi-molecule batch reports
- `nist_fetcher.py` - NIST experimental spectrum fetching
- `coverage_analysis.py` - Frequency range coverage analysis

## Top-level Modules

- `workflow.py` - Main calculation workflow orchestration
- `gm_helper.py` - ZMQ bridge script launched by Gaussian
- `dft_baseline.py` - DFT baseline management; `DFT_BASELINES` dict
- `batch.py` - Batch calculation pipeline
- `cli.py` - CLI entry point
- `slurm.py` - SLURM cluster submission
- `diagnostics.py` - Environment diagnostics

## Env Vars

- `MACE_DIPOLE_MODEL_PATH` - Override dipole model path (default: `~/mace_gaussian/mace4ir_models/pretrained_models/model_1_dipole.model`)
- `MACE_HELPER_SCRIPT_PATH` - Path to `gm_helper.py` for Gaussian External
