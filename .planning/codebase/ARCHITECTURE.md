# Architecture

**Analysis Date:** 2026-02-16

## Pattern Overview

**Overall:** Layered ML-to-Quantum-Chemistry Bridge

**Key Characteristics:**
- Real-time ML dipole injection into Gaussian via ZMQ IPC sockets
- Modular dipole calculator abstraction with pluggable backends
- Dual MACE package system (standard + dipole-enabled) with monkey-patched module loading
- Spectral analysis pipeline with mode matching via eigenvector overlap
- Factory pattern for both calculator selection and results organization

## Layers

**Entry Point & CLI:**
- Purpose: Command-line interface and workflow orchestration
- Location: `cli.py`, `gm_main.py`
- Contains: Click CLI commands, workflow state management, phase sequencing
- Depends on: All calculation and analysis modules
- Used by: End users and automated workflows

**Calculation Core (ML Integration):**
- Purpose: Geometry optimization and frequency calculations with ML potentials
- Location: `gm_main.py` (functions `geometry_optimisation()`, `run_frequency_calculation()`, `run_workflow()`)
- Contains: ASE calculator interface, Gaussian job control, ZMQ server management
- Depends on: `mace_calculators.py`, `gm_helper.py`, ASE, MACE libraries
- Used by: Workflow orchestration in `gm_main.py`

**Dipole Calculation Abstraction:**
- Purpose: Pluggable dipole calculation backends
- Location: `gm_main.py` (classes `DipoleCalculatorBase`, `EspalomaDipoleCalculator`, `MACEMLDipoleCalculator`, `XTBDipoleCalculator`, `DipoleCalculatorFactory`)
- Contains: Abstract base with `calculate_dipole()` and `calculate_dipole_derivatives()`, three implementations
- Depends on: RDKit, espaloma-charge, torch, mace-dipole-core (optional)
- Used by: Frequency calculation dispatch in `run_next_calculation()`

**MACE Module Switching:**
- Purpose: Runtime module monkey-patching to load either standard or dipole-enabled MACE
- Location: `mace_calculators.py`
- Contains: `load_standard_mace_calculator()`, `load_dipole_mace_calculator()`, `cleanup_mace_modules()`, `MACEDipoleCalculator` wrapper
- Depends on: sys, importlib, mace.calculators, mace-dipole-core
- Used by: Both `gm_main.py` and dipole calculators

**Gaussian I/O Bridge:**
- Purpose: Parse external calculation files from Gaussian, compute derivatives, write results
- Location: `gm_main.py` (functions `parse_gaussian_input()`, `update_molecule_geometry()`, `calculate_energy_and_forces()`, `run_next_calculation()`)
- Contains: File I/O, coordinate/dipole conversion (Bohr ↔ Angstrom), numpy array assembly
- Depends on: ASE atoms objects, numpy
- Used by: ZMQ message loop in frequency calculation

**Results Persistence:**
- Purpose: Organize results into molecule/calculation hierarchies, serialize to JSON
- Location: `results_manager.py`
- Contains: `ResultsManager` class with directory creation and JSON metadata writing
- Depends on: pathlib, json, datetime
- Used by: All calculation phases

**Gaussian Log Parser:**
- Purpose: Extract frequencies and IR intensities from Gaussian log files
- Location: `gaussian_parser.py`
- Contains: `GaussianLogParser` class with regex-based section extraction, handling harmonic and anharmonic sections
- Depends on: regex, pathlib
- Used by: Spectral analysis and comparison workflows

**Formatted Checkpoint Parser:**
- Purpose: Extract vibrational modes and metadata from Gaussian .fchk files
- Location: `fchk_parser.py`
- Contains: `convert_chk_to_fchk()`, `parse_fchk_section()`, `extract_modes_from_fchk()`
- Depends on: subprocess (formchk), regex, numpy
- Used by: Mode matching system

**Mode Matching System:**
- Purpose: Match vibrational modes between DFT and ML using eigenvector dot product overlap
- Location: `mode_matching.py`
- Contains: `extract_mode_data_from_checkpoint()`, `compute_reduced_masses()`, `match_modes()`, `create_alignment_matrix()`, `plot_mode_overlap_heatmap()`
- Depends on: fchk_parser, matplotlib, numpy
- Used by: Comparison workflows

**DFT Baseline Calculator:**
- Purpose: Run pure Gaussian DFT calculations as reference spectra
- Location: `dft_baseline.py`
- Contains: `run_dft_frequency_calculation()`, `run_all_dft_baselines()`, DFT method configuration
- Depends on: gaussian_parser, results_manager, subprocess, ASE
- Used by: Workflow PHASE 2

**Spectral Analysis Engine:**
- Purpose: Compare frequencies and intensities with Gaussian KDE broadening, statistical metrics
- Location: `analyze_spectra.py`
- Contains: `SpectrumAnalyzer` class with broadening, `SpectrumData` and `ComparisonMetrics` dataclasses
- Depends on: numpy, scipy.stats, matplotlib, seaborn
- Used by: Comparison workflows

**Comparison & Analysis Orchestration:**
- Purpose: Coordinate DFT baseline + ML results analysis and visualization
- Location: `comparison_workflow.py`
- Contains: `ComparisonWorkflow` class with harmonic/anharmonic mode switching, DFT finder, analysis sequencer
- Depends on: analyze_spectra, mode_matching, html_report_generator
- Used by: `run_analysis.py` and `run_analysis_harmonic.py` convenience scripts

**HTML Report Generation:**
- Purpose: Create interactive web-viewable spectral analysis reports
- Location: `html_report_generator.py`
- Contains: `HTMLReportGenerator` class with image embedding, CSS/JS, table generation
- Depends on: base64, pathlib, pandas
- Used by: Comparison workflows

**Charge Analysis (Optional):**
- Purpose: Compute and visualize partial charges from dipole calculators
- Location: `charge_analysis.py`
- Contains: `ChargeAnalyzer` class with plotting
- Depends on: numpy, matplotlib
- Used by: Optional enhanced analysis

**Environment Diagnostics:**
- Purpose: Validate availability of calculators and dependencies at startup
- Location: `diagnostics.py`
- Contains: `diagnose_python_environment()`, `test_espaloma_functionality()`, dependency checks
- Depends on: importlib, sys
- Used by: `gm_main.diagnose()` CLI command

## Data Flow

**Workflow (Phases 1-3):**

1. **PHASE 1: Geometry Optimization**
   - Load `.xyz` file → ASE Atoms object
   - Attach ML calculator (e.g., MACE-OMOL) → `geometry_optimisation()`
   - LBFGS relaxation until convergence
   - Save optimized geometry + results.json to `comparison_results/{molecule}/geometry_opt/`

2. **PHASE 2: DFT Baseline**
   - For each DFT method (e.g., B3LYP/6-31G(d,p))
   - Generate Gaussian .gjf with route `# opt freq(anharm) b3lyp/6-31G(d,p)`
   - Run pure Gaussian calculation (no external interface)
   - Parse .log for frequencies, save to `comparison_results/{molecule}/{dft_name}/`

3. **PHASE 3: ML Frequency Calculations**
   - For each energy × dipole calculator combination:
     a. Generate Gaussian .gjf with route `# freq(anharm)` + `# external="gm_helper.py"`
     b. Launch `g16 input.gjf` → Gaussian process polls ZMQ socket
     c. Server loop: `is_calc_finished()` checks socket or process exit
     d. Receive `"infile|outfile"` → Parse input file (coordinates in Bohr)
     e. Convert to Angstrom → Update ASE atoms
     f. Attach energy calculator, compute dipole derivatives
     g. Write outfile with energy + dipole derivatives (in e*Bohr)
     h. Send "ready", Gaussian continues
     i. On exit: Parse .log for frequencies, save results.json

**Analysis (Anharmonic Mode):**

1. Find DFT baseline in `comparison_results/{molecule}/` (prefer B3LYP)
2. Collect all ML results (each energy×dipole combination)
3. For each ML result:
   - Extract anharmonic frequencies/intensities from .log
   - Extract harmonic modes from .fchk for mode matching
   - Match ML modes to DFT modes via eigenvector overlap
   - Compute metrics: MAE, RMSE, R², slope/intercept (regression)
4. Generate plots (spectra, regression, mode overlap heatmaps)
5. Create HTML report with embedded plots and tables

**Analysis (Harmonic Mode):**

1. Same DFT baseline finding
2. Force extraction of harmonic fundamentals only (300+ modes if degenerate)
3. Match all modes via eigenvector overlap
4. Mode-by-mode regression (not broadened spectra)
5. Generate mode overlap heatmaps showing alignment quality

**State Management:**

- **Cached**: Geometry optimizations reused if `optimized.xyz` exists (unless `--force-optimization`)
- **Cached**: DFT baselines checked via `check_baseline_exists()`
- **Generated**: Each frequency calculation creates new timestamped directories
- **Parsed**: Results.json files consumed to populate comparison metadata

## Key Abstractions

**DipoleCalculator Abstraction:**
- Purpose: Decouple frequency calculation from dipole backend
- Examples: `EspalomaDipoleCalculator`, `MACEMLDipoleCalculator`, `XTBDipoleCalculator`
- Pattern: Abstract base class with `calculate_dipole()` and `calculate_dipole_derivatives()` methods

**SpectrumData Container:**
- Purpose: Encapsulate frequency data across different analysis modes
- Fields: `frequencies`, `intensities`, `labels` (fundamental/overtone/combination), `mode_ids` (for matching)

**ComparisonMetrics Dataclass:**
- Purpose: Bundle statistical comparison results
- Fields: MAE/RMSE/R² (frequency), slope/intercept, mode matching counts, match rates

**ResultsManager:**
- Purpose: Hide directory structure details, provide consistent APIs
- Methods: `create_molecule_directory()`, `create_frequency_directory()`, `save_optimization_results()`

**MACEDipoleCalculator Wrapper:**
- Purpose: Lazy-load and cache the dipole model, handle format conversions
- Logic: Module switching via `sys.modules` manipulation, tensor-to-numpy conversion

## Entry Points

**CLI Entry Point:**
- Location: `cli.py`
- Triggers: `python cli.py run water.xyz` (main workflow)
- Responsibilities: Click argument parsing, validation, delegate to `gm_main.run_workflow()`

**Workflow Entry Point:**
- Location: `gm_main.run_workflow()`
- Triggers: Called by `cli.py run` or directly by scripts
- Responsibilities: Orchestrate PHASE 1-3, manage results_mgr, return summary dict

**Analysis Entry Points:**
- Location: `run_analysis.py` (anharmonic), `run_analysis_harmonic.py` (harmonic-only)
- Triggers: `python run_analysis.py water`
- Responsibilities: Parse CLI args, call `comparison_workflow.analyze_molecule()` or `analyze_molecule_harmonic()`

**ZMQ Server Loop:**
- Location: `gm_main.run_frequency_calculation()` (lines 930-960)
- Triggers: Spawned as context manager during Gaussian execution
- Responsibilities: Poll socket and process, dispatch `run_next_calculation()` on incoming requests

## Error Handling

**Strategy:** Raise-early with context, log full tracebacks, clean up resources

**Patterns:**

- **ZMQ Cleanup**: `try/finally` blocks in context managers (`zmq_server()`, `zmq_client()`) ensure socket file deletion
- **Module Cleanup**: `cleanup_mace_modules()` in finally blocks after dipole calculations to restore standard MACE
- **File Not Found**: Check existence before parsing (e.g., `GaussianLogParser.__init__()` raises `FileNotFoundError`)
- **External Command Failures**: `subprocess.run(..., check=True)` propagates non-zero exits as `CalledProcessError`
- **Process Monitoring**: `is_calc_finished()` combines socket polling and process exit status to detect hard failures

## Cross-Cutting Concerns

**Logging:** Standard Python logging with module-level loggers; INFO level for workflows, DEBUG for intermediate calculations

**Validation:** Input file extension checking (cli.py), checkpoint file suffix validation (mode_matching.py), Bohr/Angstrom unit conversions with explicit constants

**Authentication:** None (local CLI tool); environment variables for paths (`MACE_DIPOLE_MODEL_PATH`, `MACE_HELPER_SCRIPT_PATH`)

**Unit Conversions:**
- Bohr ↔ Angstrom: `BOHR_TO_ANGSTROM = 0.52917721092` (gm_main.py line 417)
- e*Angstrom → e*Bohr: `dipole / 0.529177210903` (gm_main.py line 205)

**Resource Limits:** No explicit limits; relies on system memory for ASE atoms and arrays; array shapes scale with N_atoms (cubic for Hessian in theory, but Gaussian provides pre-diagonalized frequencies)

---

*Architecture analysis: 2026-02-16*
