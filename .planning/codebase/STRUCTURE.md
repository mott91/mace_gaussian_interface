# Codebase Structure

**Analysis Date:** 2026-02-16

## Directory Layout

```
mace_gaussian/
├── cli.py                           # Click CLI entry point
├── gm_main.py                       # Main workflow orchestration
├── gm_helper.py                     # ZMQ helper script called by Gaussian
├── results_manager.py               # Results directory/JSON management
├── mace_calculators.py              # MACE module loading + switching
├── gaussian_parser.py               # Gaussian .log file parser
├── fchk_parser.py                   # Gaussian .fchk file parser
├── mode_matching.py                 # Vibrational mode matching via eigenvector overlap
├── dft_baseline.py                  # DFT reference calculation runner
├── analyze_spectra.py               # Spectral analysis engine (KDE broadening, stats)
├── comparison_workflow.py           # Orchestrator for analysis (anharmonic/harmonic modes)
├── html_report_generator.py         # HTML report generation
├── charge_analysis.py               # Optional charge visualization (unused in main workflow)
├── diagnostics.py                   # Environment validation at startup
├── run_analysis.py                  # Convenience script: anharmonic analysis
├── run_analysis_harmonic.py         # Convenience script: harmonic-only analysis
├── convert_all_chk_files.py         # Batch .chk → .fchk converter
├── produce_molecules.py             # Utility for test molecule generation (unused)
│
├── comparison_results/              # Output: ML calculation results (by user)
│   └── {molecule_name}/
│       ├── geometry_opt/
│       │   ├── results.json
│       │   ├── initial.xyz
│       │   └── optimized.xyz
│       ├── {dft_method_name}/       # E.g., "b3lyp_6-31Gdp" (DFT baseline)
│       │   ├── results.json
│       │   ├── gaussian_freq.gjf
│       │   ├── gaussian_freq.log
│       │   ├── gaussian_freq.chk
│       │   └── gaussian_freq.fchk
│       └── {energy}_{dipole}/       # E.g., "mace_mp_espaloma" (ML frequency)
│           ├── results.json
│           ├── gaussian_freq.gjf
│           ├── gaussian_freq.log
│           ├── gaussian_freq.chk
│           ├── gaussian_freq.fchk
│           └── (intermediate .ipc files, removed after success)
│
├── analysis_results/                # Output: Anharmonic analysis (by user)
│   └── {molecule_name}/
│       ├── plots/
│       │   ├── spectrum_combined.png
│       │   ├── regression_combined.png
│       │   ├── mode_overlap_*.png
│       │   └── (method-specific plots)
│       ├── data/
│       │   ├── metrics_summary.json
│       │   └── (frequency tables)
│       └── report.html
│
├── analysis_results_harmonic/       # Output: Harmonic-only analysis (by user)
│   └── {molecule_name}/
│       ├── plots/
│       │   ├── mode_overlap_*.png
│       │   └── (heatmaps)
│       ├── data/
│       │   └── metrics_summary.json
│       └── report.html
│
├── docs/
│   ├── ARCHITECTURE.md              # Existing architecture notes
│   └── DEVELOPMENT.md               # Development guidelines
│
├── mace_ML_pkg/                     # Standard MACE-torch installation
│   ├── mace/
│   │   ├── calculators/
│   │   └── modules/
│   ├── scripts/
│   └── tests/
│
├── mace_dipole_pkg/                 # Custom MACE fork with dipole support
│   ├── mace_dipole_core/
│   │   ├── calculators/
│   │   │   └── mace.py              # MACECalculator_dipole class
│   │   └── modules/
│   │       └── models.py            # (monkey-patched)
│   ├── scripts/
│   └── tests/
│
├── dipole_model/
│   └── model_1.model                # Pre-trained MACE dipole moment model
│
├── .planning/
│   └── codebase/                    # GSD mapping output (this directory)
│       ├── ARCHITECTURE.md
│       └── STRUCTURE.md
│
├── .venv/                           # Python virtual environment (uv sync)
│   └── lib/python3.12/site-packages/
│       ├── mace/
│       ├── ase/
│       ├── torch/
│       └── ...
│
├── pyproject.toml                   # Project metadata, dependencies, tool config
├── uv.lock                          # Locked dependency versions
├── CLAUDE.md                        # Project guidelines and critical gotchas
└── [molecule files]                 # Input files
    ├── water.xyz
    ├── co.xyz
    ├── acoh.xyz
    └── ...
```

## Directory Purposes

**Root Level (Calculation & CLI):**
- `cli.py`: Click command-line interface with `run`, `list`, `compare`, `export`, `diagnose` commands
- `gm_main.py`: Core workflow with `run_workflow()`, `run_frequency_calculation()`, calculation loops
- `gm_helper.py`: Standalone ZMQ client script launched by Gaussian via `External=`
- Molecule `.xyz` files: Input structures for calculations

**Utilities (Calculation Support):**
- `results_manager.py`: Directory creation, results.json serialization
- `mace_calculators.py`: Module monkey-patching for dual MACE loading
- `diagnostics.py`: Startup checks for availability of calculators

**Parsers (I/O):**
- `gaussian_parser.py`: Regex-based .log file parsing (frequencies, intensities)
- `fchk_parser.py`: .fchk file parsing (modes, coordinates, masses, formchk subprocess call)

**DFT Baseline:**
- `dft_baseline.py`: Pure Gaussian runner for B3LYP/6-31G(d,p) reference spectra

**Analysis & Visualization:**
- `mode_matching.py`: Eigenvector-based mode matching, overlap heatmaps
- `analyze_spectra.py`: Gaussian KDE broadening, statistical metrics (MAE/RMSE/R²)
- `comparison_workflow.py`: Orchestration for anharmonic/harmonic analysis, DFT finder
- `html_report_generator.py`: Web report generation with embedded images and tables
- `charge_analysis.py`: (Optional) partial charge analysis/plotting

**Convenience Scripts:**
- `run_analysis.py`: One-line entry for anharmonic analysis: `python run_analysis.py water`
- `run_analysis_harmonic.py`: One-line entry for harmonic-only analysis
- `convert_all_chk_files.py`: Batch convert .chk → .fchk

**External Packages:**
- `mace_ML_pkg/`: Standard MACE-torch with energy/force calculators
- `mace_dipole_pkg/`: Custom fork with dipole prediction (via neural network)
- `dipole_model/model_1.model`: Pre-trained weights for dipole model

**Output Directories (Created at Runtime):**
- `comparison_results/{molecule}/`: ML calculation artifacts (JSON metadata + Gaussian files)
- `analysis_results/{molecule}/`: Anharmonic analysis plots, data, HTML report
- `analysis_results_harmonic/{molecule}/`: Harmonic-only analysis (mode matching emphasis)

**Documentation & Config:**
- `docs/`: Architecture notes (existing), development guidelines
- `pyproject.toml`: Python package metadata, ruff/type-check config
- `CLAUDE.md`: Critical gotchas (module cleanup, absolute paths, ZMQ sockets)

## Key File Locations

**Entry Points:**
- CLI: `cli.py` (dispatches to `gm_main.py`)
- Workflow: `gm_main.py:run_workflow()` (phases 1-3)
- Analysis: `run_analysis.py` or `comparison_workflow.py:analyze_molecule()`

**Configuration:**
- Package deps: `pyproject.toml`
- Env vars: `MACE_DIPOLE_MODEL_PATH` (default: `~/mace_gaussian/dipole_model/model_1.model`)
- Env vars: `MACE_HELPER_SCRIPT_PATH` (default: `./gm_helper.py`)
- DFT methods: `dft_baseline.py` (lines 46-67) `DFT_BASELINES` dict

**Core Logic:**
- Geometry optimization: `gm_main.py:geometry_optimisation()`
- Frequency calculation: `gm_main.py:run_frequency_calculation()` + `run_next_calculation()`
- ZMQ loop: `gm_main.py:run_frequency_calculation()` (lines 930-960) with `is_calc_finished()`
- Dipole selection: `gm_main.py:DipoleCalculatorFactory.get_calculator()`

**Results Processing:**
- Directory structure: `results_manager.py:ResultsManager`
- Log parsing: `gaussian_parser.py:GaussianLogParser`
- Mode extraction: `fchk_parser.py:extract_modes_from_fchk()`
- Mode matching: `mode_matching.py:match_modes()`

**Analysis:**
- Spectral broadening: `analyze_spectra.py:SpectrumAnalyzer.broaden_spectrum()`
- Regression metrics: `analyze_spectra.py:compute_comparison_metrics()`
- HTML generation: `html_report_generator.py:HTMLReportGenerator.generate_report()`

**Testing:**
- (No dedicated test directory; testing via CLAUDE.md gotchas)

## Naming Conventions

**Files:**

- Input molecules: `{name}.xyz` (e.g., `water.xyz`, `acoh.xyz`)
- Gaussian I/O:
  - Temporary (during workflow): `temp_{energy}_{dipole}_{timestamp}.gjf` → `.log`, `.chk`, `.fchk`
  - Final (in results dir): `gaussian_freq.gjf` → `.log`, `.chk`, `.fchk`
- DFT calculator names: `{method}_{basis}` (e.g., `b3lyp_6-31Gdp`, stripped of special chars)
- ML calculator names: `{energy}_{dipole}` (e.g., `mace_mp_espaloma`)
- Plots: `{plot_type}_{calc1}_vs_{calc2}.png` (e.g., `regression_mace_mp_espaloma_vs_b3lyp_6-31Gdp.png`)

**Directories:**

- Molecule result dirs: `comparison_results/{molecule_name}/`
- Frequency calc dirs: `{energy_calc}_{dipole_calc}/` or `{dft_method_name}/` for DFT
- Geometry opt: Always `geometry_opt/`
- Analysis output: `analysis_results/{molecule_name}/` or `analysis_results_harmonic/{molecule_name}/`
- Subfolders: `plots/`, `data/`

**Classes & Functions:**

- Calculators (dipole): `{Name}DipoleCalculator` (e.g., `EspalomaDipoleCalculator`)
- Factories: `{Name}Factory` (e.g., `DipoleCalculatorFactory`)
- Managers: `{Name}Manager` (e.g., `ResultsManager`)
- Parsers: `{Name}Parser` (e.g., `GaussianLogParser`)
- Analyzers: `{Name}Analyzer` (e.g., `SpectrumAnalyzer`)
- Private functions: `_snake_case_with_leading_underscore`
- Public functions: `snake_case`

**Variables:**

- Atoms objects: `atoms`, `mol`, `optimized_mol`
- Frequencies: `frequencies` (numpy array), `freq_cm` (single value in cm⁻¹)
- Coordinates: `coordinates`, `coords` (shape (n_atoms, 3) in Angstrom)
- Mode data: `modes` (shape (n_modes, n_atoms, 3))
- Paths: Use `Path` objects from pathlib; file paths are absolute in Gaussian .gjf files

## Where to Add New Code

**New Dipole Calculator:**
1. Create subclass of `DipoleCalculatorBase` in `gm_main.py` (lines 85-277)
2. Implement `_check_availability()` and `calculate_dipole(atoms, **kwargs)` → (dipole_vector, partial_charges)
3. Register in `DipoleCalculatorFactory._register_calculators()` (line 287)
4. Add to `preferred_order` list (line 284) if desired
5. Test via `python cli.py diagnose`

**New DFT Method:**
1. Add entry to `DFT_BASELINES` dict in `dft_baseline.py` (lines 48-67)
2. Key = optimization_calculator name (e.g., 'mace_omol')
3. Value = dict with 'method', 'basis', 'description', 'extra_keywords'
4. Code reuses `run_dft_frequency_calculation()` unchanged

**New Analysis Metric:**
1. Add field to `ComparisonMetrics` dataclass in `analyze_spectra.py` (lines 50-65)
2. Compute in `SpectrumAnalyzer.compute_comparison_metrics()` (called from comparison_workflow.py)
3. Populate in `ComparisonWorkflow.compare_ml_with_dft()` or `compare_ml_with_ml()`
4. Add to HTML table in `html_report_generator.py`

**New CLI Command:**
1. Add `@cli.command()` function in `cli.py`
2. Use Click decorators: `@click.argument()`, `@click.option()`, `@click.echo()`
3. Delegate to module functions (don't put logic in cli.py)

**New Calculation Mode (e.g., Raman):**
1. Create new dipole derivative calculator or frequency type
2. Add new frequency directory type to `ResultsManager` if needed
3. Extend `run_next_calculation()` to handle new external request type
4. Update Gaussian input generation (`ase_to_gjf()`) with new route

**Unit Tests:**
- (No existing pytest structure; add `tests/` folder with `test_*.py` files)
- Test locations: `tests/test_gaussian_parser.py`, `tests/test_mode_matching.py`, etc.
- Run: `pytest tests/` (configured in `pyproject.toml`)

## Special Directories

**Temporary Files During Frequency Calculation:**
- Location: Cwd (same as Gaussian invocation)
- Names: `temp_{energy}_{dipole}_{timestamp}.{gjf,log,chk,fchk}`
- Cleanup: Moved to results directory after success, deleted on error
- Generated: One per energy×dipole combination, sequential

**IPC Socket File:**
- Location: `.ipc_file` in cwd
- Purpose: ZMQ server listens here for incoming Gaussian external calculation requests
- Lifetime: Created at context entry, deleted on exit (in `zmq_server()` finally block)
- Critical Gotcha: If previous run crashed, `.ipc_file` persists and blocks new runs (manual deletion needed)

**.fchk Conversion Cache:**
- Location: Same directory as input `.chk` file
- Pattern: Input `file.chk` → Output `file.fchk` via `formchk` subprocess
- Cleanup: None (fchk files are small ~100 KB, preserved for mode matching)
- Dependency: Requires Gaussian `formchk` utility in PATH

**HTML Reports:**
- Location: `analysis_results/{molecule}/report.html` or `analysis_results_harmonic/{molecule}/report.html`
- Embedding: All images base64-encoded inside HTML (self-contained, can email/archive)
- Styling: Internal CSS + responsive layout (no external dependencies)

---

*Structure analysis: 2026-02-16*
