# Development Guide

## Adding a New Dipole Calculator

1. Subclass `DipoleCalculatorBase` in `gm_main.py`
2. Implement `_check_availability()` and `calculate_dipole()`
3. Register in `DipoleCalculatorFactory._register_calculators()`
4. Add to `preferred_order` list if appropriate
5. Return format: `(dipole_vector_3d, partial_charges_or_None)`

## Adding a New Analysis Metric

1. Add field to `ComparisonMetrics` dataclass in `analyze_spectra.py`
2. Calculate in `SpectrumAnalyzer.calculate_metrics()`
3. Update `HTMLReportGenerator` to display in report
4. Add to CSV export in `comparison_workflow.py`

## Test Molecules

- `water.xyz` - 3 atoms (fast testing)
- `co.xyz` - Diatomic
- `CH4_ase.xyz` - Tetrahedral symmetry, degenerate modes
- `ammonia.xyz`, `formaldehyde.xyz`
- `acoh.xyz` - Acetic acid (complex overtones/combinations)
- `gly.xyz` - Glycine (amino acid)
- `aspirin.xyz` - Larger drug molecule
- `cocaine.xyz` - Complex ring system
- `crco4.xyz` - Transition metal complex
- n-alkanes: `C2H6_ase.xyz`, `C3H8_ase.xyz`, `n_alkane_4C.xyz` through `n_alkane_8C.xyz`

Always run `python cli.py diagnose` before first calculation.

## Utility Modules

- `diagnostics.py` - Environment diagnostics
- `charge_analysis.py` - Charge visualization (`ChargeAnalyzer`)
- `dft_baseline.py` - DFT baseline management; `DFT_BASELINES` dict, `sanitize_calculator_name()`
- `produce_molecules.py` - n-alkane geometry generation (`build_n_alkane(n_carbons)`)
- `convert_all_chk_files.py` - Batch .chk → .fchk conversion
- `run_analysis_harmonic.py` - Harmonic-only analysis entry point
