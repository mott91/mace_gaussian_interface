# Development Guide

## Adding a New Dipole Calculator

1. Subclass `DipoleCalculatorBase` in `mace_gaussian/calculators/base.py`
2. Implement `_check_availability()` and `calculate_dipole()`
3. Register in `DipoleCalculatorFactory._register_calculators()` (`calculators/factory.py`)
4. Add to `preferred_order` list if appropriate
5. Return format: `(dipole_vector_3d, partial_charges_or_None)`

## Adding a New Analysis Metric

1. Add field to `ComparisonMetrics` dataclass in `analysis/analyze_spectra.py`
2. Calculate in `SpectrumAnalyzer.calculate_metrics()`
3. Update `HTMLReportGenerator` to display in report
4. Add to CSV export in `analysis/analysis_workflow.py`

## Code Quality

```bash
ruff check --fix && ruff format    # Lint + format (line length 100)
ty check                           # Type check
```
