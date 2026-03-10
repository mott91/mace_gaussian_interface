# Architecture Research

**Research Date:** 2026-02-16
**Domain:** Refactoring scientific Python code with heavy external dependencies

## Current Architecture Issues

The codebase has a monolithic main module (`gm_main.py`, ~1000 lines) that mixes:
- Workflow orchestration (phase sequencing)
- Gaussian I/O (file parsing, writing)
- ZMQ server management
- Dipole calculator classes (3 implementations + factory)
- Unit conversion utilities
- Geometry optimization
- Frequency calculation logic

This is the primary structural problem. Everything else (parsers, analysis, reports) is already reasonably modular.

## Target Module Structure

### Recommended decomposition of gm_main.py:

```
mace_gaussian/              # Package directory (rename from flat layout)
  __init__.py               # Version, package metadata
  cli.py                    # Click CLI (already separate, good)
  workflow.py               # Phase orchestration (extract from gm_main)

  calculators/
    __init__.py
    base.py                 # DipoleCalculatorBase (abstract)
    espaloma.py             # EspalomaDipoleCalculator
    mace_ml.py              # MACEMLDipoleCalculator
    xtb.py                  # XTBDipoleCalculator
    factory.py              # DipoleCalculatorFactory
    mace_loader.py          # Module loading (replace monkey-patching)

  gaussian/
    __init__.py
    io.py                   # parse_gaussian_input, write output files
    runner.py               # Subprocess management, timeout handling
    parser.py               # GaussianLogParser (move from gaussian_parser.py)
    fchk.py                 # FCHK parsing (move from fchk_parser.py)
    zmq_server.py           # ZMQ context manager and message loop

  analysis/
    __init__.py
    spectra.py              # SpectrumAnalyzer (move from analyze_spectra.py)
    mode_matching.py        # Mode matching (move from mode_matching.py)
    comparison.py           # ComparisonWorkflow (move from comparison_workflow.py)
    report.py               # HTML report generation

  utils/
    __init__.py
    units.py                # BOHR_TO_ANGSTROM, conversion functions
    validation.py           # Input validation, prerequisite checks
    results.py              # ResultsManager (move from results_manager.py)
```

### What NOT to refactor:

- **analyze_spectra.py** — already clean, self-contained. Just move.
- **mode_matching.py** — already modular. Just move.
- **html_report_generator.py** — works, self-contained. Just move.
- **comparison_workflow.py** — clean orchestrator. Just move.
- **results_manager.py** — simple, works. Just move.

The main refactoring effort is splitting `gm_main.py` and replacing `mace_calculators.py` monkey-patching.

## Replacing Module Monkey-Patching

The current approach manipulates `sys.modules` to swap MACE implementations at runtime. This is the highest-risk fragility.

### Recommended: Process isolation

**Pattern:** Run dipole-enabled MACE calculations in a subprocess with a clean Python interpreter. This completely avoids module contamination.

```python
# Instead of monkey-patching sys.modules:
import subprocess
result = subprocess.run(
    [sys.executable, "-m", "mace_gaussian.calculators.mace_ml",
     "--atoms", atoms_file, "--output", result_file],
    capture_output=True
)
```

**Tradeoff:** Adds IPC overhead (~10ms per call) but eliminates the entire class of module contamination bugs.

### Alternative: Lazy import isolation

```python
# calculators/mace_loader.py
def get_dipole_calculator():
    """Import dipole MACE only when needed, in isolated scope."""
    import importlib
    # Use importlib to load specific module without polluting sys.modules
    spec = importlib.util.spec_from_file_location(
        "mace_dipole", "/path/to/mace_dipole_pkg/..."
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.MACEDipoleCalculator
```

**Tradeoff:** Less overhead but still requires careful scope management. Less robust than process isolation.

### Recommendation: Start with lazy import isolation (simpler), fall back to process isolation if contamination issues persist.

## Testing Strategy for External Dependencies

### Tier 1: Unit tests (no external deps)
- Parser tests with committed .log/.fchk snippets
- Mode matching with synthetic eigenvector matrices
- Metric calculations with known inputs/outputs
- Unit conversion functions
- CLI argument parsing

### Tier 2: Integration tests (mocked externals)
- Gaussian subprocess calls mocked with saved outputs
- ZMQ communication with mock server/client
- Calculator factory with mocked backends

### Tier 3: System tests (real externals, CI-excluded)
- Full water molecule workflow (requires Gaussian + GPU)
- Marked with `@pytest.mark.gaussian` and `@pytest.mark.gpu`
- Run manually or in HPC CI environment

### Test data strategy:
```
tests/test_data/
  water_freq.log           # Real Gaussian output (sanitized)
  water_freq.fchk          # Real FCHK (small molecule)
  water_results.json       # Expected parsed output
  ch4_anharm.log           # Anharmonic output for overtone/combination tests
  acoh_problematic.log     # The acetic acid edge case
```

## Refactoring Order

**Critical constraint:** Each refactoring step must keep the tool functional. Never have a broken state on the branch.

### Suggested order:

1. **Add tests first** — test existing behavior before changing anything. This is the safety net.
2. **Extract utilities** — units.py, validation.py. Pure functions, zero risk.
3. **Extract calculators** — DipoleCalculator classes into calculators/ package. Factory pattern already exists.
4. **Replace monkey-patching** — mace_loader.py with lazy imports. Test thoroughly.
5. **Extract Gaussian I/O** — io.py, runner.py, zmq_server.py. These are tightly coupled, do them together.
6. **Extract workflow** — workflow.py as thin orchestrator calling modular components.
7. **Move analysis modules** — already modular, just reorganize into package structure.
8. **Update imports everywhere** — CLI, scripts, tests.

### Each step: write tests → refactor → verify tests pass → commit.

## Data Flow After Refactoring

```
cli.py
  → workflow.run_workflow()
    → Phase 1: workflow.geometry_optimization()
        → calculators.mace_loader.get_energy_calculator()
        → ASE LBFGS optimization
        → utils.results.save_optimization()

    → Phase 2: gaussian.runner.run_dft_baseline()
        → gaussian.io.write_gjf()
        → subprocess with timeout
        → gaussian.parser.parse_log()

    → Phase 3: workflow.run_frequency_calculation()
        → gaussian.io.write_gjf() (with external hook)
        → gaussian.zmq_server.serve() (context manager)
            → gaussian.io.parse_input()
            → calculators.factory.get_calculator()
            → calculator.calculate_dipole_derivatives()
            → gaussian.io.write_output()
        → gaussian.parser.parse_log()
        → utils.results.save_frequency()
```

---
*Architecture research: 2026-02-16*
