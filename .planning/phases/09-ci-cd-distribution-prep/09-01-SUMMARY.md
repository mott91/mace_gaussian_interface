---
phase: 09-ci-cd-distribution-prep
plan: 01
subsystem: code-quality
tags: [ruff, linting, type-annotations, pathlib, cleanup]
dependency_graph:
  requires: []
  provides: [clean-ruff-baseline, ci-ready-linting]
  affects: [all-mace_gaussian-modules, tests]
tech_stack:
  added: []
  patterns: [pathlib-over-os-path, builtin-type-annotations, contextlib-suppress]
key_files:
  created: []
  modified:
    - pyproject.toml
    - mace_gaussian/gaussian/parser.py
    - mace_gaussian/gaussian/zmq_server.py
    - mace_gaussian/gaussian/fchk.py
    - mace_gaussian/gaussian/io.py
    - mace_gaussian/gm_helper.py
    - mace_gaussian/cli.py
    - mace_gaussian/workflow.py
    - mace_gaussian/dft_baseline.py
    - mace_gaussian/diagnostics.py
    - mace_gaussian/calculators/base.py
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/analysis_workflow.py
    - mace_gaussian/analysis/html_report_generator.py
    - mace_gaussian/analysis/mode_matching.py
    - mace_gaussian/utils/results.py
    - tests/conftest.py
    - tests/test_calculators.py
    - tests/test_fchk_parser.py
    - tests/test_exceptions.py
decisions:
  - "Use # noqa: F401 for intentional import-for-testing in diagnostics.py (not worth removing)"
  - "Change water_dft_fchk fixture to return Path object (not str) to enable .open() in tests"
  - "B017: pytest.raises(Exception) replaced with specific exception type for test clarity"
metrics:
  duration: ~90 minutes (continued from previous session)
  completed: 2026-02-26
  tasks_completed: 2
  files_modified: 23
---

# Phase 09 Plan 01: Expand Ruff Rules and Fix All Violations Summary

**One-liner:** Expanded ruff lint rules to B, SIM, PTH, RUF, and fixed all 55 violations across 23 files, achieving zero lint errors, zero format issues, and 128 passing tests.

## What Was Built

Added four ruff rule sets (B/SIM/PTH/RUF) to `pyproject.toml` and eliminated all resulting violations through two-pass approach: auto-fixes first, then systematic manual fixes across the entire codebase.

## Tasks Completed

### Task 1: Expand ruff rules and auto-fix (commit: 10dae2a)

- Updated `pyproject.toml`: `select = ["E", "F", "W", "I", "UP"]` → `select = ["E", "F", "W", "I", "UP", "B", "SIM", "PTH", "RUF"]`
- Applied `ruff --fix` (13 auto-fixable violations: RUF022 sorted `__all__`, RUF100 unused noqa, SIM300 yoda conditions)
- Applied `ruff format` (12 files reformatted)
- Fixed stale test: `test_path_suffix` assertion updated from `dipole_model/model_1.model` to `mace4ir_models/pretrained_models/model_1_dipole.model`

### Task 2: Manual fixes for remaining 55 violations (commit: 26a98ff)

Fixed violations across 13 files:

**Type annotation modernization (UP035/UP006/UP045/F821):**
- Removed `Dict`, `List`, `Tuple` from `typing` imports in `parser.py`, `dft_baseline.py`
- Replaced all `Dict[...]`, `List[...]`, `Tuple[...]` annotations with `dict[...]`, `list[...]`, `tuple[...]`
- Replaced `Optional[np.ndarray]` with `np.ndarray | None` in `calculators/base.py`

**Pathlib migration (PTH100/PTH107/PTH110/PTH123):**
- Converted all `open(path)` → `path.open()` or `Path(path).open()` across 8 files
- Converted `os.path.abspath()` → `Path(...).resolve()` in `zmq_server.py`, `gm_helper.py`
- Converted `os.path.exists()` → `Path(...).exists()` in `workflow.py`, `results.py`
- Converted `os.remove()` → `Path.unlink()` in `zmq_server.py`
- Dropped `import os` from `gm_helper.py`, `utils/results.py`

**Simplification (SIM105):**
- `zmq_server.py` cleanup: `try: os.remove(path) except OSError: pass` → `contextlib.suppress(OSError)` with `unlink(missing_ok=True)`

**Bug-bearing patterns (B017/B028):**
- `test_exceptions.py`: `pytest.raises(Exception)` → `pytest.raises(GaussianParseError, match="bad output")`
- `test_exceptions.py`: `warnings.warn(...)` → `warnings.warn(..., stacklevel=2)`

**Unused variable cleanup (RUF059):**
- `test_fchk_parser.py`: `modes, frequencies, coords, masses, n_atoms` → `modes, frequencies, _, _, n_atoms`

**Line length (E501):**
- `dft_baseline.py`: Split two long print strings using implicit string concatenation
- `gaussian/parser.py`: Split long inline comment across two lines

**Format:**
- `analysis/mode_matching.py`: Applied `ruff format` after type annotation changes

## Verification Results

```
ruff check mace_gaussian/ tests/   -> All checks passed (0 violations)
ruff format --check mace_gaussian/ tests/ -> 41 files already formatted
pytest tests/ -q -> 128 passed, 2 skipped, 1 xfailed in 3.99s
```

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing] cli.py had 2 PTH123 violations not in original violation list**
- **Found during:** Task 2 (ruff check)
- **Issue:** `mace_gaussian/cli.py` had two `open()` calls on `Path` objects that ruff flagged
- **Fix:** Changed both to `.open()` calls
- **Files modified:** `mace_gaussian/cli.py`
- **Commit:** 26a98ff

**2. [Rule 1 - Bug] conftest.py water_dft_fchk fixture returned str not Path**
- **Found during:** Task 2, after changing test_fchk_parser.py to use `.open()`
- **Issue:** `water_dft_fchk` fixture returned `str`, but PTH123 fix required `.open()` which only works on `Path`
- **Fix:** Changed fixture to return `Path` object directly (without `str()` wrapper)
- **Files modified:** `tests/conftest.py`
- **Commit:** 26a98ff

## Self-Check: PASSED

| Item | Result |
|------|--------|
| 09-01-SUMMARY.md exists | FOUND |
| pyproject.toml exists | FOUND |
| zmq_server.py exists | FOUND |
| commit 10dae2a exists | FOUND |
| commit 26a98ff exists | FOUND |
| ruff check mace_gaussian/ tests/ | 0 violations |
| ruff format --check | 41 files already formatted |
| pytest tests/ | 128 passed, 2 skipped, 1 xfailed |
