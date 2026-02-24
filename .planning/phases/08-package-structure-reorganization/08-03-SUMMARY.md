---
phase: 08-package-structure-reorganization
plan: 03
subsystem: testing
tags: [pyproject, entry-point, imports, packaging, hatchling, ruff, pytest]

# Dependency graph
requires:
  - phase: 08-01
    provides: mace_gaussian/ package created with relative imports
  - phase: 08-02
    provides: all root files moved into mace_gaussian/, analysis/ subpackage created

provides:
  - Correct pyproject.toml: packages = ["mace_gaussian"], entry point mace-gaussian = "mace_gaussian.cli:cli"
  - All 128 tests updated to use mace_gaussian.* import paths and passing
  - mace-gaussian console script working as installed entry point
  - Analysis __init__.py fixed with lazy imports to avoid seaborn at module load time

affects: [cli, workflow, package-distribution, test-maintenance]

# Tech tracking
tech-stack:
  added: [click (installed into venv)]
  patterns: [lazy-wrapper-functions-in-__init__.py, absolute-mace_gaussian-imports-in-tests]

key-files:
  created:
    - .planning/phases/08-package-structure-reorganization/08-03-SUMMARY.md
  modified:
    - pyproject.toml
    - tests/test_gaussian_parser.py
    - tests/test_fchk_parser.py
    - tests/test_calculators.py
    - tests/test_mace_loader.py
    - tests/test_cli_validation.py
    - tests/test_mode_matching.py
    - tests/test_exceptions.py
    - tests/test_units.py
    - tests/test_validation.py
    - tests/test_regression.py
    - mace_gaussian/analysis/__init__.py
    - mace_gaussian/utils/results.py

key-decisions:
  - "Lazy wrapper functions in analysis/__init__.py instead of eager imports — avoids seaborn/pandas import at mace_gaussian.analysis package load time, enabling lightweight test imports"
  - "click installed separately into venv via uv pip install click — uv sync blocked by dgl platform incompatibility on Linux"
  - "Patch targets in tests updated to mace_gaussian.* paths to match where names are bound after reorganization"
  - "results.py lazy import of collect_version_metadata updated from utils.validation to mace_gaussian.utils.validation — was silently returning {} before"
  - "results.py Dict -> dict type annotation modernization done as part of touching the file (Rule 2 scope)"

patterns-established:
  - "All test imports use absolute mace_gaussian.* paths — no flat module shortcuts"
  - "patch() targets must match the fully-qualified module path where the name is bound"
  - "analysis/__init__.py uses lazy wrapper functions for heavy-dep submodules to enable lightweight submodule imports"

requirements-completed: [STRUCT-08, STRUCT-09, STRUCT-10]

# Metrics
duration: 7min
completed: 2026-02-25
---

# Phase 8 Plan 03: Package Wiring and Test Suite Green Summary

**pyproject.toml wired to mace_gaussian package with correct entry points, all 128 tests converted to mace_gaussian.* imports and passing green**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-24T23:22:03Z
- **Completed:** 2026-02-24T23:30:00Z
- **Tasks:** 2
- **Files modified:** 13

## Accomplishments

- pyproject.toml corrected: `mace-gaussian = "mace_gaussian.cli:cli"` entry point, `packages = ["mace_gaussian"]`, updated isort and coverage source
- All 128 tests updated from flat module paths (e.g., `from utils.exceptions import`) to `mace_gaussian.*` absolute imports — all passing
- `mace-gaussian --help` works as installed console script
- `python -c "import mace_gaussian; print(mace_gaussian.__version__)"` prints `0.2.0`

## Task Commits

Each task was committed atomically:

1. **Task 1: Update pyproject.toml and reinstall package** - `248b8d0` (chore)
2. **Task 2: Update all test imports to mace_gaussian.* and verify full suite passes** - `297ca03` (feat)

## Files Created/Modified

- `pyproject.toml` - Fixed entry points (mace-gaussian), packages, isort, coverage config
- `tests/test_gaussian_parser.py` - Updated to mace_gaussian.gaussian.parser, mace_gaussian.utils.exceptions
- `tests/test_fchk_parser.py` - Updated to mace_gaussian.gaussian.fchk
- `tests/test_calculators.py` - Updated to mace_gaussian.calculators.*, updated patch() targets and factory imports
- `tests/test_mace_loader.py` - Updated to mace_gaussian.calculators.mace_loader, updated importlib paths
- `tests/test_cli_validation.py` - Updated to mace_gaussian.cli, mace_gaussian.gaussian.runner, patch target
- `tests/test_mode_matching.py` - Updated to mace_gaussian.analysis.mode_matching, mace_gaussian.gaussian.fchk
- `tests/test_exceptions.py` - Updated to mace_gaussian.utils.exceptions
- `tests/test_units.py` - Updated to mace_gaussian.utils.units
- `tests/test_validation.py` - Updated to mace_gaussian.utils.*, patch decorators updated
- `tests/test_regression.py` - Updated to mace_gaussian.utils.results
- `mace_gaussian/analysis/__init__.py` - Converted eager imports to lazy wrapper functions (seaborn fix)
- `mace_gaussian/utils/results.py` - Fixed collect_version_metadata import path + modernized Dict->dict

## Decisions Made

- **Lazy wrappers in analysis/__init__.py:** The `analysis_workflow.py` chain imports seaborn/pandas/matplotlib at module load. Importing `mace_gaussian.analysis.mode_matching` triggered this chain, failing in test environments. Lazy wrapper functions in `__init__.py` defer the heavy import to call time, enabling lightweight submodule imports.
- **click installed separately:** `uv sync` fails because dgl 2.2.1 has no Linux wheels. Used `uv pip install click --python .venv/bin/python` to install click without touching the lockfile.
- **results.py import bug:** The `collect_version_metadata` lazy import used `from utils.validation import ...` (old flat path). `ImportError` was caught silently, returning `{}` for version_info. Fixed to `from mace_gaussian.utils.validation import ...`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed analysis/__init__.py seaborn import at package load time**
- **Found during:** Task 2 (running test suite after import updates)
- **Issue:** `mace_gaussian/analysis/__init__.py` eagerly imported `analysis_workflow` which transitively imports `seaborn`. Any test importing `from mace_gaussian.analysis.mode_matching import ...` triggered the full chain, failing with `No module named 'seaborn'`
- **Fix:** Replaced eager `from .analysis_workflow import` with lazy wrapper functions that defer the actual import to call time
- **Files modified:** `mace_gaussian/analysis/__init__.py`
- **Verification:** `test_mode_matching.py` passes, 128 tests total pass
- **Committed in:** 297ca03 (Task 2 commit)

**2. [Rule 1 - Bug] Fixed silent ImportError in results.py collect_version_metadata**
- **Found during:** Task 2 (test_frequency_results_contain_version_info failing)
- **Issue:** `results.py` used `from utils.validation import collect_version_metadata` (old flat path). The `try/except ImportError` silently swallowed the error, returning `version_info = {}`. Tests expecting `"tool_version"` in `version_info` failed.
- **Fix:** Updated import to `from mace_gaussian.utils.validation import collect_version_metadata`
- **Files modified:** `mace_gaussian/utils/results.py`
- **Verification:** `test_regression.py::TestResultsManagerMetadata` passes (3 tests)
- **Committed in:** 297ca03 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both auto-fixes essential for test suite correctness. No scope creep.

## Issues Encountered

- `uv sync` blocked by dgl 2.2.1 having no Linux wheels. Worked around by using `uv pip install -e . --python .venv/bin/python --no-deps` to reinstall just the package metadata (entry points), then installing click separately.

## Next Phase Readiness

- Phase 8 complete: mace_gaussian/ is a proper installable Python package with correct entry points and passing test suite
- STRUCT-08 (package structure), STRUCT-09 (entry points), STRUCT-10 (test suite green) all satisfied
- No blockers for Phase 9

---
*Phase: 08-package-structure-reorganization*
*Completed: 2026-02-25*
