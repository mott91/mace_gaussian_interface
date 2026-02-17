---
phase: 03-extract-utilities-conventions
plan: 01
subsystem: utils
tags: [unit-conversion, exceptions, validation, refactoring, codata]

# Dependency graph
requires:
  - phase: 02-error-handling
    provides: "exceptions.py and validation.py modules"
provides:
  - "utils/ package with units.py, exceptions.py, validation.py"
  - "Centralized CODATA 2018 physical constants"
  - "All imports updated to utils.* paths"
affects: [03-02, 04-calculator-hierarchy, 05-mace-loading, 06-gaussian-io]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Centralized physical constants in utils/units.py", "Package re-exports via __all__ in utils/__init__.py"]

key-files:
  created:
    - utils/__init__.py
    - utils/units.py
  modified:
    - utils/exceptions.py (moved from top-level)
    - utils/validation.py (moved from top-level, updated import)
    - gm_main.py
    - cli.py
    - gaussian_parser.py
    - fchk_parser.py
    - charge_analysis.py
    - dft_baseline.py
    - results_manager.py
    - tests/test_exceptions.py
    - tests/test_validation.py
    - tests/test_gaussian_parser.py
    - pyproject.toml

key-decisions:
  - "Used CODATA 2018 full precision: 0.529177210903 (Bohr) and 27.211386245988 (Hartree-to-eV)"
  - "No compatibility shims at old locations -- clean break with direct import updates"
  - "Used __all__ in utils/__init__.py for re-exports (avoids ruff F401 without redundant aliasing)"
  - "Preserved original numerical behavior by matching constant names to ~0.529 or ~27.21 values, not fixing misnamed variables"

patterns-established:
  - "All unit conversions via utils.units -- no inline constant definitions"
  - "Exception imports from utils.exceptions, validation from utils.validation"

requirements-completed: [STRUCT-01]

# Metrics
duration: 6min
completed: 2026-02-17
---

# Phase 03 Plan 01: Create Utils Package Summary

**Centralized CODATA 2018 physical constants in utils/units.py, relocated exceptions/validation into utils/ package, replaced 6 inline constant definitions across codebase**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-17T14:19:05Z
- **Completed:** 2026-02-17T14:25:20Z
- **Tasks:** 2
- **Files modified:** 16

## Accomplishments
- Created utils/ package with units.py (4 constants, 4 convenience functions), exceptions.py, validation.py
- Replaced all 6 inline unit constant definitions in gm_main.py, fchk_parser.py, charge_analysis.py, dft_baseline.py
- Updated 12+ import statements across 10 files to use utils.* paths
- All 89 tests pass with no behavioral changes

## Task Commits

Each task was committed atomically:

1. **Task 1: Create utils/ package with units.py, exceptions.py, and validation.py** - `d682c60` (feat)
2. **Task 2: Update all imports across codebase to use utils.* paths** - `8d64bf1` (refactor)

## Files Created/Modified
- `utils/__init__.py` - Package init with __all__ re-exports
- `utils/units.py` - CODATA 2018 constants (BOHR_TO_ANGSTROM, HARTREE_TO_EV, etc.) and conversion functions
- `utils/exceptions.py` - Exception hierarchy (moved from top-level, unchanged)
- `utils/validation.py` - Validation functions (moved from top-level, import path updated)
- `gm_main.py` - Replaced 6 inline constant definitions with utils.units imports
- `fchk_parser.py` - Replaced inline BOHR_TO_ANGSTROM with import
- `charge_analysis.py` - Replaced inline ANGSTROM_TO_BOHR with import
- `dft_baseline.py` - Replaced inline HARTREE_TO_EV with import
- `results_manager.py` - Updated lazy import paths for utils.validation
- `pyproject.toml` - Updated isort known-first-party (removed old modules, added utils) and coverage source

## Decisions Made
- Used CODATA 2018 full precision values matching ASE internals
- Preserved original numerical behavior: gm_main.py had variables named ANGSTROM_TO_BOHR and EV_TO_HARTREE with values ~0.529 and ~27.21 respectively (wrong names but correct math). Replaced with BOHR_TO_ANGSTROM and HARTREE_TO_EV to match the actual numerical values.
- No compatibility shims -- codebase is small enough for direct import updates
- Used __all__ for re-exports instead of redundant aliasing (cleaner, ruff-compliant)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed mock targets in test_validation.py**
- **Found during:** Task 2 (import updates)
- **Issue:** Tests patching `validation.shutil.which` failed because module moved to `utils.validation`
- **Fix:** Updated mock targets to `utils.validation.shutil.which`
- **Files modified:** tests/test_validation.py
- **Verification:** All 89 tests pass
- **Committed in:** 8d64bf1 (Task 2 commit)

**2. [Rule 1 - Bug] Updated results_manager.py lazy imports**
- **Found during:** Task 2 (import updates)
- **Issue:** results_manager.py had 2 lazy imports from `validation` that were not in the plan's explicit list
- **Fix:** Updated both to `from utils.validation import collect_version_metadata`
- **Files modified:** results_manager.py
- **Committed in:** 8d64bf1 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- utils/ package is ready for Plan 02 (conventions documentation, ResultsManager extraction)
- All imports consistently use utils.* paths
- pyproject.toml configuration updated for new package structure

---
*Phase: 03-extract-utilities-conventions*
*Completed: 2026-02-17*
