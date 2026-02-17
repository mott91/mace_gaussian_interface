---
phase: 04-extract-calculator-classes
plan: 01
subsystem: calculators
tags: [refactoring, abc, factory-pattern, dipole-calculators]

# Dependency graph
requires:
  - phase: 03-extract-utilities-conventions
    provides: "utils/ package pattern, BOHR_TO_ANGSTROM constant in utils/units.py"
provides:
  - "calculators/ package with DipoleCalculatorBase ABC, 3 implementations, factory, singleton"
  - "gm_main.py reduced by ~246 lines"
affects: [05-lazy-mace-loading, 06-extract-gaussian-io]

# Tech tracking
tech-stack:
  added: []
  patterns: [abc-base-class, factory-singleton, lazy-heavy-imports, package-extraction]

key-files:
  created:
    - calculators/__init__.py
    - calculators/base.py
    - calculators/espaloma.py
    - calculators/mace_ml.py
    - calculators/xtb.py
    - calculators/factory.py
  modified:
    - gm_main.py

key-decisions:
  - "Lazy import of MACEDipoleCalculator inside _check_availability to avoid module-level side effects from mace_calculators.py"
  - "Modernized type annotations (tuple instead of Tuple, dict instead of Dict) in new files for Python 3.12 compatibility"

patterns-established:
  - "Calculator package structure: base.py ABC, one file per implementation, factory.py with singleton"
  - "Heavy dependency imports (espaloma, xtb, mace) stay local inside methods, not at module level"

requirements-completed: [STRUCT-02]

# Metrics
duration: 4min
completed: 2026-02-17
---

# Phase 4 Plan 1: Extract Calculator Classes Summary

**Dipole calculator hierarchy (ABC + 3 implementations + factory singleton) extracted from gm_main.py into calculators/ package with 6 modules**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-17T14:55:05Z
- **Completed:** 2026-02-17T14:59:55Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments
- Created calculators/ package with clean separation: base ABC, espaloma, xtb, mace_ml implementations, factory with singleton
- Reduced gm_main.py by 246 lines (1285 -> 1039) while preserving all behavior
- All 92 existing tests pass unchanged

## Task Commits

Each task was committed atomically:

1. **Task 1: Create calculators/ package with extracted classes** - `89e9ab3` (feat)
2. **Task 2: Update gm_main.py to import from calculators package** - `d47a0e6` (refactor)

## Files Created/Modified
- `calculators/__init__.py` - Package re-exports with __all__
- `calculators/base.py` - DipoleCalculatorBase ABC with calculate_dipole and calculate_dipole_derivatives
- `calculators/espaloma.py` - EspalomaDipoleCalculator with lazy espaloma/rdkit imports
- `calculators/mace_ml.py` - MACEMLDipoleCalculator with DEFAULT_MACE_DIPOLE_MODEL constant and lazy MACEDipoleCalculator import
- `calculators/xtb.py` - XTBDipoleCalculator with lazy xtb import
- `calculators/factory.py` - DipoleCalculatorFactory class and dipole_factory singleton
- `gm_main.py` - Removed calculator hierarchy, imports from calculators package

## Decisions Made
- Moved DEFAULT_MACE_DIPOLE_MODEL constant to calculators/mace_ml.py (only consumer)
- Made MACEDipoleCalculator import lazy inside _check_availability() to avoid mace_calculators.py module-level side effects
- Modernized type annotations to lowercase tuple/dict in new files (Python 3.12 codebase)
- Fixed ruff lint issues in extracted code (unused variable test_charges, deprecated typing imports)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed lint issues in extracted code**
- **Found during:** Task 1
- **Issue:** Extracted code had ruff warnings: unused variable (test_charges in espaloma), deprecated typing imports (Tuple, Dict), unused import (XTB in xtb._check_availability)
- **Fix:** Removed unused assignment, modernized type annotations to lowercase, added noqa for XTB availability check import
- **Files modified:** calculators/base.py, calculators/espaloma.py, calculators/xtb.py, calculators/factory.py
- **Verification:** ruff check calculators/ passes clean
- **Committed in:** 89e9ab3 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor code quality improvement during extraction. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- calculators/ package ready for Phase 5 (lazy MACE loading) to modify import patterns in calculators/mace_ml.py
- gm_main.py ready for Phase 6 (Gaussian I/O extraction) to move remaining helper functions

---
*Phase: 04-extract-calculator-classes*
*Completed: 2026-02-17*
