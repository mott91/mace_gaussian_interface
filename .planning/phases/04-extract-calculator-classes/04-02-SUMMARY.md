---
phase: 04-extract-calculator-classes
plan: 02
subsystem: testing
tags: [unit-tests, mocking, abc, factory-pattern, dipole-calculators]

# Dependency graph
requires:
  - phase: 04-extract-calculator-classes
    plan: 01
    provides: "calculators/ package with DipoleCalculatorBase ABC, 3 implementations, factory"
provides:
  - "Unit tests for calculator interface compliance, factory pattern, and config constants"
  - "18 tests covering subclass checks, method presence, ABC enforcement, factory get/auto/list/error, config paths"
affects: [05-lazy-mace-loading]

# Tech tracking
tech-stack:
  added: []
  patterns: [sys-modules-mocking-for-heavy-deps, pre-import-dependency-stubbing]

key-files:
  created:
    - tests/test_calculators.py
  modified:
    - calculators/base.py
    - calculators/factory.py
    - gm_main.py

key-decisions:
  - "Pre-mock heavy dependencies via sys.modules before importing calculators package to avoid DGL/espaloma/xtb side effects"
  - "Added from __future__ import annotations to calculators/base.py, calculators/factory.py, gm_main.py for Python 3.8 compatibility with lowercase type annotations"

patterns-established:
  - "Heavy-dep test isolation: populate sys.modules with MagicMock stubs before importing packages that eagerly instantiate calculators"
  - "Factory testing: patch calculator constructors at their factory import location, not their source modules"

requirements-completed: [STRUCT-02]

# Metrics
duration: 3min
completed: 2026-02-17
---

# Phase 4 Plan 2: Calculator Unit Tests Summary

**18 unit tests for calculator hierarchy covering ABC enforcement, interface compliance, factory pattern (get/auto/list/errors), and config constant env var override**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-17T15:02:35Z
- **Completed:** 2026-02-17T15:05:43Z
- **Tasks:** 1
- **Files modified:** 4

## Accomplishments
- Created 18 tests in 3 test classes (TestCalculatorInterface, TestDipoleCalculatorFactory, TestDefaultMaceDipoleModel)
- All tests run without GPU, espaloma, xtb, or MACE model dependencies
- Fixed Python 3.8 compatibility issue with lowercase type annotations from Phase 04-01

## Task Commits

Each task was committed atomically:

1. **Task 1: Write calculator unit tests with mocked dependencies** - `5aa277c` (test)

## Files Created/Modified
- `tests/test_calculators.py` - 18 unit tests for calculator interface, factory, and config
- `calculators/base.py` - Added `from __future__ import annotations` for Python 3.8
- `calculators/factory.py` - Added `from __future__ import annotations` for Python 3.8
- `gm_main.py` - Added `from __future__ import annotations` for Python 3.8

## Decisions Made
- Pre-mock heavy dependencies (espaloma_charge, rdkit, xtb, mace_calculators) via sys.modules before any calculators import, since the package __init__.py eagerly instantiates the dipole_factory singleton
- Used `from __future__ import annotations` instead of reverting to `Tuple`/`Dict` typing, keeping the modernized annotations from 04-01 while supporting Python 3.8

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added from __future__ import annotations for Python 3.8 compatibility**
- **Found during:** Task 1 (test collection)
- **Issue:** calculators/base.py used `tuple[np.ndarray, ...]` and `dict[str, bool]` annotations (Python 3.10+), but runtime is Python 3.8.13. Test collection failed with `TypeError: 'type' object is not subscriptable`
- **Fix:** Added `from __future__ import annotations` to calculators/base.py, calculators/factory.py, and gm_main.py
- **Files modified:** calculators/base.py, calculators/factory.py, gm_main.py
- **Verification:** Full test suite (117 tests) passes
- **Committed in:** 5aa277c (Task 1 commit)

**2. [Rule 3 - Blocking] Used sys.modules pre-mocking instead of direct submodule imports**
- **Found during:** Task 1 (test collection)
- **Issue:** Importing any calculators submodule triggers `calculators/__init__.py` which imports factory.py, which instantiates the singleton, which tries to import espaloma/DGL (RuntimeError: DGL requires PyTorch >= 1.13.0)
- **Fix:** Pre-populate sys.modules with MagicMock stubs for all heavy deps before importing calculators
- **Files modified:** tests/test_calculators.py
- **Verification:** All 18 calculator tests pass
- **Committed in:** 5aa277c (Task 1 commit)

**3. [Rule 3 - Blocking] Used Python 3.8 compatible context manager syntax**
- **Found during:** Task 1 (factory tests)
- **Issue:** Parenthesized context managers `with (patch(...), patch(...)):` require Python 3.10+
- **Fix:** Used separate variable assignment + comma-separated `with p1, p2, p3:` syntax
- **Files modified:** tests/test_calculators.py
- **Verification:** All factory tests pass
- **Committed in:** 5aa277c (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (3 blocking)
**Impact on plan:** All fixes required for Python 3.8 runtime compatibility. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- calculators/ package fully tested and Python 3.8 compatible
- Ready for Phase 5 (lazy MACE loading) to modify import patterns

---
*Phase: 04-extract-calculator-classes*
*Completed: 2026-02-17*
