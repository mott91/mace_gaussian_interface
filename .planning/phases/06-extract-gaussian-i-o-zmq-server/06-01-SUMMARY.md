---
phase: 06-extract-gaussian-i-o-zmq-server
plan: "01"
subsystem: gaussian
tags: [gaussian, exceptions, io, zmq, ase, numpy]

# Dependency graph
requires:
  - phase: 03-physical-constants-and-units
    provides: BOHR_TO_ANGSTROM and HARTREE_TO_EV constants used in gaussian/io.py
  - phase: 02-project-structure-and-validation
    provides: utils/exceptions.py base exception hierarchy
provides:
  - GaussianRunError and GaussianTimeoutError exception classes in utils/exceptions.py
  - gaussian/__init__.py package stub
  - gaussian/io.py with parse_gaussian_input, write_gaussian_output, ase_to_gjf, DEFAULT_HELPER_SCRIPT
affects:
  - 06-02 (ZMQ server module will import from gaussian/)
  - 06-03 (runner module will raise GaussianRunError and GaussianTimeoutError)
  - 06-05 (caller wiring - gm_main.py will import from gaussian.io)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Extracting pure I/O functions from monolithic gm_main.py into focused gaussian/io.py module
    - gaussian/ package established as home for all Gaussian integration code

key-files:
  created:
    - gaussian/__init__.py
    - gaussian/io.py
  modified:
    - utils/exceptions.py

key-decisions:
  - "DEFAULT_HELPER_SCRIPT path in gaussian/io.py uses Path(__file__).parent.parent to reference gm_helper.py at project root (one level up from gaussian/)"
  - "gm_main.py callers not updated in this plan - import wiring deferred to plan 06-05 to minimize diff scope per task"
  - "gaussian/__init__.py is a sparse stub - full re-exports added in plan 06-05 after all submodules exist"

patterns-established:
  - "gaussian/ submodules import from utils.units rather than redefining constants"
  - "Pure I/O functions extracted without behavioral changes - body verbatim from gm_main.py"

requirements-completed: [STRUCT-04, ERR-06]

# Metrics
duration: 2min
completed: 2026-02-20
---

# Phase 6 Plan 01: gaussian/ Package Foundation Summary

**gaussian/ package created with GaussianRunError/GaussianTimeoutError exceptions and extracted I/O functions (parse_gaussian_input, write_gaussian_output, ase_to_gjf) from gm_main.py**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-20T08:59:21Z
- **Completed:** 2026-02-20T09:01:30Z
- **Tasks:** 2
- **Files modified:** 3 (1 modified, 2 created)

## Accomplishments
- Added GaussianRunError and GaussianTimeoutError as MaceGaussianError subclasses documenting intended use by gaussian.runner
- Created gaussian/__init__.py package stub enabling `import gaussian` to work
- Created gaussian/io.py extracting all 4 I/O symbols (parse_gaussian_input, write_gaussian_output, ase_to_gjf, DEFAULT_HELPER_SCRIPT) with proper imports from utils.units and ase.data

## Task Commits

Each task was committed atomically:

1. **Task 1: Add GaussianRunError and GaussianTimeoutError to utils/exceptions.py** - `e430267` (feat)
2. **Task 2: Create gaussian/ package with stub __init__.py and io.py** - `b31e7d6` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `utils/exceptions.py` - Added GaussianRunError and GaussianTimeoutError after CUDANotAvailableWarning
- `gaussian/__init__.py` - Package docstring stub
- `gaussian/io.py` - Extracted I/O layer: parse_gaussian_input, write_gaussian_output, ase_to_gjf, DEFAULT_HELPER_SCRIPT

## Decisions Made
- DEFAULT_HELPER_SCRIPT in gaussian/io.py uses `Path(__file__).parent.parent` to navigate up from `gaussian/` to project root where `gm_helper.py` lives
- gm_main.py is not modified in this plan - caller wiring is deferred to plan 06-05 to keep each plan's diff small and focused
- gaussian/__init__.py stays sparse (single docstring) until all submodules exist in later plans

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Two pre-existing test failures (test_valid_water, test_optimization_results_contain_version_info) appear intermittently due to scipy ObjSense type collision when certain test modules load in sequence. These failures are environment-sensitive and unrelated to the gaussian/ package. Both tests pass individually.

## Next Phase Readiness
- gaussian/ package is importable with core I/O functions available
- Ready for plan 06-02: ZMQ server extraction into gaussian/zmq_server.py
- gaussian/io.py imports will be wired into gm_main.py in plan 06-05

## Self-Check: PASSED

- FOUND: gaussian/__init__.py
- FOUND: gaussian/io.py
- FOUND: utils/exceptions.py
- FOUND: .planning/phases/06-extract-gaussian-i-o-zmq-server/06-01-SUMMARY.md
- FOUND commit: e430267 (Task 1)
- FOUND commit: b31e7d6 (Task 2)

---
*Phase: 06-extract-gaussian-i-o-zmq-server*
*Completed: 2026-02-20*
