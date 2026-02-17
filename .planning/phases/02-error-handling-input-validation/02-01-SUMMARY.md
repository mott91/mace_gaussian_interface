---
phase: 02-error-handling-input-validation
plan: 01
subsystem: error-handling
tags: [exceptions, validation, prerequisites, device-detection, version-metadata]

requires:
  - phase: 01-testing-infrastructure
    provides: "Test infrastructure (pytest, fixtures, conftest)"
provides:
  - "MaceGaussianError exception hierarchy (4 exceptions + 1 warning)"
  - "validation.py with 8 functions for prerequisite checks, XYZ validation, device detection, version metadata"
  - "Test suite for exceptions and validation (30 tests)"
affects: [02-02, 02-03, 03-cli-robustness, 04-gaussian-interface]

tech-stack:
  added: []
  patterns: [flat-exception-hierarchy, pure-validation-functions, mocked-system-checks]

key-files:
  created:
    - exceptions.py
    - validation.py
    - tests/test_exceptions.py
    - tests/test_validation.py
  modified:
    - pyproject.toml

key-decisions:
  - "Local torch import in detect_device for resilience when PyTorch not installed"
  - "validate_xyz_file warns on >200 atoms instead of erroring (supports large molecules)"
  - "collect_version_metadata uses try/except per dependency for graceful degradation"

patterns-established:
  - "Exception pattern: flat hierarchy under MaceGaussianError, one class per failure mode"
  - "Validation pattern: pure functions that return values or raise specific exceptions"
  - "Mock pattern: patch torch.cuda directly (not validation.torch) for local imports"

duration: 4min
completed: 2026-02-17
---

# Phase 02 Plan 01: Exception Hierarchy and Validation Module Summary

**Flat exception hierarchy with 4 domain errors and 8 pure validation functions for prerequisites, XYZ files, CUDA detection, and version metadata**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-17T12:36:22Z
- **Completed:** 2026-02-17T12:40:53Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Created MaceGaussianError hierarchy: PrerequisiteError, GaussianParseError, InputValidationError, CUDANotAvailableWarning
- Created validation.py with 8 functions covering g16/formchk checks, dipole model/helper script checks, XYZ validation, device detection, version metadata collection
- Created comprehensive test suite (30 tests) covering all exception classes and validation functions with mocked system dependencies

## Task Commits

Each task was committed atomically:

1. **Task 1: Create exception hierarchy and validation module** - `730971f` (feat) -- note: bundled with prior 02-02 commit from previous session
2. **Task 2: Create tests for exceptions and validation** - `d2ad179` (test)

## Files Created/Modified
- `exceptions.py` - Custom exception hierarchy (MaceGaussianError base + 3 subclasses + 1 warning)
- `validation.py` - 8 pure validation functions for runtime environment checks
- `tests/test_exceptions.py` - 11 tests for exception hierarchy and usage
- `tests/test_validation.py` - 19 tests for all validation functions
- `pyproject.toml` - Added exceptions and validation to ruff known-first-party

## Decisions Made
- Used local `import torch` inside `detect_device()` for resilience when PyTorch is not installed
- `validate_xyz_file` logs a warning (not error) for molecules with >200 atoms, supporting legitimate large-molecule use cases
- `collect_version_metadata` wraps each dependency version lookup in try/except for graceful degradation

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed mock paths for detect_device tests**
- **Found during:** Task 2 (test creation)
- **Issue:** `@patch("validation.torch.cuda.is_available")` failed because torch is imported locally inside detect_device, so validation module has no torch attribute
- **Fix:** Changed mock targets to `@patch("torch.cuda.is_available")` and `@patch("torch.cuda.get_device_name")` which patches the actual torch module
- **Files modified:** tests/test_validation.py
- **Verification:** All 30 tests pass
- **Committed in:** d2ad179

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minimal -- standard mock path correction for local imports.

## Issues Encountered
- Task 1 files (exceptions.py, validation.py) were already committed in a previous session's `730971f` commit that bundled 02-01 prerequisites with 02-02 work. The pyproject.toml isort config was also already updated. Task 1 was verified as complete and no duplicate commit was needed.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Exception hierarchy ready for use by Plan 02 (gaussian_parser hardening) and Plan 03 (gm_main integration)
- Validation functions ready for integration into CLI and main workflow modules
- Test patterns established for mocking system dependencies

## Self-Check: PASSED

- All 4 source/test files exist on disk
- Commit 730971f found (Task 1)
- Commit d2ad179 found (Task 2)

---
*Phase: 02-error-handling-input-validation*
*Completed: 2026-02-17*
