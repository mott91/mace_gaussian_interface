---
phase: 11-integration-wiring-fixes
plan: 01
subsystem: utils
tags: [python39, type-annotations, warnings, validation, cli, wiring]

# Dependency graph
requires:
  - phase: 02-error-handling-and-validation
    provides: exceptions.py, validation.py with validate_all_prerequisites
  - phase: 07-workflow-extraction
    provides: workflow.py run_frequency_calculation() calling save_frequency_results()
  - phase: 08-package-structure
    provides: mace_gaussian package with cli.py at module boundary

provides:
  - GaussianRunError and GaussianTimeoutError importable from mace_gaussian.utils
  - CUDANotAvailableWarning emitted via warnings.warn in detect_device() when CUDA unavailable
  - calculation_parameters dict wired into save_frequency_results() call in workflow.py
  - cli.py resolves MACE_DIPOLE_MODEL_PATH and MACE_HELPER_SCRIPT_PATH env vars before validation
  - results.py Python 3.9 compatible via from __future__ import annotations

affects: [testing, cli, workflow, audit-closure]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "warnings.warn with stacklevel=2 so warning origin points to caller of detect_device()"
    - "Module-level imports in cli.py for patchability in tests (not inside function body)"
    - "from __future__ import annotations for Python 3.9 compat with dict[str, Any] signatures"
    - "X | None instead of Optional[X] in conjunction with __future__ annotations"

key-files:
  created:
    - tests/test_validation.py (two new test classes: TestDetectDeviceWarning, TestCliEnvVarResolution)
  modified:
    - mace_gaussian/utils/results.py
    - mace_gaussian/utils/__init__.py
    - mace_gaussian/utils/validation.py
    - mace_gaussian/workflow.py
    - mace_gaussian/cli.py
    - tests/test_validation.py

key-decisions:
  - "Move validation imports to module level in cli.py — required so mace_gaussian.cli.validate_all_prerequisites is patchable in tests (local imports inside function body don't create module-level attributes)"
  - "Use # noqa: SIM117 on nested with statements in test — parenthesized with syntax (Python 3.10+) would break system Python 3.8 test runner"
  - "Convert Optional[X] to X | None in results.py — UP045 ruff violations triggered by adding from __future__ import annotations (Rule 1 auto-fix)"
  - "warnings.warn with stacklevel=2 in detect_device() — warning appears to originate from caller not from inside validation.py"
  - "Inline env var defaults in cli.py match calculator module defaults — not imported to avoid DGL/espaloma side effects at cli import time"

requirements-completed: [REPR-02, ERR-03, ERR-04]

# Metrics
duration: 4min
completed: 2026-02-27
---

# Phase 11 Plan 01: Integration Wiring Fixes Summary

**Four surgical call-site wiring fixes: GaussianRunError/GaussianTimeoutError exports, CUDANotAvailableWarning emission, calculation_parameters dict, and env var path resolution in CLI**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-27T10:34:34Z
- **Completed:** 2026-02-27T10:39:10Z
- **Tasks:** 3
- **Files modified:** 6

## Accomplishments
- GaussianRunError and GaussianTimeoutError now importable directly from `mace_gaussian.utils`
- `detect_device()` emits `CUDANotAvailableWarning` via `warnings.warn(stacklevel=2)` when CUDA unavailable, catchable with `pytest.warns` and `warnings.catch_warnings`
- `run_frequency_calculation()` in workflow.py passes `calculation_parameters={"energy_calculator": ..., "dipole_calculator": ...}` to `save_frequency_results()` — results.json always has non-empty calculator entries
- `cli.py run()` resolves `MACE_DIPOLE_MODEL_PATH` and `MACE_HELPER_SCRIPT_PATH` from env vars before calling `validate_all_prerequisites()` (was passing None silently, skipping checks)
- `results.py` gains `from __future__ import annotations` for Python 3.9 compatibility with `dict[str, Any]` signatures
- 2 new tests added for ERR-03 and ERR-04; all 131 tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix Python 3.9 compat in results.py and add GaussianRunError/GaussianTimeoutError exports** - `a7b7118` (fix)
2. **Task 2: Emit CUDANotAvailableWarning in detect_device() and wire calculation_parameters in workflow.py** - `b481d84` (fix)
3. **Task 3: Resolve env var paths in cli.py and add tests for ERR-03/ERR-04** - `7a54aba` (fix)

## Files Created/Modified
- `mace_gaussian/utils/results.py` - Added `from __future__ import annotations`, converted `Optional[X]` to `X | None`
- `mace_gaussian/utils/__init__.py` - Added `GaussianRunError` and `GaussianTimeoutError` to imports and `__all__`
- `mace_gaussian/utils/validation.py` - Added `import warnings` and `CUDANotAvailableWarning` import; emit `warnings.warn` in both CUDA-unavailable branches of `detect_device()`
- `mace_gaussian/workflow.py` - Added `calculation_parameters` dict and passed it to `save_frequency_results()`
- `mace_gaussian/cli.py` - Moved validation imports to module level; added `import os`; added env var resolution for dipole model path and helper script path
- `tests/test_validation.py` - Added `import warnings`; added `TestDetectDeviceWarning` and `TestCliEnvVarResolution` classes

## Decisions Made
- Moved validation imports to module level in `cli.py` so that `mace_gaussian.cli.validate_all_prerequisites` is patchable in tests. Imports inside a function body are local names and don't create attributes on the module, making `patch("mace_gaussian.cli.validate_all_prerequisites")` fail.
- Used `# noqa: SIM117` on nested `with` statements in tests. The parenthesized multi-context `with (...)` syntax is Python 3.10+ and would cause SyntaxError on the system Python 3.8 test runner. CI runs Python 3.10 where it is valid, but local is Python 3.8.
- Inline env var defaults in `cli.py` (not imported from calculator modules) to avoid triggering DGL/espaloma heavy dependency loads at CLI import time.
- Converted `Optional[X]` to `X | None` in results.py as a Rule 1 auto-fix — adding `from __future__ import annotations` caused ruff UP045 to flag the pre-existing `Optional` annotations.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Convert Optional[X] to X | None in results.py**
- **Found during:** Task 1 (Fix Python 3.9 compat in results.py)
- **Issue:** Adding `from __future__ import annotations` triggered ruff UP045 on pre-existing `Optional[...]` annotations in the same file — ruff now flags them as requiring conversion to `X | None` syntax
- **Fix:** Converted all 6 `Optional[...]` occurrences in method signatures to `X | None` and removed `Optional` from the `from typing import` line
- **Files modified:** `mace_gaussian/utils/results.py`
- **Verification:** `ruff check` passes, all tests pass
- **Committed in:** `7a54aba` (Task 3 commit, grouped with other ruff-clean work)

**2. [Rule 3 - Blocking] Move validation imports to module level in cli.py**
- **Found during:** Task 3 (test for ERR-04 failing with AttributeError)
- **Issue:** Test patches `mace_gaussian.cli.validate_all_prerequisites` but that attribute didn't exist in the module namespace — it was imported inside the `run()` function body as a local name
- **Fix:** Moved `from mace_gaussian.utils.validation import (...)` and `from mace_gaussian.utils.exceptions import (...)` to module level in cli.py; removed the duplicate local imports inside `run()`
- **Files modified:** `mace_gaussian/cli.py`
- **Verification:** Test passes, ruff clean, all 131 tests pass
- **Committed in:** `7a54aba` (part of Task 3 commit)

---

**Total deviations:** 2 auto-fixed (1 Rule 1 — ruff violation from annotation style change; 1 Rule 3 — blocking test failure from unpatchable local import)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep.

## Issues Encountered
- Wrong Python environment tried first (`micromamba run -n mace4ir_v2` has no pytest; system Python 3.8 has DGL incompatibility). Used `.venv/bin/python` (Python 3.12) which has all test dependencies installed. This is the correct environment for the test suite.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 4 audit gaps closed: REPR-02, ERR-03, ERR-04, and the GaussianRunError/GaussianTimeoutError export gap
- results.json now always contains non-empty `calculation_parameters` for ML runs
- CLI now validates dipole model and helper script at startup (not silently skipped)
- `CUDANotAvailableWarning` is properly a catchable Python warning
- 131 tests passing, ruff clean across all 5 modified files

---
*Phase: 11-integration-wiring-fixes*
*Completed: 2026-02-27*
