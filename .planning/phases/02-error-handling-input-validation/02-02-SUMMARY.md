---
phase: 02-error-handling-input-validation
plan: 02
subsystem: parsing
tags: [gaussian-parser, error-handling, reproducibility, version-metadata]

# Dependency graph
requires:
  - phase: 01-testing-infrastructure
    provides: "Test fixtures, conftest.py, regression test framework"
provides:
  - "GaussianParseError on empty harmonic frequencies"
  - "strict parameter for anharmonic/overtones/combination_bands parsers"
  - "Real optimization step count from LBFGS optimizer"
  - "version_info and calculation_parameters in results JSON"
affects: [02-03, 03-gaussian-interface, 04-results-analysis]

# Tech tracking
tech-stack:
  added: []
  patterns: ["strict parameter for optional error raising", "lazy import with try/except for validation module"]

key-files:
  created:
    - exceptions.py
    - validation.py
  modified:
    - gaussian_parser.py
    - gm_main.py
    - results_manager.py
    - tests/test_gaussian_parser.py
    - tests/test_regression.py
    - pyproject.toml

key-decisions:
  - "Created exceptions.py and validation.py inline as Plan 01 prerequisite (Rule 3 - blocking)"
  - "Used lazy import (try/except ImportError) for validation module in results_manager to avoid import-time issues"
  - "Default strict=False preserves backward compatibility for anharmonic/overtones/combination_bands"

patterns-established:
  - "strict parameter pattern: default False for backward compat, True to enforce presence"
  - "Lazy validation import: try/except ImportError with fallback to empty dict"

# Metrics
duration: 4min
completed: 2026-02-17
---

# Phase 02 Plan 02: Parser Error Hardening and Results Metadata Summary

**GaussianParseError on empty harmonic frequencies, strict mode for anharmonic parsers, real optimization step tracking, and version_info metadata in results JSON**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-17T12:36:32Z
- **Completed:** 2026-02-17T12:40:25Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments
- parse_harmonic_frequencies() now raises GaussianParseError when no frequencies found (always an error for freq calc logs)
- Anharmonic, overtones, and combination bands parsers accept strict parameter (default False for backward compat)
- geometry_optimisation() returns (mol, num_steps) tuple with real step count from LBFGS optimizer
- Results JSON now includes version_info (tool_version, python_version, platform, package versions) and optional calculation_parameters
- 10 new tests (7 parser error handling + 3 results metadata) all passing

## Task Commits

Each task was committed atomically:

1. **Task 1: Harden parser errors and fix optimization step tracking** - `730971f` (feat)
2. **Task 2: Add version metadata to results JSON and update tests** - `03265d2` (feat)

## Files Created/Modified
- `exceptions.py` - Custom exception hierarchy (MaceGaussianError, GaussianParseError, etc.)
- `validation.py` - Prerequisite checks, device detection, version metadata collection
- `gaussian_parser.py` - GaussianParseError import, empty harmonic raises, strict parameter on 3 methods
- `gm_main.py` - geometry_optimisation returns (mol, num_steps), caller unpacks and passes real count
- `results_manager.py` - version_info and calculation_parameters in both save methods
- `tests/test_gaussian_parser.py` - 7 new TestParserErrorHandling tests
- `tests/test_regression.py` - 3 new TestResultsManagerMetadata tests
- `pyproject.toml` - Added exceptions, validation, results_manager to known-first-party

## Decisions Made
- Created exceptions.py and validation.py as Plan 01 prerequisites (they did not exist yet but are required by this plan's imports)
- Used lazy import (try/except ImportError) for validation.collect_version_metadata in results_manager to gracefully degrade if validation module is unavailable
- Kept strict=False as default for all anharmonic-type parsers to maintain backward compatibility

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Created exceptions.py and validation.py (Plan 01 prerequisites)**
- **Found during:** Task 1 (parser error hardening)
- **Issue:** Plan 02 imports GaussianParseError from exceptions.py and collect_version_metadata from validation.py, but neither file existed (Plan 01 not yet executed)
- **Fix:** Created both foundation modules with full exception hierarchy and validation functions per Plan 01 spec
- **Files created:** exceptions.py, validation.py
- **Verification:** `from exceptions import GaussianParseError` and `from validation import collect_version_metadata` both succeed
- **Committed in:** 730971f (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Prerequisite creation was necessary to unblock Plan 02. No scope creep -- these files follow Plan 01 spec exactly.

## Issues Encountered
- Test for save_optimization_results initially used MagicMock for ASE Atoms, but ase.io.write requires real Atoms objects with proper positions array. Fixed by using real ASE Atoms constructor.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Parser error handling and version metadata are in place
- exceptions.py and validation.py foundation modules exist for Plan 01 (tests still needed) and Plan 03
- All 33 tests passing (30 passed, 2 skipped, 1 xfailed)

---
*Phase: 02-error-handling-input-validation*
*Completed: 2026-02-17*

## Self-Check: PASSED
- All 7 key files verified present on disk
- Both task commits (730971f, 03265d2) verified in git log
