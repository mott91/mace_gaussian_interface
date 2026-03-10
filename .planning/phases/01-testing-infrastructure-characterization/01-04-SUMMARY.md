---
phase: 01-testing-infrastructure-characterization
plan: 04
subsystem: testing
tags: [pytest, regression-tests, results-json, markers, coverage, water, CH4]

# Dependency graph
requires:
  - phase: 01-01
    provides: "Test fixtures (water/CH4 results.json), conftest.py with fixture accessors, pytest markers, pytest-cov config"
provides:
  - "Regression tests validating water results.json structure and frequency ranges"
  - "Regression tests validating CH4 results.json structure and frequency ranges"
  - "Marker infrastructure verification (gpu/gaussian placeholders)"
  - "Coverage reporting end-to-end verification"
affects: [02-parser-hardening]

# Tech tracking
tech-stack:
  added: []
  patterns: [class-based test grouping by molecule, JSON fixture validation via pathlib.read_text]

key-files:
  created:
    - tests/test_regression.py
  modified: []

key-decisions:
  - "Used class-based test organization (TestWaterReferenceOutput, TestCH4ReferenceOutput, TestReferenceOutputsCommitted) for clear grouping"
  - "Validated physically reasonable frequency ranges rather than exact values to avoid brittleness if ML models are updated"
  - "CH4 harmonic entries are 4 (collapsed degenerate modes) not 9 -- matched to actual JSON structure"

patterns-established:
  - "Reference output validation: check structure (keys, counts) then values (ranges, physical constraints)"
  - "Marker placeholders: pytest.skip() inside marked tests to verify marker infrastructure without real integration tests"

# Metrics
duration: 2min
completed: 2026-02-16
---

# Phase 1 Plan 4: Reference Output Regression Tests Summary

**Regression tests for water/CH4 results.json with marker verification and coverage reporting validation (TEST-05, TEST-06, TEST-07, TEST-09)**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-16T14:05:47Z
- **Completed:** 2026-02-16T14:07:52Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- 9 regression tests covering water and CH4 reference output validation
- Water tests verify 3 harmonic modes, 3 anharmonic, 3 overtones, nonzero dipole, frequencies in 400-4500 cm-1
- CH4 tests verify 4 unique harmonic entries, 9 anharmonic, bending ~1100-1500 cm-1, stretching ~2900-3300 cm-1, zero dipole
- Marker infrastructure verified: `pytest -m "not gpu and not gaussian"` deselects 2 placeholder tests
- Coverage report runs end-to-end: gaussian_parser at 84% coverage across 32 total tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Write reference output regression tests** - `4db3134` (test)

## Files Created/Modified
- `tests/test_regression.py` - 9 tests: water structure/values, CH4 structure/values, fixture existence, valid JSON, unit baseline, GPU marker, Gaussian marker

## Decisions Made
- **Class-based grouping:** Organized tests into `TestWaterReferenceOutput`, `TestCH4ReferenceOutput`, `TestReferenceOutputsCommitted` for clear scoping and pytest output readability.
- **Range-based value assertions:** Used physically reasonable frequency ranges (e.g., 400-4500 cm-1 for water) rather than exact floating-point matches, preventing test brittleness if ML models are retrained.
- **CH4 has 4 harmonic entries (not 9):** The results.json collapses degenerate modes, storing 4 unique harmonic frequencies while the 9 anharmonic fundamentals list each mode individually.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- The system Python (miniconda) does not have pytest-cov installed. Coverage commands only work with the uv venv Python (`/home/mot/mace_gaussian/.venv/bin/python`). This is expected behavior since pytest-cov is configured as a dev dependency in `[dependency-groups]`.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 4 plans in Phase 01 complete (01-01 infrastructure, 01-02/03 parser tests, 01-04 regression tests)
- 32 total tests: 29 pass, 2 skip (marker placeholders), 1 xfail (acoh bug)
- Coverage reporting functional (gaussian_parser at 84%)
- Marker filtering working for unit/GPU/Gaussian test separation
- Ready for Phase 02 (parser hardening)

## Self-Check: PASSED

All 1 created file verified present on disk. Task commit (4db3134) verified in git log.

---
*Phase: 01-testing-infrastructure-characterization*
*Completed: 2026-02-16*
