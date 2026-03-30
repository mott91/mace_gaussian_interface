---
phase: 18-spectral-quality-foundations
plan: 02
subsystem: analysis
tags: [intensity-filtering, lorentzian, html-report, regression-metrics]

# Dependency graph
requires:
  - phase: 18-01
    provides: Lorentzian broadening, SpectrumAnalyzer with bandwidth_fwhm parameter
provides:
  - Intensity filtering in ComparisonMetrics (modes with DFT < 0.1 km/mol excluded)
  - num_intensity_filtered field for tracking filtered mode count
  - Dynamic Lorentzian methodology text in HTML reports
  - FWHM threading from ComparisonWorkflow to HTMLReportGenerator
affects: [analysis-reports, benchmark-campaign]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Intensity threshold filtering with boolean mask (int_mask = dft_int >= 0.1)"
    - "Dynamic methodology text in HTML report via self.bandwidth_fwhm"

key-files:
  created:
    - tests/test_spectral_broadening.py
  modified:
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/html_report_generator.py
    - mace_gaussian/analysis/analysis_workflow.py

key-decisions:
  - "Used 0.1 km/mol as fixed intensity threshold (SPEC-02 requirement)"
  - "Frequency metrics remain completely unfiltered (D-04 decision)"
  - "SPEC-03 (stick spectrum) deferred per user decision D-01"
  - "PIPE-02 (mace_polar dipole) deferred per user request to standalone todo"

patterns-established:
  - "Intensity filtering via boolean mask before regression: int_mask = dft_int >= THRESHOLD"
  - "Edge case guard: len(filtered) > 1 before linregress, == 1 returns MAE only, == 0 returns zeros"

requirements-completed: [SPEC-02]

# Metrics
duration: 5min
completed: 2026-03-29
---

# Phase 18 Plan 02: Intensity Filtering and Report Methodology Summary

**Zero-intensity mode filtering for intensity regression (< 0.1 km/mol threshold) with dynamic Lorentzian methodology text in HTML reports**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-29T16:41:38Z
- **Completed:** 2026-03-29T16:47:26Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Intensity regression metrics (MAE, R2) now exclude modes where DFT intensity < 0.1 km/mol, removing noise from near-zero modes
- Frequency metrics remain unfiltered -- all modes contribute to MAE, RMSE, R2 for frequencies
- HTML report methodology text dynamically reflects Lorentzian broadening with the actual FWHM value
- Intensity filtering methodology documented in report text
- 3 TDD tests verify filtering behavior, frequency independence, and edge case handling

## Task Commits

Each task was committed atomically:

1. **Task 1: Add intensity filtering to calculate_metrics with tests** - `08b26a9` (feat)
2. **Task 2: Update HTML report methodology text and thread FWHM** - `fb24cba` (feat)

## Files Created/Modified
- `tests/test_spectral_broadening.py` - 3 tests for intensity filtering (threshold, freq unfiltered, all-below edge case)
- `mace_gaussian/analysis/analyze_spectra.py` - num_intensity_filtered field on ComparisonMetrics, int_mask filtering in calculate_metrics
- `mace_gaussian/analysis/html_report_generator.py` - bandwidth_fwhm parameter, dynamic Lorentzian text, intensity filter note
- `mace_gaussian/analysis/analysis_workflow.py` - self.bandwidth_fwhm stored and threaded to report generator

## Decisions Made
- Used the actual method name `calculate_metrics` (plan referenced `compare_spectra` which does not exist) -- adapted tests accordingly
- Test data adjusted: ML frequencies use +5.0 offset consistently to get exact MAE=5.0 (plan had inconsistent +10.0 for F2)
- Pre-existing type errors in analysis_workflow.py (Path vs str in load_results) logged as out-of-scope, not fixed

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test data inconsistency**
- **Found during:** Task 1 (GREEN phase)
- **Issue:** Plan's test data had F2 frequencies differing by 10.0 (1000 vs 1010) while others differed by 5.0, making MAE=6.0 not 5.0
- **Fix:** Changed ML F2 frequency from 1010.0 to 1005.0 so all pairs differ by exactly 5.0
- **Files modified:** tests/test_spectral_broadening.py
- **Verification:** mae_freq == pytest.approx(5.0) passes
- **Committed in:** 08b26a9

**2. [Rule 3 - Blocking] Adapted to actual method name**
- **Found during:** Task 1 (test writing)
- **Issue:** Plan referenced `compare_spectra()` method but actual method is `calculate_metrics()` with (ml, dft) argument order
- **Fix:** Tests call `analyzer.calculate_metrics(ml, dft)` instead of `analyzer.compare_spectra(dft, ml)`
- **Files modified:** tests/test_spectral_broadening.py
- **Verification:** All 3 tests pass
- **Committed in:** 08b26a9

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both fixes necessary for test correctness. No scope creep.

## Deferred Requirements

- **SPEC-03 (Stick spectrum overlay):** Deferred per user decision D-01. Will be addressed in Phase 23 (report overhaul) or later.
- **PIPE-02 (mace_polar dipole calculator):** Deferred per user request. Tracked in standalone todo at `.planning/todos/pending/2026-03-26-reevaluate-mace-polar-as-dipole-calculator.md`.

## Issues Encountered
- Pre-existing ty check errors in analysis_workflow.py (Path passed where str expected in load_results, mode matching type mismatch) -- these are not caused by Plan 02 changes and are out of scope

## User Setup Required
None - no external service configuration required.

## Known Stubs
None - all functionality is fully wired.

## Next Phase Readiness
- Intensity filtering and report methodology complete
- Phase 18 fully done (both plans executed)
- Ready for Phase 19 or subsequent analysis quality phases

## Self-Check: PASSED

All files exist, all commits verified, all content markers present.

---
*Phase: 18-spectral-quality-foundations*
*Completed: 2026-03-29*
