---
phase: 21-nist-experimental-overlay
plan: 01
subsystem: analysis
tags: [nistchempy, jcamp, experimental-spectra, nist-webbook, ir-spectroscopy]

# Dependency graph
requires:
  - phase: 18-spectral-quality-foundations
    provides: Lorentzian broadening infrastructure and SpectrumAnalyzer plotting methods
provides:
  - ExperimentalSpectrum dataclass and fetch_experimental_spectrum() function
  - NIST WebBook fetch with JDX caching in comparison_results/{molecule}/experimental/
  - analysis_workflow.py threading of experimental data to all plot method calls
affects: [21-02 (plot overlay rendering), batch_report, html_report_generator]

# Tech tracking
tech-stack:
  added: [nistchempy, jcamp]
  patterns: [cache-first fetch with best-effort fallback, optional parameter threading]

key-files:
  created:
    - mace_gaussian/analysis/nist_fetcher.py
  modified:
    - mace_gaussian/analysis/analysis_workflow.py
    - pyproject.toml

key-decisions:
  - "Gas-phase only filter with first-match selection from NIST compound IR spectra"
  - "Cache raw JDX text (not parsed arrays) for provenance preservation"
  - "Best-effort fetch wrapped in try/except; never blocks analysis pipeline"

patterns-established:
  - "Optional experimental data threading: add param with None default to all plot methods"
  - "Cache-first NIST fetch: check JDX file exists before network call"

requirements-completed: [NIST-01]

# Metrics
duration: 3min
completed: 2026-04-01
---

# Phase 21 Plan 01: NIST Fetcher Summary

**NIST WebBook fetcher with JDX caching, transmittance-to-absorbance conversion, and experimental data threading through analysis workflow**

## Performance

- **Duration:** 3 min
- **Started:** 2026-04-01T13:06:57Z
- **Completed:** 2026-04-01T13:10:05Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Created nist_fetcher.py module with ExperimentalSpectrum dataclass, fetch_experimental_spectrum(), and _parse_jdx_file() helper
- Added nistchempy and jcamp as project dependencies in pyproject.toml
- Wired NIST fetch into run_full_analysis() with experimental data threaded through run_single_comparison, create_combined_plots, and all downstream plot method calls

## Task Commits

Each task was committed atomically:

1. **Task 1: Create nist_fetcher.py module and add dependencies** - `618aeaf` (feat)
2. **Task 2: Wire NIST fetching into analysis_workflow.py** - `13c770e` (feat)

## Files Created/Modified
- `mace_gaussian/analysis/nist_fetcher.py` - NIST fetch, cache, JCAMP-DX parse module with ExperimentalSpectrum dataclass
- `mace_gaussian/analysis/analysis_workflow.py` - Import nist_fetcher, fetch at start of run_full_analysis, thread experimental to all plot calls
- `pyproject.toml` - Added nistchempy>=1.0.0 and jcamp>=1.2.0 dependencies

## Decisions Made
- Gas-phase only: skip non-gas spectra silently (D-03)
- First gas-phase spectrum selected when multiple exist; transmittance converted to absorbance automatically
- Raw JDX cached (not parsed arrays) for provenance; parsing is ~1ms so no caching benefit
- Entire fetch wrapped in try/except returning None; analysis never blocked by NIST failures (D-04)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Known Stubs
- `experimental=experimental` is passed to `plot_spectra_comparison`, `plot_combined_spectra`, and `plot_combined_spectra_extended` but those methods do not yet accept the parameter. Plan 21-02 will add the `experimental` parameter to the plot methods and render the overlay trace.

## Next Phase Readiness
- nist_fetcher.py ready for use; ExperimentalSpectrum flows through the entire analysis workflow
- Plan 21-02 must add `experimental` parameter to SpectrumAnalyzer plot methods to actually render the trace
- Plot methods will currently ignore the unexpected `experimental` kwarg (matplotlib/Python will raise TypeError) -- 21-02 must be executed before running analysis

---
*Phase: 21-nist-experimental-overlay*
*Completed: 2026-04-01*
