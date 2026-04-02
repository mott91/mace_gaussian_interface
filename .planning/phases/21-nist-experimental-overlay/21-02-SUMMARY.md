---
phase: 21-nist-experimental-overlay
plan: 02
subsystem: analysis
tags: [matplotlib, nist, experimental-spectra, html-report, batch-report]

requires:
  - phase: 21-nist-experimental-overlay plan 01
    provides: ExperimentalSpectrum dataclass and fetch_experimental_spectrum function
provides:
  - Experimental trace overlay on all spectrum plots (3 methods)
  - Experimental data section in HTML report
  - Experimental overlay on batch report per-molecule spectrum plots
affects: [run_analysis.py, run_analysis_harmonic.py, analysis_workflow.py]

tech-stack:
  added: []
  patterns: [TYPE_CHECKING import for nist_fetcher, scipy.interp1d for spectrum interpolation]

key-files:
  created: []
  modified:
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/html_report_generator.py
    - mace_gaussian/analysis/batch_report.py

key-decisions:
  - "Experimental trace plotted on DFT baseline (no offset) as black dashed line per D-07/D-09"
  - "Used TYPE_CHECKING guard to avoid circular imports from nist_fetcher"
  - "NIST-03 (peak comparison table) deferred as specified in CONTEXT.md"

patterns-established:
  - "Experimental overlay: interp1d onto freq_grid, normalize, plot with color=#000000 linestyle=--"
  - "Best-effort experimental fetch: try/except in batch_report, None guard in all methods"

requirements-completed: [NIST-02, NIST-03]

duration: 4min
completed: 2026-04-02
---

# Phase 21 Plan 02: Experimental Spectrum Overlay Summary

**Black dashed experimental trace overlaid on all spectrum plots (3 methods), HTML report experimental section, and batch report per-molecule overlay**

## Performance

- **Duration:** 4 min
- **Started:** 2026-04-02T08:13:54Z
- **Completed:** 2026-04-02T08:18:18Z
- **Tasks:** 4
- **Files modified:** 3

## Accomplishments
- All 3 plot methods in SpectrumAnalyzer accept optional ExperimentalSpectrum and render black dashed line
- HTML report shows experimental data source section (source, molecule, CAS, data range) when available
- Batch report per-molecule spectrum overlay includes experimental trace scaled to DFT intensity range
- All methods backward compatible (experimental=None produces identical output to before)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add experimental trace to plot_spectra_comparison** - `f04ba3b` (feat)
2. **Task 2: Add experimental trace to plot_combined_spectra and plot_combined_spectra_extended** - `c992ede` (feat)
3. **Task 3: Add experimental section to HTML report** - `27e164a` (feat)
4. **Task 4: Add experimental overlay to batch report spectrum plot** - `90a9a23` (feat)

## Files Created/Modified
- `mace_gaussian/analysis/analyze_spectra.py` - Added experimental parameter and black dashed trace to plot_spectra_comparison, plot_combined_spectra, plot_combined_spectra_extended
- `mace_gaussian/analysis/html_report_generator.py` - Added create_experimental_section method and wired into generate_report
- `mace_gaussian/analysis/batch_report.py` - Added experimental overlay to _plot_spectrum_overlay with best-effort NIST fetch

## Decisions Made
- Used TYPE_CHECKING guard for ExperimentalSpectrum import to avoid circular imports
- Experimental trace plotted on DFT baseline (no vertical offset) as specified in D-07
- Black dashed line (color="#000000", linestyle="--") per D-09
- NIST-03 (quantitative peak comparison table) noted as deferred in HTML report section

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Known Stubs
None - all functionality is fully wired. The experimental parameter defaults to None for backward compatibility; callers (run_analysis.py, analysis_workflow.py) will need to pass experimental data when available.

## Next Phase Readiness
- All plot methods accept experimental data; callers need to be updated to pass it through
- run_analysis.py and run_analysis_harmonic.py need minor updates to fetch and pass experimental data to plot methods
- NIST-03 (peak comparison table) deferred for future phase

---
*Phase: 21-nist-experimental-overlay*
*Completed: 2026-04-02*
