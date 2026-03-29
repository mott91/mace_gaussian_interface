---
phase: 18-spectral-quality-foundations
plan: 01
subsystem: analysis
tags: [lorentzian, broadening, ir-spectra, fwhm, argparse, cli]

# Dependency graph
requires: []
provides:
  - Lorentzian broadening kernel in SpectrumAnalyzer.broaden_spectrum()
  - Configurable FWHM via --fwhm CLI flag on both analysis entry points
  - Default FWHM of 10.0 cm-1 throughout the pipeline
affects: [18-02, analysis-reports, benchmark-campaign]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Lorentzian kernel: gamma_sq / ((v - v0)^2 + gamma_sq) for unnormalized peak broadening"
    - "argparse-based CLI entry points replacing raw sys.argv parsing"

key-files:
  created:
    - tests/test_spectral_broadening.py
  modified:
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/analysis_workflow.py

key-decisions:
  - "Lorentzian uses unnormalized form (peak = intensity) not area-normalized form"
  - "raise SystemExit(1) from e used instead of sys.exit(1) per ruff B904"

patterns-established:
  - "CLI entry points use argparse with descriptive help text and defaults"
  - "FWHM threaded through ComparisonWorkflow -> SpectrumAnalyzer via bandwidth_fwhm parameter"

requirements-completed: [SPEC-01]

# Metrics
duration: 4min
completed: 2026-03-29
---

# Phase 18 Plan 01: Lorentzian Broadening Summary

**Lorentzian IR line shapes replace Gaussian KDE with configurable FWHM (default 10 cm-1) wired through CLI --fwhm flag**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-29T16:35:01Z
- **Completed:** 2026-03-29T16:39:24Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Replaced Gaussian KDE broadening with physically correct Lorentzian line shapes in all IR spectrum plots
- Changed default FWHM from 8.0 to 10.0 cm-1 per computational spectroscopy standards
- Added --fwhm CLI flag to both run_analysis.py and run_analysis_harmonic.py via argparse
- Threaded bandwidth_fwhm parameter through ComparisonWorkflow and analyze_molecule functions

## Task Commits

Each task was committed atomically:

1. **Task 1: Create test scaffold and implement Lorentzian broadening** - `143a66c` (test: RED), `d50eb72` (feat: GREEN+REFACTOR)
2. **Task 2: Wire FWHM through CLI and ComparisonWorkflow** - `486160f` (feat)

## Files Created/Modified
- `tests/test_spectral_broadening.py` - Unit tests for Lorentzian peak height, FWHM verification, and default FWHM
- `mace_gaussian/analysis/analyze_spectra.py` - Lorentzian kernel in broaden_spectrum(), default FWHM 10.0, removed Gaussian broad_param
- `mace_gaussian/analysis/analysis_workflow.py` - bandwidth_fwhm parameter in ComparisonWorkflow, argparse --fwhm in both main() entry points

## Decisions Made
- Used unnormalized Lorentzian form (peak equals intensity) rather than area-normalized form, matching standard IR spectral simulation practice
- Replaced sys.exit(1) with raise SystemExit(1) from e per ruff B904 lint rule for proper exception chaining

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed stale bandwidth_fwhm=8.0 in main() example function**
- **Found during:** Task 2 (CLI wiring)
- **Issue:** analyze_spectra.py main() example function still had hardcoded bandwidth_fwhm=8.0
- **Fix:** Changed to bandwidth_fwhm=10.0 to match new default
- **Files modified:** mace_gaussian/analysis/analyze_spectra.py
- **Verification:** grep confirms no bandwidth_fwhm=8.0 remains in mace_gaussian/
- **Committed in:** 486160f (Task 2 commit)

**2. [Rule 1 - Bug] Fixed ruff B904 lint errors in new exception handlers**
- **Found during:** Task 2 (CLI wiring)
- **Issue:** raise SystemExit(1) inside except blocks needed from e for proper chaining
- **Fix:** Changed to raise SystemExit(1) from e in all 4 exception handlers
- **Files modified:** mace_gaussian/analysis/analysis_workflow.py
- **Verification:** ruff check passes clean
- **Committed in:** 486160f (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered
- Pre-existing ty type check errors in analyze_spectra.py (ComparisonMetrics floating[Any], plt.tight_layout tuple) and analysis_workflow.py (load_results path type) are unrelated to broadening changes and were left as-is per scope boundary rules.

## User Setup Required
None - no external service configuration required.

## Known Stubs
None - all functionality is fully wired.

## Next Phase Readiness
- Lorentzian broadening is live and configurable
- Plan 18-02 can build on this for zero-intensity filtering and further spectral quality improvements
- All existing analysis scripts will produce Lorentzian-broadened spectra on next run

---
*Phase: 18-spectral-quality-foundations*
*Completed: 2026-03-29*
