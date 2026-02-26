---
phase: 10-documentation
plan: 02
subsystem: documentation
tags: [quickstart, worked-example, water, IR-spectroscopy, MACE-Gaussian]

# Dependency graph
requires:
  - phase: 09-ci-distribution
    provides: README Installation section with explicit 5-step procedure
provides:
  - docs/quickstart.md — end-to-end guide clone to HTML report
  - docs/examples/water/ — self-contained worked example with committed JSON + PNG reference output
affects: [README, new-user-onboarding]

# Tech tracking
tech-stack:
  added: []
  patterns: [committed-expected-output, portable-reference-files-JSON-PNG-only]

key-files:
  created:
    - docs/quickstart.md
    - docs/examples/water/README.md
    - docs/examples/water/water.xyz
    - docs/examples/water/expected_output/geometry_opt_results.json
    - docs/examples/water/expected_output/mace_omol_mace_ml_results.json
    - docs/examples/water/expected_output/analysis_metrics_summary.json
    - docs/examples/water/expected_plots/spectrum_combined.png
    - docs/examples/water/expected_plots/regression_combined.png
  modified: []

key-decisions:
  - "quickstart.md references README#installation rather than duplicating install steps to avoid drift"
  - "expected_output/ uses JSON only (no .chk/.fchk/.log) — portable, system-independent, human-readable"
  - "Plots sourced from anharmonic analysis_results/ (not harmonic) to match the quickstart workflow path"

patterns-established:
  - "Worked examples: committed expected_output/ JSON + expected_plots/ PNG, no binary Gaussian files"
  - "Quickstart links to README#installation; does not duplicate install steps"

requirements-completed: [DOC-02, DOC-03]

# Metrics
duration: 2min
completed: 2026-02-26
---

# Phase 10 Plan 02: Quickstart and Water Worked Example Summary

**docs/quickstart.md (clone-to-report guide) and docs/examples/water/ (committed JSON + PNG reference for water H2O IR spectroscopy)**

## Performance

- **Duration:** ~2 min
- **Started:** 2026-02-26T15:55:51Z
- **Completed:** 2026-02-26T15:57:34Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- Created docs/quickstart.md covering all 5 steps (install, run, check, analyze, view report) with reference to README#installation instead of duplicating the install procedure
- Created docs/examples/water/ self-contained worked example with README, input file, 3 JSON reference outputs, and 2 PNG plots
- Committed expected output (JSON only, no binary Gaussian files) enabling inspection without running Gaussian 16

## Task Commits

Each task was committed atomically:

1. **Task 1: Create docs/quickstart.md** - `e4ad4d6` (feat)
2. **Task 2: Create docs/examples/water/ worked example** - `cbaebfd` (feat)

**Plan metadata:** (see final commit below)

## Files Created/Modified
- `docs/quickstart.md` - End-to-end quickstart: install, run water.xyz, check results, analyze, view HTML report
- `docs/examples/water/README.md` - Self-contained worked example with Expected Output and To Reproduce sections
- `docs/examples/water/water.xyz` - 3-atom water geometry (O + 2H)
- `docs/examples/water/expected_output/geometry_opt_results.json` - mace_omol geometry opt (8 steps, ~4.2s, converged)
- `docs/examples/water/expected_output/mace_omol_mace_ml_results.json` - 3 harmonic + anharmonic frequencies for best ML combination
- `docs/examples/water/expected_output/analysis_metrics_summary.json` - R² and RMSE for all 4 combinations (best: R²=0.9997, RMSE=44 cm-1)
- `docs/examples/water/expected_plots/spectrum_combined.png` - IR spectrum comparison (all 4 ML combinations vs DFT)
- `docs/examples/water/expected_plots/regression_combined.png` - Frequency regression plot

## Decisions Made
- quickstart.md references README#installation rather than duplicating the 5-step install procedure — avoids drift between quickstart and README
- expected_output/ uses JSON files only (no .chk, .fchk, .log) — portable, system-independent, human-readable reference outputs
- Plots sourced from anharmonic analysis_results/water/plots/ to match the quickstart workflow (run_analysis.py path)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- docs/quickstart.md and docs/examples/water/ are complete
- New users can follow quickstart.md to reproduce the water calculation or inspect the committed reference output without running Gaussian 16
- Phase 10 documentation continues with remaining plans (architecture diagrams, contribution guide, etc.)

---
*Phase: 10-documentation*
*Completed: 2026-02-26*
