---
phase: 19-degenerate-mode-handling
plan: 02
subsystem: analysis
tags: [heatmap-collapse, regression-markers, degenerate-groups, html-report, visualization]

# Dependency graph
requires:
  - 19-01 (DegenerateGroup, DegenerateGroupResult, build_degenerate_result, compute_subspace_overlap)
provides:
  - collapse_alignment_matrix() for heatmap collapsing
  - x_labels/y_labels override in plot_mode_overlap_heatmap()
  - Degenerate-aware regression plots with diamond markers
  - HTML report degenerate group summary table
  - Group-aware extract_mode_mapping returning 3-tuple
affects: [analysis-pipeline, html-reports, heatmaps, regression-plots]

# Tech tracking
tech-stack:
  added: []
  patterns: [collapsed-heatmap-matrix, diamond-marker-degenerate-groups, group-data-passthrough]

key-files:
  created: []
  modified:
    - mace_gaussian/analysis/mode_matching.py
    - mace_gaussian/analysis/analysis_workflow.py
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/html_report_generator.py

key-decisions:
  - "collapse_alignment_matrix uses max overlap for cross cells (group vs non-group)"
  - "Diamond markers (#D08770 warm orange) for degenerate group averages in regression plots"
  - "Degenerate group summary inserted before heatmap images in HTML report"

patterns-established:
  - "extract_mode_mapping returns 3-tuple (mapping, overlaps, DegenerateGroupResult) -- all downstream consumers updated"
  - "collapse_alignment_matrix produces custom labels for collapsed heatmap rows/columns"

requirements-completed: [MODE-05, MODE-06]

# Metrics
duration: 7min
completed: 2026-03-31
---

# Phase 19 Plan 02: Degenerate Mode Visualization Pipeline Summary

**Collapsed heatmaps, diamond-marker regression plots, and HTML report group labels wired through analysis pipeline**

## Performance

- **Duration:** 7 min 32 sec
- **Started:** 2026-03-31T13:26:48Z
- **Completed:** 2026-03-31T13:34:20Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- collapse_alignment_matrix() merges degenerate groups into single rows/columns with subspace overlap values and symmetry labels (e.g. "T @1356")
- plot_mode_overlap_heatmap() accepts optional x_labels/y_labels overrides for collapsed heatmaps
- extract_mode_mapping() now returns (mode_mapping, mode_overlaps, DegenerateGroupResult) 3-tuple
- plot_regression() uses diamond markers ("D") with warm orange (#D08770) for degenerate group averaged points
- HTMLReportGenerator shows degenerate group summary table with multiplicity and subspace overlap before heatmap images
- generate_mode_overlap_heatmaps() automatically collapses heatmap when degenerate groups are detected
- All 268 tests pass (excluding 4 pre-existing unrelated failures: 1 parser intensity mismatch, 3 SLURM config)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add collapse_alignment_matrix and wire degenerate-aware heatmaps/extract_mode_mapping** - `b94e6dc` (feat)
2. **Task 2: Degenerate-aware regression plots and HTML report labels** - `7571160` (feat)

## Files Created/Modified
- `mace_gaussian/analysis/mode_matching.py` - Added collapse_alignment_matrix(), x_labels/y_labels params to plot_mode_overlap_heatmap()
- `mace_gaussian/analysis/analysis_workflow.py` - Updated extract_mode_mapping() to 3-tuple, wired collapsed heatmaps, passed deg_result to regression and HTML report
- `mace_gaussian/analysis/analyze_spectra.py` - Added deg_result param to plot_regression(), diamond markers for degenerate groups
- `mace_gaussian/analysis/html_report_generator.py` - Added degenerate_groups param, group summary table in mode overlap section

## Decisions Made
- Cross cells in collapsed matrix use max overlap from the sub-block (conservative, shows strongest signal)
- Diamond markers use warm orange (#D08770) from Nord palette -- colorblind-safe contrast with existing blue (#5E81AC) circles
- Degenerate groups deduped across ML comparisons (same DFT ref groups) when building report data

## Deviations from Plan

None - plan executed exactly as written.

## Known Stubs

None - all functionality is fully wired.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 19 (degenerate mode handling) is complete
- Molecules with degenerate modes (e.g., CH4, BH3-NH3) will show collapsed heatmaps and diamond markers
- Molecules without degenerate modes (e.g., water) render identically to before

---
*Phase: 19-degenerate-mode-handling*
*Completed: 2026-03-31*
