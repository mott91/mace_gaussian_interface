---
phase: 19-degenerate-mode-handling
plan: 01
subsystem: analysis
tags: [numpy, linear-algebra, subspace-overlap, mode-matching, degenerate-modes]

# Dependency graph
requires: []
provides:
  - DegenerateGroup and DegenerateGroupResult dataclasses for mode matching
  - detect_degenerate_groups() function for frequency-proximity clustering
  - compute_subspace_overlap() using trace(M^T M)/k formula
  - build_degenerate_result() convenience function
  - Block-level deduplication in Gaussian parser preserving degenerate modes
affects: [19-02-visualization, analysis-workflow, regression-plots, heatmaps]

# Tech tracking
tech-stack:
  added: []
  patterns: [subspace-overlap-via-trace, block-level-deduplication, group-aware-statistics]

key-files:
  created: []
  modified:
    - mace_gaussian/analysis/mode_matching.py
    - mace_gaussian/gaussian/parser.py
    - tests/test_mode_matching.py
    - tests/test_gaussian_parser.py

key-decisions:
  - "Normalize subspace overlap by min(calc, ref) count per Pitfall 2 -- handles incomplete Hungarian matches"
  - "Block-level deduplication replaces per-frequency deduplication in parser"

patterns-established:
  - "DegenerateGroupResult wraps Hungarian matches with group awareness -- downstream consumers call .statistics() and .regression_data()"
  - "Symmetry label heuristic: 2-fold=E, 3-fold=T, else=deg-{k}"

requirements-completed: [MODE-05, MODE-06]

# Metrics
duration: 4min
completed: 2026-03-31
---

# Phase 19 Plan 01: Degenerate Mode Handling Core Summary

**Subspace overlap (trace M^T M / k) for degenerate mode groups with group-aware statistics and parser fix preserving all CH4 modes**

## Performance

- **Duration:** 4 min 28 sec
- **Started:** 2026-03-31T09:04:10Z
- **Completed:** 2026-03-31T09:08:38Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- DegenerateGroup and DegenerateGroupResult dataclasses with effective_overlap(), statistics(), and regression_data() methods
- detect_degenerate_groups() clusters DFT reference frequencies within 5 cm-1 threshold
- compute_subspace_overlap() correctly handles rotated degenerate pairs (individual dot products ~0.707, subspace overlap = 1.0)
- Gaussian parser now preserves all 9 CH4 harmonic modes instead of collapsing to 4
- 17 new tests across 4 test classes, all passing alongside 14 existing tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Degenerate group detection, subspace overlap, and group-aware result dataclasses with tests** - `4a199a1` (feat)
2. **Task 2: Fix Gaussian parser deduplication to preserve degenerate modes** - `06f15e7` (fix)

## Files Created/Modified
- `mace_gaussian/analysis/mode_matching.py` - Added DegenerateGroup, DegenerateGroupResult, detect_degenerate_groups(), compute_subspace_overlap(), build_degenerate_result()
- `mace_gaussian/gaussian/parser.py` - Replaced seen_freqs per-frequency deduplication with block-level deduplication using seen_block_keys
- `tests/test_mode_matching.py` - Added TestDegenerateDetection (5 tests), TestSubspaceOverlap (5 tests), TestGroupAwareStatistics (2 tests), TestGroupRegressionData (2 tests)
- `tests/test_gaussian_parser.py` - Updated CH4 test expectations from 4 to 9 harmonic modes

## Decisions Made
- Subspace overlap normalized by min(calc_indices, ref_indices) length to handle partial Hungarian matches (Pitfall 2)
- Block-level deduplication strategy: fingerprint each frequency block as a tuple of (freq, intensity) pairs, deduplicate identical blocks, preserve all modes within unique blocks

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Pre-existing test failure in TestParseAnharmonicFrequencies::test_water_anharmonic_count_and_values (intensity mismatch 0.6916 vs expected 0.6283) -- confirmed pre-existing, unrelated to this plan's changes. Out of scope.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- DegenerateGroupResult ready for Plan 02 to consume for visualization (heatmap collapsing, regression markers, report labels)
- build_degenerate_result() provides the entry point for analysis_workflow.py to wire into

---
*Phase: 19-degenerate-mode-handling*
*Completed: 2026-03-31*
