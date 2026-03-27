---
phase: 15-slurm-integration-batch-report
plan: "02"
subsystem: analysis
tags: [batch-report, html, visualization, cli]
dependency_graph:
  requires: ["15-01"]
  provides: ["batch_report_module", "report_cli_command"]
  affects: ["mace_gaussian/analysis/", "mace_gaussian/cli.py"]
tech_stack:
  added: [seaborn-heatmap, pandas-aggregation, base64-embedding]
  patterns: [self-contained-html, tempfile-plot-generation]
key_files:
  created:
    - mace_gaussian/analysis/batch_report.py
    - tests/test_batch_report.py
  modified:
    - mace_gaussian/cli.py
decisions:
  - "CSS extracted to _build_css() helper to keep line lengths under 100"
  - "Leaderboard highlights best (green) and worst (red) combo rows"
  - "R^2 computed from sorted harmonic frequencies only (matching plan spec)"
metrics:
  duration: "5 min"
  completed: "2026-03-27"
  tasks: 2
  files_created: 2
  files_modified: 1
requirements:
  completed: [BATCH-05]
---

# Phase 15 Plan 02: Batch Accuracy Report Summary

Multi-molecule HTML batch report aggregating R^2 and RMSE across all calculator combinations, with 4 plot types and accuracy leaderboard, accessible via `mace-gaussian report`

## What Was Done

### Task 1: Create batch report module (687846e)
- Created `mace_gaussian/analysis/batch_report.py` (799 lines)
- `aggregate_results()` walks comparison_results/ directory, loads DFT reference and ML results, computes R^2 and RMSE per molecule/combo pair
- `_plot_heatmap()` generates RMSE and R^2 heatmaps (combo x molecule) with annotated values
- `_plot_boxplots()` creates RMSE distribution box plots sorted by median, with individual data points
- `_plot_size_scaling()` creates RMSE vs atom count trend lines per combo
- `_plot_spectrum_overlay()` creates per-molecule stick spectra (DFT black, ML combos colored)
- `_generate_html()` builds self-contained HTML with base64-embedded plots, nav toolbar, and leaderboard table
- `generate_batch_report()` orchestrates everything and writes batch_report.html

### Task 2: Add report CLI command and tests (dfdb921)
- Added `mace-gaussian report` command with `--results-dir` and `--output-dir` options
- Created 4 tests validating real data aggregation, HTML generation, CLI help, and empty directory handling
- All 4 tests pass against actual comparison_results/ data (water + 13 other molecules)

## Deviations from Plan

None - plan executed exactly as written.

## Verification Results

- `ruff check` passes on all modified files
- `mace-gaussian report --help` exits 0 with expected options
- All 4 pytest tests pass
- Generated HTML contains Leaderboard, heatmap, box plots, size-scaling, and spectrum overlays
- HTML is self-contained with base64-embedded PNG plots

## Known Stubs

None - all functions are fully implemented and produce real output from existing data.

## Commits

| Task | Commit | Description |
|------|--------|-------------|
| 1 | 687846e | Batch report module with aggregation and 4 plot types |
| 2 | dfdb921 | CLI report command and integration tests |

## Self-Check: PASSED

- batch_report.py: FOUND
- test_batch_report.py: FOUND
- Commit 687846e: FOUND
- Commit dfdb921: FOUND
