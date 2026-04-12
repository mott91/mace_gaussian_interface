---
phase: 23-anharmonic-pipeline-report-overhaul
plan: 02
subsystem: analysis
tags: [plotly, css, report, presentation]
dependency_graph:
  requires: []
  provides: [plotly-builders, shared-css]
  affects: [mace_gaussian/analysis/plotly_builders.py, mace_gaussian/analysis/_shared_css.py]
tech_stack:
  added: []
  patterns: [pure-function-builders, single-source-css]
key_files:
  created:
    - mace_gaussian/analysis/plotly_builders.py
    - mace_gaussian/analysis/_shared_css.py
    - tests/test_shared_css.py
  modified:
    - tests/test_plotly_builders.py
decisions:
  - Plotly builders are pure functions (no I/O, no HTML) returning go.Figure objects
  - _shared_css.py provides build_css() consumed by both html_report_generator and batch_report
  - Changed experimental_on_grid type hint from object to Any to satisfy ty checker
metrics:
  duration: 300s
  completed: "2026-04-12T06:25:00Z"
  tasks_completed: 2
  tasks_total: 2
deviations:
  - "Rule 3: Created test files locally since 23-01 (Wave 0 test scaffold) runs in parallel worktree"
  - "Rule 1: Changed experimental_on_grid type hint from object to Any for ty compatibility"
self_check: PASSED
---

## What was built

Two presentation primitives that the rest of Phase 23 depends on:

1. **plotly_builders.py** (286 lines) — 4 pure-function Plotly figure builders:
   - `build_spectrum_figure()` — overlay ML vs DFT spectra with optional experimental
   - `build_regression_figure()` — frequency correlation with R²/RMSE annotations
   - `build_heatmap_figure()` — mode overlap matrix visualization
   - `build_timing_figure()` — bar chart comparing ML vs DFT wall-clock times

2. **_shared_css.py** — `build_css()` function producing unified CSS consumed by both
   `html_report_generator.py` and `batch_report.py` (single source of truth).

## Verification

- 11/11 tests pass (8 plotly builder tests + 3 shared CSS tests)
- ruff clean, ty clean on all new files
