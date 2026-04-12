---
phase: 23-anharmonic-pipeline-report-overhaul
plan: 04
subsystem: analysis/html_report_generator
tags: [report, html, plotly, integration, executive-summary]
dependency_graph:
  requires: [23-02, 23-03]
  provides: [HTMLReportGenerator-plotly-overhaul]
  affects: [analysis_workflow, run_analysis, run_analysis_harmonic]
tech_stack:
  added: []
  patterns: [plotly-emit-once, html-escape-helper, mode-flag-branching, explicit-stub-policy]
key_files:
  created: []
  modified:
    - mace_gaussian/analysis/html_report_generator.py
    - tests/conftest.py
decisions:
  - "Overtones: Option B (explicit stub placeholder) -- grep survey confirmed comparison dicts do not carry structured overtone records; data exists in SpectrumData labels but not as separate comparison['overtones'] dicts"
  - "Heatmaps remain as relative <img src> paths (not base64) since report.html lives alongside plots/ dir"
  - "Added _format_exp_agreement static helper to DRY experimental agreement formatting"
  - "Added _create_experimental_info_section to preserve NIST source metadata in reports"
  - "Added _create_timing_hardware_section detailed helper preserving Phase 20 timing table"
metrics:
  duration_s: 476
  completed: "2026-04-12T11:15:07Z"
  tasks_completed: 1
  tasks_total: 1
  files_modified: 2
---

# Phase 23 Plan 04: HTMLReportGenerator Plotly Overhaul Summary

Rewrote HTMLReportGenerator to Plotly-powered, mode-aware, executive-summary-first report with shared CSS from _shared_css and structured data export via report_data.

## What Was Done

### Task 1: Rewrite HTMLReportGenerator (b2f9a26)

Complete rewrite of `mace_gaussian/analysis/html_report_generator.py` (691 lines):

- **Plotly integration**: Replaced all matplotlib-base64 spectrum/regression plots with interactive Plotly figures via `build_spectrum_figure`, `build_regression_figure`, `build_combined_spectrum_figure` from plotly_builders
- **Emit-once CDN**: `_fig_to_div` tracks `_plotlyjs_emitted` boolean; first call includes CDN script, subsequent calls use `include_plotlyjs=False`. Defensive pre-write assertion raises RuntimeError if CDN appears > 1 time
- **Executive summary (D-01, D-02)**: `_create_executive_summary` renders verdict + method-cards grid from `rank_methods` / `build_verdict`. First card gets `best-method` CSS class
- **Mode flag (D-04)**: `__init__` accepts `mode: Literal["harmonic", "anharmonic"]`. Harmonic mode skips overtones section entirely
- **Shared CSS (D-04)**: Imports `build_css()` from `_shared_css`; no inline CSS. Byte-identical across modes (verified by test)
- **Report flow (D-03)**: Head > Header > Nav > Executive Summary > Combined Plots > Per-Method Sections > Experimental Info > Summary Table > Overtones (anharmonic only) > Footer
- **Per-method sections (D-11)**: Plotly spectrum div, Plotly regression div, PNG heatmap (relative `<img src>`), timing block (ML/DFT Gaussian elapsed + speedup), degenerate notes, metrics mini-table
- **Data export (D-07)**: Calls `export_report_data` at end of `generate_report` to write `report_data.json` + `summary_metrics.csv`
- **HTML escaping (T-23-01)**: `_esc` helper wraps `html.escape(str(s), quote=True)`. Applied to molecule_name, method names, experimental source/molecule_name/CAS, heatmap paths
- **Overtones section**: Option B explicit stub with `overtones-placeholder` CSS class and `TODO(phase-follow-up)` comment. Real data path (Option A) with `overtone-table` class also implemented but currently unreachable since pipeline does not emit structured overtone records

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed conftest fake_metrics_b missing num_intensity_filtered**
- **Found during:** Task 1 (test collection)
- **Issue:** `ComparisonMetrics.__init__()` requires `num_intensity_filtered` positional argument, but `fake_metrics_b` fixture in conftest.py omitted it
- **Fix:** Added `num_intensity_filtered=0` to the fixture
- **Files modified:** tests/conftest.py
- **Commit:** b2f9a26

## Test Results

All 8 `test_html_report.py` tests GREEN:
- test_per_method_section_has_plotly_div
- test_plotlyjs_emitted_once
- test_mode_flag_harmonic_skips_overtones
- test_css_shared_across_modes
- test_timing_section_renders
- test_experimental_trace_present_when_experimental_available
- test_html_well_formed
- test_html_escapes_molecule_name

All 40 Wave 1 regression tests GREEN (plotly_builders, report_data, executive_summary).

ruff check + ty check both pass clean.

## Known Stubs

| File | Location | Stub | Reason |
|------|----------|------|--------|
| html_report_generator.py | `_create_overtones_section` | overtones-placeholder div | Overtone/combination band data is embedded in SpectrumData labels but NOT surfaced as structured records in comparison dicts. TODO(phase-follow-up) marks this for future phase. |

## Self-Check: PASSED
