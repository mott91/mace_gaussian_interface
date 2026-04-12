---
phase: 23-anharmonic-pipeline-report-overhaul
plan: 01
subsystem: tests, pyproject
tags: [report, plotly, tests, pyproject, tdd, red]
dependency_graph:
  requires: []
  provides: [conftest-fixtures, red-test-suite, plotly-dependency]
  affects: [tests/conftest.py, pyproject.toml]
tech_stack:
  added: [plotly>=5.20.0]
  patterns: [pytest-importorskip-for-red-tests, fixture-composition]
key_files:
  created:
    - tests/test_plotly_builders.py
    - tests/test_report_data.py
    - tests/test_executive_summary.py
    - tests/test_html_report.py
  modified:
    - tests/conftest.py
    - tests/test_batch_report.py
    - pyproject.toml
decisions:
  - Used pytest.importorskip for modules that do not exist yet (plotly_builders, report_data, executive_summary) so collection succeeds
  - Placed Phase 23 imports (json, pathlib) at top of test_batch_report.py to satisfy ruff E402
metrics:
  duration: 252s
  completed: "2026-04-12T06:20:52Z"
  tasks_completed: 2
  tasks_total: 2
  files_created: 4
  files_modified: 3
---

# Phase 23 Plan 01: RED Test Foundation + Plotly Dependency Summary

Declared plotly>=5.20.0,<7 in pyproject.toml and created 4 RED test files (22 test functions total) plus XSS escape test appended to test_batch_report.py, with shared fixture family in conftest.py sourced from exact codebase dataclass field names.

## Task Summary

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Declare plotly dependency + conftest fixtures | 9aadfb8 | pyproject.toml, tests/conftest.py |
| 2 | Create RED test files + batch XSS test | 3e24a41 | tests/test_plotly_builders.py, tests/test_report_data.py, tests/test_executive_summary.py, tests/test_html_report.py, tests/test_batch_report.py |

## What Was Built

**Task 1 -- Dependency + Fixtures:**
- Added `plotly>=5.20.0,<7` to `[project].dependencies` in pyproject.toml
- Extended `tests/conftest.py` with 4 new fixtures: `fake_metrics` (ComparisonMetrics), `fake_spectrum` (SpectrumData), `fake_experimental` (ExperimentalSpectrum), `fake_analysis_results` (full comparison dict)
- All fixtures use exact field names verified from `analyze_spectra.py` and `nist_fetcher.py`

**Task 2 -- RED Test Files:**
- `tests/test_plotly_builders.py` (81 lines, 5 tests): spectrum figure, experimental trace, x-axis reversed, regression figure, combined multi-method figure
- `tests/test_report_data.py` (49 lines, 4 tests): JSON export write, round-trip, numpy serialization, null experimental
- `tests/test_executive_summary.py` (76 lines, 5 tests): identical spectra agreement, zero-returns-NaN, deterministic ranking, verdict with/without experimental
- `tests/test_html_report.py` (98 lines, 8 tests): plotly div, CDN emitted once, harmonic skips overtones, CSS shared, timing section, experimental trace, well-formed HTML, XSS escape (T-23-01)
- `tests/test_batch_report.py` (appended): `test_batch_report_escapes_molecule_name` for T-23-01 batch-path mitigation

**RED State:** plotly_builders, report_data, executive_summary modules do not exist -- tests skip on collection via importorskip. html_report tests fail on missing `mode`/`plotly_js` kwargs. batch XSS test fails because html.escape not yet applied. This is intentional -- drives Wave 1-3 implementation.

## Deviations from Plan

None -- plan executed exactly as written.

## Verification Results

- `pytest --collect-only` exits 0 (13 tests collected, 3 modules skipped via importorskip)
- `ruff check` passes on all 6 touched files
- `grep -c "plotly>=5.20.0,<7" pyproject.toml` = 1
- `grep -c "fake_analysis_results" tests/conftest.py` = 1
- All line count and test count acceptance criteria met

## Threat Model Coverage

| Threat | Test File | Test Function | Status |
|--------|-----------|---------------|--------|
| T-23-01 (single-molecule XSS) | tests/test_html_report.py | test_html_escapes_molecule_name | RED |
| T-23-01 (batch-path XSS) | tests/test_batch_report.py | test_batch_report_escapes_molecule_name | RED |

## Self-Check: PASSED
