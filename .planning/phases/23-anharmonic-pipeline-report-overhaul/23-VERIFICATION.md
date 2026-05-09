---
phase: 23-anharmonic-pipeline-report-overhaul
verified: 2026-04-13T14:30:00Z
status: gaps_found
score: 5/7
overrides_applied: 0
gaps:
  - truth: "batch_report.py imports _shared_css.build_css and wraps user-controlled strings with html.escape -- single source CSS + XSS mitigation (Plan 05 Task 1)"
    status: failed
    reason: "batch_report.py was never modified. No import of _shared_css, no import html, no html.escape calls. Plan 05 SUMMARY claims this was done but no commit touched batch_report.py."
    artifacts:
      - path: "mace_gaussian/analysis/batch_report.py"
        issue: "Missing `from ._shared_css import build_css` and `import html` + html.escape wrapping of user-controlled f-string interpolations"
    missing:
      - "Import and call _shared_css.build_css() instead of inline _build_css()"
      - "Import html stdlib and wrap molecule, combo, ml_hw, dft_hw, experimental.source with html.escape(..., quote=True)"
      - "test_batch_report_escapes_molecule_name must pass (currently fails)"
  - truth: "analysis_workflow.py passes mode='harmonic' or mode='anharmonic' to HTMLReportGenerator"
    status: failed
    reason: "HTMLReportGenerator is constructed without mode= kwarg in analysis_workflow.py line 894. Defaults to 'anharmonic', so anharmonic reports work but harmonic reports will incorrectly include overtones section."
    artifacts:
      - path: "mace_gaussian/analysis/analysis_workflow.py"
        issue: "Line 894: HTMLReportGenerator(...) missing mode= parameter"
    missing:
      - "Pass mode= parameter from ComparisonWorkflow context to HTMLReportGenerator constructor"
---

# Phase 23: Anharmonic Pipeline & Report Overhaul Verification Report

**Phase Goal:** The anharmonic analysis pipeline produces a thesis-quality HTML report that integrates all v1.2 features into a polished presentation
**Verified:** 2026-04-13T14:30:00Z
**Status:** gaps_found
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | The anharmonic HTML report integrates Lorentzian spectra, experimental overlay, timing breakdown, and degenerate-mode-aware mode matching into a single cohesive document | VERIFIED | html_report_generator.py (691 lines) imports plotly_builders for Lorentzian broadened spectra, experimental_on_grid for overlay, renders timing blocks per method section, heatmap PNGs, degenerate notes. 48 module tests pass. |
| 2 | The report includes per-molecule summary cards showing key metrics: R2_freq, R2_intensity, RMSE, speedup, experimental agreement | VERIFIED | executive_summary.py rank_methods + build_verdict wired into HTMLReportGenerator._create_executive_summary. Method cards rendered with CSS class "method-card", best method highlighted with "best-method". 16 executive_summary tests pass. |
| 3 | The report is visually thesis-ready: consistent styling, publication-quality interactive plots, clear narrative flow | VERIFIED | Plotly interactive figures (not static PNGs) for spectrum/regression. Shared CSS via _shared_css.build_css(). Report flow: header > nav > executive summary > combined plots > per-method sections > experimental info > summary table > overtones > footer. Emit-once CDN pattern verified by test. |
| 4 | Structured data export (report_data.json + summary_metrics.csv) accompanies every report | VERIFIED | html_report_generator.py line 122 calls export_report_data. report_data.py (188 lines) writes JSON with schema_version=1 + CSV. 16 tests pass. |
| 5 | HTMLReportGenerator is mode-aware (harmonic/anharmonic) with Plotly, shared CSS, html.escape | VERIFIED | __init__ accepts mode param, plotly_js param. Harmonic mode skips overtones (test_mode_flag_harmonic_skips_overtones passes). _esc helper wraps html.escape. CSS from _shared_css. |
| 6 | batch_report.py uses shared CSS and html.escape for XSS mitigation | FAILED | batch_report.py has no import of _shared_css, no import html, no html.escape calls. test_batch_report_escapes_molecule_name fails. Plan 05 SUMMARY claimed this was done but no code changes were committed. |
| 7 | analysis_workflow.py passes mode flag to HTMLReportGenerator | FAILED | Line 894 constructs HTMLReportGenerator without mode= parameter. Harmonic analysis will incorrectly render overtones section. |

**Score:** 5/7 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `mace_gaussian/analysis/html_report_generator.py` | Plotly-powered HTMLReportGenerator with exec summary, shared CSS, mode flag, html.escape | VERIFIED | 691 lines, all key imports wired, 8/8 tests pass |
| `mace_gaussian/analysis/plotly_builders.py` | Pure-function Plotly figure builders | VERIFIED | 286 lines, 4 functions (spectrum, regression, combined, experimental_on_grid), function-local plotly imports, 8/8 tests pass |
| `mace_gaussian/analysis/_shared_css.py` | Single source CSS for both report generators | VERIFIED | 363 lines, build_css() exports all required class names |
| `mace_gaussian/analysis/executive_summary.py` | Ranking + verdict logic | VERIFIED | 124 lines, 3 public functions, 16/16 tests pass |
| `mace_gaussian/analysis/report_data.py` | JSON + CSV exporter | VERIFIED | 188 lines, export_report_data writes JSON + CSV, 16/16 tests pass |
| `mace_gaussian/analysis/batch_report.py` | Shared CSS + html.escape applied | FAILED | No changes made -- still uses inline _build_css, no html.escape |
| `mace_gaussian/analysis/analysis_workflow.py` | Passes mode= to HTMLReportGenerator | FAILED | mode= kwarg missing from constructor call |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| html_report_generator.py | plotly_builders.py | `from .plotly_builders import` | WIRED | Line 25: imports build_spectrum_figure, build_regression_figure, build_combined_spectrum_figure, experimental_on_grid |
| html_report_generator.py | _shared_css.py | `from ._shared_css import build_css` | WIRED | Line 18 |
| html_report_generator.py | executive_summary.py | `from .executive_summary import` | WIRED | Lines 20-23: imports rank_methods, build_verdict, compute_experimental_agreement |
| html_report_generator.py | report_data.py | `from .report_data import export_report_data` | WIRED | Line 31, called at line 122 |
| batch_report.py | _shared_css.py | `from ._shared_css import build_css` | NOT_WIRED | Import does not exist |
| batch_report.py | html stdlib | `import html` | NOT_WIRED | Import does not exist |
| analysis_workflow.py | html_report_generator.py | `HTMLReportGenerator(..., mode=...)` | PARTIAL | Constructor called but mode= kwarg missing |

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| html_report_generator.py | analysis_results["comparisons"] | analysis_workflow.run_single_comparison | Yes -- DB/file-backed comparison dicts with metrics, spectra, timing | FLOWING |
| html_report_generator.py | analysis_results["experimental"] | NIST fetcher | Yes -- real NIST data when available | FLOWING |
| executive_summary.py | comparisons list | html_report_generator enrichment | Yes -- broadened spectra + agreement scores | FLOWING |
| report_data.py | analysis_results | html_report_generator | Yes -- writes real JSON/CSV | FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Phase 23 module tests pass | `pytest tests/test_plotly_builders.py tests/test_report_data.py tests/test_executive_summary.py tests/test_html_report.py` | 48 passed, 0 failed | PASS |
| Plotly builders importable | `python -c "from mace_gaussian.analysis.plotly_builders import build_spectrum_figure"` | exits 0 | PASS |
| Shared CSS importable | `python -c "from mace_gaussian.analysis._shared_css import build_css; assert len(build_css()) > 500"` | exits 0 | PASS |
| Batch report XSS test | `pytest tests/test_batch_report.py -k escapes_molecule_name` | FAILED | FAIL |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| ANAL-01 | 23-01 through 23-05 | Anharmonic analysis pipeline produces thesis-quality HTML report integrating Lorentzian spectra, experimental overlay, timing, and mode matching | SATISFIED | HTMLReportGenerator produces cohesive Plotly-powered report with all features. E2E smoke test on water confirmed (commit 38f362e). |
| ANAL-02 | 23-01 through 23-05 | Report includes per-molecule summary cards with key metrics (R-squared, RMSE, speedup, experimental agreement) | SATISFIED | executive_summary.rank_methods produces ranked cards; HTMLReportGenerator._create_executive_summary renders method-card divs with all metrics. |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| html_report_generator.py | 606 | TODO(phase-follow-up): overtone data | Info | Expected -- explicit stub policy (Option B). Overtone data not yet in comparison dicts. |
| batch_report.py | (throughout) | Raw f-string interpolation of user-controlled strings | Warning | XSS surface. Plan 05 was supposed to fix this but changes were never committed. |

### Human Verification Required

### 1. Visual Report Quality

**Test:** Open `analysis_results/water/anharmonic/report.html` in a browser
**Expected:** Thesis-quality layout with interactive Plotly figures, executive summary at top with method cards, timing breakdowns, experimental overlay on spectra
**Why human:** Visual quality and layout cannot be verified programmatically

### 2. Plotly Interactive Figures

**Test:** Hover over spectrum and regression plots in the report
**Expected:** Tooltips show wavenumber and absorbance values; zoom/pan works
**Why human:** Interactive behavior requires a browser

### Gaps Summary

Two gaps remain from Plan 05 (Wave 3), which wrote a SUMMARY claiming completion but never committed the actual code changes:

1. **batch_report.py not updated** -- The batch report still uses its own inline `_build_css()` function and has no `html.escape` wrapping of user-controlled strings. The XSS mitigation test (`test_batch_report_escapes_molecule_name`) created in Plan 01 remains RED. This is a concrete security gap (T-23-01) and a CSS consistency violation (D-04).

2. **analysis_workflow.py missing mode= kwarg** -- The `ComparisonWorkflow._generate_html_report` method constructs `HTMLReportGenerator` without passing `mode=`. While the default is `"anharmonic"` (so anharmonic reports work correctly), harmonic analyses will incorrectly include the overtones section. This is a functional gap for harmonic mode users.

Both gaps are localized to Plan 05 scope and do not affect the core report generation pipeline (Plans 01-04), which is fully functional and tested.

---

_Verified: 2026-04-13T14:30:00Z_
_Verifier: Claude (gsd-verifier)_
