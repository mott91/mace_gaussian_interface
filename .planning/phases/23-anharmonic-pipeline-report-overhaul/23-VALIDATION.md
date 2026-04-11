---
phase: 23
slug: anharmonic-pipeline-report-overhaul
status: ready
nyquist_compliant: true
wave_0_complete: false
created: 2026-04-11
updated: 2026-04-11
---

# Phase 23 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 9.0.2 (verified in mace4ir_v2 micromamba env) |
| **Config file** | `pyproject.toml [tool.pytest.ini_options]` + `tests/conftest.py` |
| **Quick run command** | `micromamba run -n mace4ir_v2 pytest tests/test_plotly_builders.py tests/test_report_data.py tests/test_executive_summary.py tests/test_html_report.py tests/test_batch_report.py -x` |
| **Full suite command** | `micromamba run -n mace4ir_v2 pytest tests/ -x` |
| **Estimated runtime** | ~10 seconds for Phase 23 subset; ~60 seconds for full suite |

---

## Sampling Rate

- **After every task commit:** Run Phase 23 quick command (the four new test files + batch_report)
- **After every plan wave:** Run quick command + `tests/test_regression.py` + `tests/test_spectral_broadening.py`
- **Before `/gsd-verify-work`:** Full suite must be green (modulo 11 known pre-existing failures per MEMORY.md)
- **Max feedback latency:** ~10 seconds (quick run)

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 23-01-01 | 01 | 0 | ANAL-01, ANAL-02 | — | pyproject declares plotly; conftest imports safely | unit (collect) | `pytest tests/conftest.py --collect-only` | N/A | ⬜ pending |
| 23-01-02 | 01 | 0 | ANAL-01, ANAL-02 | T-23-01 | RED tests for HTML escaping exist in BOTH test_html_report.py AND test_batch_report.py | unit (collect) | `pytest tests/test_html_report.py tests/test_batch_report.py --collect-only && grep -c "test_batch_report_escapes_molecule_name" tests/test_batch_report.py` | ❌ W0 creates | ⬜ pending |
| 23-02-01 | 02 | 1 | ANAL-01 | — | `_shared_css.build_css()` importable | unit | `python -c "from mace_gaussian.analysis._shared_css import build_css; assert build_css()"` | ❌ W1 creates | ⬜ pending |
| 23-02-02 | 02 | 1 | ANAL-01 | — | figure builders return go.Figure with correct traces + reversed x-axis | unit | `pytest tests/test_plotly_builders.py -x` | ❌ W1 creates | ⬜ pending |
| 23-03-01 | 03 | 1 | ANAL-02 | — | `rank_methods` deterministic; `build_verdict` format correct | unit | `pytest tests/test_executive_summary.py -x` | ❌ W1 creates | ⬜ pending |
| 23-03-02 | 03 | 1 | ANAL-01 | T-23-04 | JSON round-trip; numpy serialized | unit | `pytest tests/test_report_data.py -x` | ❌ W1 creates | ⬜ pending |
| 23-04-01 | 04 | 2 | ANAL-01, ANAL-02 | T-23-01, T-23-05 | HTML escape + Plotly emit-once + mode flag + overtones data-or-stub contract | unit + grep | `pytest tests/test_html_report.py -x && grep -cE "overtone-table\|overtones-placeholder" mace_gaussian/analysis/html_report_generator.py && grep -c "TODO(phase-follow-up)" mace_gaussian/analysis/html_report_generator.py` | ❌ W2 creates | ⬜ pending |
| 23-05-01 | 05 | 3 | ANAL-01, ANAL-02 | T-23-01 | batch_report uses shared CSS; batch_report HTML-escapes molecule/combo/hardware strings; workflow passes mode flag | unit + regression + grep | `pytest tests/test_batch_report.py tests/test_html_report.py tests/test_plotly_builders.py tests/test_report_data.py tests/test_executive_summary.py -x && grep -c "html.escape" mace_gaussian/analysis/batch_report.py && grep -c "_esc(" mace_gaussian/analysis/batch_report.py` | N/A | ⬜ pending |
| 23-05-01b | 05 | 3 | ANAL-01 | T-23-01 | batch_report XSS escape test passes (RED→GREEN transition) | unit | `pytest tests/test_batch_report.py::test_batch_report_escapes_molecule_name -x` | N/A | ⬜ pending |
| 23-05-02 | 05 | 3 | ANAL-01, ANAL-02 | — | End-to-end smoke on water fixture + overtones section content check | integration | `micromamba run -n mace4ir_v2 python run_analysis_harmonic.py water && test -f analysis_results_harmonic/water/report_data.json` | ❌ produces artifact | ⬜ pending |
| 23-05-03 | 05 | 3 | ANAL-01, ANAL-02 | — | Human visual verification in browser (no JS errors, overtones section has real data OR explicit placeholder) | checkpoint:human-verify | (manual) `xdg-open analysis_results_harmonic/water/report.html` | N/A | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test_plotly_builders.py` — RED tests for `build_spectrum_figure`, `build_regression_figure`, `build_combined_spectrum_figure`, `experimental_on_grid` (≥5 test functions)
- [ ] `tests/test_report_data.py` — RED tests for `export_report_data` JSON round-trip + numpy serialization + `experimental=None` (≥4 test functions)
- [ ] `tests/test_executive_summary.py` — RED tests for `compute_experimental_agreement`, `rank_methods`, `build_verdict` (≥5 test functions)
- [ ] `tests/test_html_report.py` — RED tests for Plotly integration, mode flag, shared CSS, HTML well-formed, HTML escape (T-23-01) (≥7 test functions)
- [ ] `tests/test_batch_report.py` — RED test `test_batch_report_escapes_molecule_name` for T-23-01 batch-path mitigation (≥1 new test function, appended to existing file)
- [ ] `tests/conftest.py` — append `fake_metrics`, `fake_spectrum`, `fake_experimental`, `fake_analysis_results` fixtures
- [ ] `pyproject.toml` — declare `plotly>=5.20.0,<7` in `[project].dependencies`

All of the above are created by Plan 23-01 (Wave 0). Wave 1+ tests MUST turn green as modules are built. The batch_report XSS escape test stays RED until Plan 23-05 Task 1 applies `html.escape` wrapping in `batch_report.py`.

---

## Overtones Section Data-or-Stub Contract (Plan 23-04 Task 1)

Per the Plan 04 revision, `_create_overtones_section` in `html_report_generator.py` MUST emit one of two shapes — NEVER a silent empty section:

**Option A (real data path)** — triggered when `analysis_results["overtones"]` or `comparison["overtones"]` is populated:
- Renders `<table class="overtone-table data-table">` with real frequency/intensity/type rows
- Grep target: `overtone-table` appears in source file

**Option B (explicit stub path)** — triggered when no overtone data is threaded through the pipeline:
- Renders `<div class="overtones-placeholder">` with visible stub text "Placeholder — populated in a follow-up phase"
- Grep target: `overtones-placeholder` appears in source file

**Always required in source regardless of runtime path:**
- `# TODO(phase-follow-up)` comment in `_create_overtones_section` or module-level

**Grep verification:**
```
grep -c "def _create_overtones_section" mace_gaussian/analysis/html_report_generator.py  # == 1
grep -cE "overtone-table|overtones-placeholder" mace_gaussian/analysis/html_report_generator.py  # >= 1
grep -c "TODO(phase-follow-up)" mace_gaussian/analysis/html_report_generator.py  # >= 1
```

---

## T-23-01 XSS Mitigation Coverage Map

| Path | Component | Test | Plan |
|------|-----------|------|------|
| Single-molecule report | `HTMLReportGenerator` in `html_report_generator.py` | `tests/test_html_report.py::test_html_escapes_molecule_name` | 23-01 (RED) → 23-04 (GREEN) |
| Batch report | `batch_report._generate_html`, `_build_timing_html`, `_build_leaderboard_html` | `tests/test_batch_report.py::test_batch_report_escapes_molecule_name` | 23-01 (RED) → 23-05 (GREEN) |

**Grep verification for batch path after Plan 05:**
```
grep -c "^import html$" mace_gaussian/analysis/batch_report.py  # == 1
grep -c "html.escape" mace_gaussian/analysis/batch_report.py  # >= 1 (inside _esc helper)
grep -c "_esc(" mace_gaussian/analysis/batch_report.py  # >= 5 (helper def + 4+ call sites)
grep -cE "\{mol_name\}|\{r\['molecule'\]\}|\{row\['combo'\]\}" mace_gaussian/analysis/batch_report.py  # == 0 (no raw interpolations)
```

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Browser JS errors | ANAL-01 | No headless browser in env | Open `analysis_results_harmonic/water/report.html` in firefox/chrome, check DevTools console for red errors |
| Visual regression vs design expectation | ANAL-01, ANAL-02 | Aesthetic judgment | Open report, confirm executive summary at top, method cards, Plotly interactivity, experimental trace visible |
| Plotly interactive zoom/hover works | ANAL-01 (D-06 interactivity motivation) | Requires user gestures | Hover over spectrum trace to see tooltip; click-drag to zoom; double-click to reset |
| Harmonic vs anharmonic visual equivalence | ANAL-01 (D-04 same product) | Side-by-side visual comparison | Open both `analysis_results_harmonic/water/report.html` and `analysis_results/water/report.html`; confirm identical CSS, layout, header; anharmonic has overtones section, harmonic does not |
| Overtones section content (anharmonic only) | ANAL-01 (D-04) | Aesthetic + content inspection | Confirm overtones section shows EITHER a real data table OR the explicit `overtones-placeholder` div with visible stub text — NEVER a silent empty section under the `<h2>` heading |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: every task runs full Phase 23 suite after commit
- [x] Wave 0 covers all MISSING references (four new test files + conftest extension + batch_report XSS RED test)
- [x] No watch-mode flags
- [x] Feedback latency < 10s (quick run)
- [x] `nyquist_compliant: true` set in frontmatter
- [x] T-23-01 mitigation has end-to-end coverage: single-molecule (23-04) + batch (23-05), both with RED→GREEN tests
- [x] Overtones section has explicit data-or-stub contract (not a word-match heuristic)

**Approval:** plan-review pending (ready for /gsd-execute-phase 23)
</content>
</invoke>
