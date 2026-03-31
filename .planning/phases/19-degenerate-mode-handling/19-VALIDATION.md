---
phase: 19
slug: degenerate-mode-handling
status: draft
nyquist_compliant: true
wave_0_complete: true
created: 2026-03-30
---

# Phase 19 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest |
| **Config file** | pyproject.toml `[tool.pytest.ini_options]` |
| **Quick run command** | `pytest tests/test_mode_matching.py -x` |
| **Full suite command** | `pytest tests/ --strict-markers -ra` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `pytest tests/test_mode_matching.py tests/test_gaussian_parser.py -x`
- **After every plan wave:** Run `pytest tests/ --strict-markers -ra`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 19-01-01 | 01 | 1 | MODE-05 | unit | `pytest tests/test_mode_matching.py::TestDegenerateDetection -x` | ❌ W0 | ⬜ pending |
| 19-01-02 | 01 | 1 | MODE-05 | unit | `pytest tests/test_mode_matching.py::TestSubspaceOverlap -x` | ❌ W0 | ⬜ pending |
| 19-01-03 | 01 | 1 | MODE-05 | integration | `pytest tests/test_mode_matching.py::TestCH4DegenerateOverlap -x` | ❌ W0 | ⬜ pending |
| 19-01-04 | 01 | 1 | MODE-06 | unit | `pytest tests/test_mode_matching.py::TestGroupAwareStatistics -x` | ❌ W0 | ⬜ pending |
| 19-01-05 | 01 | 1 | MODE-06 | unit | `pytest tests/test_mode_matching.py::TestGroupRegressionData -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

*Wave 0 is handled within Plan 01 Task 1 (TDD task): tests are written first (RED phase) then implementation follows (GREEN phase) in a single self-contained task. This is compliant because the test-first artifact is a discrete step within the task, and the verify command validates all tests pass before the task is marked done.*

*Existing infrastructure: pytest configured, CH4 fixtures (ch4_dft_fchk, ch4_ml_fchk, ch4_dft_log) already available in conftest.py*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Diamond markers visible in regression plots | MODE-06 (D-10) | Visual verification of plot aesthetics | Run `python run_analysis_harmonic.py methane`, open regression plot PNG, confirm diamond markers for degenerate groups |
| Heatmap collapsed rows/columns readable | MODE-05 (D-09) | Visual verification of label layout | Run harmonic analysis on CH4, open heatmap PNG, confirm group labels like "T2 @1356" |
| HTML report shows group labels | MODE-05 (D-08) | Visual verification of report rendering | Open HTML report, check mode matching section for "T2 (3-fold)" labels |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
