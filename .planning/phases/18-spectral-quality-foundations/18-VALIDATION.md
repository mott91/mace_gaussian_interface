---
phase: 18
slug: spectral-quality-foundations
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-28
---

# Phase 18 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest >= 7.0.0 |
| **Config file** | `pyproject.toml` [tool.pytest.ini_options] |
| **Quick run command** | `pytest tests/test_spectral_broadening.py -x -q --no-header` |
| **Full suite command** | `pytest tests/ -ra` |
| **Estimated runtime** | ~5 seconds |

---

## Sampling Rate

- **After every task commit:** Run `pytest tests/test_spectral_broadening.py -x -q --no-header`
- **After every plan wave:** Run `pytest tests/ -ra`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 5 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 18-01-01 | 01 | 0 | SPEC-01, SPEC-02 | unit stubs | `pytest tests/test_spectral_broadening.py --collect-only` | ❌ W0 | ⬜ pending |
| 18-02-01 | 02 | 1 | SPEC-01 | unit | `pytest tests/test_spectral_broadening.py::test_lorentzian_fwhm -x` | ❌ W0 | ⬜ pending |
| 18-02-02 | 02 | 1 | SPEC-01 | unit | `pytest tests/test_spectral_broadening.py::test_lorentzian_peak_height -x` | ❌ W0 | ⬜ pending |
| 18-02-03 | 02 | 1 | SPEC-01 | unit | `pytest tests/test_spectral_broadening.py::test_default_fwhm -x` | ❌ W0 | ⬜ pending |
| 18-03-01 | 03 | 1 | SPEC-02 | unit | `pytest tests/test_spectral_broadening.py::test_intensity_filter_threshold -x` | ❌ W0 | ⬜ pending |
| 18-03-02 | 03 | 1 | SPEC-02 | unit | `pytest tests/test_spectral_broadening.py::test_frequency_metrics_unfiltered -x` | ❌ W0 | ⬜ pending |
| 18-03-03 | 03 | 1 | SPEC-02 | unit | `pytest tests/test_spectral_broadening.py::test_intensity_filter_all_below -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test_spectral_broadening.py` — stubs for SPEC-01, SPEC-02 (new file)
- [ ] No framework install needed — pytest already configured in pyproject.toml

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Lorentzian plot looks correct visually | SPEC-01 | Visual inspection of peak shapes | Run `python run_analysis.py water`, check spectrum plot has Lorentzian peaks |
| HTML report shows correct FWHM text | SPEC-01 | HTML rendering check | Open HTML report, verify FWHM parameter displayed correctly |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 5s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
