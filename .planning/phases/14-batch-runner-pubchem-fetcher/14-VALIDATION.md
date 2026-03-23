---
phase: 14
slug: batch-runner-pubchem-fetcher
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-23
---

# Phase 14 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.0+ |
| **Config file** | pyproject.toml `[tool.pytest.ini_options]` |
| **Quick run command** | `pytest tests/test_pubchem.py tests/test_batch.py -x -q` |
| **Full suite command** | `pytest tests/ -ra` |
| **Estimated runtime** | ~5 seconds |

---

## Sampling Rate

- **After every task commit:** Run `pytest tests/test_pubchem.py tests/test_batch.py -x -q`
- **After every plan wave:** Run `pytest tests/ -ra`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 10 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 14-01-01 | 01 | 1 | BATCH-01 | unit | `pytest tests/test_pubchem.py -x` | ❌ W0 | ⬜ pending |
| 14-02-01 | 02 | 1 | BATCH-02 | unit | `pytest tests/test_batch.py::test_batch_sequential -x` | ❌ W0 | ⬜ pending |
| 14-02-02 | 02 | 1 | BATCH-03 | unit | `pytest tests/test_batch.py::test_manifest_restart -x` | ❌ W0 | ⬜ pending |
| 14-02-03 | 02 | 1 | BATCH-04 | unit | `pytest tests/test_batch.py::test_skip_dft_flag -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test_pubchem.py` — stubs for BATCH-01 (fetch command, PubChem API mocking, error cases)
- [ ] `tests/test_batch.py` — stubs for BATCH-02, BATCH-03, BATCH-04 (batch loop, manifest CRUD, restart logic, flag passthrough)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| PubChem live API fetch | BATCH-01 | Requires network access to PubChem | Run `mace-gaussian fetch aspirin` and verify `molecules/aspirin.xyz` has valid 3D coordinates |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 10s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
