---
phase: 2
slug: content
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-07
---

# Phase 2 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest (existing, 131 tests passing) |
| **Config file** | `pyproject.toml` |
| **Quick run command** | `python presentation/generate_pptx_v2.py && python -c "from pptx import Presentation; prs = Presentation('presentation/presentation_v2.pptx'); print(f'{len(prs.slides)} slides OK')"` |
| **Full suite command** | `pytest` |
| **Estimated runtime** | ~5 seconds |

---

## Sampling Rate

- **After every task commit:** Run quick run command above
- **After every plan wave:** Run quick command + manual visual spot-check of output pptx
- **Before `/gsd:verify-work`:** Generator runs clean, slide count is 13, content verified
- **Max feedback latency:** ~5 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| CONT-01 | 01 | 1 | CONT-01 | smoke | `python presentation/generate_pptx_v2.py` (no crash = pass) | ✅ | ⬜ pending |
| CONT-02 | 01 | 1 | CONT-02 | smoke | Count slides: expect 13 | ✅ via pptx parse | ⬜ pending |
| CONT-03 | 01 | 1 | CONT-03 | smoke | Search slide text for `water/report.html` and `aspirin/report.html` | ✅ via pptx parse | ⬜ pending |
| CONT-04 | 01 | 1 | CONT-04 | smoke | Search slide text for `zmq_server.py`, absence of `gm_main.py` | ✅ via pptx parse | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

None — existing infrastructure covers all phase requirements. No new test files or framework changes needed.

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Results table visually readable and sorted correctly | CONT-01 | Font size and column width require visual inspection | Open pptx, verify water table shows 8 rows sorted by MAE ascending, aspirin shows 4 rows |
| HTML report links and ASCII art legible at slide scale | CONT-03 | Layout quality requires visual check | Open pptx, verify ASCII art and report paths are readable on the overview slide |
| ZMQ flow diagram arrows and labels correct | CONT-04 | Diagram layout requires visual inspection | Open pptx, verify flow diagram shows hook point and zmq_server.py label |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 5s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
