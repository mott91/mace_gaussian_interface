# Phase 1: Narrative & Structure — Validation Strategy

**Phase:** 1
**Phase slug:** narrative-structure
**Date:** 2026-03-06

---

## Test Framework

| Property | Value |
|----------|-------|
| Framework | Manual visual inspection (no automated test for .pptx content) |
| Config file | none |
| Quick run command | `python presentation/generate_pptx_v2.py` |
| Full suite command | `python presentation/generate_pptx_v2.py && echo "Saved OK"` |

---

## Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command |
|--------|----------|-----------|-------------------|
| NARR-01 | 11 slides in correct order | smoke | `python -c "from pptx import Presentation; p=Presentation('presentation/presentation_v2.pptx'); print(len(p.slides), 'slides')"` |
| NARR-02 | Status slide contains research journey text (no "CI", no "131 tests") | smoke | Manual inspection of generated .pptx |
| NARR-03 | VPT2 slide absent; "anharmonic" appears only in IR slide footer | smoke | `python -c "from pptx import Presentation; p=Presentation('presentation/presentation_v2.pptx'); [print(s.shapes[0].text if s.shapes else '') for s in p.slides]"` |

---

## Sampling Rate

- **Per task commit:** `python presentation/generate_pptx_v2.py` (must exit 0)
- **Per wave merge:** Visual review of generated .pptx slide count and slide order
- **Phase gate:** All 11 slides present in correct order, no NameError on generation

---

## Wave 0 Gaps (Prerequisites)

- [ ] Fix output path in `main()` before first run — current path `/home/mot/mace_gaussian/presentation/presentation_v2.pptx` will cause `FileNotFoundError`
- [ ] No pytest fixtures needed — this is a pure generation script
