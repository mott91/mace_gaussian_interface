---
phase: 01-narrative-structure
plan: 01
subsystem: presentation
tags: [python-pptx, slides, generate_pptx_v2]

requires: []
provides:
  - generate_pptx_v2.py restructured to 11-slide layout (VPT2 slide removed)
  - slide_results_overview function: terminal ls listing of 8 molecule HTML reports
  - slide_results_table function: harmonic benchmark table sorted by frequency MAE
  - presentation/presentation_v2.pptx generated and verified (11 slides)
affects:
  - 01-02 (rewrites existing slide functions — depends on clean 11-slide structure)

tech-stack:
  added: [python-pptx (installed via pip3 for verification)]
  patterns:
    - content_slide helper for all new results slides
    - Unicode literals for special characters (subscripts, arrows, checkmarks)

key-files:
  created: []
  modified:
    - presentation/generate_pptx_v2.py

key-decisions:
  - "Delete slide_vpt2 entirely (function + call) — harmonic-first narrative, VPT2 is future work"
  - "Replace per-molecule slides (water, aspirin) with cross-molecule overview + table slides"
  - "Output path fixed to relative presentation/presentation_v2.pptx (was hardcoded absolute path)"
  - "slide_results_table data sourced directly from analysis_results_harmonic/water/data/metrics_summary.json values"

patterns-established:
  - "New results slides use content_slide (not two_col_slide) for density"
  - "Table rows: GREEN for mace_ml combos, YELLOW for mace_mp rows"

requirements-completed: [NARR-01, NARR-03]

duration: 15min
completed: 2026-03-07
---

# Phase 01 Plan 01: Narrative Structure SUMMARY

**11-slide generate_pptx_v2.py with VPT2 slide removed, relative output path fixed, and two new results slides (HTML report listing + harmonic benchmark table) replacing the old per-molecule slides**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-03-07T07:10:00Z
- **Completed:** 2026-03-07T07:26:59Z
- **Tasks:** 2
- **Files modified:** 2 (generate_pptx_v2.py, presentation_v2.pptx)

## Accomplishments

- Deleted `slide_vpt2` function and call — 12-slide script is now 11-slide
- Fixed hardcoded absolute output path (`/home/mot/mace_gaussian/...`) to relative `presentation/presentation_v2.pptx`
- Added `slide_results_overview`: terminal `ls` idiom listing 8 molecules with mode/combo counts
- Added `slide_results_table`: harmonic benchmark table with real water metrics, 8 combos sorted by frequency MAE (mace_anicc best at 19 cm⁻¹, mace_mp worst at 106 cm⁻¹)
- Removed `slide_results_water` and `slide_results_aspirin` function definitions
- Script runs clean, generates `presentation_v2.pptx` with exactly 11 slides

## Task Commits

1. **Task 1: Fix output path and delete slide_vpt2 entirely** - `e3ec5f7` (feat)
2. **Task 2: Add slide_results_overview and slide_results_table, verify 11 slides** - `73a71f1` (feat)

## Files Created/Modified

- `presentation/generate_pptx_v2.py` - Restructured: VPT2 slide removed, path fixed, two new results slide functions, main() updated
- `presentation/presentation_v2.pptx` - Generated output (proof script runs, 11 slides verified)

## Decisions Made

- Deleted `slide_vpt2` entirely rather than repurposing it — harmonic-first narrative means VPT2 theory is out of scope for this presentation
- Used `content_slide` for both new results slides (not `two_col_slide`) — more lines needed for the benchmark table
- Table data sourced from the real `analysis_results_harmonic/water/data/metrics_summary.json` values
- `mace_anicc`/`mace_off`/`mace_omol` rows with `mace_ml` dipole marked GREEN (top performers); `mace_mp` rows YELLOW (outlier)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Installed python-pptx into system Python**
- **Found during:** Task 2 verification
- **Issue:** `.venv` Python lacked `pptx` module; `uv add` blocked by dgl wheel incompatibility on this platform
- **Fix:** Ran `pip3 install python-pptx` into system Python 3.8 for verification purposes only
- **Files modified:** None (system package install)
- **Verification:** `python3 presentation/generate_pptx_v2.py` exits 0; slide count confirmed 11
- **Committed in:** not committed (system install, not project code)

---

**Total deviations:** 1 auto-fixed (1 blocking — environment gap)
**Impact on plan:** No scope creep. python-pptx already listed as a dependency; gap was the venv vs system Python distinction on this machine.

## Issues Encountered

- `.venv` did not have `python-pptx` installed; `uv add python-pptx` failed because `dgl==2.2.1` has no wheel for `manylinux_2_31_x86_64`. Resolved by installing into system Python for verification. The project dependency itself is unchanged.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `generate_pptx_v2.py` has a clean 11-slide structure with correct function names in `main()`
- Plan 02 can rewrite individual slide functions (motivation, IR theory, architecture, etc.) without touching the structural wiring
- No blockers

---
*Phase: 01-narrative-structure*
*Completed: 2026-03-07*
