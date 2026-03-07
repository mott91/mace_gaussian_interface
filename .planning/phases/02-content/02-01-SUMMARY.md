---
phase: 02-content
plan: 01
subsystem: presentation
tags: [pptx, python-pptx, data-visualization, json, ir-spectroscopy]

# Dependency graph
requires:
  - phase: 01-narrative-structure
    provides: generate_pptx_v2.py with correct slide order and narrative arc

provides:
  - JSON-sourced results tables for water (8 combos) and aspirin (4 combos) sorted by MAE ascending
  - slide_scaling() with ASCII bar chart showing water/glycine/aspirin speedups
  - Rewritten slide_zmq() using single-column content_slide() with zmq_server.py
  - Cleaned slide_architecture() with no Python function name annotations
  - load_combo_metrics() and parse_combo_name() data helpers
  - 13-slide presentation_v2.pptx with all content gaps filled

affects: [03-polish, presentation delivery]

# Tech tracking
tech-stack:
  added: [json (stdlib), os (stdlib)]
  patterns:
    - Load metrics from analysis_results_harmonic/{molecule}/data/metrics_summary.json
    - Sort combos by mae_freq ascending (lowest error = best = shown first)
    - Uniform TEXT color for all data rows (color-coding deferred to Phase 3)
    - content_slide() for all pedagogical narrative slides (not two_col_slide())

key-files:
  created: []
  modified:
    - presentation/generate_pptx_v2.py

key-decisions:
  - "Results table sorted by MAE ascending (not R² descending as originally discussed)"
  - "Data loaded from analysis_results_harmonic/ metrics_summary.json (not comparison_results/ directory walk)"
  - "Uniform TEXT color for data rows — per-row color-coding deferred to Phase 3 polish"
  - "slide_zmq() uses content_slide() single-column layout (not two_col_slide()) for cleaner narrative"
  - "slide_results_table() produces two slides (water + aspirin) adding one slide to total"

patterns-established:
  - "load_combo_metrics(molecule): canonical way to load sorted combo metrics from JSON"
  - "parse_combo_name(name): splits combo name like mace_omol_mace_ml into (energy, dipole) tuple"

requirements-completed: [CONT-01, CONT-02, CONT-03, CONT-04]

# Metrics
duration: 4min
completed: 2026-03-07
---

# Phase 2 Plan 01: Content Fill Summary

**13-slide PPTX generator with JSON-sourced water/aspirin results tables, ASCII speedup bar chart, and rewritten ZMQ slide using zmq_server.py in single-column layout**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-07T12:10:37Z
- **Completed:** 2026-03-07T12:14:39Z
- **Tasks:** 4 (Tasks 1, 2a, 2b, 3)
- **Files modified:** 1

## Accomplishments

- Replaced hardcoded water-only results table with JSON-sourced two-slide table (water 8 combos + aspirin 4 combos), sorted by MAE ascending
- Added slide_scaling() with ASCII bar chart explicitly showing the O(N) vs O(N³) speedup scaling argument
- Rewrote slide_zmq() from two_col_slide() to content_slide() with four-beat narrative using zmq_server.py (not stale gm_main.py)
- Stripped Python function name annotations from slide_architecture() boxes, replaced with component-level descriptions
- Added load_combo_metrics() and parse_combo_name() helpers for clean data loading

## Task Commits

Each task was committed atomically:

1. **Task 1: Rewrite slide_zmq() and clean slide_architecture()** - `e3607b8` (feat)
2. **Task 2a: Add load_combo_metrics() and parse_combo_name() helpers** - `1a08e72` (feat)
3. **Task 2b: Replace slide_results_table() and update slide_results_overview()** - `22664a3` (feat)
4. **Task 3: Add slide_scaling() and wire main()** - `a6c59a2` (feat)

**Plan metadata:** (docs commit — see below)

## Files Created/Modified

- `presentation/generate_pptx_v2.py` - Added helpers, rewrote 4 slide functions, added slide_scaling(), updated main() for 13-slide sequence

## Decisions Made

- Results table sorted by MAE ascending (lowest error = best shown first), per CONTEXT.md — overrides original R² descending from discuss-phase
- Data loaded from `analysis_results_harmonic/{molecule}/data/metrics_summary.json` — overrides original comparison_results/ directory walk approach
- Uniform TEXT color for all data rows — per-row color-coding (green/yellow/red thresholds) deferred to Phase 3 polish
- slide_zmq() uses single-column content_slide() for a cleaner pedagogical flow
- slide_results_table() now generates 2 slides (water + aspirin), bumping total from 11 to 13 once scaling slide added

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- The project .venv does not have python-pptx installed; python-pptx is available in the `mace4ir_v2` conda env. Used `/home/mot/micromamba/envs/mace4ir_v2/bin/python` for verification commands. This is a pre-existing environment configuration issue, not introduced by this plan.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All 6 automated checks pass: 13-slide count, zmq_server.py present, gm_main.py absent, water/report.html, aspirin/report.html, benchmark.py --scaling
- Phase 3 (Polish & Speaking Notes) can proceed: color-code data rows green/yellow/red by threshold, add slide numbers, add speaking notes

---
*Phase: 02-content*
*Completed: 2026-03-07*
