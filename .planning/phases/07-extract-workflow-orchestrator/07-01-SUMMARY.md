---
phase: 07-extract-workflow-orchestrator
plan: 01
subsystem: pipeline
tags: [workflow, orchestrator, pipeline, mace, gaussian, zmq, ase]

# Dependency graph
requires:
  - phase: 06-extract-gaussian-io
    provides: gaussian.* package (runner, io, parser, fchk, zmq_server)
provides:
  - workflow.py with run_pipeline(), run_geometry_optimization(), run_dft_baselines(), run_ml_calculations(), run_frequency_calculation()
  - All low-level ZMQ callback helpers (update_molecule_geometry, calculate_energy_and_forces, calculate_hessian, calculate_dipole_properties, run_next_calculation)
  - geometry_optimisation() and calculator() helpers in workflow.py
affects:
  - 07-extract-workflow-orchestrator (plans 02+: cli.py wiring, gm_main.py deletion)
  - cli.py will call run_pipeline() directly

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Lazy import of dft_baseline inside run_dft_baselines() body to avoid heavy DGL/espaloma at import time"
    - "Three named stage functions (run_geometry_optimization, run_dft_baselines, run_ml_calculations) plus top-level run_pipeline() coordinator"
    - "Cluster seam pattern: stages 2 and 3 depend only on stage 1 geometry output, enabling independent execution"

key-files:
  created:
    - workflow.py
  modified: []

key-decisions:
  - "Optional[np.ndarray] replaced with np.ndarray | None (modern Python type syntax, ruff UP045)"
  - "detect_device() called without assignment (return value unused — side effect is logging only)"
  - "json.load() result discarded in opt_file cache-hit path (opt_metadata unused, charge/spin set to defaults)"
  - "Dead code excluded: setup_output_directory (no callers), analyze_molecular_charges (not in pipeline path), print_diagnostics (CLI concern), CHARGE_ANALYSIS_AVAILABLE (dead import block)"
  - "dft_baseline lazy import stays inside run_dft_baselines() body — prevents DGL/espaloma side effects at module load"

patterns-established:
  - "Stage function pattern: each stage is a named public function callable independently or via run_pipeline()"
  - "Lazy heavy-dep import: dft_baseline (and by extension DGL/espaloma) only loaded when stage 2 actually runs"

requirements-completed: [STRUCT-06]

# Metrics
duration: 3min
completed: 2026-02-24
---

# Phase 7 Plan 1: Extract Workflow Orchestrator Summary

**workflow.py created with 718 lines: run_pipeline() replacing run_workflow() monolith, three named stage functions, all ZMQ callback helpers, lazy dft_baseline import, dead code excluded**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-02-24T16:57:59Z
- **Completed:** 2026-02-24T17:00:30Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Created workflow.py (718 lines) as single authoritative pipeline orchestrator
- Extracted run_workflow() monolith into four named stage functions: run_geometry_optimization, run_dft_baselines, run_ml_calculations, plus run_pipeline() coordinator
- Preserved all 7 low-level helpers verbatim from gm_main.py (update_molecule_geometry, calculate_energy_and_forces, calculate_hessian, calculate_dipole_properties, run_next_calculation, geometry_optimisation, calculator)
- Kept dft_baseline import lazy inside run_dft_baselines() — module imports cleanly without DGL/espaloma
- Excluded all dead code: setup_output_directory, analyze_molecular_charges, CHARGE_ANALYSIS_AVAILABLE, print_diagnostics, __main__ block

## Task Commits

1. **Task 1: Create workflow.py with all extracted pipeline functions** - `936bdb5` (feat)

**Plan metadata:** (docs commit — see below)

## Files Created/Modified
- `/home/mot/mace_gaussian/workflow.py` - Complete pipeline orchestrator with stage functions, coordinator, and ZMQ callback helpers (718 lines)

## Decisions Made
- Replaced `Optional[np.ndarray]` with `np.ndarray | None` to satisfy ruff UP045 modern Python annotation rule
- `detect_device()` called for its side effect (logging) without storing return value — `device` variable was unused
- `json.load()` result discarded in existing-geometry cache path — metadata was loaded but charge/spin always set to defaults (0.0/1.0)
- Dead code (`setup_output_directory`, `analyze_molecular_charges`, `print_diagnostics`) excluded as specified — `setup_output_directory` had no callers, `analyze_molecular_charges` is not in the pipeline path, `print_diagnostics` belongs in CLI layer

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed ruff lint errors: Optional -> X|None, unused variables**
- **Found during:** Task 1 (after initial write, pre-commit ruff check)
- **Issue:** `Optional[np.ndarray]` triggered UP045; `device = detect_device()` and `opt_metadata = json.load(f)` triggered F841 unused variable warnings
- **Fix:** Replaced `Optional[np.ndarray]` with `np.ndarray | None` in two function signatures; removed `device =` assignment; removed `opt_metadata =` assignment (json.load result discarded)
- **Files modified:** workflow.py
- **Verification:** `ruff check workflow.py` passes with "All checks passed!"
- **Committed in:** 936bdb5 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - ruff lint compliance)
**Impact on plan:** Necessary for ruff clean requirement. No scope creep.

## Issues Encountered
None beyond the ruff fixes above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- workflow.py is ready for cli.py to import run_pipeline() in plan 07-02
- gm_main.py still exists (deletion deferred to later plan per phase design)
- All stage functions testable independently via their public interfaces

---
*Phase: 07-extract-workflow-orchestrator*
*Completed: 2026-02-24*

## Self-Check: PASSED

- FOUND: /home/mot/mace_gaussian/workflow.py
- FOUND: 936bdb5
