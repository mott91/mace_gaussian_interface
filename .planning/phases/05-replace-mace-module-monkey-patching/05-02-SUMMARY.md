---
phase: 05-replace-mace-module-monkey-patching
plan: 02
subsystem: calculators
tags: [mace, cleanup, monkey-patching, pickle-module, migration]

# Dependency graph
requires:
  - phase: 05-01
    provides: "Safe MACE dipole loading via pickle_module remapping in calculators/mace_loader.py"
provides:
  - "Complete removal of mace_calculators.py and all monkey-patching references"
  - "Updated documentation reflecting new safe loading approach"
  - "Clean pyproject.toml and test configuration for calculators package"
affects: [future phases, documentation, new contributors]

# Tech tracking
tech-stack:
  added: []
  patterns: ["pickle_module remapping replaces sys.modules monkey-patching"]

key-files:
  created: []
  modified:
    - tests/test_calculators.py
    - pyproject.toml
    - CLAUDE.md
    - docs/ARCHITECTURE.md

key-decisions:
  - "Replaced mace_calculators mock with full mace_dipole_core hierarchy in test pre-mock list for reliable import resolution"

patterns-established:
  - "Mock mace_dipole_core.* hierarchy (not mace_calculators) when testing calculators package"

requirements-completed: [STRUCT-03]

# Metrics
duration: 2min
completed: 2026-02-19
---

# Phase 05 Plan 02: Delete mace_calculators.py Summary

**Removed legacy monkey-patching module and updated all references across docs, config, and tests to reflect pickle_module remapping approach**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-19T16:50:51Z
- **Completed:** 2026-02-19T16:52:44Z
- **Tasks:** 1
- **Files modified:** 5

## Accomplishments
- Deleted mace_calculators.py (106 lines of monkey-patching code removed)
- Updated test pre-mock list from single mace_calculators entry to full mace_dipole_core hierarchy
- Updated pyproject.toml known-first-party and coverage source to reference calculators package
- Updated CLAUDE.md and docs/ARCHITECTURE.md to describe the new safe loading mechanism

## Task Commits

Each task was committed atomically:

1. **Task 1: Delete mace_calculators.py and update all references** - `ba9033d` (chore)

## Files Created/Modified
- `mace_calculators.py` - Deleted (was the legacy monkey-patching module)
- `tests/test_calculators.py` - Updated mock targets from mace_calculators to mace_dipole_core hierarchy
- `pyproject.toml` - Updated known-first-party and coverage source references
- `CLAUDE.md` - Replaced monkey-patching gotcha with pickle_module remapping description
- `docs/ARCHITECTURE.md` - Updated Dual MACE Package Architecture description

## Decisions Made
- Replaced single `mace_calculators` mock with full `mace_dipole_core.*` hierarchy (5 entries) in test pre-mock list for reliable import resolution with Python's import machinery

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 05 (Replace MACE Module Monkey-Patching) is fully complete
- All monkey-patching code removed, safe pickle_module approach in place
- 128 tests passing, no regressions
- Ready for next phase in roadmap

---
*Phase: 05-replace-mace-module-monkey-patching*
*Completed: 2026-02-19*
