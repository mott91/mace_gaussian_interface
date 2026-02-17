# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** Phase 4 - Extract Calculator Classes (in progress)

## Current Position

Phase: 4 of 10 (Extract Calculator Classes) -- COMPLETE
Plan: 2 of 2 in current phase
Status: Phase 04 complete, ready for Phase 05
Last activity: 2026-02-17 — Completed 04-02: Calculator unit tests with mocked dependencies

Progress: [████░░░░░░] 40%

## Performance Metrics

**Velocity:**
- Total plans completed: 11
- Average duration: 4.5 min
- Total execution time: 0.82 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | 22 min | 5.5 min |
| 02 | 3 | 13 min | 4.3 min |
| 03 | 2 | 9 min | 4.5 min |
| 04 | 2 | 7 min | 3.5 min |

**Recent Trend:**
- Last 5 plans: 02-03 (5 min), 03-01 (6 min), 03-02 (3 min), 04-01 (4 min), 04-02 (3 min)
- Trend: Consistent fast execution on well-scoped plans

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Focus on refactoring before new features (code quality and reproducibility are prerequisites for distribution)
- Keep all 4 ML model combinations (comparison across combinations is part of thesis methodology)
- Use quality model profile for GSD agents (important project, want thorough analysis and planning)
- Reproducibility over pip-installability (thesis needs clone-and-reproduce, not package distribution yet)
- [Phase 01]: Moved dev dependencies from orphaned [project] key to [dependency-groups] section for proper uv lockfile inclusion
- [Phase 01-03]: Vib-E2 section has 14 values per mode (freq + thermodynamic properties), not just frequencies; test first N entries
- [Phase 01-03]: Import only pure numerical functions from mode_matching to avoid formchk dependency in tests
- [Phase 01-04]: Used range-based frequency assertions (not exact values) for regression test robustness
- [Phase 01-04]: CH4 harmonic has 4 entries (collapsed degenerate modes), not 9
- [Phase 02-02]: Created exceptions.py and validation.py inline as Plan 01 prerequisites (Rule 3 blocking)
- [Phase 02-02]: Used lazy import (try/except ImportError) for validation module in results_manager
- [Phase 02-02]: Default strict=False preserves backward compatibility for anharmonic-type parsers
- [Phase 02-03]: Timeout tests gracefully handle gm_main import failures (heavy deps like DGL/espaloma)
- [Phase 02-03]: Validation only on run command, not list/diagnose/compare/export
- [Phase 03-01]: Used CODATA 2018 full precision constants (0.529177210903, 27.211386245988)
- [Phase 03-01]: Preserved original numerical behavior when replacing misnamed constants (ANGSTROM_TO_BOHR=0.529 -> BOHR_TO_ANGSTROM)
- [Phase 03-01]: No compatibility shims -- clean break with direct import updates across codebase
- [Phase 03-02]: Corrected angstrom_to_bohr test expected value to match CODATA 2018 derived inverse
- [Phase 04-01]: Lazy import of MACEDipoleCalculator inside _check_availability to avoid module-level side effects
- [Phase 04-01]: Modernized type annotations (tuple/dict) in new calculator files for Python 3.12
- [Phase 04-02]: Pre-mock heavy deps via sys.modules before importing calculators package to avoid DGL/espaloma/xtb side effects in tests
- [Phase 04-02]: Added from __future__ import annotations to calculators/base.py, factory.py, gm_main.py for Python 3.8 compatibility

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 1 considerations:**
- Acetic acid parser bug (commit a4384c4) needs test case to capture expected vs actual behavior
- Coverage targets should be pragmatic (focus on parsers, mode matching, calculators — not necessarily 100%)

**Phase 5 known complexity:**
- MACE module loading requires deeper research on importlib.util patterns for isolated loading
- CUDA initialization state management needs careful handling
- Rollback strategy needed if lazy imports fail

## Session Continuity

Last session: 2026-02-17 (phase execution)
Stopped at: Completed 04-02-PLAN.md (Phase 04 complete)
Resume file: None
