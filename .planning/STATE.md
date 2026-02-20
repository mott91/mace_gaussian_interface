# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** Phase 6 - Extract Gaussian I/O and ZMQ Server

## Current Position

Phase: 6 of 10 (Extract Gaussian I/O and ZMQ Server)
Plan: 3 of 5 in current phase -- COMPLETE
Status: Phase 06-03 complete, ready for Phase 06-04
Last activity: 2026-02-20 — Completed 06-03: Create gaussian/zmq_server.py with GaussianZMQServer class

Progress: [██████░░░░] 59%

## Performance Metrics

**Velocity:**
- Total plans completed: 14
- Average duration: 4.1 min
- Total execution time: 0.95 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | 22 min | 5.5 min |
| 02 | 3 | 13 min | 4.3 min |
| 03 | 2 | 9 min | 4.5 min |
| 04 | 2 | 7 min | 3.5 min |
| 05 | 2 | 6 min | 3.0 min |
| 06 | 3 | 9 min | 3.0 min |

**Recent Trend:**
- Last 5 plans: 05-01 (4 min), 05-02 (2 min), 06-01 (2 min), 06-02 (4 min), 06-03 (3 min)
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
- [Phase 05-01]: Used consistent mock package hierarchy in tests for reliable import resolution with Python's import machinery
- [Phase 05-01]: Falls-through test verifies real mace.modules.models class returned when dipole module lacks attribute
- [Phase 05-01]: Docstring avoids literal cleanup/sys.modules strings for source-grep test reliability
- [Phase 05-02]: Replaced mace_calculators mock with full mace_dipole_core hierarchy in test pre-mock list
- [Phase 06-01]: DEFAULT_HELPER_SCRIPT in gaussian/io.py uses Path(__file__).parent.parent to reach project root (one level up from gaussian/)
- [Phase 06-01]: gm_main.py callers not updated in plan 06-01 - import wiring deferred to plan 06-05 to minimize per-plan diff scope
- [Phase 06-01]: gaussian/__init__.py stays sparse stub until all submodules exist in later plans
- [Phase 06-02]: No behavioral changes — gaussian/parser.py and gaussian/fchk.py are verbatim copies with added docstrings; originals preserved until plan 05 cleans up callers
- [Phase 06-03]: LINGER=0 applied after socket() and before bind() to prevent socket.close() blocking forever on Gaussian crash (STRUCT-07)
- [Phase 06-03]: No open() placeholder file creation in __enter__ — socket.bind() creates IPC file itself; original zmq_server() had this bug
- [Phase 06-03]: __exit__ returns False (not suppress); nested try/finally ensures cleanup even if socket.close() raises

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 1 considerations:**
- Acetic acid parser bug (commit a4384c4) needs test case to capture expected vs actual behavior
- Coverage targets should be pragmatic (focus on parsers, mode matching, calculators — not necessarily 100%)

**Phase 5 status:** COMPLETE
- Safe loading mechanism (pickle_module remapping) implemented and tested in 05-01
- Legacy mace_calculators.py deleted, all references cleaned up in 05-02

**Phase 6 in progress:**
- 06-01 complete: gaussian/ package foundation (exceptions, io.py)
- 06-02 complete: gaussian/parser.py and gaussian/fchk.py created (verbatim copies)
- 06-03 complete: gaussian/zmq_server.py with GaussianZMQServer class (LINGER=0 fix, no placeholder)
- 06-04 next: next plan in phase

## Session Continuity

Last session: 2026-02-20 (phase execution)
Stopped at: Completed 06-03-PLAN.md
Resume file: None
