# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** Phase 7 - Extract Workflow Orchestrator

## Current Position

Phase: 7 of 10 (Extract Workflow Orchestrator)
Plan: 1 of 4 in current phase -- COMPLETE
Status: Phase 07 in progress, plan 01 done
Last activity: 2026-02-24 — Completed 07-01: workflow.py created with run_pipeline() and stage functions extracted from gm_main.py

Progress: [████████░░] 72%

## Performance Metrics

**Velocity:**
- Total plans completed: 16
- Average duration: 3.9 min
- Total execution time: 0.97 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | 22 min | 5.5 min |
| 02 | 3 | 13 min | 4.3 min |
| 03 | 2 | 9 min | 4.5 min |
| 04 | 2 | 7 min | 3.5 min |
| 05 | 2 | 6 min | 3.0 min |
| 06 | 5 | 15 min | 3.0 min |
| 07 | 1+ | 3 min | 3.0 min |

**Recent Trend:**
- Last 5 plans: 06-02 (4 min), 06-03 (3 min), 06-04 (2 min), 06-05 (4 min), 07-01 (3 min)
- Trend: Consistent fast execution on well-scoped plans

*Updated after each plan completion*
| Phase 06 P04 | 2 | 1 tasks | 1 files |
| Phase 06 P05 | 4 | 2 tasks | 12 files |
| Phase 07 P01 | 3 | 1 tasks | 1 files |

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
- [Phase 06-04]: stdout=PIPE and stderr=PIPE always captured so both GaussianTimeoutError and GaussianRunError include full Gaussian output for diagnostics
- [Phase 06-04]: proc.kill() (SIGKILL) on timeout, not proc.terminate() (SIGTERM) — Gaussian ignores SIGTERM
- [Phase 06-04]: runner imports only from gaussian.zmq_server and utils.exceptions — clean dependency boundary, no gaussian.io or gaussian.parser
- [Phase 06]: ruff auto-fix removed unused ase.data.chemical_symbols from gm_main.py after parse_gaussian_input moved to gaussian.io
- [Phase 06]: test_cli_validation.py GAUSSIAN_TIMEOUT_SECONDS replaced with DEFAULT_TIMEOUT_SECONDS from gaussian.runner throughout
- [Phase 07-01]: Optional[np.ndarray] replaced with np.ndarray | None (modern Python type syntax, ruff UP045)
- [Phase 07-01]: detect_device() called for side effect (logging) without storing return value
- [Phase 07-01]: dft_baseline lazy import stays inside run_dft_baselines() body — prevents DGL/espaloma side effects at module load
- [Phase 07-01]: Dead code excluded: setup_output_directory (no callers), analyze_molecular_charges (not in pipeline path), print_diagnostics (CLI concern)

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 1 considerations:**
- Acetic acid parser bug (commit a4384c4) needs test case to capture expected vs actual behavior
- Coverage targets should be pragmatic (focus on parsers, mode matching, calculators — not necessarily 100%)

**Phase 5 status:** COMPLETE
- Safe loading mechanism (pickle_module remapping) implemented and tested in 05-01
- Legacy mace_calculators.py deleted, all references cleaned up in 05-02

**Phase 6 status:** COMPLETE
- 06-01 complete: gaussian/ package foundation (exceptions, io.py)
- 06-02 complete: gaussian/parser.py and gaussian/fchk.py created (verbatim copies)
- 06-03 complete: gaussian/zmq_server.py with GaussianZMQServer class (LINGER=0 fix, no placeholder)
- 06-04 complete: gaussian/runner.py with run_gaussian_with_zmq (SIGKILL timeout, stdout/stderr capture)
- 06-05 complete: Final wiring - all callers use gaussian.* imports, gaussian_parser.py/fchk_parser.py deleted, gaussian/ is single authoritative source

**Phase 7 status:** IN PROGRESS
- 07-01 complete: workflow.py created with run_pipeline() + stage functions extracted from gm_main.py (718 lines, ruff clean)

## Session Continuity

Last session: 2026-02-24 (phase execution)
Stopped at: Completed 07-01-PLAN.md
Resume file: None
