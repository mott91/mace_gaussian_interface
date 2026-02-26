---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: unknown
last_updated: "2026-02-26T16:12:33.622Z"
progress:
  total_phases: 10
  completed_phases: 10
  total_plans: 30
  completed_plans: 30
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** Phase 10 - Documentation

## Current Position

Phase: 10 of 10 (Documentation)
Plan: 4 of 4 in current phase -- COMPLETE
Status: Phase 10 COMPLETE — 4/4 plans done
Last activity: 2026-02-26 — Completed 10-04: compare/export runtime body fix (COMING SOON banners replaced with honest 'Not yet implemented' messages, DOC-04 satisfied)

Progress: [██████████] 100%

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
| 07 | 2 | 8 min | 4.0 min |
| 08 | 1+ | 2 min | 2.0 min |

**Recent Trend:**
- Last 5 plans: 06-04 (2 min), 06-05 (4 min), 07-01 (3 min), 07-02 (5 min), 08-01 (2 min)
- Trend: Consistent fast execution on well-scoped plans

*Updated after each plan completion*
| Phase 06 P04 | 2 | 1 tasks | 1 files |
| Phase 06 P05 | 4 | 2 tasks | 12 files |
| Phase 07 P01 | 3 | 1 tasks | 1 files |
| Phase 07 P02 | 5 | 2 tasks | 3 files |
| Phase 08 P01 | 2 | 2 tasks | 13 files |
| Phase 08 P02 | 4 | 2 tasks | 15 files |
| Phase 08 P03 | 7 | 2 tasks | 13 files |
| Phase 09 P01 | 90 | 2 tasks | 23 files |
| Phase 09 P02 | 5 | 1 tasks | 1 files |
| Phase 09 P03 | 2 | 1 tasks | 1 files |
| Phase 10 P02 | 2 | 2 tasks | 8 files |
| Phase 10 P01 | 2 | 2 tasks | 2 files |
| Phase 10 P03 | 2 | 1 tasks | 1 files |
| Phase 10 P04 | 1 | 2 tasks | 1 files |

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
- [Phase 07-02]: No try/except around 'from workflow import run_pipeline' — workflow.py internal imports don't have heavy-dep failure risk unlike gm_main
- [Phase 07-02]: testingStuff/test_refactoring.py gm_main import left as-is — legacy scratch, acceptable dead reference per plan spec
- [Phase 08-01]: workflow.run_pipeline import commented out in mace_gaussian/__init__.py — workflow.py not yet inside package (Plan 02 activates it)
- [Phase 08-01]: ruff --select I,F401 --fix applied to auto-sort 2 import blocks after manual relative import conversion (espaloma.py, runner.py)
- [Phase 08-01]: mace_loader.py lazy import inside _check_availability converted to from .mace_loader — same lazy semantics, now relative
- [Phase 08-02]: cli.py uses absolute mace_gaussian.* imports (not relative) — entry boundary, must work standalone
- [Phase 08-02]: gm_helper.py uses only stdlib + zmq — runs as standalone subprocess invoked by Gaussian, relative imports would fail
- [Phase 08-02]: run_analysis_main() and run_analysis_harmonic_main() added to analysis_workflow.py — shims delegate to package
- [Phase 08-02]: comparison_workflow.py renamed to analysis_workflow.py during git mv to match its orchestrator role
- [Phase 08-03]: Lazy wrapper functions in analysis/__init__.py to defer seaborn/pandas import — avoids module-level side effects when importing lightweight submodules like mode_matching
- [Phase 08-03]: results.py lazy import of collect_version_metadata was silently failing (flat path), fixed to mace_gaussian.utils.validation
- [Phase 08-03]: click installed separately into venv — uv sync blocked by dgl 2.2.1 having no Linux wheels
- [Phase 09-01]: Use # noqa: F401 for intentional import-for-testing in diagnostics.py (not worth removing the test intent)
- [Phase 09-01]: water_dft_fchk fixture changed to return Path (not str) to enable .open() — all callers accept both str and Path
- [Phase 09-01]: B017 pytest.raises(Exception) replaced with specific exception type for improved test clarity
- [Phase 09-02]: pip install -e . --no-deps avoids dgl==2.2.1 Windows-only wheels that block uv sync on ubuntu-latest
- [Phase 09-02]: ruff==0.15.1 pinned to match uv.lock for CI reproducibility
- [Phase 09-02]: lint and test jobs run in parallel (no needs: field) for fastest CI feedback
- [Phase 09-02]: No --cov-fail-under threshold — coverage is informational only (locked decision)
- [Phase 09]: Manual pip install steps shown explicitly before mentioning install_mace_packages.sh convenience script
- [Phase 10-01]: compare and export docstrings reference run_analysis.py/run_analysis_harmonic.py as existing alternatives — avoids implying CLI will gain compare/export features in near term
- [Phase 10-01]: Quickstart section added before Usage section (not duplicating 5-step Installation from Phase 9)
- [Phase 10-03]: methods.md written in passive voice past tense throughout — ZMQ injection described as novel software contribution with four paragraphs; placeholder citations in [Author, Year] format for easy journal reformatting
- [Phase 10-02]: quickstart.md references README#installation rather than duplicating install steps to avoid drift
- [Phase 10-02]: expected_output/ uses JSON only (no .chk/.fchk/.log) — portable, system-independent, human-readable reference outputs
- [Phase 10-02]: Plots sourced from anharmonic analysis_results/ to match the quickstart workflow path
- [Phase 10-04]: String literals in compare() message split across lines to stay within 100-char ruff limit while preserving output text

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

**Phase 7 status:** COMPLETE
- 07-01 complete: workflow.py created with run_pipeline() + stage functions extracted from gm_main.py (718 lines, ruff clean)
- 07-02 complete: cli.py rewired to workflow.run_pipeline(), print_diagnostics inlined using dipole_factory, gm_main.py deleted — STRUCT-06 fully satisfied

**Phase 8 status:** COMPLETE
- 08-01 complete: mace_gaussian/ package created via git mv (17 renames, history preserved), all internal imports converted to relative dot-notation — STRUCT-08 satisfied
- 08-02 complete: all root files moved into mace_gaussian/, analysis/ subpackage created, root shims rewritten as 3-line delegators, run_pipeline re-export activated — STRUCT-08, STRUCT-10 satisfied
- 08-03 complete: pyproject.toml wired (mace-gaussian entry point, packages = ["mace_gaussian"]), all 128 tests updated to mace_gaussian.* imports and passing — STRUCT-08, STRUCT-09, STRUCT-10 satisfied

## Session Continuity

Last session: 2026-02-26 (phase execution)
Stopped at: Completed 10-04-PLAN.md (compare/export runtime body fix — Phase 10 all plans done)
Resume file: None
