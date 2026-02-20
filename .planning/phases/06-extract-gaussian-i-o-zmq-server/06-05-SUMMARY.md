---
phase: 06-extract-gaussian-i-o-zmq-server
plan: "05"
subsystem: refactoring
tags: [gaussian, zmq, import-cleanup, package-structure]

# Dependency graph
requires:
  - phase: 06-01
    provides: gaussian/io.py with DEFAULT_HELPER_SCRIPT, ase_to_gjf, parse_gaussian_input, write_gaussian_output
  - phase: 06-02
    provides: gaussian/parser.py and gaussian/fchk.py (verbatim copies)
  - phase: 06-03
    provides: gaussian/zmq_server.py with GaussianZMQServer
  - phase: 06-04
    provides: gaussian/runner.py with run_gaussian_with_zmq and DEFAULT_TIMEOUT_SECONDS
provides:
  - gaussian/__init__.py with 13 public re-exports completing the package API
  - All callers (gm_main, dft_baseline, mode_matching, convert_all_chk_files) using gaussian.* paths
  - gm_main.run_frequency_calculation() using run_gaussian_with_zmq from gaussian.runner
  - Deleted gaussian_parser.py and fchk_parser.py (now authoritative in gaussian/ package)
  - All test files importing from gaussian.* paths
affects: [future phases using gaussian package, any caller of gm_main, dft_baseline, mode_matching]

# Tech tracking
tech-stack:
  added: []
  patterns: [clean break import updates - no backwards compat shims, gaussian package as single authoritative source]

key-files:
  created: []
  modified:
    - gaussian/__init__.py
    - pyproject.toml
    - gm_main.py
    - dft_baseline.py
    - mode_matching.py
    - convert_all_chk_files.py
    - tests/test_gaussian_parser.py
    - tests/test_fchk_parser.py
    - tests/test_mode_matching.py
    - tests/test_cli_validation.py

key-decisions:
  - "ruff auto-fix removed unused ase.data.chemical_symbols import from gm_main.py (parse_gaussian_input moved to gaussian.io)"
  - "Inline from fchk_parser import convert_chk_to_fchk in run_frequency_calculation() removed since convert_chk_to_fchk is now hoisted to top-level imports"
  - "test_cli_validation.py timeout tests updated to use DEFAULT_TIMEOUT_SECONDS from gaussian.runner instead of GAUSSIAN_TIMEOUT_SECONDS from gm_main"

patterns-established:
  - "gaussian/ package is now the single authoritative source for all Gaussian I/O, parsing, ZMQ, and runner code"
  - "No backwards-compat shims - clean break with direct import updates per project pattern"

requirements-completed: [STRUCT-04, STRUCT-05, STRUCT-07]

# Metrics
duration: 4min
completed: 2026-02-20
---

# Phase 6 Plan 05: Final Wiring - gaussian Package Complete Summary

**gaussian/ package fully wired as single authoritative source: all callers updated to gaussian.* imports, run_frequency_calculation() uses run_gaussian_with_zmq(), gaussian_parser.py and fchk_parser.py deleted**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-20T09:15:28Z
- **Completed:** 2026-02-20T09:19:14Z
- **Tasks:** 2
- **Files modified:** 10 (+ 2 deleted)

## Accomplishments

- gaussian/__init__.py now exposes all 13 public symbols from the package (GaussianLogParser, GaussianZMQServer, convert_chk_to_fchk, ase_to_gjf, DEFAULT_TIMEOUT_SECONDS, extract_modes_from_fchk, get_fchk_from_chk, is_calc_finished, parse_fchk_section, parse_gaussian_input, parse_gaussian_log, run_gaussian_with_zmq, write_gaussian_output)
- gm_main.py refactored: zmq_server(), is_calc_finished(), parse_gaussian_input(), write_gaussian_output(), ase_to_gjf(), DEFAULT_HELPER_SCRIPT, GAUSSIAN_TIMEOUT_SECONDS all removed; run_frequency_calculation() now delegates to run_gaussian_with_zmq from gaussian.runner
- gaussian_parser.py and fchk_parser.py deleted; gaussian/ package is now the sole source
- pyproject.toml updated: known-first-party replaces gaussian_parser/fchk_parser with gaussian; coverage source updated accordingly

## Task Commits

Each task was committed atomically:

1. **Task 1: Update gaussian/__init__.py with public re-exports and update pyproject.toml** - `7327a0c` (feat)
2. **Task 2: Update all callers, refactor gm_main.run_frequency_calculation, delete old files** - `14c2ded` (feat)

**Plan metadata:** (docs commit below)

## Files Created/Modified

- `gaussian/__init__.py` - Full public re-exports for all 13 symbols in the gaussian package
- `pyproject.toml` - Updated isort known-first-party and coverage source
- `gm_main.py` - Removed 7 function defs + 2 constants; added gaussian.* imports; refactored run_frequency_calculation() to use run_gaussian_with_zmq
- `dft_baseline.py` - gaussian_parser -> gaussian.parser, fchk_parser -> gaussian.fchk
- `mode_matching.py` - fchk_parser -> gaussian.fchk
- `convert_all_chk_files.py` - fchk_parser -> gaussian.fchk
- `tests/test_gaussian_parser.py` - gaussian_parser -> gaussian.parser
- `tests/test_fchk_parser.py` - fchk_parser -> gaussian.fchk
- `tests/test_mode_matching.py` - fchk_parser -> gaussian.fchk (2 inline imports)
- `tests/test_cli_validation.py` - gm_main.GAUSSIAN_TIMEOUT_SECONDS -> gaussian.runner.DEFAULT_TIMEOUT_SECONDS
- `gaussian_parser.py` - DELETED
- `fchk_parser.py` - DELETED

## Decisions Made

- ruff auto-fix (run as part of verification) removed unused `from ase.data import chemical_symbols` from gm_main.py since parse_gaussian_input was moved to gaussian.io — no behavior change, import was dead code
- Inline `from fchk_parser import convert_chk_to_fchk` inside run_frequency_calculation was removed (no longer needed — convert_chk_to_fchk is now a top-level import from gaussian.fchk)
- test_cli_validation.py timeout tests updated to import DEFAULT_TIMEOUT_SECONDS from gaussian.runner; comments updated from "gm_main has heavy deps" to "gaussian.runner may have deps" to reflect the new import path

## Deviations from Plan

None - plan executed exactly as written. The ruff auto-fix of unused chemical_symbols import is a natural consequence of moving parse_gaussian_input out of gm_main, and ruff --fix is prescribed in the plan's verification step.

## Issues Encountered

- `tests/test_cli_validation.py` could not be run in isolation due to missing `click` module in the test environment (pre-existing issue, not caused by this plan). All other 42 tests in the refactored test files pass cleanly.
- Pre-existing scipy.optimize `ObjSense` registration conflict affects test_regression and test_validation tests (unrelated to this plan's changes).

## Next Phase Readiness

Phase 6 is now complete:
- gaussian/ package is the single authoritative source for all Gaussian I/O, parsing, ZMQ, and runner code
- All top-level files (gaussian_parser.py, fchk_parser.py) deleted
- All callers use gaussian.* imports
- 42 tests pass cleanly

Ready for Phase 7.

## Self-Check: PASSED

- SUMMARY.md: FOUND
- gaussian/__init__.py: FOUND
- gaussian_parser.py deleted: CONFIRMED
- fchk_parser.py deleted: CONFIRMED
- Commit 7327a0c: FOUND
- Commit 14c2ded: FOUND

---
*Phase: 06-extract-gaussian-i-o-zmq-server*
*Completed: 2026-02-20*
