---
phase: 08-package-structure-reorganization
plan: 02
subsystem: infra
tags: [python, packaging, imports, refactor, analysis]

# Dependency graph
requires:
  - phase: 08-01
    provides: mace_gaussian/ package created via git mv with relative imports on all subpackages
provides:
  - mace_gaussian/cli.py with absolute mace_gaussian.* imports
  - mace_gaussian/workflow.py with relative imports (formerly root workflow.py)
  - mace_gaussian/gm_helper.py (ZMQ bridge, stdlib-only, safe for Gaussian standalone invocation)
  - mace_gaussian/dft_baseline.py with relative imports
  - mace_gaussian/diagnostics.py (stdlib-only)
  - mace_gaussian/analysis/ subpackage with analyze_spectra, mode_matching, html_report_generator, analysis_workflow
  - mace_gaussian/analysis/__init__.py exporting analyze_molecule, analyze_molecule_harmonic, run_analysis_main, run_analysis_harmonic_main
  - run_analysis.py and run_analysis_harmonic.py as 3-line shims (no sys.path)
  - scripts/ containing produce_molecules.py, convert_all_chk_files.py, charge_analysis.py, comparison_workflow.py
  - mace_gaussian/__init__.py re-exporting run_pipeline from .workflow
affects:
  - 08-03 (cleanup and test verification)
  - any future package distribution or CLI entry points

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Package submodule pattern: library code in mace_gaussian/, root has only thin shims"
    - "Relative imports inside package, absolute mace_gaussian.* imports at entry boundaries (cli.py, shims)"
    - "analysis/ subpackage groups analysis+reporting modules with __init__.py public API"
    - "Cross-subpackage relative: analysis/mode_matching uses ..gaussian.fchk (parent-relative)"
    - "Entry functions co-located in analysis_workflow.py (run_analysis_main, run_analysis_harmonic_main)"

key-files:
  created:
    - mace_gaussian/analysis/__init__.py
    - scripts/comparison_workflow.py
  modified:
    - mace_gaussian/cli.py
    - mace_gaussian/workflow.py
    - mace_gaussian/gm_helper.py
    - mace_gaussian/dft_baseline.py
    - mace_gaussian/diagnostics.py
    - mace_gaussian/analysis/analyze_spectra.py
    - mace_gaussian/analysis/mode_matching.py
    - mace_gaussian/analysis/html_report_generator.py
    - mace_gaussian/analysis/analysis_workflow.py
    - mace_gaussian/__init__.py
    - run_analysis.py
    - run_analysis_harmonic.py
    - scripts/produce_molecules.py
    - scripts/convert_all_chk_files.py
    - scripts/charge_analysis.py

key-decisions:
  - "cli.py uses absolute mace_gaussian.* imports (not relative) — it is the entry boundary, must be importable standalone"
  - "gm_helper.py uses only stdlib + zmq (no internal imports) — runs as standalone script invoked by Gaussian, relative imports would fail"
  - "run_analysis_main() and run_analysis_harmonic_main() added to analysis_workflow.py — shims delegate to package, not inline"
  - "seaborn and click missing from venv are pre-existing env issues, not caused by this plan"
  - "dft_baseline.py pre-existing ruff issues (typing.Dict, line length) are out of scope per SCOPE BOUNDARY rule"

patterns-established:
  - "All library source files live inside mace_gaussian/ package — no stragglers at root"
  - "Root shims (run_analysis.py, run_analysis_harmonic.py) are 3 lines each, no sys.path manipulation"

requirements-completed: [STRUCT-08, STRUCT-10]

# Metrics
duration: 4min
completed: 2026-02-25
---

# Phase 8 Plan 02: Package Consolidation with analysis/ Subpackage Summary

**All pipeline source files moved into mace_gaussian/, analysis/ subpackage created with relative imports, root shims rewritten as 3-line delegators**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-25T~16:07Z
- **Completed:** 2026-02-25T~16:11Z
- **Tasks:** 2
- **Files modified:** 15

## Accomplishments
- All flat root source files (cli.py, workflow.py, gm_helper.py, dft_baseline.py, diagnostics.py) moved into mace_gaussian/ via git mv (history preserved)
- mace_gaussian/analysis/ subpackage created with all four analysis modules (analyze_spectra, mode_matching, html_report_generator, analysis_workflow)
- comparison_workflow.py renamed to analysis_workflow.py during move; entry functions added
- mace_gaussian/__init__.py now re-exports run_pipeline from .workflow (activated after 08-01)
- root run_analysis.py and run_analysis_harmonic.py rewritten as 3-line shims with no sys.path manipulation
- scripts/ populated with produce_molecules.py, convert_all_chk_files.py, charge_analysis.py, comparison_workflow.py (thin delegator)

## Task Commits

Each task was committed atomically:

1. **Task 1: Move flat root source files into mace_gaussian/ and create analysis/ subpackage** - `c277278` (feat)
2. **Task 2: Update all imports in moved files and rewrite root shims** - `7e033ee` (feat)

**Plan metadata:** (docs commit created next)

## Files Created/Modified
- `mace_gaussian/analysis/__init__.py` - Analysis subpackage public API (analyze_molecule, analyze_molecule_harmonic, run_analysis_main, run_analysis_harmonic_main)
- `mace_gaussian/analysis/analysis_workflow.py` - Renamed from comparison_workflow.py; imports updated to relative; run_analysis_main() and run_analysis_harmonic_main() added
- `mace_gaussian/analysis/mode_matching.py` - Cross-subpackage import updated: gaussian.fchk → ..gaussian.fchk
- `mace_gaussian/workflow.py` - All imports converted to relative (.calculators, .gaussian.*, .utils.*, .dft_baseline)
- `mace_gaussian/cli.py` - All imports converted to absolute mace_gaussian.* (entry boundary)
- `mace_gaussian/dft_baseline.py` - Imports converted to relative (.gaussian.parser, .utils.*, lazy .gaussian.fchk)
- `mace_gaussian/__init__.py` - Activated run_pipeline re-export from .workflow
- `run_analysis.py` - Rewritten as 3-line shim delegating to run_analysis_main
- `run_analysis_harmonic.py` - Rewritten as 3-line shim delegating to run_analysis_harmonic_main
- `scripts/comparison_workflow.py` - New thin delegator for researchers

## Decisions Made
- cli.py uses absolute mace_gaussian.* imports (not relative) because it is the CLI entry boundary — must work when invoked as `python mace_gaussian/cli.py` or via entry_points
- gm_helper.py is stdlib+zmq only with no internal imports — safe for Gaussian's standalone subprocess invocation where relative imports would fail
- Entry point functions (run_analysis_main, run_analysis_harmonic_main) added to analysis_workflow.py so root shims can be 3 lines without inlining logic
- comparison_workflow.py renamed to analysis_workflow.py during the git mv to match its new role as the analysis package orchestrator

## Deviations from Plan

None — plan executed exactly as written. The seaborn/click missing-from-venv issue is a pre-existing environment gap not introduced by this plan.

## Issues Encountered

- seaborn not installed in venv: `from mace_gaussian.analysis import analyze_molecule` fails at module load time because analyze_spectra.py imports seaborn at the top level. This is pre-existing (comparison_workflow.py had the same dependency). Import structure is confirmed correct via AST parsing and pattern-matching.
- click not installed in venv: cli import can't be tested at runtime. Same pre-existing status.

## Next Phase Readiness
- All library code is now inside mace_gaussian/ package
- Relative import chain is verified via pattern checks and AST parsing
- mace_gaussian.__version__ = "0.2.0" and run_pipeline import confirmed working
- Phase 08-03 can proceed with tests and final cleanup

## Self-Check: PASSED

All 17 expected files confirmed present. Both task commits (c277278, 7e033ee) confirmed in git log.

---
*Phase: 08-package-structure-reorganization*
*Completed: 2026-02-25*
