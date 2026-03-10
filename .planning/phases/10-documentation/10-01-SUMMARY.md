---
phase: 10-documentation
plan: 01
subsystem: cli
tags: [click, docstrings, help-text, readme]

# Dependency graph
requires:
  - phase: 09-ci-cd
    provides: README Installation section with 5-step procedure
provides:
  - Accurate CLI --help text for all five subcommands (no misleading placeholder content)
  - README link to docs/quickstart.md
affects: [docs/quickstart.md (Plan 02 will create this target)]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Honest placeholder pattern: 'Not yet implemented' with actual alternatives, not 'Coming soon' feature lists"

key-files:
  created: []
  modified:
    - mace_gaussian/cli.py
    - README.md

key-decisions:
  - "compare and export docstrings reference run_analysis.py/run_analysis_harmonic.py as existing alternatives — avoids implying the CLI will gain compare/export features in the near term"
  - "Quickstart link added before Usage section (not duplicating Installation) — quickstart.md to be created in Plan 02"

patterns-established:
  - "Unimplemented CLI commands say 'Not yet implemented' and point to existing script alternatives"

requirements-completed: [DOC-01, DOC-04]

# Metrics
duration: 2min
completed: 2026-02-26
---

# Phase 10 Plan 01: CLI Docstring Fixes and README Quickstart Link Summary

**Replaced misleading CLI placeholder docstrings with honest "Not yet implemented" text pointing to run_analysis.py alternatives, and added docs/quickstart.md link to README**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-26T15:55:51Z
- **Completed:** 2026-02-26T15:58:00Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- `compare --help` and `export --help` now say "Not yet implemented" with concrete alternatives (run_analysis.py, run_analysis_harmonic.py and JSON file locations)
- `diagnose --help` now explicitly lists Gaussian 16, formchk, and CUDA in its Checks section
- `run --help` now shows `--output-dir my_results` in the Examples block
- README.md `## Quickstart` section added linking to `docs/quickstart.md` (to be created in Plan 02)
- ruff passes on mace_gaussian/cli.py with line length 100

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix compare and export CLI docstrings** - `20bb2c0` (fix)
2. **Task 2: Verify README Installation and add quickstart link** - `0b796b7` (feat)

**Plan metadata:** (final docs commit follows)

## Files Created/Modified
- `mace_gaussian/cli.py` - Updated docstrings for compare, export, diagnose, and run commands
- `README.md` - Added ## Quickstart section with docs/quickstart.md link

## Decisions Made
- `compare` and `export` docstrings reference `run_analysis.py`/`run_analysis_harmonic.py` as the actual alternatives — makes the existing workflow discoverable without implying CLI will grow these features soon
- Quickstart section added before the existing "Usage" section (not duplicating the 5-step Installation procedure already written in Phase 9)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed two E501 line-too-long ruff violations in compare docstring**
- **Found during:** Task 1 (Fix compare and export CLI docstrings)
- **Issue:** Initial docstring text had two lines exceeding 100-character ruff limit (103 and 102 chars)
- **Fix:** Wrapped long comment text at column 100 within docstring strings
- **Files modified:** mace_gaussian/cli.py
- **Verification:** `ruff check mace_gaussian/cli.py` exits 0
- **Committed in:** 20bb2c0 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - line length)
**Impact on plan:** Minor formatting fix only, no scope creep.

## Issues Encountered
None beyond the ruff line-length issue auto-fixed above.

## Next Phase Readiness
- Plan 02 (quickstart guide) can now be executed; the README link to docs/quickstart.md is already in place
- All five CLI subcommands have accurate --help text ready for user-facing documentation

## Self-Check: PASSED

All files exist and all task commits verified:
- mace_gaussian/cli.py: FOUND
- README.md: FOUND
- 10-01-SUMMARY.md: FOUND
- Commit 20bb2c0: FOUND
- Commit 0b796b7: FOUND

---
*Phase: 10-documentation*
*Completed: 2026-02-26*
