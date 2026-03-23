---
phase: 14-batch-runner-pubchem-fetcher
plan: 01
subsystem: cli
tags: [pubchem, requests, ase, sdf, xyz, click]

requires:
  - phase: 13-calculator-expansion-acoh-bug-fix
    provides: CLI validation callbacks and calculator constants
provides:
  - "fetch_3d_structure() PubChem API client in mace_gaussian/pubchem.py"
  - "mace-gaussian fetch <name> CLI command with --force and --output-dir"
  - "requests>=2.27.0 as declared dependency"
affects: [14-02-batch-runner, documentation]

tech-stack:
  added: [requests]
  patterns: [PubChem PUG REST API name-to-SDF-to-XYZ pipeline]

key-files:
  created: [mace_gaussian/pubchem.py, tests/test_pubchem.py]
  modified: [mace_gaussian/cli.py, pyproject.toml]

key-decisions:
  - "FileExistsError exits 0 (skip with warning), ValueError exits 1 (not found) per D-03/D-04"
  - "SDF-to-XYZ conversion via ASE read/write (not RDKit) for consistency with project patterns"

patterns-established:
  - "PubChem fetch: requests.get -> SDF -> ase.io.read -> ase.io.write XYZ"
  - "CLI error handling: import at function scope to avoid heavy deps at CLI load time"

requirements-completed: [BATCH-01]

duration: 17min
completed: 2026-03-23
---

# Phase 14 Plan 01: PubChem Fetcher Summary

**PubChem 3D structure fetcher via PUG REST API with CLI command `mace-gaussian fetch <name>` saving XYZ to molecules/ directory**

## Performance

- **Duration:** 17 min
- **Started:** 2026-03-23T11:27:52Z
- **Completed:** 2026-03-23T11:44:26Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- `fetch_3d_structure()` function handles PubChem name lookup, SDF download, ASE conversion to XYZ
- 6 unit tests with mocked HTTP responses covering success, unknown molecule, no 3D conformer, file exists, force overwrite, and network timeout
- `mace-gaussian fetch` CLI command with `--force` and `--output-dir` options
- `requests>=2.27.0` added to pyproject.toml dependencies

## Task Commits

Each task was committed atomically:

1. **Task 1: Create pubchem.py module and test_pubchem.py (TDD)** - `a377b6f` (test: RED), `bc7a4e4` (feat: GREEN)
2. **Task 2: Add fetch CLI command and requests dependency** - `dc8e3bb` (feat)

## Files Created/Modified
- `mace_gaussian/pubchem.py` - PubChem PUG REST API client with fetch_3d_structure()
- `tests/test_pubchem.py` - 6 unit tests with mocked HTTP for all fetch scenarios
- `mace_gaussian/cli.py` - Added fetch command after run, before list
- `pyproject.toml` - Added requests>=2.27.0 to dependencies

## Decisions Made
- FileExistsError exits with code 0 (skip is not an error, per D-04)
- ValueError exits with code 1 (unknown molecule is a real error, per D-03)
- Used ASE for SDF-to-XYZ conversion (already a core dependency, consistent with project)
- Lazy import of pubchem module inside fetch command to avoid loading heavy deps at CLI startup

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- pubchem.py is available for import by batch command (Plan 02)
- fetch command registered in CLI group, ready for user invocation
- All 24 tests pass (6 pubchem + 18 CLI validation)

---
*Phase: 14-batch-runner-pubchem-fetcher*
*Completed: 2026-03-23*

## Self-Check: PASSED

All files exist, all commits verified.
