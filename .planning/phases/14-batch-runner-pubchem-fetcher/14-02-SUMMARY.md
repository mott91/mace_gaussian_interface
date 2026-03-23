---
phase: 14-batch-runner-pubchem-fetcher
plan: 02
subsystem: cli
tags: [batch, manifest, json, click, restart-safety, atomic-write]

requires:
  - phase: 14-01
    provides: PubChem fetcher CLI (independent, no runtime dependency)
provides:
  - "mace-gaussian batch CLI command for multi-molecule pipeline runs"
  - "Manifest-driven per-calculator restart (batch_manifest.json)"
  - "Per-calculator failure isolation with error recording"
  - "Atomic manifest writes via tempfile + os.replace"
affects: [15-slurm-dft, 16-benchmark-campaign, 17-batch-report]

tech-stack:
  added: []
  patterns: [manifest-driven-restart, atomic-json-write, per-calculator-granularity]

key-files:
  created:
    - mace_gaussian/batch.py
    - tests/test_batch.py
  modified:
    - mace_gaussian/cli.py

key-decisions:
  - "Batch calls three stage functions directly (not run_pipeline) for per-calculator restart granularity"
  - "Atomic manifest write via tempfile.mkstemp + os.replace prevents corruption on interrupt"
  - "Omitted cleanup_stale_scratch_dirs from batch CLI (not available in worktree; non-blocking)"

patterns-established:
  - "Manifest-based restart: load manifest -> check status -> skip complete -> run pending -> save after each"
  - "Per-calculator granularity: energy_calc_dipole_calc keys under each molecule in manifest"

requirements-completed: [BATCH-02, BATCH-03, BATCH-04]

duration: 6min
completed: 2026-03-23
---

# Phase 14 Plan 02: Batch Runner Summary

**Manifest-driven batch runner with per-calculator restart, failure isolation, and atomic JSON writes for crash-safe multi-molecule pipeline execution**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-23T11:47:50Z
- **Completed:** 2026-03-23T11:54:03Z
- **Tasks:** 2 (Task 1 was TDD: RED + GREEN)
- **Files modified:** 3

## Accomplishments
- `mace_gaussian/batch.py` (310 lines) with load_manifest, save_manifest, parse_batch_file, run_batch
- Per-calculator granularity: each energy x dipole combination tracked independently in manifest
- Failure isolation: failed combos recorded with error details, batch continues to next
- Restart safety: combinations marked "complete" skipped on re-run
- `mace-gaussian batch molecules.txt` CLI command with same options as `run`
- 16 unit tests covering manifest CRUD, batch file parsing, restart, failure isolation

## Task Commits

Each task was committed atomically:

1. **Task 1 RED: Failing tests** - `34b0454` (test)
2. **Task 1 GREEN: batch.py implementation + test fixes** - `90eb789` (feat)
3. **Task 2: Add batch CLI command** - `912e290` (feat)

## Files Created/Modified
- `mace_gaussian/batch.py` - Batch runner with manifest-driven per-calculator restart
- `tests/test_batch.py` - 16 unit tests for batch module
- `mace_gaussian/cli.py` - Added `batch` command with same options as `run`

## Decisions Made
- Batch calls `run_geometry_optimization`, `run_dft_baselines`, `run_frequency_calculation` directly (not `run_pipeline`) for per-calculator restart granularity
- Atomic manifest write via `tempfile.mkstemp` + `os.replace` prevents corruption on Ctrl+C
- `os.replace` kept with noqa comment rather than `Path.replace` to satisfy acceptance criteria and maintain explicit POSIX atomic rename semantics

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Omitted cleanup_stale_scratch_dirs from batch CLI**
- **Found during:** Task 2 (CLI batch command)
- **Issue:** Plan references `from mace_gaussian.utils.scratch import cleanup_stale_scratch_dirs` but this module does not exist in the worktree (was added in Phase 13.2 which is not merged into this branch)
- **Fix:** Omitted the import and call; batch CLI works without scratch cleanup. Will be available after branch merge.
- **Files modified:** mace_gaussian/cli.py
- **Verification:** CLI imports and all tests pass

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Scratch cleanup is a convenience feature, not critical for batch runner correctness. Will be available after merge with main branch containing Phase 13.2 changes.

## Issues Encountered
- Worktree is based on an older branch (pre-Phase 13.2) missing `utils/scratch.py` and `keep_scratch` parameter on `run_dft_baselines`. Implemented against the actual API available in the worktree.

## Known Stubs
None - all functions are fully implemented.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Batch runner ready for multi-molecule campaigns
- SLURM DFT offload (Phase 15) can integrate with batch by wrapping DFT baseline stage
- Benchmark campaign (Phase 16) can use `mace-gaussian batch` directly

---
*Phase: 14-batch-runner-pubchem-fetcher*
*Completed: 2026-03-23*
