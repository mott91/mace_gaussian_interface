---
phase: 15-slurm-integration-batch-report
plan: 01
subsystem: hpc
tags: [slurm, ssh, scp, batch, dft, subprocess]

# Dependency graph
requires:
  - phase: 14-pubchem-fetcher-batch-runner
    provides: batch.py manifest system, batch CLI command
provides:
  - SLURM DFT job submission via SSH/SCP
  - sacct-based polling with exponential backoff
  - Result retrieval with local formchk fallback
  - --dft-on-cluster and --slurm-template CLI options
affects: [16-benchmark-campaign, 15-02-batch-report]

# Tech tracking
tech-stack:
  added: []
  patterns: [SSH subprocess wrappers with backoff, SLURM template-based submission]

key-files:
  created:
    - mace_gaussian/slurm.py
    - templates/slurm_dft.sh
    - tests/test_slurm.py
  modified:
    - mace_gaussian/batch.py
    - mace_gaussian/cli.py

key-decisions:
  - "formchk included in SLURM template (cluster has g16) with local fallback if it fails"
  - "sacct polling with 60s sleep increments for Ctrl-C interruptibility"
  - "Lazy imports for slurm module inside dft_on_cluster conditional branch"

patterns-established:
  - "SSH/SCP wrappers with ConnectTimeout=30, BatchMode=yes, StrictHostKeyChecking=accept-new"
  - "Exponential backoff: min(2^attempt * 30, 600) seconds, max 5 retries"
  - "SLURM template placeholders: {molecule}, {remote_dir}, {gjf_filename}, {chk_filename}, {fchk_filename}"

requirements-completed: [HPC-01, HPC-02]

# Metrics
duration: 5min
completed: 2026-03-27
---

# Phase 15 Plan 01: SLURM Integration Summary

**SLURM DFT offloading via SSH/SCP with sacct polling, exponential backoff retry, and manifest-tracked job state**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-27T19:44:10Z
- **Completed:** 2026-03-27T19:49:39Z
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments
- Complete SLURM module with submit, poll, retrieve workflow and SSH backoff retry
- User-editable SLURM template with formchk step and configurable resources
- Batch runner branches to SLURM when --dft-on-cluster is set, with dft_failed manifest marking
- 17 unit tests covering command construction, sacct parsing, backoff logic, and poll behavior

## Task Commits

Each task was committed atomically:

1. **Task 1: Create SLURM module and job template** - `72fe254` (feat)
2. **Task 2: Integrate SLURM into batch runner and CLI** - `23e1896` (feat)
3. **Task 3: Unit tests for SLURM module** - `25df55e` (test)

## Files Created/Modified
- `mace_gaussian/slurm.py` - SSH/SCP wrappers, submit_dft_jobs, poll_jobs, retrieve_results
- `templates/slurm_dft.sh` - User-editable SLURM job template with formchk step
- `mace_gaussian/batch.py` - Extended run_batch with dft_on_cluster/slurm_template params
- `mace_gaussian/cli.py` - Added --dft-on-cluster and --slurm-template options to batch command
- `tests/test_slurm.py` - 17 tests with mocked SSH/SCP subprocess calls

## Decisions Made
- formchk included in SLURM template since the cluster runs Gaussian 16; local fallback covers cases where it fails on cluster
- sacct polling uses 60-second sleep increments instead of one long sleep for Ctrl-C interruptibility
- SLURM imports are lazy (inside conditional) to avoid loading the module when --dft-on-cluster is not used

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Default Python environment (3.8) has incompatible DGL/PyTorch; tests run correctly under mace4ir_v2 conda environment as documented in project memory

## User Setup Required

None - no external service configuration required.

## Known Stubs

None - all functions are fully implemented. The SLURM workflow requires a real SSH target to exercise end-to-end but all logic paths are unit-tested with mocks.

## Next Phase Readiness
- SLURM integration complete, ready for batch report (plan 15-02)
- Batch runner can now offload DFT to HPC cluster for the ~25-molecule benchmark campaign (Phase 16)

---
*Phase: 15-slurm-integration-batch-report*
*Completed: 2026-03-27*
