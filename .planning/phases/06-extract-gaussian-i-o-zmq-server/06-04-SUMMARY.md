---
phase: 06-extract-gaussian-i-o-zmq-server
plan: "04"
subsystem: gaussian
tags: [zmq, subprocess, timeout, sigkill, exceptions, runner]

# Dependency graph
requires:
  - phase: 06-extract-gaussian-i-o-zmq-server
    provides: GaussianZMQServer and is_calc_finished from gaussian/zmq_server.py
  - phase: 02-add-typed-exceptions
    provides: GaussianRunError and GaussianTimeoutError from utils/exceptions.py
provides:
  - run_gaussian_with_zmq function in gaussian/runner.py
  - DEFAULT_TIMEOUT_SECONDS constant reading GAUSSIAN_TIMEOUT_SECONDS env var
affects:
  - 06-05 (gm_main.py caller wiring will import from gaussian.runner)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - stdout=PIPE / stderr=PIPE on Popen for post-mortem diagnostics
    - proc.kill() (SIGKILL) for hard timeout enforcement (Gaussian ignores SIGTERM)
    - Decode captured bytes with errors="replace" for binary-safe output inclusion in exceptions
    - on_request callback pattern separates ML logic from subprocess lifecycle

key-files:
  created:
    - gaussian/runner.py
  modified: []

key-decisions:
  - "stdout=PIPE and stderr=PIPE always captured so both GaussianTimeoutError and GaussianRunError include full Gaussian output for diagnostics"
  - "proc.kill() (SIGKILL) on timeout, not proc.terminate() (SIGTERM) — Gaussian ignores SIGTERM"
  - "DEFAULT_TIMEOUT_SECONDS reads GAUSSIAN_TIMEOUT_SECONDS env var with 86400 (24h) fallback, matching prior gm_main.py behavior"
  - "on_request callback exceptions: send 'error' to Gaussian first, then re-raise (no silent swallowing)"
  - "runner depends only on gaussian.zmq_server and utils.exceptions — no dependency on gaussian.io or gaussian.parser"

patterns-established:
  - "Subprocess runner pattern: Popen with PIPE, service loop via is_calc_finished, timeout check with kill/wait, read captured streams after exit"

requirements-completed: [STRUCT-05, STRUCT-07, ERR-06]

# Metrics
duration: 2min
completed: 2026-02-20
---

# Phase 6 Plan 04: gaussian/runner.py Summary

**Subprocess runner with ZMQ callback loop, SIGKILL timeout, and stdout/stderr capture for typed Gaussian exceptions**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-20T09:11:39Z
- **Completed:** 2026-02-20T09:13:05Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Created `gaussian/runner.py` with `run_gaussian_with_zmq` consolidating timeout/kill/ZMQ logic from `gm_main.py`
- Enforces hard timeout via `proc.kill()` (SIGKILL) — Gaussian ignores SIGTERM
- Always captures stdout/stderr via PIPE; includes decoded output in both `GaussianTimeoutError` and `GaussianRunError` messages
- `DEFAULT_TIMEOUT_SECONDS` reads `GAUSSIAN_TIMEOUT_SECONDS` env var with 86400 (24h) fallback

## Task Commits

Each task was committed atomically:

1. **Task 1: Create gaussian/runner.py with run_gaussian_with_zmq** - `3a14ef4` (feat)

**Plan metadata:** (to be added with docs commit)

## Files Created/Modified
- `gaussian/runner.py` - Subprocess runner for Gaussian with ZMQ external interface, timeout, and typed exceptions

## Decisions Made
- stdout=PIPE and stderr=PIPE always captured so both GaussianTimeoutError and GaussianRunError include full Gaussian output for diagnostics
- proc.kill() (SIGKILL) on timeout, not proc.terminate() (SIGTERM) — Gaussian ignores SIGTERM
- DEFAULT_TIMEOUT_SECONDS reads GAUSSIAN_TIMEOUT_SECONDS env var with 86400 (24h) fallback, matching prior gm_main.py behavior
- on_request callback exceptions: send "error" to Gaussian first, then re-raise (no silent swallowing)
- runner imports only from gaussian.zmq_server and utils.exceptions — clean dependency boundary, no gaussian.io or gaussian.parser

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Pre-existing test failures (test_cli_validation.py missing `click`, test ordering issue in test_regression.py) confirmed not introduced by this plan via stash test.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- gaussian/runner.py is importable and provides run_gaussian_with_zmq
- Ready for plan 06-05: wire gm_main.py callers to use gaussian.runner instead of inline subprocess logic

---
*Phase: 06-extract-gaussian-i-o-zmq-server*
*Completed: 2026-02-20*

## Self-Check: PASSED

- FOUND: gaussian/runner.py
- FOUND: 06-04-SUMMARY.md
- FOUND commit: 3a14ef4 (feat(06-04): create gaussian/runner.py with run_gaussian_with_zmq)
