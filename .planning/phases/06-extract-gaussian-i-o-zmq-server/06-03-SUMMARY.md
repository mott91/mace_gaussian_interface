---
phase: 06-extract-gaussian-i-o-zmq-server
plan: "03"
subsystem: infra
tags: [zmq, ipc, gaussian, context-manager, socket]

# Dependency graph
requires:
  - phase: 06-extract-gaussian-i-o-zmq-server
    provides: gaussian/ package with __init__.py foundation (06-01)
provides:
  - GaussianZMQServer class-based context manager in gaussian/zmq_server.py
  - is_calc_finished function in gaussian/zmq_server.py
  - STRUCT-07 fix: LINGER=0 on ZMQ socket before bind
  - Correct stale IPC file removal without placeholder file creation
affects:
  - 06-05-PLAN (wires gm_main.py callers to new gaussian.zmq_server)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Class-based context manager with nested try/finally for guaranteed resource cleanup"
    - "LINGER=0 set before socket.bind() to prevent hang on crash"

key-files:
  created:
    - gaussian/zmq_server.py
  modified: []

key-decisions:
  - "LINGER=0 applied after socket() and before bind() to prevent socket.close() from blocking forever on Gaussian crash (STRUCT-07)"
  - "No open() placeholder file creation in __enter__ — socket.bind() creates IPC file itself; this was the original bug"
  - "__exit__ returns False to not suppress exceptions; nested try/finally ensures cleanup even if socket.close() raises"
  - "Instance attributes (running, socket, socket_path) set in __init__; verify script tests on instance not class"

patterns-established:
  - "ZMQ context manager pattern: setsockopt(LINGER, 0) before bind(), nested try/finally for cleanup ordering"

requirements-completed: [STRUCT-05, STRUCT-07]

# Metrics
duration: 3min
completed: 2026-02-20
---

# Phase 6 Plan 03: Extract ZMQ Server Summary

**GaussianZMQServer class-based context manager with LINGER=0 race-condition fix and correct stale IPC cleanup, extracted into gaussian/zmq_server.py**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-02-20T09:08:33Z
- **Completed:** 2026-02-20T09:11:00Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Created `gaussian/zmq_server.py` with `GaussianZMQServer` class implementing STRUCT-05 (dedicated ZMQ module)
- Fixed STRUCT-07 race condition: `setsockopt(zmq.LINGER, 0)` called after `socket()` and before `bind()`, preventing indefinite hang when Gaussian crashes with pending reply in outgoing buffer
- Fixed original `zmq_server()` bug: removed `open(addr, "x")` placeholder file creation; `socket.bind()` creates the IPC file itself
- Nested `try/finally` in `__exit__` guarantees `socket.close()` then `ctx.term()` then IPC file removal run regardless of exceptions

## Task Commits

Each task was committed atomically:

1. **Task 1: Create gaussian/zmq_server.py with GaussianZMQServer class** - `67666a6` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `/home/mot/mace_gaussian/gaussian/zmq_server.py` - GaussianZMQServer context manager and is_calc_finished function

## Decisions Made
- LINGER=0 must be set after `ctx.socket(zmq.REP)` and before `socket.bind()` — this is the correct ordering for pyzmq
- No placeholder file creation: the original `open(addr, "x")` call in `gm_main.zmq_server()` was a bug that created a regular file at the IPC path before `bind()`, interfering with socket creation
- `__exit__` returns `False` (not suppress) — documented in code with comment

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Plan verify script used `hasattr(GaussianZMQServer, 'running')` on the class, but `running` is an instance attribute set in `__init__`. Verified correctness using an instance instead — this is correct Python design; the class attribute check in the verify script was the anomaly, not the implementation.

## Next Phase Readiness
- `gaussian/zmq_server.py` is importable and fully functional
- `gm_main.py` callers NOT updated yet — deferred to plan 06-05 per plan scope boundaries
- Ready for plan 06-04 (next plan in phase)

---
*Phase: 06-extract-gaussian-i-o-zmq-server*
*Completed: 2026-02-20*
