# Phase 6: Extract Gaussian I/O & ZMQ Server - Context

**Gathered:** 2026-02-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Extract Gaussian-specific code from `gm_main.py` into a focused `gaussian/` package with four modules: `io.py` (I/O functions), `parser.py` + `fchk.py` (log/FCHK parsers), `zmq_server.py` (real-time dipole injection server), and `runner.py` (Gaussian subprocess management). The full pipeline must remain functional throughout. Analysis modules and workflow orchestration are separate phases.

</domain>

<decisions>
## Implementation Decisions

### ZMQ server lifecycle
- Implement as a Python context manager (`with GaussianZMQServer(...) as server:`)
- On startup: auto-detect and delete stale `.ipc_file` from previous crashes before binding (known gotcha from CLAUDE.md)
- Cleanup (socket close + IPC file delete) runs in `__exit__` whether exiting cleanly or via exception — guaranteed via try/finally
- Expose minimal state to caller: socket path and a `running` flag; no internal message loop details leaked

### Subprocess runner
- Hard kill (`SIGKILL`) on timeout — Gaussian doesn't respond to SIGTERM reliably
- Timeout duration passed as parameter (no hardcoded magic number); default from existing behavior
- On failure: raise a custom exception (from the existing exception hierarchy) with stdout/stderr captured
- Log timeout and exit-code failures at ERROR level before raising

### Module split boundaries
- Follow success criteria exactly: `gaussian/io.py` (parse_gaussian_input, write output files), `gaussian/parser.py` (GaussianLogParser), `gaussian/fchk.py` (FCHK parser), `gaussian/zmq_server.py` (ZMQ server context manager), `gaussian/runner.py` (subprocess runner with timeout)
- Functions that are pure I/O (read/write files, no parsing logic) go to `io.py`
- Parser classes that interpret Gaussian output formats go to `parser.py`/`fchk.py`

### Migration strategy
- Update all callers in `gm_main.py` to import from `gaussian/` — no backwards-compat shims (internal code, clean break follows established project pattern)
- No re-exports from `gm_main.py` after migration
- Full pipeline smoke-tested on water molecule after each sub-plan

### Claude's Discretion
- All implementation details above (chosen by Claude, not user-specified)
- User explicitly deferred all decisions: "you'll handle it well, I trust you, make reasonable choices"
- `gaussian/__init__.py` contents (what to re-export publicly)
- Exact class/function signatures within each module

</decisions>

<specifics>
## Specific Ideas

- ZMQ stale IPC file cleanup is a documented gotcha in CLAUDE.md — must be handled automatically, not left to user
- Follow the no-shims convention established in Phase 3 (clean break, direct import updates)
- Exception hierarchy already exists from Phase 2 — use it for subprocess/timeout errors

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 06-extract-gaussian-i-o-zmq-server*
*Context gathered: 2026-02-19*
