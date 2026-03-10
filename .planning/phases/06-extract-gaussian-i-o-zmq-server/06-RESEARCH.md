# Phase 6: Extract Gaussian I/O & ZMQ Server - Research

**Researched:** 2026-02-19
**Domain:** Python module extraction, pyzmq IPC, subprocess management, Python context managers
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**ZMQ server lifecycle:**
- Implement as a Python context manager (`with GaussianZMQServer(...) as server:`)
- On startup: auto-detect and delete stale `.ipc_file` from previous crashes before binding (known gotcha from CLAUDE.md)
- Cleanup (socket close + IPC file delete) runs in `__exit__` whether exiting cleanly or via exception — guaranteed via try/finally
- Expose minimal state to caller: socket path and a `running` flag; no internal message loop details leaked

**Subprocess runner:**
- Hard kill (`SIGKILL`) on timeout — Gaussian doesn't respond to SIGTERM reliably
- Timeout duration passed as parameter (no hardcoded magic number); default from existing behavior
- On failure: raise a custom exception (from the existing exception hierarchy) with stdout/stderr captured
- Log timeout and exit-code failures at ERROR level before raising

**Module split boundaries:**
- Follow success criteria exactly: `gaussian/io.py` (parse_gaussian_input, write output files), `gaussian/parser.py` (GaussianLogParser), `gaussian/fchk.py` (FCHK parser), `gaussian/zmq_server.py` (ZMQ server context manager), `gaussian/runner.py` (subprocess runner with timeout)
- Functions that are pure I/O (read/write files, no parsing logic) go to `io.py`
- Parser classes that interpret Gaussian output formats go to `parser.py`/`fchk.py`

**Migration strategy:**
- Update all callers in `gm_main.py` to import from `gaussian/` — no backwards-compat shims (internal code, clean break follows established project pattern)
- No re-exports from `gm_main.py` after migration
- Full pipeline smoke-tested on water molecule after each sub-plan

### Claude's Discretion

- All implementation details above (chosen by Claude, not user-specified)
- User explicitly deferred all decisions: "you'll handle it well, I trust you, make reasonable choices"
- `gaussian/__init__.py` contents (what to re-export publicly)
- Exact class/function signatures within each module

### Deferred Ideas (OUT OF SCOPE)

None — discussion stayed within phase scope.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-04 | Gaussian I/O functions extracted from gm_main.py into gaussian/ package | All I/O functions catalogued: `parse_gaussian_input` (line 134), `write_gaussian_output` (line 281), `ase_to_gjf` (line 403) go to `gaussian/io.py`. Caller map shows only `gm_main.py` and `testingStuff/` import them. |
| STRUCT-05 | ZMQ server context manager extracted into dedicated module | `zmq_server()` at line 83 and `is_calc_finished()` at line 113 extracted to `gaussian/zmq_server.py`. Class-based context manager design documented. |
| STRUCT-07 | ZMQ socket cleanup race condition fixed (proper LINGER settings) | Verified: default LINGER=-1 causes `socket.close()` to block forever if pending outgoing message not received. Fix: `socket.setsockopt(zmq.LINGER, 0)` before bind. Verified with pyzmq 27.0.2. |
</phase_requirements>

## Summary

Phase 6 is a targeted extraction of Gaussian-specific code from `gm_main.py` into a new `gaussian/` package with five modules. The extraction is well-bounded: all functions being moved live in `gm_main.py`, their callers are confined to that same file (plus a stale `testingStuff/` script), and no other modules in the production codebase import them directly. The parser files `gaussian_parser.py` and `fchk_parser.py` already exist as top-level modules and move wholesale into the package with updated imports.

The LINGER race condition (STRUCT-07) is the only bug fix in this phase. The current `zmq_server()` context manager uses pyzmq's default LINGER=-1, which causes `socket.close()` to block indefinitely if there is a pending unsent "ready" reply when Gaussian crashes. Verified on pyzmq 27.0.2: setting `socket.setsockopt(zmq.LINGER, 0)` before bind causes close to return immediately and discard pending messages. This is correct for a fire-and-forget reply model.

The subprocess runner (`gaussian/runner.py`) consolidates two divergent implementations: `gm_main.py` uses `subprocess.Popen` with a manual time-check loop and hard `proc.kill()`, while `dft_baseline.py` uses `proc.wait(timeout=...)` with `subprocess.TimeoutExpired`. Phase 6 creates a single authoritative runner for ML frequency calculations; `dft_baseline.py`'s runner is out of scope. Two new exception subclasses (`GaussianRunError` and `GaussianTimeoutError`) must be added to `utils/exceptions.py` since the existing hierarchy has no subprocess-specific types.

**Primary recommendation:** Extract functions directly — no behavioral changes, no new abstractions beyond the context manager class for ZMQ. Follow the no-shims convention from Phases 3-5 (clean import updates, delete nothing lingering).

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pyzmq | 27.0.2 (installed) | ZMQ IPC transport for Gaussian external interface | Already used; `zmq.REP`/`zmq.REQ` pair is the established protocol |
| subprocess | stdlib | Run Gaussian (g16) process | Already used; `Popen` for async launch with manual timeout |
| contextlib | stdlib | `contextmanager` decorator or ABC for context managers | Already used in `gm_main.py` via `@contextmanager` |
| pathlib | stdlib | Path manipulation | Already used across project; use `Path` not string ops |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| signal | stdlib | SIGKILL on Gaussian timeout | `proc.kill()` sends SIGKILL on POSIX; no extra import needed (Popen.kill() is SIGKILL) |
| os | stdlib | IPC file path resolution, stale file cleanup | `os.path.abspath()` for IPC path, `os.path.exists()` / `os.remove()` |
| logging | stdlib | ERROR-level logging on failures | Already configured project-wide |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `@contextmanager` decorator | Class-based `__enter__`/`__exit__` | Class is preferred: exposes `socket_path` and `running` attributes clearly; decorator makes it harder to access state. Locked decision: class-based. |
| `proc.terminate()` (SIGTERM) | `proc.kill()` (SIGKILL) | Gaussian ignores SIGTERM. Locked decision: SIGKILL. |
| `proc.wait(timeout=N)` + `TimeoutExpired` | Manual time loop | Manual loop is required because ZMQ polling must run concurrently — can't block on `proc.wait()`. |

## Architecture Patterns

### Target Project Structure

```
mace_gaussian/
├── gaussian/
│   ├── __init__.py          # Re-exports public API: GaussianZMQServer, parse_gaussian_input, etc.
│   ├── io.py                # parse_gaussian_input, write_gaussian_output, ase_to_gjf
│   ├── parser.py            # GaussianLogParser (moved from gaussian_parser.py), parse_gaussian_log
│   ├── fchk.py              # All fchk_parser.py contents: convert_chk_to_fchk, parse_fchk_section, etc.
│   ├── zmq_server.py        # GaussianZMQServer class + is_calc_finished()
│   └── runner.py            # run_gaussian_ml (replaces Popen+ZMQ code in run_frequency_calculation)
├── gaussian_parser.py       # DELETED — contents moved to gaussian/parser.py
├── fchk_parser.py           # DELETED — contents moved to gaussian/fchk.py
├── gm_main.py               # Updated: imports from gaussian/, keeps run_workflow etc.
├── utils/
│   └── exceptions.py        # ADD: GaussianRunError, GaussianTimeoutError subclasses
├── tests/
│   ├── test_gaussian_parser.py  # Update import: from gaussian.parser import GaussianLogParser
│   ├── test_fchk_parser.py      # Update import: from gaussian.fchk import ...
│   └── test_mode_matching.py    # Update import: from gaussian.fchk import ...
```

### Pattern 1: Class-Based Context Manager for ZMQ Server

**What:** `GaussianZMQServer` exposes `socket_path`, `running`, and `socket` as attributes. The caller accesses the socket directly from the context manager object.
**When to use:** Any time Gaussian needs to call back via ZMQ IPC.

```python
# gaussian/zmq_server.py
import logging
import os
import zmq
from utils.exceptions import MaceGaussianError

logger = logging.getLogger(__name__)


class GaussianZMQServer:
    """Context manager for Gaussian ZMQ IPC server.

    Handles stale IPC file cleanup on entry, LINGER=0 on socket,
    and guaranteed cleanup of both socket and IPC file on exit.

    Usage:
        with GaussianZMQServer(".ipc_file") as server:
            while not is_calc_finished(proc, server.socket):
                msg = server.socket.recv_string()
                ...
                server.socket.send_string("ready")
    """

    def __init__(self, ipc_file: str):
        self.socket_path = os.path.abspath(ipc_file)
        self.running = False
        self.socket = None
        self._ctx = None

    def __enter__(self) -> "GaussianZMQServer":
        # Remove stale IPC file from previous crash (documented CLAUDE.md gotcha)
        if os.path.exists(self.socket_path):
            try:
                os.remove(self.socket_path)
                logger.warning(f"Removed stale IPC file: {self.socket_path}")
            except OSError as e:
                logger.warning(f"Could not remove stale IPC file {self.socket_path}: {e}")

        self._ctx = zmq.Context()
        self.socket = self._ctx.socket(zmq.REP)
        # STRUCT-07: Set LINGER=0 before bind to prevent socket.close() blocking
        # if Gaussian crashes while a 'ready' reply is pending delivery.
        self.socket.setsockopt(zmq.LINGER, 0)
        self.socket.bind(f"ipc://{self.socket_path}")
        self.running = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Guaranteed cleanup: runs whether exiting cleanly or via exception
        self.running = False
        try:
            if self.socket is not None:
                self.socket.close()
        finally:
            try:
                if self._ctx is not None:
                    self._ctx.term()
            finally:
                if os.path.exists(self.socket_path):
                    try:
                        os.remove(self.socket_path)
                    except OSError:
                        pass
        return False  # Do not suppress exceptions


def is_calc_finished(proc, socket) -> bool:
    """Check if Gaussian calculation finished or next ML step requested.

    Polls ZMQ socket for incoming message (ML step request) and process
    exit status. Returns True when Gaussian process has exited.
    """
    import time
    while True:
        if socket.poll(timeout=10) != 0:
            return False   # New message: another ML step requested
        elif proc.poll() is not None:
            return True    # Process exited: calculation finished
        else:
            time.sleep(1)
```

### Pattern 2: Subprocess Runner with SIGKILL Timeout

**What:** `run_gaussian_subprocess` wraps `Popen`, manages the ZMQ callback loop, enforces hard timeout via `proc.kill()`, raises typed exceptions.
**When to use:** All ML frequency calculations that need the ZMQ external interface.

```python
# gaussian/runner.py
import logging
import os
import subprocess
import time
from pathlib import Path

from gaussian.zmq_server import GaussianZMQServer, is_calc_finished
from utils.exceptions import GaussianRunError, GaussianTimeoutError

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT_SECONDS = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))


def run_gaussian_with_zmq(
    gjf_file: str,
    on_request,             # callable(msg: str) -> str  (processes one ML step)
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    ipc_file: str = ".ipc_file",
) -> None:
    """Run Gaussian with ZMQ external interface.

    Raises:
        GaussianTimeoutError: if elapsed > timeout_seconds (hard kill via SIGKILL)
        GaussianRunError: if Gaussian exits with non-zero return code
    """
    proc = subprocess.Popen(["g16", gjf_file])
    start = time.time()

    with GaussianZMQServer(ipc_file) as server:
        while not is_calc_finished(proc, server.socket):
            elapsed = time.time() - start
            if elapsed > timeout_seconds:
                proc.kill()
                proc.wait()
                logger.error(
                    f"Gaussian timed out after {elapsed / 3600:.1f}h "
                    f"(limit: {timeout_seconds}s)"
                )
                raise GaussianTimeoutError(
                    f"Gaussian timed out after {elapsed:.0f}s "
                    f"(GAUSSIAN_TIMEOUT_SECONDS={timeout_seconds})"
                )
            msg = server.socket.recv_string()
            try:
                reply = on_request(msg)
                server.socket.send_string(reply)
            except Exception as e:
                server.socket.send_string("error")
                raise

    proc.wait()
    if proc.returncode != 0:
        logger.error(f"Gaussian exited with code {proc.returncode}")
        raise GaussianRunError(
            f"Gaussian (g16) exited with code {proc.returncode} for {gjf_file}"
        )
```

### Pattern 3: New Exception Subclasses in utils/exceptions.py

```python
# Add to utils/exceptions.py

class GaussianRunError(MaceGaussianError):
    """Gaussian subprocess exited with non-zero return code.

    Raised by gaussian.runner when g16 exits with an error.
    Callers can catch this to handle Gaussian-specific failures.
    """


class GaussianTimeoutError(MaceGaussianError):
    """Gaussian subprocess exceeded the configured timeout.

    Raised by gaussian.runner when elapsed time > timeout_seconds.
    The Gaussian process is hard-killed (SIGKILL) before this is raised.
    """
```

### Pattern 4: Module Moves (gaussian_parser.py -> gaussian/parser.py, fchk_parser.py -> gaussian/fchk.py)

Both files move wholesale. Only the import paths change.

```python
# gaussian/parser.py  (was: gaussian_parser.py)
# One change: internal import stays the same (from utils.exceptions import GaussianParseError)
# External callers updated from:
#   from gaussian_parser import GaussianLogParser
# To:
#   from gaussian.parser import GaussianLogParser

# gaussian/fchk.py  (was: fchk_parser.py)
# One change: internal import stays the same (from utils.units import BOHR_TO_ANGSTROM)
# External callers updated from:
#   from fchk_parser import extract_modes_from_fchk, ...
# To:
#   from gaussian.fchk import extract_modes_from_fchk, ...
```

### Pattern 5: gaussian/io.py — Pure Gaussian I/O Functions

```python
# gaussian/io.py
# Extract from gm_main.py:
#   - parse_gaussian_input() [lines 134-177]
#   - write_gaussian_output() [lines 281-348]
#   - ase_to_gjf() [lines 403-426]
# Internal deps: from utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV
# External deps: numpy, ase.data.chemical_symbols
```

### Pattern 6: gaussian/__init__.py — Public API Re-exports

```python
# gaussian/__init__.py
"""Gaussian integration package: I/O, parsing, ZMQ server, subprocess runner."""

from gaussian.fchk import (
    convert_chk_to_fchk,
    extract_modes_from_fchk,
    get_fchk_from_chk,
    parse_fchk_section,
)
from gaussian.io import ase_to_gjf, parse_gaussian_input, write_gaussian_output
from gaussian.parser import GaussianLogParser, parse_gaussian_log
from gaussian.runner import run_gaussian_with_zmq
from gaussian.zmq_server import GaussianZMQServer, is_calc_finished
```

### Anti-Patterns to Avoid

- **Leaving gaussian_parser.py and fchk_parser.py as stubs:** Delete them. The project pattern (Phases 3-5) is clean break, not compatibility shims.
- **Using `proc.terminate()` instead of `proc.kill()`:** Gaussian ignores SIGTERM. Always use `proc.kill()` (SIGKILL).
- **Omitting `socket.setsockopt(zmq.LINGER, 0)` before bind:** This is the STRUCT-07 bug. Without it, `socket.close()` may block forever when Gaussian crashes mid-reply.
- **Creating IPC file before bind:** The current `zmq_server()` creates a placeholder file via `open(addr, "x")` before `socket.bind()`. This is unnecessary — `socket.bind()` creates the IPC socket file itself. Remove the file-creation step.
- **Leaving `GAUSSIAN_TIMEOUT_SECONDS` module-level constant in gm_main.py:** Move the env-var read into `gaussian/runner.py` as `DEFAULT_TIMEOUT_SECONDS`. Update `test_cli_validation.py` accordingly.
- **Circular imports:** `gaussian/zmq_server.py` must not import from `gaussian/runner.py` and vice versa. Keep them separate layers.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| ZMQ context lifetime | Custom open/close tracking | `zmq.Context()` as context manager (already supports `__enter__`/`__exit__`) | pyzmq handles context termination and socket cleanup |
| IPC file locking | Custom lock files | OS-level IPC file presence check + `os.remove()` | Simple and sufficient; ZMQ IPC is local-only |
| Subprocess timeout | `threading.Timer` | `proc.kill()` in time-check loop | Already works; don't replace with threading |
| Custom binary format for Gaussian output | Custom byte packing | Existing `write_gaussian_output()` logic (Fortran D-exponent format) | Gaussian expects exactly this format; don't change it |

**Key insight:** The ZMQ + subprocess pattern is bespoke to the Gaussian external interface. Extract it cleanly but don't redesign it — it works. The only change is LINGER=0.

## Common Pitfalls

### Pitfall 1: ZMQ LINGER=-1 Blocks socket.close() Forever

**What goes wrong:** When Gaussian crashes mid-calculation, there may be a pending "ready" reply in the socket's outgoing buffer. With default LINGER=-1, `socket.close()` waits indefinitely for the message to be delivered. The process hangs.
**Why it happens:** pyzmq default LINGER is -1 (wait forever). When the context manager exits, `socket.__exit__` calls `socket.close()`.
**How to avoid:** Call `socket.setsockopt(zmq.LINGER, 0)` immediately after `ctx.socket(zmq.REP)` and before `socket.bind()`. Verified: with pyzmq 27.0.2, this takes effect immediately.
**Warning signs:** Process hangs after Gaussian crash; `ps aux` shows Python still running; `lsof` shows IPC socket still open.

### Pitfall 2: Stale IPC File Prevents Bind

**What goes wrong:** If a previous run crashed and left `.ipc_file` on disk, the new `socket.bind(f"ipc://{addr}")` raises `zmq.error.ZMQError: Address already in use`.
**Why it happens:** ZMQ IPC sockets create a filesystem socket file. If the previous context manager's `finally` block did not run (e.g., SIGKILL to the Python process), the file remains.
**How to avoid:** `GaussianZMQServer.__enter__` must check for and delete any existing file at `socket_path` before calling `socket.bind()`. Already implemented in current `zmq_server()` — preserve this logic.
**Warning signs:** `ZMQError: Address already in use` on startup; `.ipc_file` present in the working directory.

### Pitfall 3: Unnecessary IPC Placeholder File Creation

**What goes wrong:** The current `zmq_server()` creates a placeholder file with `open(addr, "x")` before calling `socket.bind()`. This is unnecessary — `socket.bind()` creates the socket file itself and will fail if the placeholder exists.
**Why it happens:** Historical artifact, possibly from debugging.
**How to avoid:** Remove the `open(addr, "x")` block entirely from the new `GaussianZMQServer.__enter__`. Only delete the stale file at entry; do not pre-create it.
**Warning signs:** `FileExistsError` during bind if the placeholder creation is kept but the stale-file-removal is also kept (they interfere).

### Pitfall 4: test_cli_validation.py Imports GAUSSIAN_TIMEOUT_SECONDS from gm_main

**What goes wrong:** `tests/test_cli_validation.py` at lines 68 and 82 imports `GAUSSIAN_TIMEOUT_SECONDS` from `gm_main`. After extraction, this constant moves to `gaussian/runner.py`.
**Why it happens:** Test tests the constant's value; import path changes with extraction.
**How to avoid:** Update `test_cli_validation.py` to import from `gaussian.runner` (or `os.getenv` fallback already present handles it gracefully since it wraps in try/except ImportError).
**Warning signs:** ImportError in test collection for `test_cli_validation.py`.

### Pitfall 5: callers in testingStuff/ and mode_matching.py Need Updating

**What goes wrong:** `testingStuff/test_refactoring.py` imports `parse_gaussian_input` and `write_gaussian_output` from `gm_main`. `mode_matching.py` imports `extract_modes_from_fchk` and `get_fchk_from_chk` from `fchk_parser`. `convert_all_chk_files.py` imports `convert_chk_to_fchk` from `fchk_parser`. `dft_baseline.py` imports `convert_chk_to_fchk` from `fchk_parser` (inline, line 377) and `parse_gaussian_log` from `gaussian_parser` (line 17).
**Why it happens:** Non-obvious callers outside `gm_main.py`.
**How to avoid:** Full caller sweep before deleting old files.

**Complete caller map:**

| Symbol | Current Source | New Source | Files to Update |
|--------|---------------|------------|-----------------|
| `GaussianLogParser` | `gaussian_parser` | `gaussian.parser` | `tests/test_gaussian_parser.py` |
| `parse_gaussian_log` | `gaussian_parser` | `gaussian.parser` | `gm_main.py` (line 35), `dft_baseline.py` (line 17) |
| `extract_modes_from_fchk` | `fchk_parser` | `gaussian.fchk` | `mode_matching.py` (line 21), `tests/test_mode_matching.py` (lines 166, 180), `tests/test_fchk_parser.py` (line 12) |
| `parse_fchk_section` | `fchk_parser` | `gaussian.fchk` | `tests/test_fchk_parser.py` (line 12) |
| `get_fchk_from_chk` | `fchk_parser` | `gaussian.fchk` | `mode_matching.py` (line 21) |
| `convert_chk_to_fchk` | `fchk_parser` | `gaussian.fchk` | `gm_main.py` (line 725), `dft_baseline.py` (line 377), `convert_all_chk_files.py` (line 11) |
| `parse_gaussian_input` | `gm_main` | `gaussian.io` | `gm_main.py` (line 378), `testingStuff/test_refactoring.py` |
| `write_gaussian_output` | `gm_main` | `gaussian.io` | `gm_main.py` (line 398), `testingStuff/test_refactoring.py` |
| `ase_to_gjf` | `gm_main` | `gaussian.io` | `gm_main.py` (line 662) |
| `zmq_server` (function) | `gm_main` | `gaussian.zmq_server.GaussianZMQServer` | `gm_main.py` (line 673) |
| `is_calc_finished` | `gm_main` | `gaussian.zmq_server` | `gm_main.py` (line 674) |
| `GAUSSIAN_TIMEOUT_SECONDS` | `gm_main` | `gaussian.runner` | `gm_main.py`, `tests/test_cli_validation.py` |

### Pitfall 6: pyproject.toml isort known-first-party Needs Update

**What goes wrong:** `pyproject.toml` lists `"gaussian_parser"` and `"fchk_parser"` as known first-party. After migration these become `"gaussian"`. Ruff will incorrectly sort `from gaussian.parser import ...` as third-party.
**Why it happens:** isort config not updated alongside code moves.
**How to avoid:** Add `"gaussian"` to `known-first-party`, remove `"gaussian_parser"` and `"fchk_parser"`.
**Warning signs:** Ruff import ordering violations after migration.

### Pitfall 7: Coverage Configuration Must Include gaussian Package

**What goes wrong:** `[tool.coverage.run] source` in `pyproject.toml` lists `"gaussian_parser"` and `"fchk_parser"` separately. After migration, both are submodules of `"gaussian"`.
**Why it happens:** Coverage source list not updated.
**How to avoid:** Add `"gaussian"` to source list, remove `"gaussian_parser"` and `"fchk_parser"`.

## Code Examples

### Current zmq_server() (to be replaced)

```python
# gm_main.py lines 83-110 — current implementation
@contextmanager
def zmq_server(file):
    addr = os.path.abspath(file)
    if os.path.exists(addr):
        try:
            os.remove(addr)
        except Exception as e:
            print(f"Warning: could not remove old IPC file {addr}: {e}")
    try:
        with zmq.Context() as ctx:
            with ctx.socket(zmq.REP) as socket:
                with open(addr, "x"):   # BUG: unnecessary placeholder, interferes with bind
                    pass
                socket.bind(f"ipc://{addr}")
                yield socket
    finally:
        if os.path.exists(addr):
            os.remove(addr)
# NOTE: Missing LINGER=0 — causes hang on Gaussian crash (STRUCT-07)
# NOTE: open(addr, "x") placeholder is an anti-pattern — remove it
```

### Correct GaussianZMQServer.__enter__ (removes both bugs)

```python
def __enter__(self) -> "GaussianZMQServer":
    # Remove stale IPC file (no open() placeholder — bind creates the socket file)
    if os.path.exists(self.socket_path):
        try:
            os.remove(self.socket_path)
            logger.warning(f"Removed stale IPC file: {self.socket_path}")
        except OSError as e:
            logger.warning(f"Could not remove stale IPC file {self.socket_path}: {e}")

    self._ctx = zmq.Context()
    self.socket = self._ctx.socket(zmq.REP)
    self.socket.setsockopt(zmq.LINGER, 0)  # STRUCT-07 fix
    self.socket.bind(f"ipc://{self.socket_path}")
    self.running = True
    return self
```

### parse_gaussian_input Signature (unchanged, moved to gaussian/io.py)

```python
# gaussian/io.py
# Source: gm_main.py lines 134-177
def parse_gaussian_input(infile: str) -> tuple[int, int, int, int, np.ndarray, list]:
    """Parse Gaussian external calculation input file.
    Returns: (natoms, deriv, charge, spin, coordinates_angstrom, atomnames)
    """
    # ... body unchanged ...
```

### write_gaussian_output Signature (unchanged, moved to gaussian/io.py)

```python
# gaussian/io.py
# Source: gm_main.py lines 281-348
def write_gaussian_output(
    outfile: str,
    natoms: int,
    energy: float,           # eV
    gradient: np.ndarray,    # eV/Angstrom, shape (natoms, 3)
    dipole: np.ndarray,      # e*Bohr, shape (3,)
    dipole_derivatives: np.ndarray,  # e, shape (3*natoms, 3)
    hessian: Optional[np.ndarray],   # Hartree/Bohr^2, shape (3*natoms, 3*natoms)
    deriv: int,
) -> None:
    # Fortran D-exponent format: "1.234D+05" not "1.234E+05"
    # DO NOT CHANGE: Gaussian requires exactly this format
    # ... body unchanged ...
```

### run_frequency_calculation Refactoring (gm_main.py caller site)

```python
# gm_main.py — after Phase 6
# Before:
#   proc = subprocess.Popen(["g16", gjf_temp])
#   with zmq_server(".ipc_file") as socket:
#       while not is_calc_finished(proc, socket):
#           ...  (timeout check + recv + run_next_calculation + send)

# After:
from gaussian.runner import run_gaussian_with_zmq

def _make_ml_handler(mol, calc, dipole_method):
    def on_request(msg: str) -> str:
        run_next_calculation(mol, msg, calc, dipole_method=dipole_method)
        return "ready"
    return on_request

run_gaussian_with_zmq(
    gjf_temp,
    on_request=_make_ml_handler(mol, calc, dipole_calculator_name),
    timeout_seconds=GAUSSIAN_TIMEOUT_SECONDS,  # or read from gaussian.runner.DEFAULT_TIMEOUT_SECONDS
)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Functions in `gm_main.py` | `gaussian/` package | Phase 6 | Clean module boundary |
| `@contextmanager` decorator | Class-based `GaussianZMQServer` | Phase 6 | Exposes `socket_path`, `running` state |
| LINGER=-1 (default) | LINGER=0 explicit | Phase 6 | Prevents hang on Gaussian crash |
| `open(addr, "x")` placeholder | Removed | Phase 6 | Eliminates interference with bind |
| Inline timeout logic in `run_frequency_calculation` | `run_gaussian_with_zmq` in `gaussian/runner.py` | Phase 6 | Single place to maintain timeout/kill logic |
| `gaussian_parser.py` top-level | `gaussian/parser.py` | Phase 6 | Package-scoped |
| `fchk_parser.py` top-level | `gaussian/fchk.py` | Phase 6 | Package-scoped |

**Deprecated/outdated after Phase 6:**
- `gaussian_parser.py` top-level file: deleted
- `fchk_parser.py` top-level file: deleted
- `zmq_server()` function in `gm_main.py`: replaced by `GaussianZMQServer` class
- `is_calc_finished()` function in `gm_main.py`: moved to `gaussian/zmq_server.py`
- `parse_gaussian_input()`, `write_gaussian_output()`, `ase_to_gjf()` in `gm_main.py`: moved to `gaussian/io.py`

## Open Questions

1. **Should `run_next_calculation()` stay in gm_main.py or move to gaussian/runner.py?**
   - What we know: `run_next_calculation` orchestrates: parse input → update geometry → calculate energy/forces → calculate hessian → calculate dipole → write output. It calls `dipole_factory.get_calculator()` from `calculators/`. This creates a cross-package dependency between `gaussian/` and `calculators/`.
   - What's unclear: Is `run_next_calculation` "Gaussian I/O" or "workflow orchestration"? The CONTEXT.md module boundaries put I/O in `gaussian/io.py` but pure parsing/writing only. The orchestration function arguably belongs to `gm_main.py`.
   - Recommendation: Leave `run_next_calculation()` in `gm_main.py`. It is the glue between the Gaussian I/O layer (`gaussian/io.py`) and the calculator layer (`calculators/`). Moving it to `gaussian/` would introduce an import from `gaussian` back to `calculators`, which is a cross-package coupling anti-pattern. The runner passes it as a callback (`on_request`), maintaining clean separation.

2. **Where does `GAUSSIAN_TIMEOUT_SECONDS` live after migration?**
   - What we know: It is defined in `gm_main.py` and imported by `tests/test_cli_validation.py`. It needs to be accessible from `gaussian/runner.py` as the default.
   - Recommendation: Define `DEFAULT_TIMEOUT_SECONDS = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))` in `gaussian/runner.py`. Update `test_cli_validation.py` to import from `gaussian.runner` (or use the env-var fallback that is already in the test). Remove from `gm_main.py`.

3. **Should `ase_to_gjf()` move to gaussian/io.py or stay in gm_main.py?**
   - What we know: `ase_to_gjf()` references `DEFAULT_HELPER_SCRIPT` (env var for the path to `gm_helper.py`). It is only called from `gm_main.py` (line 662) in the test of production code. The function writes `.gjf` files — pure I/O.
   - Recommendation: Move to `gaussian/io.py` with `DEFAULT_HELPER_SCRIPT` also moving to `gaussian/io.py` (it's logically a Gaussian configuration constant, not a gm_main workflow setting). The caller in `gm_main.py` imports `ase_to_gjf` from `gaussian.io`.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis: `gm_main.py` (1042 lines, all functions catalogued with line numbers), `gaussian_parser.py` (377 lines), `fchk_parser.py` (322 lines), `gm_helper.py` (38 lines), `utils/exceptions.py`
- `pyproject.toml` ruff and coverage configuration
- pyzmq 27.0.2 installed locally — verified LINGER behavior interactively

### Secondary (MEDIUM confidence)
- Phase 3 RESEARCH.md — established no-shims convention and import update patterns
- Phase 4 RESEARCH.md — established `calculators/` package structure as the model for `gaussian/`
- Phase 5 RESEARCH.md — confirmed clean-break migration strategy

### Tertiary (LOW confidence)
- ZMQ LINGER semantics: verified empirically with pyzmq 27.0.2 (default LINGER=-1 confirmed; setsockopt LINGER=0 verified working) — no external docs consulted

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — pyzmq already installed and verified; stdlib only beyond that
- Architecture: HIGH — all source functions catalogued with line numbers, all callers identified
- Pitfalls: HIGH — LINGER bug reproduced empirically; caller map built from direct grep analysis
- Caller map: HIGH — grep across all .py files, verified against test files

**Research date:** 2026-02-19
**Valid until:** 2026-03-19 (stable — no external API dependencies, pure internal refactoring)
