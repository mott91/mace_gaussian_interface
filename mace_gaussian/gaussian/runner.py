"""Gaussian subprocess runner with ZMQ external interface.

Provides run_gaussian_with_zmq, which launches Gaussian (g16), manages
the ZMQ callback loop, enforces a hard timeout via SIGKILL, captures
stdout/stderr for diagnostics, and raises typed exceptions on failure.
"""
from __future__ import annotations

import logging
import os
import subprocess
import time

from gaussian.zmq_server import GaussianZMQServer, is_calc_finished
from utils.exceptions import GaussianRunError, GaussianTimeoutError

logger = logging.getLogger(__name__)

# Default timeout: 24 hours. Override with GAUSSIAN_TIMEOUT_SECONDS env var.
# This constant replaces GAUSSIAN_TIMEOUT_SECONDS in gm_main.py.
DEFAULT_TIMEOUT_SECONDS: int = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))


def run_gaussian_with_zmq(
    gjf_file: str,
    on_request: object,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    ipc_file: str = ".ipc_file",
) -> None:
    """Run Gaussian (g16) with ZMQ external interface callback loop.

    Launches Gaussian, then services ZMQ messages in a loop until the
    process exits. Each message contains "infile|outfile" paths; the
    caller-supplied on_request callback processes one ML step and returns
    the reply string (typically "ready").

    Enforces a hard timeout: if elapsed time exceeds timeout_seconds, the
    Gaussian process is killed via SIGKILL (proc.kill()) and
    GaussianTimeoutError is raised. Gaussian ignores SIGTERM, so SIGKILL
    is required.

    stdout and stderr are always captured (stdout=PIPE, stderr=PIPE) so
    that GaussianTimeoutError and GaussianRunError messages contain the
    Gaussian output for diagnostics.

    Args:
        gjf_file: Path to the Gaussian input .gjf file.
        on_request: Callable(msg: str) -> str. Called for each ZMQ message.
            Should process one ML calculation step and return the reply
            string. Exceptions from on_request are re-raised after sending
            "error" to Gaussian.
        timeout_seconds: Hard timeout in seconds. Defaults to
            DEFAULT_TIMEOUT_SECONDS (GAUSSIAN_TIMEOUT_SECONDS env var or
            86400 = 24h).
        ipc_file: Path to ZMQ IPC socket file. Default: ".ipc_file".

    Raises:
        GaussianTimeoutError: elapsed > timeout_seconds. Gaussian is
            hard-killed (SIGKILL) before this is raised. Exception message
            includes captured stdout/stderr.
        GaussianRunError: Gaussian exits with non-zero return code.
            Exception message includes captured stdout/stderr.
    """
    proc = subprocess.Popen(
        ["g16", gjf_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    start = time.time()

    with GaussianZMQServer(ipc_file) as server:
        while not is_calc_finished(proc, server.socket):
            elapsed = time.time() - start
            if elapsed > timeout_seconds:
                proc.kill()
                proc.wait()
                stdout_data = proc.stdout.read().decode(errors="replace")
                stderr_data = proc.stderr.read().decode(errors="replace")
                logger.error(
                    "Gaussian timed out after %.1fh (limit: %ds, gjf: %s)",
                    elapsed / 3600,
                    timeout_seconds,
                    gjf_file,
                )
                raise GaussianTimeoutError(
                    f"Gaussian timed out after {elapsed:.0f}s "
                    f"(GAUSSIAN_TIMEOUT_SECONDS={timeout_seconds}, gjf={gjf_file})\n"
                    f"stdout: {stdout_data}\nstderr: {stderr_data}"
                )

            msg = server.socket.recv_string()
            try:
                reply = on_request(msg)
                server.socket.send_string(reply)
            except Exception:
                server.socket.send_string("error")
                raise

    proc.wait()
    if proc.returncode != 0:
        stdout_data = proc.stdout.read().decode(errors="replace")
        stderr_data = proc.stderr.read().decode(errors="replace")
        logger.error(
            "Gaussian exited with code %d (gjf: %s)",
            proc.returncode,
            gjf_file,
        )
        raise GaussianRunError(
            f"Gaussian (g16) exited with code {proc.returncode} for {gjf_file}\n"
            f"stdout: {stdout_data}\nstderr: {stderr_data}"
        )
