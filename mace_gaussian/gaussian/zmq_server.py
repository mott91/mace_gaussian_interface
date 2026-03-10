"""ZMQ IPC server for Gaussian external interface.

Provides GaussianZMQServer, a context manager that binds a ZMQ REP socket
for Gaussian's external interface protocol, and is_calc_finished() for
monitoring calculation completion.
"""

from __future__ import annotations

import contextlib
import logging
import time
from pathlib import Path

import zmq

logger = logging.getLogger(__name__)


class GaussianZMQServer:
    """Context manager for Gaussian ZMQ IPC server.

    Handles stale IPC file cleanup on entry, LINGER=0 on socket,
    and guaranteed cleanup of both socket and IPC file on exit.

    The LINGER=0 setting (STRUCT-07 fix) prevents socket.close() from
    blocking forever when Gaussian crashes with a pending 'ready' reply
    in the outgoing buffer. Without it, pyzmq default LINGER=-1 causes
    the process to hang indefinitely.

    Usage:
        with GaussianZMQServer(".ipc_file") as server:
            while not is_calc_finished(proc, server.socket):
                msg = server.socket.recv_string()
                ...
                server.socket.send_string("ready")

    Attributes:
        socket_path: Absolute path to the IPC socket file.
        running: True while the context manager is active.
        socket: The bound zmq.REP socket (available inside the with block).
    """

    def __init__(self, ipc_file: str) -> None:
        self.socket_path: Path = Path(ipc_file).resolve()
        self.running: bool = False
        self.socket: zmq.Socket | None = None
        self._ctx: zmq.Context | None = None

    def __enter__(self) -> GaussianZMQServer:
        # Remove stale IPC file from a previous crash (documented CLAUDE.md gotcha).
        # Do NOT create a placeholder — socket.bind() creates the IPC socket file itself.
        if self.socket_path.exists():
            try:
                self.socket_path.unlink()
                logger.warning("Removed stale IPC file: %s", self.socket_path)
            except OSError as e:
                logger.warning("Could not remove stale IPC file %s: %s", self.socket_path, e)

        self._ctx = zmq.Context()
        self.socket = self._ctx.socket(zmq.REP)
        # STRUCT-07: Set LINGER=0 before bind so socket.close() returns immediately
        # if Gaussian crashes while a 'ready' reply is pending delivery.
        self.socket.setsockopt(zmq.LINGER, 0)
        self.socket.bind(f"ipc://{self.socket_path}")
        self.running = True
        return self

    def __exit__(
        self,
        exc_type: type | None,
        exc_val: BaseException | None,
        exc_tb: object,
    ) -> bool:
        # Guaranteed cleanup: runs whether exiting cleanly or via exception.
        self.running = False
        try:
            if self.socket is not None:
                self.socket.close()
        finally:
            try:
                if self._ctx is not None:
                    self._ctx.term()
            finally:
                with contextlib.suppress(OSError):
                    self.socket_path.unlink(missing_ok=True)
        return False  # Do not suppress exceptions


def is_calc_finished(proc: object, socket: zmq.Socket) -> bool:
    """Check if Gaussian calculation is finished or next ML step requested.

    Polls the ZMQ socket and the process exit status. Returns True when the
    Gaussian process has exited (calculation complete or crashed), False when
    a new message is available (ML step requested).

    Args:
        proc: subprocess.Popen instance for the running Gaussian process.
        socket: Bound zmq.REP socket from GaussianZMQServer.

    Returns:
        True if Gaussian process has exited, False if a new message arrived.
    """
    while True:
        # Poll socket for messages with 10 ms timeout; returns 0 if no message.
        if socket.poll(timeout=10) != 0:
            return False  # New message: another ML step requested
        elif proc.poll() is not None:
            return True  # Process exited: calculation finished or crashed
        else:
            time.sleep(1)
