"""Scratch directory utility for Gaussian intermediate files.

Provides a context manager that creates per-run scratch directories
under .scratch/ and cleans them up on exit. Includes utilities for
detecting and removing stale scratch directories from crashed runs.
"""

from __future__ import annotations

import logging
import shutil
import time
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)

SCRATCH_ROOT = ".scratch"
STALE_THRESHOLD_SECONDS = 24 * 3600  # 24 hours


@contextmanager
def scratch_dir(
    run_name: str,
    base_dir: str | Path = ".",
    keep_on_failure: bool = False,
) -> Generator[Path, None, None]:
    """Create a temporary scratch directory for a single run.

    Creates `{base_dir}/.scratch/{run_name}/` on entry and removes it
    on exit (both normal and exceptional). If ``keep_on_failure`` is
    True, the directory is preserved when an exception occurs.

    Args:
        run_name: Name for this run's scratch subdirectory.
        base_dir: Parent directory where .scratch/ lives.
        keep_on_failure: If True, keep the directory when an exception
            is raised (useful for debugging failed runs).

    Yields:
        Path to the created scratch directory.
    """
    base = Path(base_dir)
    scratch_path = base / SCRATCH_ROOT / run_name
    scratch_path.mkdir(parents=True, exist_ok=True)
    logger.debug("Created scratch directory: %s", scratch_path)

    try:
        yield scratch_path
    except BaseException:
        if keep_on_failure:
            logger.warning("Keeping scratch directory for debugging: %s", scratch_path)
        else:
            shutil.rmtree(scratch_path, ignore_errors=True)
            logger.debug("Cleaned up scratch directory after failure: %s", scratch_path)
        raise
    else:
        shutil.rmtree(scratch_path, ignore_errors=True)
        logger.debug("Cleaned up scratch directory: %s", scratch_path)


def cleanup_stale_scratch_dirs(base_dir: str | Path = ".") -> list[str]:
    """Remove scratch directories older than the stale threshold.

    Scans `{base_dir}/.scratch/` for subdirectories whose mtime is
    older than ``STALE_THRESHOLD_SECONDS`` and removes them. If the
    .scratch/ parent is empty after cleanup, it is also removed.

    Args:
        base_dir: Parent directory where .scratch/ lives.

    Returns:
        List of removed directory names.
    """
    scratch_root = Path(base_dir) / SCRATCH_ROOT
    if not scratch_root.exists():
        return []

    removed: list[str] = []
    now = time.time()

    for child in sorted(scratch_root.iterdir()):
        if not child.is_dir():
            continue
        age = now - child.stat().st_mtime
        if age > STALE_THRESHOLD_SECONDS:
            shutil.rmtree(child, ignore_errors=True)
            removed.append(child.name)
            logger.info("Removed stale scratch directory: %s (%.1fh old)", child.name, age / 3600)

    # Remove empty .scratch/ parent
    if scratch_root.exists() and not any(scratch_root.iterdir()):
        scratch_root.rmdir()
        logger.debug("Removed empty scratch root: %s", scratch_root)

    return removed


def list_stale_scratch_dirs(base_dir: str | Path = ".") -> list[tuple[str, float]]:
    """List scratch directories older than the stale threshold.

    Same scan as ``cleanup_stale_scratch_dirs`` but returns info
    without deleting anything.

    Args:
        base_dir: Parent directory where .scratch/ lives.

    Returns:
        List of (name, age_hours) tuples for stale directories.
    """
    scratch_root = Path(base_dir) / SCRATCH_ROOT
    if not scratch_root.exists():
        return []

    stale: list[tuple[str, float]] = []
    now = time.time()

    for child in sorted(scratch_root.iterdir()):
        if not child.is_dir():
            continue
        age = now - child.stat().st_mtime
        if age > STALE_THRESHOLD_SECONDS:
            stale.append((child.name, age / 3600))

    return stale
