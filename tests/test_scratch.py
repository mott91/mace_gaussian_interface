"""Tests for scratch directory utility."""

import os
import time

import pytest

from mace_gaussian.utils.scratch import (
    SCRATCH_ROOT,
    STALE_THRESHOLD_SECONDS,
    cleanup_stale_scratch_dirs,
    list_stale_scratch_dirs,
    scratch_dir,
)


class TestScratchDir:
    """Tests for the scratch_dir context manager."""

    def test_creates_subdirectory_under_scratch(self, tmp_path):
        """scratch_dir() creates subdirectory under .scratch/ and yields its path."""
        with scratch_dir("test_run", base_dir=tmp_path) as sd:
            assert sd.exists()
            assert sd.parent.name == SCRATCH_ROOT
            assert sd.name == "test_run"
            assert sd.parent == tmp_path / SCRATCH_ROOT

    def test_deletes_directory_on_normal_exit(self, tmp_path):
        """scratch_dir() deletes directory on normal exit (success path)."""
        with scratch_dir("test_run", base_dir=tmp_path) as sd:
            # Create a file inside to ensure rmtree works
            (sd / "test.txt").write_text("hello")
            assert sd.exists()
        # After exiting, directory should be gone
        assert not sd.exists()

    def test_deletes_directory_on_exception(self, tmp_path):
        """scratch_dir() deletes directory on exception (failure path)."""
        with pytest.raises(ValueError, match="test error"):
            with scratch_dir("test_run", base_dir=tmp_path) as sd:
                (sd / "test.txt").write_text("hello")
                raise ValueError("test error")
        assert not sd.exists()

    def test_keeps_directory_on_failure_when_flag_set(self, tmp_path):
        """scratch_dir(keep_on_failure=True) preserves directory on exception."""
        with pytest.raises(ValueError, match="test error"):
            with scratch_dir("test_run", base_dir=tmp_path, keep_on_failure=True) as sd:
                (sd / "test.txt").write_text("hello")
                raise ValueError("test error")
        assert sd.exists()
        assert (sd / "test.txt").read_text() == "hello"

    def test_catches_base_exception(self, tmp_path):
        """scratch_dir() catches BaseException (KeyboardInterrupt), not just Exception."""
        with pytest.raises(KeyboardInterrupt):
            with scratch_dir("test_run", base_dir=tmp_path) as sd:
                (sd / "test.txt").write_text("hello")
                raise KeyboardInterrupt()
        # Directory should still be cleaned up
        assert not sd.exists()


class TestCleanupStaleScratchDirs:
    """Tests for cleanup_stale_scratch_dirs."""

    def test_removes_stale_dirs_keeps_recent(self, tmp_path):
        """cleanup_stale_scratch_dirs() removes dirs older than threshold, keeps recent ones."""
        scratch_root = tmp_path / SCRATCH_ROOT
        scratch_root.mkdir()

        # Create a stale directory (backdate mtime)
        stale_dir = scratch_root / "stale_run"
        stale_dir.mkdir()
        old_time = time.time() - STALE_THRESHOLD_SECONDS - 3600  # 1 hour past threshold
        os.utime(stale_dir, (old_time, old_time))

        # Create a recent directory
        recent_dir = scratch_root / "recent_run"
        recent_dir.mkdir()

        removed = cleanup_stale_scratch_dirs(base_dir=tmp_path)

        assert "stale_run" in removed
        assert not stale_dir.exists()
        assert recent_dir.exists()

    def test_removes_empty_scratch_parent(self, tmp_path):
        """cleanup_stale_scratch_dirs() removes empty .scratch/ parent after cleaning."""
        scratch_root = tmp_path / SCRATCH_ROOT
        scratch_root.mkdir()

        # Create only a stale directory
        stale_dir = scratch_root / "stale_run"
        stale_dir.mkdir()
        old_time = time.time() - STALE_THRESHOLD_SECONDS - 3600
        os.utime(stale_dir, (old_time, old_time))

        cleanup_stale_scratch_dirs(base_dir=tmp_path)

        # Parent .scratch/ should be removed since it's now empty
        assert not scratch_root.exists()


class TestListStaleScratchDirs:
    """Tests for list_stale_scratch_dirs."""

    def test_returns_stale_dirs_only(self, tmp_path):
        """list_stale_scratch_dirs() returns (name, age_hours) tuples for stale dirs only."""
        scratch_root = tmp_path / SCRATCH_ROOT
        scratch_root.mkdir()

        # Create a stale directory (backdate mtime by 30 hours)
        stale_dir = scratch_root / "stale_run"
        stale_dir.mkdir()
        old_time = time.time() - 30 * 3600  # 30 hours old
        os.utime(stale_dir, (old_time, old_time))

        # Create a recent directory
        recent_dir = scratch_root / "recent_run"
        recent_dir.mkdir()

        result = list_stale_scratch_dirs(base_dir=tmp_path)

        assert len(result) == 1
        name, age_hours = result[0]
        assert name == "stale_run"
        assert age_hours >= 29.0  # At least ~30 hours, with some tolerance
        # Should NOT include the recent dir
        assert all(name != "recent_run" for name, _ in result)
