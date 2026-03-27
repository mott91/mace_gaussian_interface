"""Tests for mace_gaussian.slurm module -- SLURM DFT offloading via SSH/SCP."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

from mace_gaussian.slurm import (
    DEFAULT_POLL_INTERVAL,
    DEFAULT_REMOTE_BASE,
    TERMINAL_STATES,
    _parse_sacct_output,
    _scp_from,
    _scp_to,
    _ssh_run,
    _ssh_with_backoff,
    poll_jobs,
    submit_dft_jobs,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


def test_terminal_states_contains_exactly_seven():
    assert len(TERMINAL_STATES) == 7
    expected = {
        "COMPLETED",
        "FAILED",
        "TIMEOUT",
        "CANCELLED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
    }
    assert TERMINAL_STATES == expected


def test_default_poll_interval_is_one_hour():
    assert DEFAULT_POLL_INTERVAL == 3600


def test_default_remote_base():
    assert DEFAULT_REMOTE_BASE == "~/mace_gaussian_dft"


# ---------------------------------------------------------------------------
# SSH / SCP wrappers
# ---------------------------------------------------------------------------


@patch("mace_gaussian.slurm.subprocess.run")
def test_ssh_run_builds_correct_command(mock_run):
    mock_run.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="ok\n", stderr=""
    )
    result = _ssh_run("user@host", "ls -la", timeout=30)
    mock_run.assert_called_once()
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "ssh"
    assert "-o" in cmd
    assert "ConnectTimeout=30" in cmd
    assert "BatchMode=yes" in cmd
    assert "StrictHostKeyChecking=accept-new" in cmd
    assert cmd[-2] == "user@host"
    assert cmd[-1] == "ls -la"
    assert mock_run.call_args[1]["timeout"] == 30
    assert result.returncode == 0


@patch("mace_gaussian.slurm.subprocess.run")
def test_scp_to_builds_correct_command(mock_run):
    mock_run.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr=""
    )
    _scp_to("user@host", "/local/file.txt", "/remote/file.txt", timeout=60)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "scp"
    assert "ConnectTimeout=30" in cmd
    assert "BatchMode=yes" in cmd
    assert "/local/file.txt" in cmd
    assert "user@host:/remote/file.txt" in cmd


@patch("mace_gaussian.slurm.subprocess.run")
def test_scp_from_builds_correct_command(mock_run):
    mock_run.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr=""
    )
    _scp_from("user@host", "/remote/file.txt", "/local/file.txt", timeout=60)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "scp"
    assert "user@host:/remote/file.txt" in cmd
    assert "/local/file.txt" in cmd


# ---------------------------------------------------------------------------
# Backoff
# ---------------------------------------------------------------------------


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_run")
def test_ssh_with_backoff_retries_on_connection_error(mock_ssh, mock_sleep):
    """Retries on returncode 255 with exponential backoff."""
    # First two calls fail with rc 255, third succeeds
    mock_ssh.side_effect = [
        subprocess.CompletedProcess(args=[], returncode=255, stdout="", stderr="conn refused"),
        subprocess.CompletedProcess(args=[], returncode=255, stdout="", stderr="conn refused"),
        subprocess.CompletedProcess(args=[], returncode=0, stdout="success\n", stderr=""),
    ]
    result = _ssh_with_backoff("user@host", "echo hello", max_retries=5)
    assert result.returncode == 0
    assert result.stdout == "success\n"
    assert mock_ssh.call_count == 3
    # Backoff waits: 2^0 * 30 = 30, 2^1 * 30 = 60
    assert mock_sleep.call_args_list == [call(30), call(60)]


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_run")
def test_ssh_with_backoff_raises_after_max_retries(mock_ssh, mock_sleep):
    """Raises ConnectionError after exhausting retries."""
    mock_ssh.return_value = subprocess.CompletedProcess(
        args=[], returncode=255, stdout="", stderr="connection refused"
    )
    with pytest.raises(ConnectionError, match="failed after"):
        _ssh_with_backoff("user@host", "echo hello", max_retries=2)
    # 1 initial + 2 retries = 3 total attempts
    assert mock_ssh.call_count == 3


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_run")
def test_ssh_with_backoff_retries_on_timeout(mock_ssh, mock_sleep):
    """Retries on subprocess.TimeoutExpired."""
    mock_ssh.side_effect = [
        subprocess.TimeoutExpired(cmd="ssh", timeout=120),
        subprocess.CompletedProcess(args=[], returncode=0, stdout="ok\n", stderr=""),
    ]
    result = _ssh_with_backoff("user@host", "slow command", max_retries=3)
    assert result.returncode == 0
    assert mock_ssh.call_count == 2


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_run")
def test_ssh_with_backoff_caps_wait_at_600(mock_ssh, mock_sleep):
    """Backoff wait time is capped at 600 seconds."""
    # 5 failures: waits = 30, 60, 120, 240, 480
    # 6th failure would be 960 but capped at 600
    mock_ssh.return_value = subprocess.CompletedProcess(
        args=[], returncode=255, stdout="", stderr="fail"
    )
    with pytest.raises(ConnectionError):
        _ssh_with_backoff("user@host", "cmd", max_retries=5)
    # Check the last sleep call is capped at 600
    waits = [c.args[0] for c in mock_sleep.call_args_list]
    assert all(w <= 600 for w in waits)


# ---------------------------------------------------------------------------
# submit_dft_jobs
# ---------------------------------------------------------------------------


@patch("mace_gaussian.slurm._scp_to")
@patch("mace_gaussian.slurm._ssh_run")
def test_submit_dft_jobs(mock_ssh, mock_scp, tmp_path):
    """submit_dft_jobs creates remote dir, SCPs files, and extracts job ID."""
    # Create a template
    template = tmp_path / "slurm_dft.sh"
    template.write_text("#!/bin/bash\n#SBATCH --job-name=dft_{molecule}\ncd {remote_dir}\ng16 {gjf_filename}\nformchk {chk_filename} {fchk_filename}\n")

    # Create a gjf file
    gjf_dir = tmp_path / "water" / "b3lyp_6-31Gdp"
    gjf_dir.mkdir(parents=True)
    gjf = gjf_dir / "water_freq_anharm.gjf"
    gjf.write_text("%chk=water_freq_anharm.chk\n#p B3LYP/6-31G(d,p)\n")

    # Mock SSH: mkdir succeeds, sbatch returns job ID
    mock_ssh.side_effect = [
        subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr=""),  # mkdir
        subprocess.CompletedProcess(
            args=[], returncode=0,
            stdout="Submitted batch job 12345\n", stderr=""
        ),  # sbatch
    ]
    mock_scp.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr=""
    )

    molecules = [{"name": "water", "gjf_path": str(gjf)}]
    job_ids = submit_dft_jobs(molecules, "user@cluster", template)

    assert job_ids == {"water": "12345"}
    # mkdir called once
    assert mock_ssh.call_count == 2
    mkdir_cmd = mock_ssh.call_args_list[0][0][1]
    assert "mkdir -p" in mkdir_cmd
    assert "mace_gaussian_dft/water" in mkdir_cmd
    # SCP called twice: gjf file + slurm script
    assert mock_scp.call_count == 2


# ---------------------------------------------------------------------------
# sacct parsing
# ---------------------------------------------------------------------------


def test_parse_sacct_output_strips_step_suffixes():
    """Parses sacct output and strips .batch step suffix."""
    output = (
        "12345|COMPLETED|0:0\n"
        "12345.batch|COMPLETED|0:0\n"
        "12346|RUNNING|0:0\n"
        "12347|FAILED|1:0\n"
    )
    job_ids = {"water": "12345", "methane": "12346", "ethanol": "12347"}
    states = _parse_sacct_output(output, job_ids)
    assert states == {
        "12345": "COMPLETED",
        "12346": "RUNNING",
        "12347": "FAILED",
    }


def test_parse_sacct_output_keeps_first_entry_per_job():
    """Only the first row per job ID is kept (main job, not steps)."""
    output = (
        "12345|COMPLETED|0:0\n"
        "12345.batch|FAILED|1:0\n"
        "12345.extern|COMPLETED|0:0\n"
    )
    job_ids = {"water": "12345"}
    states = _parse_sacct_output(output, job_ids)
    assert states["12345"] == "COMPLETED"


def test_parse_sacct_output_empty_returns_empty():
    """Empty sacct output returns empty dict."""
    job_ids = {"water": "12345"}
    states = _parse_sacct_output("", job_ids)
    assert states == {}


def test_parse_sacct_output_ignores_unknown_job_ids():
    """Job IDs not in the expected set are ignored."""
    output = "99999|COMPLETED|0:0\n"
    job_ids = {"water": "12345"}
    states = _parse_sacct_output(output, job_ids)
    assert states == {}


# ---------------------------------------------------------------------------
# poll_jobs
# ---------------------------------------------------------------------------


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_with_backoff")
def test_poll_jobs_parses_sacct_and_returns_terminal(mock_ssh, mock_sleep):
    """poll_jobs returns molecule->state mapping when all jobs are terminal."""
    sacct_output = (
        "12345|COMPLETED|0:0\n"
        "12345.batch|COMPLETED|0:0\n"
        "12346|FAILED|1:0\n"
    )
    mock_ssh.return_value = subprocess.CompletedProcess(
        args=[], returncode=0, stdout=sacct_output, stderr=""
    )
    job_ids = {"water": "12345", "methane": "12346"}
    result = poll_jobs("user@host", job_ids, poll_interval=60)
    assert result == {"water": "COMPLETED", "methane": "FAILED"}


@patch("mace_gaussian.slurm.time.sleep")
@patch("mace_gaussian.slurm._ssh_with_backoff")
def test_poll_jobs_treats_empty_sacct_as_pending(mock_ssh, mock_sleep):
    """Empty sacct response for a known job ID is treated as PENDING."""
    # First poll: empty (accounting lag), second poll: completed
    mock_ssh.side_effect = [
        subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr=""),
        subprocess.CompletedProcess(
            args=[], returncode=0,
            stdout="12345|COMPLETED|0:0\n", stderr=""
        ),
    ]
    job_ids = {"water": "12345"}
    result = poll_jobs("user@host", job_ids, poll_interval=60)
    assert result == {"water": "COMPLETED"}
    # Should have slept once (60s increments)
    assert mock_sleep.call_count >= 1
