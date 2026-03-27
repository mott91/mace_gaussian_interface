"""SLURM DFT offloading via SSH/SCP.

Submits DFT baseline jobs to a SLURM cluster, polls for completion
via sacct, and retrieves results via SCP. Designed for the batch
workflow: submit all jobs upfront, then poll periodically.

Requires passwordless SSH to the target host (key-based auth).
"""

from __future__ import annotations

import json
import logging
import subprocess
import time
from pathlib import Path

logger = logging.getLogger(__name__)

TERMINAL_STATES = frozenset(
    {
        "COMPLETED",
        "FAILED",
        "TIMEOUT",
        "CANCELLED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
    }
)

DEFAULT_REMOTE_BASE = "~/mace_gaussian_dft"
DEFAULT_POLL_INTERVAL = 3600  # 1 hour per user request

# Common SSH options: no password prompts, auto-accept new hosts
_SSH_OPTS = [
    "-o", "ConnectTimeout=30",
    "-o", "BatchMode=yes",
    "-o", "StrictHostKeyChecking=accept-new",
]


# ---------------------------------------------------------------------------
# SSH / SCP wrappers
# ---------------------------------------------------------------------------


def _ssh_run(
    host: str, command: str, timeout: int = 60
) -> subprocess.CompletedProcess:
    """Run a command on a remote host via SSH.

    Parameters
    ----------
    host : str
        SSH target, e.g. ``user@hostname``
    command : str
        Shell command to execute remotely
    timeout : int
        Local subprocess timeout in seconds

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd = ["ssh", *_SSH_OPTS, host, command]
    logger.debug("SSH: %s", " ".join(cmd))
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout
    )


def _scp_to(
    host: str, local: str, remote: str, timeout: int = 120
) -> subprocess.CompletedProcess:
    """Copy a local file to a remote host via SCP.

    Parameters
    ----------
    host : str
        SSH target
    local : str
        Local file path
    remote : str
        Remote destination path
    timeout : int
        Subprocess timeout in seconds

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd = ["scp", *_SSH_OPTS, local, f"{host}:{remote}"]
    logger.debug("SCP ->: %s", " ".join(cmd))
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout
    )


def _scp_from(
    host: str, remote: str, local: str, timeout: int = 120
) -> subprocess.CompletedProcess:
    """Copy a remote file to the local machine via SCP.

    Parameters
    ----------
    host : str
        SSH target
    remote : str
        Remote file path
    local : str
        Local destination path
    timeout : int
        Subprocess timeout in seconds

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd = ["scp", *_SSH_OPTS, f"{host}:{remote}", local]
    logger.debug("SCP <-: %s", " ".join(cmd))
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout
    )


def _ssh_with_backoff(
    host: str, command: str, max_retries: int = 5
) -> subprocess.CompletedProcess:
    """Run an SSH command with exponential backoff on connection failure.

    Retries on SSH connection errors (returncode 255) or subprocess
    timeout. Wait times: ``min(2 ** attempt * 30, 600)`` seconds.

    Parameters
    ----------
    host : str
        SSH target
    command : str
        Shell command to execute remotely
    max_retries : int
        Maximum number of retry attempts

    Returns
    -------
    subprocess.CompletedProcess

    Raises
    ------
    ConnectionError
        After exhausting all retries
    """
    last_error = None
    for attempt in range(max_retries + 1):
        try:
            result = _ssh_run(host, command, timeout=120)
            if result.returncode != 255:
                return result
            last_error = f"SSH connection failed (rc=255): {result.stderr.strip()}"
        except subprocess.TimeoutExpired as exc:
            last_error = f"SSH timed out: {exc}"

        if attempt < max_retries:
            wait = min(2 ** attempt * 30, 600)
            logger.warning(
                "SSH attempt %d/%d failed (%s), retrying in %ds...",
                attempt + 1,
                max_retries,
                last_error,
                wait,
            )
            time.sleep(wait)

    raise ConnectionError(
        f"SSH connection to {host} failed after {max_retries + 1} attempts: {last_error}"
    )


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def submit_dft_jobs(
    molecules: list[dict],
    host: str,
    template_path: Path,
    results_dir: str = "comparison_results",
) -> dict[str, str]:
    """Submit DFT jobs to a SLURM cluster via SSH.

    All jobs are submitted upfront before returning (per D-02).

    Parameters
    ----------
    molecules : list[dict]
        Each dict has keys ``name`` and ``gjf_path``.
    host : str
        SSH target, e.g. ``user@hostname``
    template_path : Path
        Path to SLURM job template with placeholders
    results_dir : str
        Local results directory (used for naming convention)

    Returns
    -------
    dict[str, str]
        Mapping of molecule name to SLURM job ID
    """
    template_text = template_path.read_text()
    job_ids: dict[str, str] = {}

    for mol in molecules:
        name = mol["name"]
        gjf_path = mol["gjf_path"]
        remote_dir = f"{DEFAULT_REMOTE_BASE}/{name}"
        gjf_filename = f"{name}_freq_anharm.gjf"
        chk_filename = f"{name}_freq_anharm.chk"
        fchk_filename = f"{name}_freq_anharm.fchk"

        # Fill template placeholders
        script = template_text.format(
            molecule=name,
            remote_dir=remote_dir,
            gjf_filename=gjf_filename,
            chk_filename=chk_filename,
            fchk_filename=fchk_filename,
        )

        # Create remote directory
        result = _ssh_run(host, f"mkdir -p {remote_dir}")
        if result.returncode != 0:
            logger.error(
                "Failed to create remote dir for %s: %s", name, result.stderr
            )
            continue

        # SCP the .gjf file
        result = _scp_to(host, gjf_path, f"{remote_dir}/{gjf_filename}")
        if result.returncode != 0:
            logger.error("Failed to SCP gjf for %s: %s", name, result.stderr)
            continue

        # Write filled SLURM script to a temp file and SCP it
        local_script = Path(gjf_path).parent / "slurm_dft.sh"
        local_script.write_text(script)
        result = _scp_to(
            host, str(local_script), f"{remote_dir}/slurm_dft.sh"
        )
        if result.returncode != 0:
            logger.error(
                "Failed to SCP slurm script for %s: %s", name, result.stderr
            )
            continue

        # Submit via sbatch
        result = _ssh_run(host, f"sbatch {remote_dir}/slurm_dft.sh")
        if result.returncode != 0:
            logger.error(
                "sbatch failed for %s: %s", name, result.stderr
            )
            continue

        # Extract job ID from "Submitted batch job NNNNN"
        job_id = result.stdout.strip().split()[-1]
        job_ids[name] = job_id
        logger.info("Submitted %s as SLURM job %s", name, job_id)

    return job_ids


def _parse_sacct_output(
    output: str, job_ids: dict[str, str]
) -> dict[str, str]:
    """Parse sacct output into a job-ID-to-state mapping.

    Handles step suffixes (e.g. ``12345.batch`` -> ``12345``) and
    keeps only the first entry per job ID (the main job row).

    Parameters
    ----------
    output : str
        Raw sacct ``--parsable2`` output
    job_ids : dict[str, str]
        Mapping of molecule name to expected job ID

    Returns
    -------
    dict[str, str]
        Mapping of job ID to SLURM state string
    """
    # Build reverse map: job_id -> molecule_name
    id_set = set(job_ids.values())
    states: dict[str, str] = {}

    for line in output.strip().splitlines():
        parts = line.split("|")
        if len(parts) < 2:
            continue
        raw_id = parts[0].split(".")[0]  # strip step suffix
        state = parts[1].strip()
        if raw_id in id_set and raw_id not in states:
            states[raw_id] = state

    return states


def poll_jobs(
    host: str,
    job_ids: dict[str, str],
    poll_interval: int = DEFAULT_POLL_INTERVAL,
) -> dict[str, str]:
    """Block until all SLURM jobs reach a terminal state.

    Polls via ``sacct`` over SSH with hourly status updates.
    Uses 60-second sleep increments for Ctrl-C interruptibility.

    Parameters
    ----------
    host : str
        SSH target
    job_ids : dict[str, str]
        Mapping of molecule name to SLURM job ID
    poll_interval : int
        Seconds between sacct polls (default: 3600)

    Returns
    -------
    dict[str, str]
        Mapping of molecule name to final SLURM state
    """
    # Reverse map: job_id -> molecule_name
    id_to_name = {jid: name for name, jid in job_ids.items()}
    all_ids = ",".join(job_ids.values())

    while True:
        sacct_cmd = (
            f"sacct -j {all_ids} "
            f"--format=JobID,State,ExitCode --noheader --parsable2"
        )
        result = _ssh_with_backoff(host, sacct_cmd)
        states = _parse_sacct_output(result.stdout, job_ids)

        # For job IDs not yet in sacct, treat as PENDING (accounting lag)
        for jid in job_ids.values():
            if jid not in states:
                states[jid] = "PENDING"

        # Count states
        n_complete = sum(1 for s in states.values() if s == "COMPLETED")
        n_failed = sum(1 for s in states.values() if s in TERMINAL_STATES and s != "COMPLETED")
        n_running = sum(1 for s in states.values() if s == "RUNNING")
        n_pending = sum(1 for s in states.values() if s not in TERMINAL_STATES and s != "RUNNING")

        # Check if all terminal
        all_terminal = all(s in TERMINAL_STATES for s in states.values())

        logger.info(
            "SLURM status: %d completed, %d running, %d pending, %d failed",
            n_complete,
            n_running,
            n_pending,
            n_failed,
        )
        print(
            f"SLURM status: {n_complete} completed, {n_running} running, "
            f"{n_pending} pending, {n_failed} failed"
        )

        if all_terminal:
            # Build final molecule_name -> state mapping
            return {
                id_to_name[jid]: state
                for jid, state in states.items()
                if jid in id_to_name
            }

        # Sleep in 60-second increments for Ctrl-C interruptibility
        slept = 0
        while slept < poll_interval:
            time.sleep(60)
            slept += 60


def retrieve_results(
    host: str,
    molecules: list[str],
    results_dir: str = "comparison_results",
) -> dict[str, bool]:
    """Retrieve DFT results from the cluster via SCP.

    Downloads ``.log``, ``.chk``, and ``.fchk`` files. If ``.fchk``
    was not produced on the cluster, converts locally using formchk.
    Parses the ``.log`` file and saves ``results.json``.

    Parameters
    ----------
    host : str
        SSH target
    molecules : list[str]
        Molecule names to retrieve
    results_dir : str
        Local base directory for results

    Returns
    -------
    dict[str, bool]
        Mapping of molecule name to retrieval success
    """
    from .gaussian.fchk import convert_chk_to_fchk
    from .gaussian.parser import parse_gaussian_log

    success_map: dict[str, bool] = {}

    for name in molecules:
        remote_dir = f"{DEFAULT_REMOTE_BASE}/{name}"
        local_dir = Path(results_dir) / name / "b3lyp_6-31Gdp"
        local_dir.mkdir(parents=True, exist_ok=True)

        base = f"{name}_freq_anharm"
        ok = True

        # Retrieve .log file (required)
        log_remote = f"{remote_dir}/{base}.log"
        log_local = str(local_dir / f"{base}.log")
        result = _scp_from(host, log_remote, log_local)
        if result.returncode != 0:
            logger.error("Failed to retrieve .log for %s: %s", name, result.stderr)
            success_map[name] = False
            continue

        # Retrieve .chk file
        chk_remote = f"{remote_dir}/{base}.chk"
        chk_local = str(local_dir / f"{base}.chk")
        result = _scp_from(host, chk_remote, chk_local)
        if result.returncode != 0:
            logger.warning("Failed to retrieve .chk for %s: %s", name, result.stderr)
            ok = False

        # Try to retrieve .fchk (may have been produced by formchk on cluster)
        fchk_remote = f"{remote_dir}/{base}.fchk"
        fchk_local = str(local_dir / f"{base}.fchk")
        result = _scp_from(host, fchk_remote, fchk_local)
        if result.returncode != 0:
            # formchk failed on cluster -- convert locally as fallback
            if Path(chk_local).exists():
                try:
                    convert_chk_to_fchk(chk_local, fchk_local)
                    logger.info("Converted .chk to .fchk locally for %s", name)
                except Exception as exc:
                    logger.error(
                        "Local formchk failed for %s: %s", name, exc
                    )
                    ok = False
            else:
                logger.warning("No .chk available for local formchk for %s", name)
                ok = False

        # Parse .log and save results.json
        try:
            parsed = parse_gaussian_log(log_local)
            results_json_path = local_dir / "results.json"
            results_data = {
                "molecule": name,
                "calculator_type": "dft",
                "method": "B3LYP/6-31G(d,p)",
                "source": "slurm_cluster",
                "host": host,
                "frequencies": parsed.get("frequencies", {}),
                "ir_intensities": parsed.get("ir_intensities", {}),
            }
            with results_json_path.open("w") as f:
                json.dump(results_data, f, indent=2)
            logger.info("Saved results.json for %s", name)
        except Exception as exc:
            logger.error("Failed to parse .log for %s: %s", name, exc)
            ok = False

        success_map[name] = ok

    return success_map
