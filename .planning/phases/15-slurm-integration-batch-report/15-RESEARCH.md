# Phase 15: SLURM Integration & Batch Report - Research

**Researched:** 2026-03-27
**Domain:** HPC job submission (SLURM via SSH/SCP), multi-molecule HTML report generation
**Confidence:** HIGH

## Summary

This phase has two independent deliverables: (1) SLURM DFT offloading via SSH/SCP with sacct-based polling, and (2) a multi-molecule HTML batch report aggregating accuracy metrics from existing `comparison_results/` JSON files.

The SLURM integration extends the existing `batch.py` manifest system. The batch runner already calls `run_dft_baselines()` as Stage 2 -- when `--dft-on-cluster` is set, this stage is replaced by: generate `.gjf` files locally, SCP them to the cluster, sbatch a templated SLURM script, poll via `sacct` over SSH, SCP results back, and convert `.chk` to `.fchk` locally. The batch report is a new `mace-gaussian report` CLI command that reads `comparison_results/{molecule}/{combo}/results.json` files and produces a self-contained HTML report with embedded plots (same pattern as `coverage_analysis/coverage_report.html`).

**Primary recommendation:** Implement SLURM and report as separate modules (`mace_gaussian/slurm.py` and `mace_gaussian/analysis/batch_report.py`), each with its own CLI command, sharing no state except the `comparison_results/` directory convention.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Template-based SLURM submission -- `slurm_dft.sh` template in repo, filled by code, user-editable
- **D-02:** All DFT jobs submitted upfront (batch sbatch), then poll periodically
- **D-03:** Fixed remote directory: `~/mace_gaussian_dft/{molecule}/` on cluster
- **D-04:** Passwordless SSH assumed
- **D-05:** SCP for file transfer (no shared filesystem)
- **D-06:** SSH + `sacct` polling for job state
- **D-07:** Blocking poll loop with hourly status updates
- **D-08:** SCP results back; pull `.chk` files and convert to `.fchk` locally (formchk not on cluster)
- **D-09:** DFT results land in existing `comparison_results/` structure
- **D-10:** Primary view is accuracy leaderboard table (calculator combo x molecule, R^2/RMSE)
- **D-11:** Four plot types: summary heatmap, box plots, size-scaling trend, per-molecule spectra
- **D-12:** Navigation toolbar with section shortcuts
- **D-13:** Report reads existing `comparison_results/` JSON files, no re-computation
- **D-14:** New CLI command: `mace-gaussian report [--output-dir ...]`
- **D-15:** SLURM invoked via `--dft-on-cluster user@host` flag on batch command
- **D-16:** Default behavior is blocking poll
- **D-17:** Failed jobs: log error, mark `dft_failed` in manifest, continue
- **D-18:** SSH failures: exponential backoff retry

### Claude's Discretion
- SLURM template content (module loads, resource defaults, Gaussian-specific settings)
- sacct query format and job state parsing
- Exact polling interval and backoff parameters
- Report HTML template and CSS styling
- Box plot and heatmap implementation details
- How to detect molecule size (atom count) for size-scaling plot
- Manifest schema extensions for SLURM job IDs and remote paths

### Deferred Ideas (OUT OF SCOPE)
- `--no-wait` flag for non-blocking submit-then-exit workflow
- Separate `slurm submit/status/retrieve` subcommand group
- Auto-retry failed SLURM jobs
- Per-molecule config for different calculator options in batch
- Interactive HTML spectrum viewer with Plotly
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| HPC-01 | `mace-gaussian batch molecules.txt --dft-on-cluster <host>` submits DFT jobs to SLURM via SSH, polls, retrieves | SLURM module with SSH subprocess calls; extends batch.py Stage 2; manifest tracks job IDs |
| HPC-02 | SLURM job script includes `formchk` so `.fchk` produced on cluster | **CONTRADICTS CONTEXT D-08** which says formchk NOT on cluster. Research finding: formchk IS available locally (`/usr/local/g16-c.02/g16/formchk`). D-08 says "pull .chk and convert locally". The SLURM template should include `formchk` IF the cluster has Gaussian (which it does -- the whole point is running DFT there). Recommend: include `formchk` in SLURM template (cluster has g16), also convert locally as fallback |
| BATCH-05 | Multi-molecule HTML report with aggregated R^2 and RMSE per calculator combination | batch_report.py reads comparison_results JSON, computes stats, generates self-contained HTML with embedded matplotlib/seaborn plots |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- **Code quality:** `ruff check --fix && ruff format` (line length 100), `ty check`
- **Plots:** PNG format, DPI=300, seaborn "colorblind" palette, Arial/Helvetica font, R^2 and RMSE in regression plots
- **File naming:** DFT dirs `b3lyp_6-31Gdp`, ML dirs `mace_mp_espaloma`
- **Results structure:** `comparison_results/{molecule}/{energy}_{dipole}/results.json`
- **Environment:** `micromamba activate mace4ir_v2` (not `.venv`)
- **Gaussian absolute paths:** Never use relative paths in `.gjf` files

## Standard Stack

### Core
| Library | Purpose | Why Standard |
|---------|---------|--------------|
| subprocess | SSH/SCP/sacct calls | stdlib, no external SSH libraries needed for simple commands |
| shutil | File operations | stdlib |
| json | Manifest and results I/O | Already used throughout |
| pathlib | Path handling | Already used throughout |
| time | Poll loop timing | stdlib |
| click | CLI extension | Already used for all commands |

### Supporting (for report)
| Library | Purpose | When to Use |
|---------|---------|-------------|
| matplotlib | Heatmap, box plots, size-scaling chart | Already in project deps |
| seaborn | Styling, colorblind palette | Already in project deps |
| pandas | Data aggregation across molecules | Already in project deps |
| base64 | Embed PNG plots in HTML | Already used in coverage_analysis |
| numpy | Statistics (R^2, RMSE computation) | Already in project deps |

### Not Needed
| Instead of | Don't Use | Reason |
|------------|-----------|--------|
| paramiko | subprocess ssh/scp | Adds dependency; passwordless SSH + subprocess is simpler and already decided |
| fabric | subprocess ssh/scp | Same reason |
| jinja2 | f-string HTML templates | Project already uses inline HTML generation (coverage_analysis pattern) |
| plotly | matplotlib | Interactive plots deferred; matplotlib matches existing report pattern |

## Architecture Patterns

### Recommended Project Structure
```
mace_gaussian/
  slurm.py                    # SLURM submission, polling, retrieval
  analysis/
    batch_report.py            # Multi-molecule report generator
    coverage_analysis.py       # Existing (untouched)
    html_report_generator.py   # Existing (untouched)
templates/
  slurm_dft.sh                 # SLURM job template (user-editable)
```

### Pattern 1: SLURM Module (`slurm.py`)

**What:** A module with functions for SSH-based SLURM job management.
**Key functions:**
- `submit_dft_jobs(molecules, host, template_path, results_dir)` -- generates .gjf locally, SCPs to cluster, sbatch's each
- `poll_jobs(host, job_ids, interval=3600)` -- blocking loop, SSH + sacct to check states
- `retrieve_results(host, molecules, results_dir)` -- SCP results back, run formchk locally
- `_ssh_run(host, command)` -- wrapper around `subprocess.run(["ssh", host, command])`
- `_scp_to(host, local, remote)` -- wrapper around `subprocess.run(["scp", local, f"{host}:{remote}"])`
- `_scp_from(host, remote, local)` -- wrapper around `subprocess.run(["scp", f"{host}:{remote}", local])`

**Integration with batch.py:**
```python
# In batch.py Stage 2, when dft_on_cluster is set:
if dft_on_cluster:
    from .slurm import submit_dft_jobs, poll_jobs, retrieve_results
    job_ids = submit_dft_jobs(molecules, host, ...)
    poll_jobs(host, job_ids)
    retrieve_results(host, molecules, results_dir)
else:
    run_dft_baselines(...)  # existing local path
```

### Pattern 2: SLURM Template (`slurm_dft.sh`)

**What:** A bash script template with placeholders that the code fills in.
**Example:**
```bash
#!/bin/bash
#SBATCH --job-name=dft_{molecule}
#SBATCH --output={remote_dir}/slurm_%j.out
#SBATCH --error={remote_dir}/slurm_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=batch

# Load Gaussian
module load gaussian/g16

cd {remote_dir}

# Run Gaussian
g16 {gjf_filename}

# Convert checkpoint to formatted checkpoint
formchk {chk_filename} {fchk_filename}
```

**Placeholders:** `{molecule}`, `{remote_dir}`, `{gjf_filename}`, `{chk_filename}`, `{fchk_filename}`

**User customization:** User edits `templates/slurm_dft.sh` for their cluster (partition name, module load command, walltime, memory).

### Pattern 3: sacct Job State Parsing

**What:** Parse SLURM job states via sacct over SSH.
**Command:**
```bash
ssh user@host "sacct -j <job_id1>,<job_id2> --format=JobID,State,ExitCode --noheader --parsable2"
```
**Output format:** `job_id|state|exit_code` (pipe-delimited with `--parsable2`)

**SLURM job states to handle:**
| State | Meaning | Action |
|-------|---------|--------|
| PENDING | Queued | Keep polling |
| RUNNING | Executing | Keep polling |
| COMPLETED | Finished successfully | Retrieve results |
| FAILED | Job error | Mark dft_failed |
| TIMEOUT | Walltime exceeded | Mark dft_failed |
| CANCELLED | User cancelled | Mark dft_failed |
| NODE_FAIL | Node failure | Mark dft_failed |
| OUT_OF_MEMORY | OOM killed | Mark dft_failed |

**Terminal states:** COMPLETED, FAILED, TIMEOUT, CANCELLED, NODE_FAIL, OUT_OF_MEMORY, PREEMPTED

### Pattern 4: Batch Report (`batch_report.py`)

**What:** Reads all `comparison_results/{molecule}/{combo}/results.json`, computes R^2 and RMSE per combo across molecules, generates HTML.
**Data flow:**
1. Walk `comparison_results/` to find all `results.json` files
2. For each molecule+combo: extract harmonic frequencies from results.json (ML) and from `b3lyp_6-31Gdp/results.json` (DFT reference)
3. Use existing harmonic analysis CSVs if available (`analysis_results_harmonic/{mol}/data/comparison_{combo}.csv`) -- these already have matched DFT vs ML frequency pairs with R^2 and RMSE
4. Aggregate into a DataFrame: rows = molecules, columns = combos, cells = R^2/RMSE
5. Generate plots: heatmap, box plots, size-scaling, per-molecule spectra
6. Build self-contained HTML with embedded base64 plots

**Atom count for size-scaling:** Read the optimized `.xyz` file from `comparison_results/{molecule}/geometry_opt/optimized.xyz` and count atoms with `ase.io.read()`. Or read `results.json` and count from there if geometry is not available.

### Pattern 5: Manifest Schema Extension

**Current manifest:**
```json
{
  "molecules": {
    "water": {
      "xyz_path": "/path/to/water.xyz",
      "geometry_opt": "complete",
      "dft_baseline": "complete",
      "combinations": { "mace_mp_espaloma": { "status": "complete" } }
    }
  }
}
```

**Extended for SLURM:**
```json
{
  "molecules": {
    "water": {
      "xyz_path": "/path/to/water.xyz",
      "geometry_opt": "complete",
      "dft_baseline": "complete",
      "slurm": {
        "job_id": "12345",
        "host": "user@rune03",
        "remote_dir": "~/mace_gaussian_dft/water",
        "status": "COMPLETED",
        "submitted_at": "2026-03-27T10:00:00",
        "completed_at": "2026-03-27T14:30:00"
      },
      "combinations": { ... }
    }
  }
}
```

### Anti-Patterns to Avoid
- **Do not use paramiko/fabric:** Subprocess SSH is simpler, already decided, and passwordless SSH is confirmed.
- **Do not poll more often than hourly:** Cluster etiquette -- user explicitly requested hourly intervals.
- **Do not re-run analysis in the report:** Report reads existing JSON/CSV data; no frequency calculations or mode matching.
- **Do not modify existing analysis scripts:** `run_analysis.py`, `run_analysis_harmonic.py` stay untouched.
- **Do not use `os.system()`:** Always use `subprocess.run()` with proper error handling and timeouts.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SSH/SCP wrappers | Paramiko bindings | `subprocess.run(["ssh", ...])` | Simpler, no new dep, passwordless SSH already works |
| HTML templating | Jinja2 engine | f-string + base64 embed pattern | Matches existing `coverage_analysis/coverage_report.html` |
| R^2 computation | Manual formula | `scipy.stats.linregress` or existing analysis code | Already used in run_analysis_harmonic.py |
| RMSE computation | Manual formula | `np.sqrt(np.mean((a-b)**2))` | One-liner, no library needed |
| Mode matching | Custom matcher | Existing `analysis_results_harmonic/` CSVs | Pre-computed by harmonic analysis |

## Common Pitfalls

### Pitfall 1: SSH Hangs on Unreachable Host
**What goes wrong:** `subprocess.run(["ssh", ...])` blocks indefinitely if host is unreachable.
**Why it happens:** SSH default timeout is very long or infinite.
**How to avoid:** Always use `-o ConnectTimeout=30 -o ServerAliveInterval=60 -o ServerAliveCountMax=3` flags.
**Warning signs:** Poll loop appears frozen with no output.

### Pitfall 2: sacct Returns Empty for Recent Jobs
**What goes wrong:** `sacct` may not show very recently submitted jobs (accounting database lag).
**Why it happens:** SLURM accounting daemon writes in batches, not real-time.
**How to avoid:** Store job IDs from sbatch output; if sacct returns empty for a known job ID, treat as PENDING (not failed). Add a short grace period (5-10 min) after submission before first poll.
**Warning signs:** Jobs reported as "not found" immediately after submission.

### Pitfall 3: SCP Glob Expansion
**What goes wrong:** `scp user@host:~/dir/*.chk .` may fail if remote shell does not expand globs.
**Why it happens:** Glob expansion depends on remote shell configuration.
**How to avoid:** Transfer specific files by name, not globs. Build the exact file list from the manifest.
**Warning signs:** "No such file" errors from SCP.

### Pitfall 4: formchk on Cluster vs Local
**What goes wrong:** Requirement HPC-02 says include formchk in SLURM script, but CONTEXT D-08 says convert locally.
**Why it happens:** Information collected at different times during discussion.
**How to avoid:** Include `formchk` in the SLURM template (the cluster has Gaussian, so formchk is available). Also convert locally as fallback if .fchk was not produced on cluster. This satisfies both HPC-02 and D-08.
**Warning signs:** `.fchk` missing after retrieval.

### Pitfall 5: Incomplete Results Directory After SCP
**What goes wrong:** DFT results land in wrong directory structure, not matching `comparison_results/{molecule}/b3lyp_6-31Gdp/`.
**Why it happens:** Remote directory convention differs from local convention.
**How to avoid:** After SCP retrieval, explicitly place files in the correct local directory using `ResultsManager.create_frequency_directory()` and save metadata via `results_mgr.save_frequency_results()`.
**Warning signs:** `mace-gaussian list` does not show the DFT results.

### Pitfall 6: Report Fails When No Harmonic Analysis Exists
**What goes wrong:** Report expects `analysis_results_harmonic/{mol}/data/comparison_*.csv` but harmonic analysis was never run for some molecules.
**Why it happens:** Batch runner runs harmonic analysis per-molecule (Stage 4), but it may fail or be skipped.
**How to avoid:** Report should compute R^2/RMSE directly from `comparison_results/` JSON data when CSVs are unavailable. Fall back to comparing harmonic frequency lists from DFT and ML results.json.
**Warning signs:** Report shows empty cells or crashes on missing CSVs.

### Pitfall 7: Hourly Poll Loop Blocks Signal Handling
**What goes wrong:** `time.sleep(3600)` blocks Ctrl-C handling.
**Why it happens:** Long sleep is not interruptible on some systems.
**How to avoid:** Use a loop with shorter sleeps (e.g., 60s intervals) that checks total elapsed time, or use `signal.signal(SIGINT, handler)` to catch interrupts.
**Warning signs:** User cannot Ctrl-C out of the poll loop.

## Code Examples

### SSH Command Wrapper
```python
import subprocess

def ssh_run(host: str, command: str, timeout: int = 60) -> subprocess.CompletedProcess:
    """Run a command on remote host via SSH."""
    return subprocess.run(
        [
            "ssh",
            "-o", "ConnectTimeout=30",
            "-o", "BatchMode=yes",
            "-o", "StrictHostKeyChecking=no",
            host,
            command,
        ],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
```

### sacct Parsing
```python
def parse_sacct_output(output: str) -> dict[str, str]:
    """Parse sacct --parsable2 output into {job_id: state} mapping."""
    states = {}
    for line in output.strip().splitlines():
        parts = line.split("|")
        if len(parts) >= 2:
            job_id = parts[0].split(".")[0]  # Strip step suffix (e.g., "12345.batch")
            state = parts[1]
            if job_id not in states:  # Keep first (main job) entry
                states[job_id] = state
    return states

TERMINAL_STATES = frozenset({
    "COMPLETED", "FAILED", "TIMEOUT", "CANCELLED",
    "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED",
})
```

### sbatch Submission and Job ID Extraction
```python
def submit_slurm_job(host: str, script_path: str) -> str:
    """Submit SLURM job and return job ID."""
    result = ssh_run(host, f"sbatch {script_path}", timeout=30)
    if result.returncode != 0:
        raise RuntimeError(f"sbatch failed: {result.stderr}")
    # Output: "Submitted batch job 12345"
    job_id = result.stdout.strip().split()[-1]
    return job_id
```

### Report Data Aggregation
```python
import json
from pathlib import Path
import numpy as np

def aggregate_results(results_dir: str = "comparison_results") -> dict:
    """Aggregate R^2 and RMSE across molecules and calculator combos."""
    base = Path(results_dir)
    data = {}  # {molecule: {combo: {r2: float, rmse: float}}}

    for mol_dir in sorted(base.iterdir()):
        if not mol_dir.is_dir():
            continue
        molecule = mol_dir.name
        data[molecule] = {}

        # Find DFT reference frequencies
        dft_dir = mol_dir / "b3lyp_6-31Gdp"
        if not dft_dir.exists():
            continue
        dft_json = dft_dir / "results.json"
        if not dft_json.exists():
            continue
        with dft_json.open() as f:
            dft_data = json.load(f)
        dft_freqs = [h["freq_cm"] for h in dft_data["frequencies"].get("harmonic", [])]

        # Compare each ML combo
        for combo_dir in sorted(mol_dir.iterdir()):
            if combo_dir.name in ("geometry_opt", "b3lyp_6-31Gdp"):
                continue
            ml_json = combo_dir / "results.json"
            if not ml_json.exists():
                continue
            with ml_json.open() as f:
                ml_data = json.load(f)
            ml_freqs = [h["freq_cm"] for h in ml_data["frequencies"].get("harmonic", [])]

            if len(dft_freqs) == len(ml_freqs) and len(dft_freqs) > 0:
                dft_arr = np.array(sorted(dft_freqs))
                ml_arr = np.array(sorted(ml_freqs))
                rmse = float(np.sqrt(np.mean((dft_arr - ml_arr) ** 2)))
                ss_res = np.sum((dft_arr - ml_arr) ** 2)
                ss_tot = np.sum((dft_arr - np.mean(dft_arr)) ** 2)
                r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
                data[molecule][combo_dir.name] = {"r2": r2, "rmse": rmse}

    return data
```

### Exponential Backoff for SSH
```python
import time

def ssh_with_backoff(host: str, command: str, max_retries: int = 5) -> subprocess.CompletedProcess:
    """SSH with exponential backoff on connection failure."""
    for attempt in range(max_retries):
        try:
            result = ssh_run(host, command, timeout=60)
            if result.returncode == 0 or result.returncode != 255:
                return result  # 255 = SSH connection error
            raise ConnectionError(f"SSH error: {result.stderr}")
        except (subprocess.TimeoutExpired, ConnectionError) as e:
            wait = min(2 ** attempt * 30, 600)  # 30s, 60s, 120s, 240s, 480s
            click.echo(f"  SSH failed (attempt {attempt + 1}/{max_retries}): {e}")
            click.echo(f"  Retrying in {wait}s...")
            time.sleep(wait)
    raise ConnectionError(f"SSH to {host} failed after {max_retries} attempts")
```

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| Manual DFT on cluster | Automated via --dft-on-cluster | Zero manual intervention |
| Per-molecule analysis only | Aggregated batch report | Thesis-ready cross-molecule comparison |
| Separate run + analyze steps | batch command handles both | Single-command workflow |

## Open Questions

1. **SLURM template `module load` command**
   - What we know: Cluster has Gaussian 16 (it runs DFT calculations)
   - What's unclear: Exact module name on the cluster (e.g., `module load gaussian/g16` vs `module load gaussian`)
   - Recommendation: Use `module load gaussian/g16` as default in template; user edits for their cluster. Add comment in template explaining this.

2. **formchk on cluster (HPC-02 vs D-08 conflict)**
   - What we know: HPC-02 requires formchk in SLURM script. D-08 says convert locally.
   - Recommendation: Do both. Include formchk in SLURM template. After SCP retrieval, check if .fchk exists; if not, convert locally with formchk. This satisfies both requirements.

3. **R^2/RMSE from raw results vs pre-computed CSVs**
   - What we know: Harmonic analysis CSVs have pre-computed matched pairs. Raw results.json has frequency lists.
   - What's unclear: Whether harmonic analysis will have been run for all molecules before report.
   - Recommendation: Primary path reads existing CSVs. Fallback computes from results.json frequency lists using simple sorted-frequency matching (less accurate but functional).

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| ssh | SLURM submission | Yes | OpenSSH 8.4p1 | -- |
| scp | File transfer | Yes | (bundled with OpenSSH) | -- |
| formchk | .chk to .fchk conversion | Yes | Gaussian 16 C.02 | -- |
| matplotlib | Report plots | Yes | (in project deps) | -- |
| seaborn | Report styling | Yes | (in project deps) | -- |
| pandas | Data aggregation | Yes | (in project deps) | -- |
| g16 | DFT on cluster | Yes (local) | C.02 | Cluster has it too |

**Missing dependencies with no fallback:** None.

## Sources

### Primary (HIGH confidence)
- Existing codebase: `batch.py`, `cli.py`, `workflow.py`, `dft_baseline.py` -- direct inspection of integration points
- Existing codebase: `coverage_analysis/coverage_report.html` -- HTML report pattern with embedded base64 plots
- Existing codebase: `mace_gaussian/analysis/coverage_analysis.py` -- plot generation with matplotlib/seaborn
- Existing codebase: `comparison_results/*/results.json` -- data schema for report aggregation

### Secondary (MEDIUM confidence)
- SLURM documentation for sacct output format and job states -- well-established, stable API
- SSH/SCP subprocess patterns -- standard Unix tools, widely documented

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already in project, no new dependencies
- Architecture: HIGH -- extends existing patterns (manifest, CLI, HTML report)
- Pitfalls: HIGH -- based on direct codebase inspection and known SLURM behaviors

**Research date:** 2026-03-27
**Valid until:** 2026-04-27 (stable domain, no rapidly changing dependencies)
