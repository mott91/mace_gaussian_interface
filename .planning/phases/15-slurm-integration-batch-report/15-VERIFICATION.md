---
phase: 15-slurm-integration-batch-report
verified: 2026-03-27T21:05:00Z
status: gaps_found
score: 9/11 must-haves verified
gaps:
  - truth: "ruff check passes on all modified files"
    status: failed
    reason: "Two E501 line-too-long violations introduced by phase 15 code"
    artifacts:
      - path: "mace_gaussian/batch.py"
        issue: "Line 284 is 102 chars (limit 100): generator expression with STATUS_COMPLETE in Stage 4"
      - path: "mace_gaussian/cli.py"
        issue: "Line 540 is 109 chars: --energy-calculators help text extended with mace_polar in phase 15"
    missing:
      - "Wrap batch.py line 284 generator across two lines"
      - "Split cli.py line 540 help string across continuation lines"
human_verification:
  - test: "End-to-end SLURM submission to real cluster"
    expected: "DFT jobs submitted via sbatch, polled via sacct, results SCP'd back to comparison_results/"
    why_human: "Requires real SSH target with SLURM; cannot verify without live cluster"
  - test: "Local formchk fallback triggers when cluster formchk fails"
    expected: "retrieve_results converts .chk to .fchk locally when cluster SCP of .fchk fails"
    why_human: "Requires simulated partial cluster failure; not testable without infrastructure"
---

# Phase 15: SLURM Integration and Batch Report — Verification Report

**Phase Goal:** Users can offload DFT baseline calculations to a SLURM cluster automatically, and a multi-molecule HTML report aggregates accuracy across all molecules and calculator combinations.
**Verified:** 2026-03-27T21:05:00Z
**Status:** gaps_found (2 ruff E501 violations; all functional goals achieved)
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User can run `mace-gaussian batch molecules.txt --dft-on-cluster user@host` and DFT jobs are submitted to SLURM via SSH | VERIFIED | `--dft-on-cluster` option present in batch CLI; `run_batch` branches to `submit_dft_jobs` when set |
| 2 | SLURM job script includes `formchk` so `.fchk` is produced on the cluster | VERIFIED | `templates/slurm_dft.sh` line 31: `formchk {chk_filename} {fchk_filename}` |
| 3 | Poll loop checks `sacct` hourly and prints status updates | VERIFIED | `poll_jobs` calls `sacct` via `_ssh_with_backoff`, prints status line each iteration, sleeps in 60s increments |
| 4 | Results are SCP'd back to `comparison_results/{molecule}/b3lyp_6-31Gdp/` | VERIFIED | `retrieve_results` SCPs `.log`, `.chk`, `.fchk` to `{results_dir}/{name}/b3lyp_6-31Gdp/` |
| 5 | A failed SLURM job marks that molecule `dft_failed` in manifest without halting | VERIFIED | `batch.py` line 388: `manifest["molecules"][mol_name]["dft_baseline"] = "dft_failed"` for non-COMPLETED states |
| 6 | SSH connection failures use exponential backoff retry | VERIFIED | `_ssh_with_backoff`: `min(2 ** attempt * 30, 600)`, max 5 retries, raises `ConnectionError` after exhaustion |
| 7 | User can run `mace-gaussian report` and get a multi-molecule HTML report | VERIFIED | `report` command in CLI, delegates to `generate_batch_report` |
| 8 | Report contains accuracy leaderboard, heatmap, box plots, size-scaling, spectrum overlays | VERIFIED | All 5 sections present in `_generate_html`; tests confirm HTML contains "Leaderboard" and base64 plots |
| 9 | Report reads existing `comparison_results/` JSON files without re-computation | VERIFIED | `aggregate_results` walks directory, reads `results.json` files; produces 79 rows across 14 molecules in live test |
| 10 | Report HTML is self-contained with embedded base64 plots | VERIFIED | `_encode_plot` uses `base64.b64encode`; test confirms `data:image/png;base64,` in output |
| 11 | `ruff check` passes on all phase 15 files | FAILED | Two E501 violations: `batch.py:284` (102 chars) and `cli.py:540` (109 chars) |

**Score:** 10/11 truths verified (11th is a code quality issue, not a functional failure)

---

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `mace_gaussian/slurm.py` | SLURM submission, polling, retrieval | VERIFIED | 482 lines; exports `submit_dft_jobs`, `poll_jobs`, `retrieve_results`, `TERMINAL_STATES` |
| `templates/slurm_dft.sh` | User-editable SLURM template with formchk | VERIFIED | Contains `#SBATCH`, `module load gaussian/g16`, `formchk {chk_filename} {fchk_filename}` |
| `mace_gaussian/batch.py` | Extended batch runner with SLURM integration | VERIFIED | `run_batch` signature has `dft_on_cluster: str | None`, `slurm_template: str | None`; SLURM block at lines 320-402 |
| `mace_gaussian/cli.py` | CLI with `--dft-on-cluster` on batch command | VERIFIED | Options present; `dft_on_cluster` and `slurm_template` passed through to `run_batch` |
| `mace_gaussian/analysis/batch_report.py` | Multi-molecule report generator | VERIFIED | 799 lines; exports `generate_batch_report`, `aggregate_results`, 5 plot functions |
| `tests/test_slurm.py` | Unit tests for SLURM module | VERIFIED | 13 tests; all pass; SSH/SCP mocked via `subprocess.run` patch |
| `tests/test_batch_report.py` | Integration tests for batch report | VERIFIED | 4 tests; all pass; uses real `comparison_results/` data |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `mace_gaussian/cli.py` | `mace_gaussian/batch.py` | `--dft-on-cluster` passed to `run_batch` | VERIFIED | `cli.py:620-621`: `dft_on_cluster=dft_on_cluster, slurm_template=slurm_template` |
| `mace_gaussian/batch.py` | `mace_gaussian/slurm.py` | `from .slurm import submit_dft_jobs, poll_jobs, retrieve_results` | VERIFIED | Lazy import at `batch.py:323-325` inside `if dft_on_cluster` branch |
| `mace_gaussian/slurm.py` | `templates/slurm_dft.sh` | reads template and fills placeholders | VERIFIED | `submit_dft_jobs` reads `template_path.read_text()` and calls `.format(molecule=..., remote_dir=..., ...)` |
| `mace_gaussian/analysis/batch_report.py` | `comparison_results/{molecule}/{combo}/results.json` | walks directory, loads JSON | VERIFIED | `aggregate_results` iterates `mol_dir.iterdir()`, loads `results.json`; confirmed with 79 rows of real data |
| `mace_gaussian/cli.py` | `mace_gaussian/analysis/batch_report.py` | `from mace_gaussian.analysis.batch_report import generate_batch_report` | VERIFIED | `cli.py:474`: lazy import inside `report` command function |

---

## Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|--------------|--------|--------------------|--------|
| `batch_report.py:_generate_html` | `df` (DataFrame) | `aggregate_results()` reads `results.json` files | Yes — 79 rows from 14 molecules in live test | FLOWING |
| `batch_report.py:_build_leaderboard_html` | `leaderboard` | `df.groupby("combo").agg(...)` | Yes — real R2/RMSE from JSON files | FLOWING |
| `batch_report.py:_plot_heatmap` | `pivot` | `df.pivot_table(...)` | Yes — real metric values per combo x molecule | FLOWING |

---

## Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| `slurm` module importable with all functions | `python -c "from mace_gaussian.slurm import submit_dft_jobs, poll_jobs, retrieve_results, TERMINAL_STATES; assert len(TERMINAL_STATES) == 7"` | Exit 0 | PASS |
| `batch` CLI shows `--dft-on-cluster` and `--slurm-template` | `mace-gaussian batch --help` | Both options present in output | PASS |
| `report` CLI shows `--results-dir` and `--output-dir` | `mace-gaussian report --help` | Both options present, exit 0 | PASS |
| `aggregate_results` produces real data | `aggregate_results("comparison_results")` → 79 rows, 14 molecules | R2 range 0.975–1.0, RMSE 0–129 cm⁻¹ | PASS |
| `generate_batch_report` creates valid HTML with embedded plots | `test_generate_batch_report_creates_html` | HTML >1000 bytes, contains "Leaderboard" and base64 data URIs | PASS |
| All 21 unit/integration tests pass | `pytest tests/test_slurm.py tests/test_batch_report.py -v` | 21 passed, 0 failed in 14.9s | PASS |
| `ruff check` on phase 15 files | `ruff check mace_gaussian/slurm.py mace_gaussian/batch.py mace_gaussian/cli.py mace_gaussian/analysis/batch_report.py` | 2 E501 violations found | FAIL |

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| HPC-01 | 15-01 | User can run `mace-gaussian batch molecules.txt --dft-on-cluster <host>` to submit DFT baseline jobs to SLURM cluster via SSH, poll for completion, and retrieve results automatically | SATISFIED | `--dft-on-cluster` on batch CLI; full submit/poll/retrieve workflow in `batch.py:320-402` |
| HPC-02 | 15-01 | SLURM submission includes `formchk` in the job script so `.fchk` is produced on the cluster without requiring local conversion | SATISFIED | `templates/slurm_dft.sh:31`: `formchk {chk_filename} {fchk_filename}`; local fallback in `retrieve_results` if cluster formchk fails |
| BATCH-05 | 15-02 | Batch run produces a multi-molecule HTML report with aggregated R² and RMSE per calculator combination across all molecules | SATISFIED | `mace-gaussian report` generates HTML with leaderboard (R²+RMSE), heatmap, box plots, size-scaling, spectrum overlays; reads existing JSON without re-computation |

All 3 required requirements are SATISFIED. No orphaned requirements found for phase 15.

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `mace_gaussian/batch.py` | 284 | E501: line too long (102 > 100) | Warning | Style violation only; logic is correct |
| `mace_gaussian/cli.py` | 540 | E501: line too long (109 > 100) | Warning | Style violation only; help text still displays correctly |

No stub patterns, empty implementations, hardcoded data, or TODO/FIXME markers found in any phase 15 file.

**Anti-pattern classification:** Both violations are style-only (line length). Neither prevents the goal. Classified as Warning.

---

## Human Verification Required

### 1. End-to-End SLURM Job Submission

**Test:** Run `mace-gaussian batch molecules.txt --dft-on-cluster user@rune03` against a real cluster node with Gaussian 16 and SLURM installed.
**Expected:** Jobs submitted via `sbatch`, job IDs printed, sacct polled hourly, `.log`/`.chk`/`.fchk` files SCP'd back to `comparison_results/{molecule}/b3lyp_6-31Gdp/`, `results.json` created.
**Why human:** Requires a live SLURM cluster with SSH access; all logic is unit-tested with mocks but cannot be validated end-to-end without infrastructure.

### 2. Local formchk Fallback Path

**Test:** Submit a job where `formchk` fails on the cluster (e.g., comment it out of the template) and verify `retrieve_results` falls back to local conversion.
**Expected:** `.fchk` file produced locally from `.chk` via the `convert_chk_to_fchk` call in `slurm.py:448`.
**Why human:** Requires simulating a partial cluster failure during SCP retrieval; not reproducible without a real cluster environment.

---

## Gaps Summary

One gap found — two `ruff` line-length violations (E501) in files modified by phase 15:

1. `batch.py:284` — the Stage 4 `complete_count` generator expression runs 102 chars. Fix: split the generator across the existing multi-line `sum(...)` structure.
2. `cli.py:540` — the `--energy-calculators` help text was extended in phase 15 to include `mace_polar`, pushing it to 109 chars. Fix: use a parenthesized string continuation.

These are code-quality violations only. All functional goals (HPC-01, HPC-02, BATCH-05) are fully achieved: the SLURM workflow is implemented and wired, the SLURM template contains `formchk`, the CLI exposes both new flags, the batch runner branches correctly, manifest tracking and `dft_failed` marking work, backoff retry is implemented, the batch report generates real HTML with embedded plots from live data, and all 21 tests pass.

---

_Verified: 2026-03-27T21:05:00Z_
_Verifier: Claude (gsd-verifier)_
