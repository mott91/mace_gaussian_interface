# Phase 15: SLURM Integration & Batch Report - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-03-27
**Phase:** 15-slurm-integration-batch-report
**Areas discussed:** SLURM job lifecycle, Batch report content, CLI integration, Error recovery

---

## SLURM Job Lifecycle

### Submission Method

| Option | Description | Selected |
|--------|-------------|----------|
| Template .sh on disk | SLURM job script template in repo, filled in by mace-gaussian | ✓ |
| Inline generated script | Build SLURM script from CLI flags (--partition, --nodes, etc.) | |
| Config file approach | Cluster config YAML/TOML, script generated from config | |

**User's choice:** Template .sh on disk
**Notes:** Recommended approach — user edits template for their cluster.

### Job Completion Detection

| Option | Description | Selected |
|--------|-------------|----------|
| SSH + sacct polling | Periodically SSH and run sacct/squeue, poll every 30-60s | ✓ (modified) |
| Watch output file | Poll for sentinel file via SSH | |
| Blocking SSH wait | SSH stays open with srun --wait | |

**User's choice:** SSH + sacct polling, but with much less frequent polling (~hourly)
**Notes:** User explicitly wants to avoid hammering the cluster with SSH calls. Poll hourly, not every 30-60s.

### Job Parallelism

| Option | Description | Selected |
|--------|-------------|----------|
| All at once | Submit all DFT jobs upfront, poll until all complete | ✓ |
| One at a time | Submit one, wait, then submit next | |
| Configurable batch size | --slurm-concurrent N | |

**User's choice:** All at once
**Notes:** SLURM handles scheduling. Maximizes cluster utilization.

### Remote File Location

| Option | Description | Selected |
|--------|-------------|----------|
| User-specified remote dir | --remote-dir flag | |
| Fixed convention | ~/mace_gaussian_dft/{molecule}/ | ✓ |
| Tempdir on cluster | Unique temp dir per batch run | |

**User's choice:** Fixed convention
**Notes:** Predictable location, no extra flag needed.

---

## Batch Report Content

### Primary View

| Option | Description | Selected |
|--------|-------------|----------|
| Accuracy leaderboard | Calculator combo x molecule table with R²/RMSE, sortable | ✓ |
| Per-molecule deep dive | One section per molecule with its own plots | |
| Both views | Leaderboard summary + per-molecule sections | |

**User's choice:** Accuracy leaderboard
**Notes:** Directly answers the thesis question: which combo wins?

### Plot Types

| Option | Description | Selected |
|--------|-------------|----------|
| Summary heatmap | Calculator combo x molecule, colored by RMSE/R² | ✓ |
| Box plots per combo | Distribution across molecules per calculator combo | ✓ |
| Size-scaling trend | RMSE vs molecule size (atom count) | ✓ |
| Per-molecule spectra | Embedded ML vs DFT spectrum overlays | ✓ |

**User's choice:** All four types
**Notes:** User also requested navigation toolbar with shortcuts to each section, matching existing report style.

### Data Source

| Option | Description | Selected |
|--------|-------------|----------|
| Read existing results | Parse comparison_results/ JSON files | ✓ |
| Re-run analysis pipeline | Re-compute from raw .log/.fchk files | |

**User's choice:** Read existing results
**Notes:** User confirmed existing single-molecule reports stay untouched — batch report is an aggregation layer.

### Entry Point

| Option | Description | Selected |
|--------|-------------|----------|
| CLI command | mace-gaussian report [--output-dir ...] | ✓ |
| Standalone script | run_batch_report.py | |
| Auto-generated after batch | Batch auto-generates report at end | |

**User's choice:** CLI command
**Notes:** Consistent with fetch/batch/run CLI-first pattern.

---

## CLI Integration

### SLURM Invocation

| Option | Description | Selected |
|--------|-------------|----------|
| Flag on batch command | --dft-on-cluster user@host on existing batch command | ✓ |
| Separate submit/retrieve | Two separate commands | |
| Subcommand group | mace-gaussian slurm submit/status/retrieve | |

**User's choice:** Flag on batch command
**Notes:** Extends existing batch command naturally.

### Blocking Behavior

| Option | Description | Selected |
|--------|-------------|----------|
| Block with polling | Submit, poll hourly, retrieve, continue. Single command. | ✓ |
| Submit-then-exit | Submit and exit, retrieve later | |
| Both modes | --no-wait flag for non-blocking | |

**User's choice:** Block with polling
**Notes:** Walk-away workflow preferred.

---

## Error Recovery

### SLURM Job Failure

| Option | Description | Selected |
|--------|-------------|----------|
| Log and continue | Mark failed in manifest, continue polling others | ✓ |
| Abort all | Cancel remaining on first failure | |
| Auto-retry once | Resubmit failed jobs once | |

**User's choice:** Log and continue
**Notes:** Failed molecules retried on next batch run via manifest restart.

### SSH Connection Failure

| Option | Description | Selected |
|--------|-------------|----------|
| Retry with backoff | Exponential backoff, warn after N failures | ✓ |
| Fail immediately | Exit with error on SSH failure | |
| Switch to manual | Print manual retrieval instructions | |

**User's choice:** Retry with backoff
**Notes:** Jobs still running on cluster even if SSH is temporarily down.

---

## Claude's Discretion

- SLURM template content (module loads, resource defaults)
- sacct query format and job state parsing
- Exact polling interval and backoff parameters
- Report HTML template and CSS styling
- Box plot and heatmap implementation
- Manifest schema extensions for SLURM tracking

## Deferred Ideas

- --no-wait flag for non-blocking submit-then-exit
- Separate slurm subcommand group
- Auto-retry failed SLURM jobs
- Per-molecule config for different calculator options
- Interactive Plotly spectrum viewer
