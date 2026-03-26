# Phase 15: SLURM Integration & Batch Report - Context

**Gathered:** 2026-03-27
**Status:** Ready for planning

<domain>
## Phase Boundary

Two capabilities:
1. **SLURM DFT offload** — `mace-gaussian batch molecules.txt --dft-on-cluster user@host` submits DFT baseline jobs to a SLURM cluster via SSH, polls for completion, retrieves results, and resumes the ML pipeline automatically.
2. **Batch report** — `mace-gaussian report` generates a multi-molecule HTML report with aggregated accuracy metrics (R², RMSE) across all molecules and calculator combinations.

These are independent: the batch report reads existing `comparison_results/` and works regardless of whether DFT was run locally or on a cluster. Existing single-molecule analysis scripts (`run_analysis.py`, `run_analysis_harmonic.py`, coverage analysis) remain untouched.

</domain>

<decisions>
## Implementation Decisions

### SLURM Job Submission
- **D-01:** Template-based submission — a SLURM job script template (e.g. `slurm_dft.sh`) lives in the repo. `mace-gaussian` fills in molecule name, paths, and resource requests, then SCPs and sbatch's it. User edits the template for their cluster.
- **D-02:** All DFT jobs submitted upfront (batch sbatch), then poll periodically. SLURM handles scheduling. Maximizes cluster utilization.
- **D-03:** Fixed remote directory convention: `~/mace_gaussian_dft/{molecule}/` on the cluster. No flag needed.
- **D-04:** Passwordless SSH assumed (already set up per STATE.md).
- **D-05:** SCP for file transfer (no shared filesystem per STATE.md).

### SLURM Job Polling
- **D-06:** SSH + `sacct` polling to check job state. Poll interval ~hourly to avoid hammering the cluster with SSH connections.
- **D-07:** Blocking poll loop — submit all jobs, then enter poll loop with hourly status updates. Single command, walk away. When all done, retrieve results and continue pipeline.

### SLURM Result Retrieval
- **D-08:** SCP results back after job completion. Pull `.chk` files and convert to `.fchk` locally (formchk not available on cluster per STATE.md).
- **D-09:** DFT results land in existing `comparison_results/` structure (same as local DFT runs per STATE.md).

### Batch Report — Content & Structure
- **D-10:** Primary view is an accuracy leaderboard — big table of calculator combo x molecule, cells show R²/RMSE. Sortable. Highlights best/worst combos. Directly answers the thesis question.
- **D-11:** Four plot types included:
  - Summary heatmap (calculator combo x molecule, colored by RMSE/R²)
  - Box plots per calculator combo (distribution of RMSE/R² across molecules)
  - Size-scaling trend (RMSE vs number of atoms)
  - Per-molecule spectrum overlays (ML vs DFT)
- **D-12:** Navigation toolbar with shortcuts to each section (matching existing report style).
- **D-13:** Report reads existing `comparison_results/{molecule}/{combo}/results.json` files. No re-computation. Fast, decoupled from pipeline execution.

### Batch Report — Entry Point
- **D-14:** New CLI command: `mace-gaussian report [--output-dir ...]`. Independent of batch command — can be run any time results exist.

### CLI Integration
- **D-15:** SLURM invoked via flag on batch command: `mace-gaussian batch molecules.txt --dft-on-cluster user@host`. Extends existing batch command.
- **D-16:** Default behavior is blocking poll (submit → poll hourly → retrieve → continue). Single command workflow.

### Error Recovery
- **D-17:** Failed SLURM jobs: log error, mark molecule as failed in manifest, continue polling remaining jobs. Failed molecules retried on next `batch` run via manifest-based restart.
- **D-18:** SSH connection failures: exponential backoff retry. After N consecutive failures, warn user but keep trying (jobs still running on cluster).

### Claude's Discretion
- SLURM template content (module loads, resource defaults, Gaussian-specific settings)
- sacct query format and job state parsing
- Exact polling interval and backoff parameters
- Report HTML template and CSS styling
- Box plot and heatmap implementation details
- How to detect molecule size (atom count) for size-scaling plot
- manifest schema extensions for SLURM job IDs and remote paths

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Existing Pipeline & Batch Infrastructure
- `mace_gaussian/workflow.py` — 3-stage pipeline orchestrator; `run_dft_baselines()` is stage 2 (the "cluster seam" documented in module docstring)
- `mace_gaussian/batch.py` — Manifest-driven batch runner with per-calculator restart; SLURM extends this
- `mace_gaussian/cli.py` — Click CLI entry point; `batch` command to extend with `--dft-on-cluster`

### Results & Analysis
- `mace_gaussian/utils/results.py` — `ResultsManager` for per-molecule directory structure and results JSON
- `mace_gaussian/gaussian/parser.py` — Gaussian log parser for frequency/intensity data
- `mace_gaussian/gaussian/fchk.py` — `.chk` to `.fchk` conversion (needed after SCP retrieval)

### Existing Report Patterns
- `coverage_analysis/` — Frequency range coverage HTML report with embedded plots (Phase 13.4 pattern to follow)
- `run_analysis.py` — Single-molecule anharmonic analysis (existing, not replaced)
- `run_analysis_harmonic.py` — Single-molecule harmonic analysis (existing, not replaced)

### Requirements
- `.planning/REQUIREMENTS.md` §Batch Workflow — BATCH-05 (multi-molecule HTML report)
- `.planning/REQUIREMENTS.md` §HPC / SLURM — HPC-01, HPC-02 (SLURM submission + formchk)

### Prior Phase Context
- `.planning/phases/14-batch-runner-pubchem-fetcher/14-CONTEXT.md` — Batch runner decisions (D-08 through D-14)
- `.planning/phases/13.4-frequency-range-coverage-analysis/` — HTML report generation pattern

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `batch.py` manifest system (`load_manifest`, `save_manifest`, `parse_batch_file`) — extend for SLURM job tracking
- `ResultsManager` — reads/writes per-molecule results JSON; report aggregates these
- `coverage_analysis/` HTML report generation — template for self-contained HTML with embedded plots
- `scratch_dir()` context manager — pattern for temporary file management
- `convert_chk_to_fchk()` in `gaussian/fchk.py` — needed after SCP retrieval from cluster

### Established Patterns
- Click `@cli.command()` with `@click.option()` for CLI commands
- Atomic JSON writes via `tempfile.mkstemp` + `os.replace` (batch manifest pattern)
- `comparison_results/{molecule}/{energy}_{dipole}/results.json` directory structure
- Seaborn "colorblind" palette, DPI=300 PNG plots, Arial/Helvetica font

### Integration Points
- `cli.py` — new `report` command; `--dft-on-cluster` option added to existing `batch` command
- `batch.py` — SLURM submission/polling logic integrates with existing per-calculator manifest loop
- `comparison_results/` — report reads from same directory structure batch writes to
- `workflow.py` — `run_dft_baselines()` stage 2 replaced by SLURM offload when `--dft-on-cluster` is set

</code_context>

<specifics>
## Specific Ideas

- User explicitly wants hourly polling to avoid hammering the cluster with SSH calls — respect cluster etiquette
- Navigation toolbar in report should match existing report style (shortcuts to sections)
- All four plot types selected: heatmap, box plots, size-scaling, per-molecule spectra
- Report is thesis-ready: accuracy leaderboard answers "which calculator combo wins?"
- Existing single-molecule reports stay untouched — batch report is a new aggregation layer on top

</specifics>

<deferred>
## Deferred Ideas

- `--no-wait` flag for non-blocking submit-then-exit workflow (submit and retrieve as separate steps)
- Separate `slurm submit/status/retrieve` subcommand group
- Auto-retry failed SLURM jobs
- Per-molecule config for different calculator options in batch
- Interactive HTML spectrum viewer with Plotly (separate todo)

</deferred>

---

*Phase: 15-slurm-integration-batch-report*
*Context gathered: 2026-03-27*
