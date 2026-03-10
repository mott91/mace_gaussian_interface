# Phase 7: Extract Workflow Orchestrator - Context

**Gathered:** 2026-02-22
**Status:** Ready for planning

<domain>
## Phase Boundary

Extract the pipeline coordination logic from `gm_main.py` into a new `workflow.py` module. The orchestrator sequences three independent-until-comparison stages: geometry optimization, DFT baseline, and ML frequency calculations. `gm_main.py` is deleted. The CLI delegates to `workflow.py` functions. Pipeline behavior is unchanged — this is a structural extraction only.

</domain>

<decisions>
## Implementation Decisions

### workflow.py interface shape
- Expose **separate per-stage functions** plus a top-level coordinator:
  - `run_geometry_optimization()` — stage 1, already exists in gm_main.py
  - `run_dft_baselines()` — stage 2, wraps `dft_baseline.run_all_dft_baselines()`
  - `run_ml_calculations()` — stage 3, loops over energy/dipole combinations
  - `run_pipeline()` — top-level coordinator, calls the three stages in order
- CLI calls `run_pipeline()` (or individual stage functions) directly — no subprocess
- `workflow.py` is a pure Python module, not a CLI script

### DFT / ML independence (cluster seam)
- DFT baseline and ML calculations are **independent stages** — both depend only on the optimized geometry output, not on each other
- Each stage reads inputs from and writes results to disk (via ResultsManager) so stages can be run independently
- This creates a natural split point for future cluster execution: run DFT on cluster → copy results → continue from `run_ml_calculations()` or `run_pipeline()` with `include_dft_baselines=False`
- Full cluster automation (SLURM job submission, SSH, file sync) is **deferred** — design the seam now, implement cluster mechanics later

### gm_main.py removal
- **Hard delete** on day one — no deprecation shim, no redirect, clean break
- Low-level helpers in gm_main.py (`update_molecule_geometry`, `calculate_energy_and_forces`, `calculate_hessian`, `calculate_dipole_properties`, `run_next_calculation`, `geometry_optimisation`, `calculator`, `analyze_molecular_charges`, `setup_output_directory`) either move to `workflow.py` or are inlined — researcher to determine best home
- `print_diagnostics()` stays with CLI or moves to `diagnostics.py` if it doesn't already exist there

### Stage control / skip flags
- `run_pipeline()` accepts `include_dft_baselines: bool` (already in `run_workflow()`) to skip DFT
- Calculator selection (`energy_calculators`, `dipole_calculators`) remain as parameters
- `force_optimization: bool` remains to bypass cached geometry
- CLI flags map 1:1 to `run_pipeline()` parameters — no logic in CLI layer

### Claude's Discretion
- Where exactly each low-level helper ends up (workflow.py vs. another module)
- Whether `run_geometry_optimization` and `run_frequency_calculation` remain in workflow.py or become internal helpers
- Exact module-level docstrings and logging output format

</decisions>

<specifics>
## Specific Ideas

- `run_workflow()` in gm_main.py (line 595) is the direct source — nearly 1:1 becomes `run_pipeline()` in workflow.py
- The three-phase comment structure already in `run_workflow()` (PHASE 1 / PHASE 2 / PHASE 3) maps cleanly to the three stage functions
- The existing `--skip-dft-baseline` CLI flag should continue to work after refactoring

</specifics>

<deferred>
## Deferred Ideas

- **Cluster execution**: Automatic SLURM job submission for DFT baseline with SSH file sync — future phase or milestone
- **True parallelism**: Running DFT and ML simultaneously in the same invocation — deferred pending cluster design

</deferred>

---

*Phase: 07-extract-workflow-orchestrator*
*Context gathered: 2026-02-22*
