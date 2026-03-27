---
gsd_state_version: 1.0
milestone: v1.1
milestone_name: — Batch Benchmarking & Calculator Expansion
status: verifying
stopped_at: Completed 15-02-PLAN.md (batch accuracy report)
last_updated: "2026-03-27T20:00:17.545Z"
progress:
  total_phases: 10
  completed_phases: 7
  total_plans: 16
  completed_plans: 16
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-03)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** Phase 15 — SLURM Integration & Batch Report

## Current Position

Phase: 15 of 17 (SLURM Integration & Batch Report)
Plan: 2 of 2 complete
Status: All plans complete — ready for verification
Last activity: 2026-03-27 -- Plan 15-02 complete (batch accuracy report)

## Performance Metrics

**Velocity (v1.0):**

- Total plans completed: 32
- Average duration: ~3.9 min
- Total execution time: ~2.1 hours

**Recent Trend:**

- Last plans were fast (2–7 min each) on well-scoped work
- v1.1 plans expected similar cadence for Phase 13; Phase 15–16 may be longer (HPC/campaign work)

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting v1.1 work:

- [v1.0 complete]: mace_off and mace_anicc already in workflow.py::calculator() — CLI wiring only needed
- [v1.0 complete]: xTB dipole in calculators/xtb.py — unit bug (e*Bohr vs e*Angstrom) must be verified before production use
- [v1.1 start]: ORCA/VPT2 confirmed out of scope — no external hook in ORCA (research finding)
- [v1.1 start]: PubChem fetcher uses requests directly against PUG REST API (not pubchempy)
- [v1.1 start]: Batch SLURM uses SSH/SCP subprocess; qsub.sh on disk is the template
- [Phase 13-calculator-expansion-acoh-bug-fix]: Format B entries return ir_intensity=0.0 — ML external log format does not emit IR intensities; downstream plotting code handles this correctly
- [Phase 13-calculator-expansion-acoh-bug-fix]: Negative frequencies captured in Format B regex (imaginary modes from unconverged geometry are valid data points)
- [13-02]: mace_anicc uses mace_anicc(device="cuda") only — no model= or default_dtype= (TypeError if passed; different API from mace_mp/mace_off/mace_omol)
- [13-02]: Element guard placed OUTSIDE try/except in run_frequency_calculation so ValueError propagates instead of being silently caught
- [13-01]: Use callback= not type=click.Choice for comma-separated options (click.Choice validates atomically, rejecting "mace_mp,mace_omol" as a whole string)
- [13-01]: VALID_ENERGY_CALCULATORS and VALID_DIPOLE_CALCULATORS exported at module level for test assertions
- [13.1-01]: No try/except around get_dielectric_derivatives -- loud failure, no silent fallback to finite differences
- [13.1-01]: dmu_dr reshape via transpose(1,2,0).reshape(3N,3) matches base-class layout
- [13.1-01]: calculate_polarizability uses direct model forward pass since ASE calculate() path has commented-out polarizability extraction
- [13.1-02]: Unit conversion (Angstrom^3 -> Bohr^3) in workflow.py, not io.py -- io.py receives pre-converted values
- [13.1-02]: dalpha_dr stored in mol.info but NOT written to Gaussian output file (Python-only pipeline data)
- [13.1-02]: Voigt packing order [axx,axy,ayy,axz,ayz,azz] matches Gaussian convention
- [13.2-01]: scratch_dir uses @contextmanager with try/except BaseException/else pattern for cleanup
- [13.2-01]: MACE_IPC_PATH env var with .ipc_file fallback in gm_helper.py
- [13.2-01]: output_dir param uses basename-only %chk paths for Gaussian scratch isolation
- [13.2-02]: IPC socket path set via MACE_IPC_PATH in env dict passed to Popen, not global os.environ
- [13.2-02]: keep_scratch resolved from CLI flag first, then MACE_KEEP_SCRATCH env var as fallback
- [13.2-02]: Stale scratch cleanup runs on every mace-gaussian run invocation before pipeline starts
- [13.3-01]: Hungarian algorithm (linear_sum_assignment) replaces greedy loop in match_modes() for bijective 1-to-1 mode pairing
- [13.3-01]: Unmatched calc modes (n_calc > n_ref) get (None, 0.0) in matches dict; downstream filters with `if dft_idx is not None`
- [13.3-02]: Confidence threshold 0.5 for solid/dashed borders and filled/open markers — matches match_modes logging threshold
- [13.3-02]: Border color #4a4a4a and POINT_COLOR #5E81AC reused from existing styling for visual coherence
- [13.4-01]: pd.cut with right=False for left-inclusive region intervals [low, high); NaN metrics for empty regions (mode_count=0)
- [13.4-02]: YlOrRd colormap for heatmaps; empty cells masked then annotated "--" in gray; NaN-to-null JSON via recursive converter
- [Phase 13.5]: mace-torch upgraded to main (0.3.15) for mace_polar import; model='polar-1-l', no dispersion kwarg, opt-in only
- [Phase 14]: FileExistsError exits 0 (skip), ValueError exits 1 (not found) per D-03/D-04
- [Phase 14]: Batch calls three stage functions directly (not run_pipeline) for per-calculator restart granularity
- [Phase 14]: Atomic manifest write via tempfile.mkstemp + os.replace for crash-safe JSON
- [15-01]: formchk in SLURM template (cluster has g16) with local fallback if cluster formchk fails
- [15-01]: sacct polling uses 60s sleep increments for Ctrl-C interruptibility
- [15-01]: SLURM imports lazy inside dft_on_cluster conditional to avoid loading when not needed
- [Phase 15]: CSS extracted to _build_css() helper for line-length compliance; leaderboard uses green/red row highlighting

### Roadmap Evolution

- Phase 13.1 inserted after Phase 13: Calculator Acceleration & Polarizability Passthrough (INSERTED) — autograd dipole derivatives + polarizability passthrough discovered during code exploration 2026-03-10
- Phase 13.2 inserted after Phase 13: Temp file cleanup — scratch directory for intermediate Gaussian files (URGENT)
- Phase 13.3 inserted after Phase 13.2: Hungarian Optimal Mode Matching (INSERTED) — fix greedy mode matching correctness bug before benchmark
- Phase 13.4 inserted after Phase 13.3: Frequency Range Coverage Analysis (INSERTED) — thesis diagnostic for ML training set gaps by frequency region
- Phase 13.5 inserted after Phase 13.4: MACE-POLAR-1 Energy Calculator (INSERTED) — new energy calculator with long-range electrostatics

### Pending Todos

None yet.

### Blockers/Concerns

- **xTB dipole unit bug** (Phase 13): verify whether xtb.py divides by BOHR_TO_ANGSTROM before declaring xtb usable as dipole calculator; risk of ~1.89x factor error
- **mace_anicc model file**: confirm `ani500k_large_CC.model` exists on disk before Phase 13 planning
- **SLURM answers confirmed** (Phase 15): SCP (no shared filesystem); formchk NOT on cluster — pull `.chk` and convert locally; DFT results land in `comparison_results/` (existing structure, not a new dir); passwordless SSH is set up

## Session Continuity

Last session: 2026-03-27T20:00:17.541Z
Stopped at: Completed 15-02-PLAN.md (batch accuracy report)
