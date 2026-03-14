---
gsd_state_version: 1.0
milestone: v1.1
milestone_name: — Batch Benchmarking & Calculator Expansion
status: completed
stopped_at: Completed 13.2-02-PLAN.md (Phase 13.2 complete)
last_updated: "2026-03-14T19:11:25.003Z"
last_activity: 2026-03-14 — Phase 13.2 complete (Scratch directory workflow integration + CLI controls)
progress:
  total_phases: 7
  completed_phases: 3
  total_plans: 7
  completed_plans: 7
  percent: 100
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-03)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** v1.1 — Phase 13: Calculator Expansion & Acoh Bug Fix

## Current Position

Phase: 13.2 of 17 (Temp File Cleanup — Scratch Directory for Intermediate Gaussian Files)
Plan: 2 of 2 in current phase (COMPLETE)
Status: Phase 13.2 complete (all plans done)
Last activity: 2026-03-14 — Phase 13.2 complete (Scratch directory workflow integration + CLI controls)

Progress: [████████████████████] 100% (7/7 plans in v1.1 complete)

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

### Roadmap Evolution

- Phase 13.1 inserted after Phase 13: Calculator Acceleration & Polarizability Passthrough (INSERTED) — autograd dipole derivatives + polarizability passthrough discovered during code exploration 2026-03-10
- Phase 13.2 inserted after Phase 13: Temp file cleanup — scratch directory for intermediate Gaussian files (URGENT)

### Pending Todos

None yet.

### Blockers/Concerns

- **xTB dipole unit bug** (Phase 13): verify whether xtb.py divides by BOHR_TO_ANGSTROM before declaring xtb usable as dipole calculator; risk of ~1.89x factor error
- **mace_anicc model file**: confirm `ani500k_large_CC.model` exists on disk before Phase 13 planning
- **SLURM answers confirmed** (Phase 15): SCP (no shared filesystem); formchk NOT on cluster — pull `.chk` and convert locally; DFT results land in `comparison_results/` (existing structure, not a new dir); passwordless SSH is set up

## Session Continuity

Last session: 2026-03-14T19:06:28Z
Stopped at: Completed 13.2-02-PLAN.md (Phase 13.2 complete)
