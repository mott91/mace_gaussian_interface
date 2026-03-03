---
gsd_state_version: 1.0
milestone: v1.1
milestone_name: — Batch Benchmarking & Calculator Expansion
status: planning
stopped_at: Completed 13-03-PLAN.md (acoh anharmonic parser fix)
last_updated: "2026-03-03T17:35:01.953Z"
last_activity: 2026-03-03 — v1.1 roadmap created (Phases 13–17 defined, 17/17 requirements mapped)
progress:
  total_phases: 5
  completed_phases: 0
  total_plans: 3
  completed_plans: 1
  percent: 35
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-03)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** v1.1 — Phase 13: Calculator Expansion & Acoh Bug Fix

## Current Position

Phase: 13 of 17 (Calculator Expansion & Acoh Bug Fix)
Plan: 0 of TBD in current phase
Status: Ready to plan Phase 13
Last activity: 2026-03-03 — v1.1 roadmap created (Phases 13–17 defined, 17/17 requirements mapped)

Progress: [███████░░░░░░░░░░░░░] 35% (12/17 phases complete; v1.1 not started)

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

### Pending Todos

None yet.

### Blockers/Concerns

- **xTB dipole unit bug** (Phase 13): verify whether xtb.py divides by BOHR_TO_ANGSTROM before declaring xtb usable as dipole calculator; risk of ~1.89x factor error
- **mace_anicc model file**: confirm `ani500k_large_CC.model` exists on disk before Phase 13 planning
- **SLURM answers confirmed** (Phase 15): SCP (no shared filesystem); formchk NOT on cluster — pull `.chk` and convert locally; DFT results land in `comparison_results/` (existing structure, not a new dir); passwordless SSH is set up

## Session Continuity

Last session: 2026-03-03T17:35:01.951Z
Stopped at: Completed 13-03-PLAN.md (acoh anharmonic parser fix)
Resume file: None
