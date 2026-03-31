---
gsd_state_version: 1.0
milestone: v1.2
milestone_name: -- Analysis Quality Overhaul
status: planning
stopped_at: Completed 19-01-PLAN.md
last_updated: "2026-03-31T09:09:39.573Z"
last_activity: 2026-03-28 -- v1.2 roadmap created (7 phases, 15 requirements)
progress:
  total_phases: 7
  completed_phases: 0
  total_plans: 2
  completed_plans: 1
  percent: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-28)

**Core value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.
**Current focus:** v1.2 -- Analysis Quality Overhaul, Phase 18 ready to plan

## Current Position

Phase: 18 of 24 (Spectral Quality Foundations)
Plan: --
Status: Ready to plan
Last activity: 2026-03-28 -- v1.2 roadmap created (7 phases, 15 requirements)

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity (v1.0 + v1.1):**

- Total plans completed: 48 (32 v1.0 + 16 v1.1)
- v1.0 average duration: ~3.9 min
- v1.1 average duration: ~5 min (more complex phases)
- Total execution time: ~3.4 hours

**Recent Trend:**

- v1.1 plans ranged 2-10 min; SLURM/batch plans took longer
- v1.2 phases are analysis-layer focused; expect 3-7 min per plan

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting v1.2 work:

- [v1.1 complete]: Phases 16-17 (benchmark campaign, docs) deferred to v1.3
- [v1.2 start]: Lorentzian broadening and zero-intensity filtering are upstream of all visual features
- [v1.2 start]: NIST overlay uses nistchempy + jcamp (two new deps); best-effort, never blocks analysis
- [v1.2 start]: VPT2 spike uses Psience (McCoy Group), not PyVPT2 (Psi4 hard dep blocks it)
- [v1.2 start]: Phase numbering starts at 18 (16-17 reserved for deferred v1.3 work)
- [Phase 19]: Subspace overlap normalized by min(calc, ref) count for partial Hungarian matches
- [Phase 19]: Block-level deduplication in parser replaces per-frequency deduplication

### Pending Todos

7 pending todos -- see `.planning/todos/pending/`

### Blockers/Concerns

- **mace_polar dipole failure**: Phase 18 includes investigation of mace_polar as dipole calculator (PIPE-02)
- **NIST coverage**: Not all 25 benchmark molecules will have gas-phase IR in NIST; Phase 21 must handle missing data gracefully
- **Psience feasibility unknown**: Phase 24 is explicitly a time-boxed spike; may conclude "not feasible"

## Session Continuity

Last session: 2026-03-31T09:09:39.570Z
Stopped at: Completed 19-01-PLAN.md
Resume file: None
