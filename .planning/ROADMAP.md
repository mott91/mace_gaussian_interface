# Roadmap: MACE-Gaussian Interface

## Milestones

- ✅ **v1.0 — Refactoring & Distribution** — Phases 1–12 (shipped 2026-02-28)
- 🚧 **v1.1 — Batch Benchmarking & Calculator Expansion** — Phases 13–17 (in progress)

## Phases

<details>
<summary>✅ v1.0 — Refactoring & Distribution (Phases 1–12) — SHIPPED 2026-02-28</summary>

- [x] **Phase 1: Testing Infrastructure & Characterization** — 4/4 plans — completed 2026-02-27
- [x] **Phase 2: Error Handling & Input Validation** — 3/3 plans — completed 2026-02-27
- [x] **Phase 3: Extract Utilities & Conventions** — 2/2 plans — completed 2026-02-27
- [x] **Phase 4: Extract Calculator Classes** — 2/2 plans — completed 2026-02-27
- [x] **Phase 5: Replace MACE Module Monkey-Patching** — 2/2 plans — completed 2026-02-19
- [x] **Phase 6: Extract Gaussian I/O & ZMQ Server** — 5/5 plans — completed 2026-02-20
- [x] **Phase 7: Extract Workflow Orchestrator** — 2/2 plans — completed 2026-02-24
- [x] **Phase 8: Package Structure & Reorganization** — 3/3 plans — completed 2026-02-24
- [x] **Phase 9: CI/CD & Distribution Prep** — 3/3 plans — completed 2026-02-26
- [x] **Phase 10: Documentation** — 4/4 plans — completed 2026-02-26
- [x] **Phase 11: Integration Wiring Fixes (Gap Closure)** — 1/1 plan — completed 2026-02-27
- [x] **Phase 12: Distribution Polish (Gap Closure)** — 1/1 plan — completed 2026-02-27

Full archive: `.planning/milestones/v1.0-ROADMAP.md`

</details>

---

### 🚧 v1.1 — Batch Benchmarking & Calculator Expansion (In Progress)

**Milestone Goal:** Expand from 4 to 7+ ML calculator combinations, add batch and HPC tooling, and run a systematic ~25-molecule benchmark campaign to answer the thesis question: does energy surface quality or dipole model quality dominate IR accuracy?

## Phase Details

### Phase 13: Calculator Expansion & Acoh Bug Fix
**Goal**: Users can run `mace-gaussian` with two new energy calculators (mace_off, mace_anicc) and the acetic acid frequency parser no longer fails on regression plots. (xTB deferred pending supervisor discussion.)
**Depends on**: Nothing (all code already partially exists; no new infrastructure required)
**Requirements**: CALC-01, CALC-02, FIX-01
**Plans**: 3 plans

Plans:
- [x] 13-01-PLAN.md — CLI validation: add mace_off/mace_anicc to --energy-calculators callback + --optimization-calculator Choice
- [x] 13-02-PLAN.md — workflow.py: add mace_anicc branch with correct API + element guard at call sites
- [x] 13-03-PLAN.md — Parser fix: dual-format anharmonic section detection + H/L-prefix regex + xfail removal

### Phase 13.1: Calculator Acceleration & Polarizability Passthrough (INSERTED)

**Goal**: MACEDipoleCalculator uses autograd (single model call) instead of 2×3N finite-difference passes for dipole derivatives; real polarizability is threaded through the pipeline and written to Gaussian output; dalpha_dr is available in the Python pipeline for future Raman analysis.
**Depends on:** Phase 13
**Requirements**: (implementation-driven, no formal requirement IDs)
**Plans**: 2 plans

Plans:
- [ ] 13.1-01-PLAN.md — mace_loader.py: add calculate_dipole_derivatives() autograd override + calculate_polarizability() method + _last_dalpha_dr storage
- [ ] 13.1-02-PLAN.md — workflow.py + io.py: expand calculate_dipole_properties() to 4-tuple, thread polarizability/dalpha_dr, write real polarizability to Gaussian output

### Phase 14: Batch Runner & PubChem Fetcher
**Goal**: Users can fetch 3D structures by name and run the full pipeline over a list of molecules with per-molecule failure isolation and restart safety.
**Depends on**: Phase 13 (calculator expansion ensures full set of combinations is available for batch)
**Requirements**: BATCH-01, BATCH-02, BATCH-03, BATCH-04
**Success Criteria** (what must be TRUE):
  1. User runs `mace-gaussian fetch aspirin` and receives `aspirin.xyz` with valid 3D coordinates in the current directory.
  2. User runs `mace-gaussian batch molecules.txt` and the pipeline executes sequentially for each molecule, with per-molecule results written to `batch_results/`.
  3. If the batch run is interrupted and restarted, already-complete molecules are skipped and the run resumes from the first incomplete molecule.
  4. `batch_manifest.json` exists after any batch run and records the status of every molecule (complete, failed, pending).
  5. User can run `mace-gaussian batch molecules.txt --skip-dft-baseline` to execute ML-only calculations without triggering DFT baseline jobs.
**Plans**: TBD

### Phase 15: SLURM Integration & Batch Report
**Goal**: Users can offload DFT baseline calculations to a SLURM cluster automatically, and a multi-molecule HTML report aggregates accuracy across all molecules and calculator combinations.
**Depends on**: Phase 14 (batch runner must exist before SLURM can extend it; multi-molecule data must exist before batch report has content)
**Requirements**: HPC-01, HPC-02, BATCH-05
**Success Criteria** (what must be TRUE):
  1. User runs `mace-gaussian batch molecules.txt --dft-on-cluster rune03` and DFT jobs are submitted to the SLURM cluster via SSH, polled until complete, and results retrieved without manual intervention.
  2. The SLURM job script includes `formchk` so `.fchk` files are produced on the cluster without requiring local conversion.
  3. A SLURM job failure for one molecule (walltime exceeded, Gaussian error) marks that molecule `dft_failed` in the manifest and does not halt the remaining batch.
  4. After a completed batch run, `batch_report.html` exists and contains aggregated R² and RMSE per calculator combination across all molecules that completed successfully.
**Plans**: TBD

### Phase 16: Benchmark Campaign Execution
**Goal**: A systematic ~25-molecule benchmark is run through the full pipeline with all 7+ calculator combinations, producing aggregated results that answer the thesis question about energy vs. dipole model dominance.
**Depends on**: Phase 15 (SLURM offload and batch report required to run benchmark at scale and aggregate results)
**Requirements**: BENCH-01, BENCH-02, BENCH-03
**Success Criteria** (what must be TRUE):
  1. The benchmark dataset of ~25 molecules (size-scaling series CH4→C10H22 + functional group diversity series) has been run through the full pipeline and results exist in `batch_results/`.
  2. All 7+ calculator combinations (existing 4 + mace_off×espaloma, mace_off×mace_ml, mace_anicc×espaloma, mace_anicc×mace_ml, xtb×xtb) appear in the batch report with populated R² and RMSE values.
  3. The batch report presents aggregated ML vs. DFT accuracy stratified by calculator combination and by molecule size/functional-group class, sufficient to support thesis conclusions.
**Plans**: TBD

### Phase 17: Architecture & Development Docs Update
**Goal**: The two developer-facing documentation files accurately reflect the current `mace_gaussian/` package layout, module names, and workflows — eliminating stale references to pre-refactor structure.
**Depends on**: Phase 16 (benchmark campaign may surface final architectural decisions; docs should reflect the settled state of the codebase)
**Requirements**: FIX-02, FIX-03
**Success Criteria** (what must be TRUE):
  1. `docs/ARCHITECTURE.md` references only module paths that exist in the current `mace_gaussian/` package (no pre-refactor names remain).
  2. `docs/DEVELOPMENT.md` correctly describes how to add a new calculator, run tests, and reproduce the benchmark workflow using the current package layout.
  3. A new contributor could follow `docs/DEVELOPMENT.md` step-by-step to add a new energy calculator without consulting source code directly.
**Plans**: TBD

---

## Progress

**Execution Order:** 13 → 13.1 → 14 → 15 → 16 → 17

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Testing Infrastructure | v1.0 | 4/4 | Complete | 2026-02-27 |
| 2. Error Handling | v1.0 | 3/3 | Complete | 2026-02-27 |
| 3. Extract Utilities | v1.0 | 2/2 | Complete | 2026-02-27 |
| 4. Extract Calculators | v1.0 | 2/2 | Complete | 2026-02-27 |
| 5. MACE Model Loading | v1.0 | 2/2 | Complete | 2026-02-19 |
| 6. Gaussian I/O & ZMQ | v1.0 | 5/5 | Complete | 2026-02-20 |
| 7. Workflow Orchestrator | v1.0 | 2/2 | Complete | 2026-02-24 |
| 8. Package Structure | v1.0 | 3/3 | Complete | 2026-02-24 |
| 9. CI/CD & Distribution | v1.0 | 3/3 | Complete | 2026-02-26 |
| 10. Documentation | v1.0 | 4/4 | Complete | 2026-02-26 |
| 11. Integration Wiring | v1.0 | 1/1 | Complete | 2026-02-27 |
| 12. Distribution Polish | v1.0 | 1/1 | Complete | 2026-02-27 |
| 13. Calculator Expansion & Acoh Fix | v1.1 | Complete    | 2026-03-03 | 2026-03-03 |
| 13.1. Calculator Acceleration & Polarizability | v1.1 | 0/2 | In planning | - |
| 14. Batch Runner & PubChem Fetcher | v1.1 | 0/TBD | Not started | - |
| 15. SLURM Integration & Batch Report | v1.1 | 0/TBD | Not started | - |
| 16. Benchmark Campaign | v1.1 | 0/TBD | Not started | - |
| 17. Docs Update | v1.1 | 0/TBD | Not started | - |
