# Roadmap: MACE-Gaussian Interface

## Milestones

- ✅ **v1.0 -- Refactoring & Distribution** -- Phases 1-12 (shipped 2026-02-28)
- ✅ **v1.1 -- Batch Benchmarking & Calculator Expansion** -- Phases 13-15 (shipped 2026-03-28)
- 🚧 **v1.2 -- Analysis Quality Overhaul** -- Phases 18-24 (in progress)
- 📋 **v1.3 -- Benchmark Campaign & Docs** -- Phases 16-17 (planned, deferred from v1.1)

## Phases

<details>
<summary>v1.0 -- Refactoring & Distribution (Phases 1-12) -- SHIPPED 2026-02-28</summary>

- [x] **Phase 1: Testing Infrastructure & Characterization** -- 4/4 plans -- completed 2026-02-27
- [x] **Phase 2: Error Handling & Input Validation** -- 3/3 plans -- completed 2026-02-27
- [x] **Phase 3: Extract Utilities & Conventions** -- 2/2 plans -- completed 2026-02-27
- [x] **Phase 4: Extract Calculator Classes** -- 2/2 plans -- completed 2026-02-27
- [x] **Phase 5: Replace MACE Module Monkey-Patching** -- 2/2 plans -- completed 2026-02-19
- [x] **Phase 6: Extract Gaussian I/O & ZMQ Server** -- 5/5 plans -- completed 2026-02-20
- [x] **Phase 7: Extract Workflow Orchestrator** -- 2/2 plans -- completed 2026-02-24
- [x] **Phase 8: Package Structure & Reorganization** -- 3/3 plans -- completed 2026-02-24
- [x] **Phase 9: CI/CD & Distribution Prep** -- 3/3 plans -- completed 2026-02-26
- [x] **Phase 10: Documentation** -- 4/4 plans -- completed 2026-02-26
- [x] **Phase 11: Integration Wiring Fixes (Gap Closure)** -- 1/1 plan -- completed 2026-02-27
- [x] **Phase 12: Distribution Polish (Gap Closure)** -- 1/1 plan -- completed 2026-02-27

Full archive: `.planning/milestones/v1.0-ROADMAP.md`

</details>

<details>
<summary>v1.1 -- Batch Benchmarking & Calculator Expansion (Phases 13-15) -- SHIPPED 2026-03-28</summary>

- [x] **Phase 13: Calculator Expansion & Acoh Bug Fix** -- 3/3 plans -- completed 2026-03-03
- [x] **Phase 13.1: Calculator Acceleration & Polarizability Passthrough** -- 2/2 plans -- completed 2026-03-11
- [x] **Phase 13.2: Temp File Cleanup / Scratch Dir** -- 2/2 plans -- completed 2026-03-14
- [x] **Phase 13.3: Hungarian Optimal Mode Matching** -- 2/2 plans -- completed 2026-03-15
- [x] **Phase 13.4: Frequency Range Coverage Analysis** -- 2/2 plans -- completed 2026-03-16
- [x] **Phase 13.5: MACE-POLAR-1 Energy Calculator** -- 1/1 plans -- completed 2026-03-23
- [x] **Phase 14: Batch Runner & PubChem Fetcher** -- 2/2 plans -- completed 2026-03-23
- [x] **Phase 15: SLURM Integration & Batch Report** -- 2/2 plans -- completed 2026-03-27

Full archive: `.planning/milestones/v1.1-ROADMAP.md`, `.planning/milestones/v1.1-REQUIREMENTS.md`

</details>

---

### v1.2 -- Analysis Quality Overhaul (In Progress)

**Milestone Goal:** Overhaul the anharmonic analysis pipeline, add spectral comparison capabilities, and fix data quality issues -- so that the benchmark campaign (v1.3) produces thesis-ready results.

- [x] **Phase 18: Spectral Quality Foundations** - Lorentzian broadening, zero-intensity filtering, stick+broadened spectra, mace_polar dipole investigation (completed 2026-03-30)
- [x] **Phase 19: Degenerate Mode Handling** - Subspace overlap for degenerate modes, correct statistics without double-counting (completed 2026-03-31)
- [x] **Phase 20: Wall-Clock Timing** - Per-molecule per-calculator timing instrumentation and batch report timing table (completed 2026-04-02)
- [x] **Phase 21: NIST Experimental Overlay** - Fetch, cache, and overlay experimental IR spectra from NIST WebBook (completed 2026-04-02)
- [ ] **Phase 22: Early SLURM Submission** - Submit DFT jobs per-molecule immediately after geometry optimization
- [ ] **Phase 23: Anharmonic Pipeline & Report Overhaul** - Thesis-quality HTML report integrating all v1.2 features
- [ ] **Phase 24: VPT2 Research Spike** - Time-boxed Psience/VPT2 feasibility evaluation on water

## Phase Details

### Phase 18: Spectral Quality Foundations
**Goal**: Simulated IR spectra are physically correct (Lorentzian line shapes) and intensity statistics are meaningful (zero-intensity modes excluded from regression)
**Depends on**: Nothing (first phase of v1.2; all code changes are in analysis layer)
**Requirements**: SPEC-01, SPEC-02, SPEC-03, PIPE-02
**Success Criteria** (what must be TRUE):
  1. Running analysis on any molecule produces a Lorentzian-broadened IR spectrum plot with configurable FWHM (default 10 cm-1), replacing the previous Gaussian KDE broadening
  2. Analysis report displays both a stick spectrum (discrete lines at each mode frequency) and a broadened Lorentzian spectrum side by side or overlaid
  3. Intensity regression metrics (R-squared, RMSE) exclude modes with DFT IR intensity below 0.1 km/mol, while frequency metrics still include all modes
  4. The mace_polar dipole calculator failure mode is either fixed (produces correct dipole derivatives) or documented with explicit skip logic and a clear error message explaining why it fails
**Plans**: 2 plans
Plans:
- [x] 18-01-PLAN.md -- Lorentzian broadening, FWHM default change, CLI --fwhm wiring
- [x] 18-02-PLAN.md -- Zero-intensity filtering for intensity regression, HTML report methodology update
**UI hint**: yes
**Note**: SPEC-03 (stick spectrum) deferred per D-01; PIPE-02 (mace_polar dipole) deferred per user request

### Phase 19: Degenerate Mode Handling
**Goal**: Mode matching correctly handles degenerate modes (e.g., methane T2, ammonia E) without systematically pessimistic overlap scores
**Depends on**: Phase 18 (spectral foundations should be stable before refining mode matching)
**Requirements**: MODE-05, MODE-06
**Success Criteria** (what must be TRUE):
  1. Degenerate modes (frequencies within 5 cm-1 threshold) are automatically detected and grouped in mode matching output
  2. Degenerate groups use subspace overlap (trace of M^T M / k) instead of individual eigenvector dot products, producing overlap scores that reflect true subspace alignment
  3. Mode matching statistics (average overlap, count of confident matches) correctly account for degenerate groups without double-counting individual modes within a group
**Plans**: 2 plans
Plans:
- [x] 19-01-PLAN.md -- Degenerate group detection, subspace overlap, group-aware statistics, parser fix
- [x] 19-02-PLAN.md -- Collapsed heatmaps, diamond-marker regression plots, HTML report labels

### Phase 20: Wall-Clock Timing
**Goal**: Every pipeline run records wall-clock timing per stage, and batch reports include ML vs DFT cost comparison
**Depends on**: Phase 18 (timing instrumentation is independent but should integrate with the analysis layer after spectral foundations are stable)
**Requirements**: TIME-01, TIME-02
**Success Criteria** (what must be TRUE):
  1. After running `mace-gaussian run molecule.xyz`, the results JSON contains wall-clock timing for geometry optimization, frequency calculation, and total pipeline time, broken down per calculator combination
  2. After a batch run, the batch report HTML includes a timing comparison table showing ML wall-clock time, DFT wall-clock time, and speedup factor (DFT/ML) per molecule
  3. Timing values are reproducible (not dominated by GPU warmup noise) -- warmup calls and CUDA synchronization are applied before timing starts
**Plans**: TBD

### Phase 21: NIST Experimental Overlay
**Goal**: Users can compare computed spectra against real experimental data from the NIST WebBook
**Depends on**: Phase 18 (Lorentzian broadening must exist for visually meaningful comparison with experimental spectra)
**Requirements**: NIST-01, NIST-02, NIST-03
**Success Criteria** (what must be TRUE):
  1. User can fetch an experimental IR spectrum from NIST WebBook by molecule name, and the spectrum is cached locally (no re-download on subsequent runs)
  2. When experimental data is available, the analysis report overlays the experimental spectrum on the computed ML and DFT spectra plots
  3. When experimental data is available, a peak position comparison table shows experimental vs computed peak positions with error metrics (MAE, RMSE in cm-1)
  4. When experimental data is unavailable for a molecule, the analysis runs normally without error -- the experimental overlay section is simply absent
**Plans**: TBD
**UI hint**: yes

### Phase 22: Early SLURM Submission
**Goal**: Batch runs submit DFT jobs to the cluster as early as possible, overlapping DFT queue time with local ML calculations
**Depends on**: Phase 20 (timing instrumentation should be in place so early submission benefits are measurable)
**Requirements**: PIPE-01
**Success Criteria** (what must be TRUE):
  1. During `mace-gaussian batch molecules.txt --dft-on-cluster host`, each molecule's DFT job is submitted to SLURM immediately after its geometry optimization completes, not after all molecules' ML calculations finish
  2. ML frequency calculations for subsequent molecules proceed in parallel with DFT queue time on the cluster
  3. The batch manifest correctly tracks both ML and DFT status per molecule without race conditions (concurrent manifest access is safe)
**Plans**: TBD

### Phase 23: Anharmonic Pipeline & Report Overhaul
**Goal**: The anharmonic analysis pipeline produces a thesis-quality HTML report that integrates all v1.2 features into a polished presentation
**Depends on**: Phase 18, Phase 19, Phase 20, Phase 21 (all analysis features must be stable before the capstone integration)
**Requirements**: ANAL-01, ANAL-02
**Success Criteria** (what must be TRUE):
  1. The anharmonic analysis HTML report integrates Lorentzian spectra, experimental overlay (when available), timing breakdown, and degenerate-mode-aware mode matching into a single cohesive document
  2. The report includes per-molecule summary cards showing key metrics: frequency R-squared, intensity R-squared, RMSE, speedup factor, and experimental agreement (when available)
  3. The report is visually thesis-ready: consistent styling, publication-quality plots (300 DPI, proper axis labels), and a clear narrative flow from molecule overview to detailed comparison
**Plans**: TBD
**UI hint**: yes

### Phase 24: VPT2 Research Spike
**Goal**: A time-boxed investigation determines whether the Psience library can provide an alternative VPT2 engine using MACE force constants
**Depends on**: Nothing (isolated research spike, can run any time but sequenced last to avoid blocking core features)
**Requirements**: VPT2-01
**Success Criteria** (what must be TRUE):
  1. A 2-hour spike has been completed that attempts to run Psience VPT2 on water using MACE-computed force constants from a Gaussian .fchk file
  2. The spike produces a written assessment: either a proof-of-concept showing VPT2 anharmonic frequencies from Psience that can be compared against Gaussian VPT2 results, or a documented explanation of why integration is blocked (missing API, Fermi resonance handling, force constant format incompatibility)
**Plans**: TBD

---

### v1.3 -- Benchmark Campaign & Docs (Planned, deferred from v1.1)

### Phase 16: Benchmark Campaign Execution
**Goal**: A systematic ~25-molecule benchmark is run through the full pipeline with all 7+ calculator combinations, producing aggregated results that answer the thesis question about energy vs. dipole model dominance.
**Depends on**: Phase 23 (analysis quality overhaul must be complete so benchmark produces thesis-ready output)
**Requirements**: BENCH-01, BENCH-02, BENCH-03
**Success Criteria** (what must be TRUE):
  1. The benchmark dataset of ~25 molecules (size-scaling series CH4->C10H22 + functional group diversity series) has been run through the full pipeline and results exist in `batch_results/`.
  2. All 7+ calculator combinations (existing 4 + mace_off x espaloma, mace_off x mace_ml, mace_anicc x espaloma, mace_anicc x mace_ml, xtb x xtb) appear in the batch report with populated R-squared and RMSE values.
  3. The batch report presents aggregated ML vs. DFT accuracy stratified by calculator combination and by molecule size/functional-group class, sufficient to support thesis conclusions.
**Plans**: TBD

### Phase 17: Architecture & Development Docs Update
**Goal**: The two developer-facing documentation files accurately reflect the current `mace_gaussian/` package layout, module names, and workflows.
**Depends on**: Phase 16 (benchmark campaign may surface final architectural decisions; docs should reflect the settled state)
**Requirements**: FIX-02, FIX-03
**Success Criteria** (what must be TRUE):
  1. `docs/ARCHITECTURE.md` references only module paths that exist in the current `mace_gaussian/` package (no pre-refactor names remain).
  2. `docs/DEVELOPMENT.md` correctly describes how to add a new calculator, run tests, and reproduce the benchmark workflow using the current package layout.
  3. A new contributor could follow `docs/DEVELOPMENT.md` step-by-step to add a new energy calculator without consulting source code directly.
**Plans**: TBD

## Progress

**Execution Order:** 18 -> 19 -> 21 -> 20 -> 22 -> 23 -> 24 (Phase 20 deferred pending prof discussion on timing methodology)

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
| 13. Calculator Expansion & Acoh Fix | v1.1 | 3/3 | Complete | 2026-03-03 |
| 13.1. Calculator Acceleration | v1.1 | 2/2 | Complete | 2026-03-11 |
| 13.2. Scratch Dir | v1.1 | 2/2 | Complete | 2026-03-14 |
| 13.3. Hungarian Mode Matching | v1.1 | 2/2 | Complete | 2026-03-15 |
| 13.4. Coverage Analysis | v1.1 | 2/2 | Complete | 2026-03-16 |
| 13.5. MACE-POLAR-1 | v1.1 | 1/1 | Complete | 2026-03-23 |
| 14. Batch Runner & PubChem | v1.1 | 2/2 | Complete | 2026-03-23 |
| 15. SLURM & Batch Report | v1.1 | 2/2 | Complete | 2026-03-27 |
| 18. Spectral Quality Foundations | v1.2 | 2/2 | Complete    | 2026-03-30 |
| 19. Degenerate Mode Handling | v1.2 | 2/2 | Complete    | 2026-03-31 |
| 20. Wall-Clock Timing | v1.2 | 1/1 | Complete    | 2026-04-02 |
| 21. NIST Experimental Overlay | v1.2 | 2/2 | Complete   | 2026-04-02 |
| 22. Early SLURM Submission | v1.2 | 0/TBD | Not started | - |
| 23. Anharmonic Pipeline & Report | v1.2 | 0/TBD | Not started | - |
| 24. VPT2 Research Spike | v1.2 | 0/TBD | Not started | - |
| 16. Benchmark Campaign | v1.3 | 0/TBD | Not started | - |
| 17. Docs Update | v1.3 | 1/1 | Complete    | 2026-04-02 |
