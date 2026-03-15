# Requirements: MACE-Gaussian Interface

**Defined:** 2026-03-03
**Milestone:** v1.1 — Batch Benchmarking & Calculator Expansion
**Core Value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.

## v1.1 Requirements

### Calculator Expansion

- [x] **CALC-01**: User can specify `mace_off` as `--energy-calculators` choice in CLI (organic-optimized PES, H/C/N/O/F/P/S/Cl/Br/I)
- [x] **CALC-02**: User can specify `mace_anicc` as `--energy-calculators` choice in CLI (CC-quality PES, CHNO only)
- [ ] **CALC-03**: User can specify `xtb` as `--energy-calculators` choice in CLI (GFN2-xTB semiempirical)
- [ ] **CALC-04**: xTB dipole calculator produces correct dipole derivatives (unit bug verified and fixed; `xtb` usable as `--dipole-calculators` choice)

### Scratch Directory Isolation

- [x] **SCRATCH-01**: All intermediate Gaussian files (.gjf, .log, .chk, .ipc_file) are created inside `.scratch/` subdirectory, not the project root — both ML frequency and DFT baseline runs
- [x] **SCRATCH-02**: Scratch directory is auto-deleted on success (after moving files to results) and on failure (unless `--keep-scratch` flag or `MACE_KEEP_SCRATCH=1` env var is set)
- [x] **SCRATCH-03**: `mace-gaussian run` auto-cleans stale scratch dirs (>24h) on startup; `mace-gaussian diagnose` reports stale scratch dirs

### Mode Matching

- [x] **MATCH-01**: `match_modes()` uses Hungarian algorithm (`linear_sum_assignment`) for globally optimal bijective 1-to-1 ML-to-DFT mode pairing — no two ML modes can claim the same DFT mode
- [x] **MATCH-02**: Unmatched modes (when mode counts differ) appear as `(None, 0.0)` in match results; low-overlap pairs are logged but kept (not dropped)
- [x] **MATCH-03**: Mode overlap heatmaps annotate Hungarian-matched cells with borders (solid for confident >= 0.5, dashed for uncertain < 0.5) and include a legend
- [x] **MATCH-04**: Regression plots differentiate confident matches (filled markers) from low-overlap matches (open circles)

### Frequency Range Coverage Analysis

- [ ] **COV-01**: Per-region error metrics (RMSE, MAE, mean % error, mode count) computed for 7 chemically meaningful frequency regions (400-700, 700-1000, 1000-1300, 1300-1500, 1500-1800, 1800-2800, 2800-4000 cm-1) from harmonic comparison data
- [ ] **COV-02**: Cross-molecule heatmap (molecules x regions matrix, colored by RMSE) generated per calculator combination, with gray "--" annotation for empty cells
- [ ] **COV-03**: Cross-molecule bar chart (mean RMSE by region with std dev error bars and mode count annotations) generated per calculator combination
- [ ] **COV-04**: Standalone self-contained HTML report with embedded heatmap and bar chart plots, summary tables per calculator combination
- [ ] **COV-05**: Entry-point script `run_coverage_analysis.py` accepting molecule names as CLI arguments, producing PNG plots, JSON metrics, and HTML report in `coverage_analysis/` directory

### Batch Workflow

- [ ] **BATCH-01**: User can run `mace-gaussian fetch <molecule-name>` to download a 3D XYZ structure from PubChem
- [ ] **BATCH-02**: User can run `mace-gaussian batch molecules.txt` to process multiple molecules sequentially through the full pipeline
- [ ] **BATCH-03**: Batch run produces a per-molecule status manifest (`batch_manifest.json`) that survives interruption — restarting skips already-complete molecules
- [ ] **BATCH-04**: User can run `mace-gaussian batch molecules.txt --skip-dft-baseline` to run ML calculations only
- [ ] **BATCH-05**: Batch run produces a multi-molecule HTML report with aggregated R² and RMSE per calculator combination across all molecules

### HPC / SLURM

- [ ] **HPC-01**: User can run `mace-gaussian batch molecules.txt --dft-on-cluster <host>` to submit DFT baseline jobs to a SLURM cluster via SSH, poll for completion, and retrieve results automatically
- [ ] **HPC-02**: SLURM submission includes `formchk` in the job script so `.fchk` is produced on the cluster without requiring local conversion

### Benchmark Campaign

- [ ] **BENCH-01**: Benchmark dataset of ~25 molecules is run through the full pipeline (size-scaling series: CH4 -> C10H22; functional group diversity series: water, ammonia, formaldehyde, formic acid, methanol, acetaldehyde, acetic acid, dimethylamine, ethanol, acetone, glycine, aspirin, cocaine + 2 drug-like)
- [ ] **BENCH-02**: Benchmark results include all 7+ calculator combinations (existing 4 + mace_off x espaloma, mace_off x mace_ml, mace_anicc x espaloma, mace_anicc x mace_ml, xtb x xtb)
- [ ] **BENCH-03**: Aggregated batch report summarises ML vs DFT accuracy by calculator combination and molecule size/class

### Bug Fixes & Docs

- [x] **FIX-01**: Acetic acid (acoh) frequency parser bug fixed — regression plot frequency matching works correctly and xfail test is promoted to passing
- [ ] **FIX-02**: `docs/ARCHITECTURE.md` updated to reflect current `mace_gaussian/` package layout and module names
- [ ] **FIX-03**: `docs/DEVELOPMENT.md` updated to reflect current package layout (adding calculators, running tests, etc.)

## v2 Requirements

### Open-Source QM Backend

- **OQM-01**: Investigate Psi4 Python API as candidate for driving VPT2 with ML-supplied dipole derivatives (feasibility prototype)
- **OQM-02**: Pure Python VPT2 solver integration to eliminate Gaussian dependency entirely

### Advanced Calculators

- **CALC-05**: TorchANI/ANI-2x wired as energy calculator option
- **CALC-06**: Autograd-based dipole derivatives for MACE-ML (~10x speedup over finite differences)

### Experimental Validation

- **EXP-01**: NIST WebBook experimental IR spectrum overlay in analysis report (JCAMP-DX parser + HTTP downloader)

## Out of Scope

| Feature | Reason |
|---------|--------|
| ORCA/VPT2 integration | No external hook in ORCA — technically blocked; requires ORCA developer cooperation |
| Periodic systems (surfaces, crystals) | Different project — QM engine, dipole model, and analysis pipeline all need replacing simultaneously |
| Born effective charge ML models | Research frontier, not packaged; 6+ month research contribution |
| PyPI distribution | Clone-and-reproduce is the target for thesis |
| Web interface or GUI | CLI appropriate for research use |
| `compare`/`export` CLI commands | Stubs, deferred to v2 |
| AIMNet2 integration | Installation friction high, marginal benefit over mace_off |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| CALC-01 | Phase 13 | Complete |
| CALC-02 | Phase 13 | Complete |
| CALC-03 | Phase 13 | Pending |
| CALC-04 | Phase 13 | Pending |
| SCRATCH-01 | Phase 13.2 | Complete |
| SCRATCH-02 | Phase 13.2 | Complete |
| SCRATCH-03 | Phase 13.2 | Complete |
| MATCH-01 | Phase 13.3 | Complete |
| MATCH-02 | Phase 13.3 | Complete |
| MATCH-03 | Phase 13.3 | Complete |
| MATCH-04 | Phase 13.3 | Complete |
| COV-01 | Phase 13.4 | Pending |
| COV-02 | Phase 13.4 | Pending |
| COV-03 | Phase 13.4 | Pending |
| COV-04 | Phase 13.4 | Pending |
| COV-05 | Phase 13.4 | Pending |
| BATCH-01 | Phase 14 | Pending |
| BATCH-02 | Phase 14 | Pending |
| BATCH-03 | Phase 14 | Pending |
| BATCH-04 | Phase 14 | Pending |
| BATCH-05 | Phase 15 | Pending |
| HPC-01 | Phase 15 | Pending |
| HPC-02 | Phase 15 | Pending |
| BENCH-01 | Phase 16 | Pending |
| BENCH-02 | Phase 16 | Pending |
| BENCH-03 | Phase 16 | Pending |
| FIX-01 | Phase 13 | Complete |
| FIX-02 | Phase 17 | Pending |
| FIX-03 | Phase 17 | Pending |

**Coverage:**
- v1.1 requirements: 29 total
- Mapped to phases: 29
- Unmapped: 0

---
*Requirements defined: 2026-03-03*
*Last updated: 2026-03-15 after Phase 13.4 planning (COV-01 through COV-05 added)*
