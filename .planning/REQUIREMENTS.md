# Requirements: MACE-Gaussian Interface

**Defined:** 2026-03-28
**Milestone:** v1.2 -- Analysis Quality Overhaul
**Core Value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.

## v1.1 Requirements

### Calculator Expansion

- [x] **CALC-01**: User can specify `mace_off` as `--energy-calculators` choice in CLI (organic-optimized PES, H/C/N/O/F/P/S/Cl/Br/I)
- [x] **CALC-02**: User can specify `mace_anicc` as `--energy-calculators` choice in CLI (CC-quality PES, CHNO only)
- [ ] **CALC-03**: User can specify `xtb` as `--energy-calculators` choice in CLI (GFN2-xTB semiempirical)
- [ ] **CALC-04**: xTB dipole calculator produces correct dipole derivatives (unit bug verified and fixed; `xtb` usable as `--dipole-calculators` choice)

### MACE-POLAR-1 Integration

- [x] **POLAR-01**: User can specify `mace_polar` as `--energy-calculators` choice in CLI (MACE-POLAR-1 with long-range electrostatics, 83 elements, trained on OMol25)
- [x] **POLAR-02**: User can specify `mace_polar` as `--optimization-calculator` choice in CLI
- [x] **POLAR-03**: `calculator("mace_polar")` in workflow.py calls `mace_polar(model="polar-1-l", device="cuda", default_dtype="float64")` with no dispersion kwarg, and DFT baseline config exists

### Scratch Directory Isolation

- [x] **SCRATCH-01**: All intermediate Gaussian files (.gjf, .log, .chk, .ipc_file) are created inside `.scratch/` subdirectory, not the project root -- both ML frequency and DFT baseline runs
- [x] **SCRATCH-02**: Scratch directory is auto-deleted on success (after moving files to results) and on failure (unless `--keep-scratch` flag or `MACE_KEEP_SCRATCH=1` env var is set)
- [x] **SCRATCH-03**: `mace-gaussian run` auto-cleans stale scratch dirs (>24h) on startup; `mace-gaussian diagnose` reports stale scratch dirs

### Mode Matching

- [x] **MATCH-01**: `match_modes()` uses Hungarian algorithm (`linear_sum_assignment`) for globally optimal bijective 1-to-1 ML-to-DFT mode pairing -- no two ML modes can claim the same DFT mode
- [x] **MATCH-02**: Unmatched modes (when mode counts differ) appear as `(None, 0.0)` in match results; low-overlap pairs are logged but kept (not dropped)
- [x] **MATCH-03**: Mode overlap heatmaps annotate Hungarian-matched cells with borders (solid for confident >= 0.5, dashed for uncertain < 0.5) and include a legend
- [x] **MATCH-04**: Regression plots differentiate confident matches (filled markers) from low-overlap matches (open circles)

### Frequency Range Coverage Analysis

- [x] **COV-01**: Per-region error metrics (RMSE, MAE, mean % error, mode count) computed for 7 chemically meaningful frequency regions (400-700, 700-1000, 1000-1300, 1300-1500, 1500-1800, 1800-2800, 2800-4000 cm-1) from harmonic comparison data
- [x] **COV-02**: Cross-molecule heatmap (molecules x regions matrix, colored by RMSE) generated per calculator combination, with gray "--" annotation for empty cells
- [x] **COV-03**: Cross-molecule bar chart (mean RMSE by region with std dev error bars and mode count annotations) generated per calculator combination
- [x] **COV-04**: Standalone self-contained HTML report with embedded heatmap and bar chart plots, summary tables per calculator combination
- [x] **COV-05**: Entry-point script `run_coverage_analysis.py` accepting molecule names as CLI arguments, producing PNG plots, JSON metrics, and HTML report in `coverage_analysis/` directory

### Batch Workflow

- [x] **BATCH-01**: User can run `mace-gaussian fetch <molecule-name>` to download a 3D XYZ structure from PubChem
- [x] **BATCH-02**: User can run `mace-gaussian batch molecules.txt` to process multiple molecules sequentially through the full pipeline
- [x] **BATCH-03**: Batch run produces a per-molecule status manifest (`batch_manifest.json`) that survives interruption -- restarting skips already-complete molecules
- [x] **BATCH-04**: User can run `mace-gaussian batch molecules.txt --skip-dft-baseline` to run ML calculations only
- [x] **BATCH-05**: Batch run produces a multi-molecule HTML report with aggregated R-squared and RMSE per calculator combination across all molecules

### HPC / SLURM

- [x] **HPC-01**: User can run `mace-gaussian batch molecules.txt --dft-on-cluster <host>` to submit DFT baseline jobs to a SLURM cluster via SSH, poll for completion, and retrieve results automatically
- [x] **HPC-02**: SLURM submission includes `formchk` in the job script so `.fchk` is produced on the cluster without requiring local conversion

### Benchmark Campaign

- [ ] **BENCH-01**: Benchmark dataset of ~25 molecules is run through the full pipeline (size-scaling series: CH4 -> C10H22; functional group diversity series: water, ammonia, formaldehyde, formic acid, methanol, acetaldehyde, acetic acid, dimethylamine, ethanol, acetone, glycine, aspirin, cocaine + 2 drug-like)
- [ ] **BENCH-02**: Benchmark results include all 7+ calculator combinations (existing 4 + mace_off x espaloma, mace_off x mace_ml, mace_anicc x espaloma, mace_anicc x mace_ml, xtb x xtb)
- [ ] **BENCH-03**: Aggregated batch report summarises ML vs DFT accuracy by calculator combination and molecule size/class

### Bug Fixes & Docs

- [x] **FIX-01**: Acetic acid (acoh) frequency parser bug fixed -- regression plot frequency matching works correctly and xfail test is promoted to passing
- [ ] **FIX-02**: `docs/ARCHITECTURE.md` updated to reflect current `mace_gaussian/` package layout and module names
- [ ] **FIX-03**: `docs/DEVELOPMENT.md` updated to reflect current package layout (adding calculators, running tests, etc.)

## v1.2 Requirements

### Spectral Quality

- [x] **SPEC-01**: Simulated IR spectra use Lorentzian line shapes instead of Gaussian KDE, with configurable FWHM (default 10 cm-1)
- [ ] **SPEC-02**: Modes with IR intensity below 0.1 km/mol are filtered from intensity regression metrics (but retained in frequency metrics)
- [ ] **SPEC-03**: Analysis reports display both stick spectrum and broadened Lorentzian spectrum

### Mode Matching (v1.2)

- [ ] **MODE-05**: Degenerate modes (within 5 cm-1 threshold) are detected and grouped; subspace overlap (trace of M^T M / k) is used for matching quality
- [ ] **MODE-06**: Mode matching statistics correctly handle degenerate groups without double-counting

### Timing

- [ ] **TIME-01**: Wall-clock timing is recorded per molecule per calculator combination (geometry opt, frequency calc, total)
- [ ] **TIME-02**: Batch report includes ML vs DFT timing comparison table and speedup factors

### Experimental Comparison

- [ ] **NIST-01**: User can fetch experimental IR spectrum from NIST WebBook by molecule name (cached locally)
- [ ] **NIST-02**: Analysis report overlays experimental spectrum on computed spectra plot when available
- [ ] **NIST-03**: Quantitative peak position comparison (experimental vs computed) with error metrics

### Pipeline Improvements

- [ ] **PIPE-01**: Batch runner submits DFT jobs to SLURM immediately after each molecule's geometry optimization, not after all ML calcs complete
- [ ] **PIPE-02**: Mace_polar dipole calculator failure mode investigated and either fixed or documented with skip logic

### Anharmonic Analysis

- [ ] **ANAL-01**: Anharmonic analysis pipeline produces thesis-quality HTML report integrating Lorentzian spectra, experimental overlay, timing, and mode matching
- [ ] **ANAL-02**: Report includes per-molecule summary cards with key metrics (R-squared, RMSE, speedup, experimental agreement)

### VPT2 Research

- [ ] **VPT2-01**: Time-boxed (2hr) spike evaluates Psience library for VPT2 anharmonic calculations using MACE force constants on water molecule

## v2 Requirements

### Open-Source QM Backend

- **OQM-01**: Investigate Psi4 Python API as candidate for driving VPT2 with ML-supplied dipole derivatives (feasibility prototype)
- **OQM-02**: Pure Python VPT2 solver integration to eliminate Gaussian dependency entirely

### Advanced Calculators

- **CALC-05**: TorchANI/ANI-2x wired as energy calculator option
- **CALC-06**: Autograd-based dipole derivatives for MACE-ML (~10x speedup over finite differences)


## Out of Scope

| Feature | Reason |
|---------|--------|
| ORCA/VPT2 integration | No external hook in ORCA -- technically blocked; requires ORCA developer cooperation |
| Periodic systems (surfaces, crystals) | Different project -- QM engine, dipole model, and analysis pipeline all need replacing simultaneously |
| Born effective charge ML models | Research frontier, not packaged; 6+ month research contribution |
| PyPI distribution | Clone-and-reproduce is the target for thesis |
| Web interface or GUI | CLI appropriate for research use |
| `compare`/`export` CLI commands | Stubs, deferred to v2 |
| AIMNet2 integration | Installation friction high, marginal benefit over mace_off |

## Traceability

### v1.1 Requirements

| Requirement | Phase | Status |
|-------------|-------|--------|
| CALC-01 | Phase 13 | Complete |
| CALC-02 | Phase 13 | Complete |
| CALC-03 | Phase 13 | Pending |
| CALC-04 | Phase 13 | Pending |
| POLAR-01 | Phase 13.5 | Complete |
| POLAR-02 | Phase 13.5 | Complete |
| POLAR-03 | Phase 13.5 | Complete |
| SCRATCH-01 | Phase 13.2 | Complete |
| SCRATCH-02 | Phase 13.2 | Complete |
| SCRATCH-03 | Phase 13.2 | Complete |
| MATCH-01 | Phase 13.3 | Complete |
| MATCH-02 | Phase 13.3 | Complete |
| MATCH-03 | Phase 13.3 | Complete |
| MATCH-04 | Phase 13.3 | Complete |
| COV-01 | Phase 13.4 | Complete |
| COV-02 | Phase 13.4 | Complete |
| COV-03 | Phase 13.4 | Complete |
| COV-04 | Phase 13.4 | Complete |
| COV-05 | Phase 13.4 | Complete |
| BATCH-01 | Phase 14 | Complete |
| BATCH-02 | Phase 14 | Complete |
| BATCH-03 | Phase 14 | Complete |
| BATCH-04 | Phase 14 | Complete |
| BATCH-05 | Phase 15 | Complete |
| HPC-01 | Phase 15 | Complete |
| HPC-02 | Phase 15 | Complete |
| BENCH-01 | Phase 16 | Pending |
| BENCH-02 | Phase 16 | Pending |
| BENCH-03 | Phase 16 | Pending |
| FIX-01 | Phase 13 | Complete |
| FIX-02 | Phase 17 | Pending |
| FIX-03 | Phase 17 | Pending |

### v1.2 Requirements

| Requirement | Phase | Status |
|-------------|-------|--------|
| SPEC-01 | Phase 18 | Complete |
| SPEC-02 | Phase 18 | Pending |
| SPEC-03 | Phase 18 | Pending |
| PIPE-02 | Phase 18 | Pending |
| MODE-05 | Phase 19 | Pending |
| MODE-06 | Phase 19 | Pending |
| TIME-01 | Phase 20 | Pending |
| TIME-02 | Phase 20 | Pending |
| NIST-01 | Phase 21 | Pending |
| NIST-02 | Phase 21 | Pending |
| NIST-03 | Phase 21 | Pending |
| PIPE-01 | Phase 22 | Pending |
| ANAL-01 | Phase 23 | Pending |
| ANAL-02 | Phase 23 | Pending |
| VPT2-01 | Phase 24 | Pending |

**Coverage:**
- v1.1 requirements: 32 total, 32 mapped, 0 unmapped
- v1.2 requirements: 15 total, 15 mapped, 0 unmapped

---
*Requirements defined: 2026-03-03*
*Last updated: 2026-03-28 -- v1.2 roadmap created with phase mappings*
