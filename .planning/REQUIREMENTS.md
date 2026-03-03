# Requirements: MACE-Gaussian Interface

**Defined:** 2026-03-03
**Milestone:** v1.1 — Batch Benchmarking & Calculator Expansion
**Core Value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.

## v1.1 Requirements

### Calculator Expansion

- [ ] **CALC-01**: User can specify `mace_off` as `--energy-calculators` choice in CLI (organic-optimized PES, H/C/N/O/F/P/S/Cl/Br/I)
- [ ] **CALC-02**: User can specify `mace_anicc` as `--energy-calculators` choice in CLI (CC-quality PES, CHNO only)
- [ ] **CALC-03**: User can specify `xtb` as `--energy-calculators` choice in CLI (GFN2-xTB semiempirical)
- [ ] **CALC-04**: xTB dipole calculator produces correct dipole derivatives (unit bug verified and fixed; `xtb` usable as `--dipole-calculators` choice)

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

- [ ] **BENCH-01**: Benchmark dataset of ~25 molecules is run through the full pipeline (size-scaling series: CH₄ → C₁₀H₂₂; functional group diversity series: water, ammonia, formaldehyde, formic acid, methanol, acetaldehyde, acetic acid, dimethylamine, ethanol, acetone, glycine, aspirin, cocaine + 2 drug-like)
- [ ] **BENCH-02**: Benchmark results include all 7+ calculator combinations (existing 4 + mace_off×espaloma, mace_off×mace_ml, mace_anicc×espaloma, mace_anicc×mace_ml, xtb×xtb)
- [ ] **BENCH-03**: Aggregated batch report summarises ML vs DFT accuracy by calculator combination and molecule size/class

### Bug Fixes & Docs

- [ ] **FIX-01**: Acetic acid (acoh) frequency parser bug fixed — regression plot frequency matching works correctly and xfail test is promoted to passing
- [ ] **FIX-02**: `docs/ARCHITECTURE.md` updated to reflect current `mace_gaussian/` package layout and module names
- [ ] **FIX-03**: `docs/DEVELOPMENT.md` updated to reflect current package layout (adding calculators, running tests, etc.)

## v2 Requirements

### Open-Source QM Backend

- **OQM-01**: Investigate Psi4 Python API as candidate for driving VPT2 with ML-supplied dipole derivatives (feasibility prototype)
- **OQM-02**: Pure Python VPT2 solver integration to eliminate Gaussian dependency entirely

### Advanced Calculators

- **CALC-05**: TorchANI/ANI-2x wired as energy calculator option
- **CALC-06**: Autograd-based dipole derivatives for MACE-ML (~10× speedup over finite differences)

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
| CALC-01 | Phase 13 | Pending |
| CALC-02 | Phase 13 | Pending |
| CALC-03 | Phase 13 | Pending |
| CALC-04 | Phase 13 | Pending |
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
| FIX-01 | Phase 13 | Pending |
| FIX-02 | Phase 17 | Pending |
| FIX-03 | Phase 17 | Pending |

**Coverage:**
- v1.1 requirements: 17 total
- Mapped to phases: 17
- Unmapped: 0 ✓

---
*Requirements defined: 2026-03-03*
*Last updated: 2026-03-03 after roadmap creation (ROADMAP.md Phases 13–17 written)*
