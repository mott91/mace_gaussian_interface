# Milestones

## v1.0 — Refactoring & Distribution (Shipped: 2026-02-28)

**Phases completed:** 12 phases, 32 plans
**Requirements satisfied:** 38/38
**Test suite:** 131 passing, 1 xfailed (known acoh parser edge case)
**Package:** 7,648 LOC (mace_gaussian/) + 2,147 LOC (tests/)
**Timeline:** Sep 2025 → Feb 2026 (~108 min GSD execution)

**Key accomplishments:**
1. Built 131-test pytest suite with fixtures, GPU/Gaussian markers, and regression baselines for water and CH4 molecules
2. Reorganized monolithic scripts into proper `mace_gaussian/` package with 29 modules and `mace-gaussian` CLI entry point
3. Replaced MACE module monkey-patching (`sys.modules` cleanup) with safe `mace_loader.py` using pickle_module remapping
4. Extracted `gaussian/` subpackage: LINGER=0 ZMQ race-condition fix, SIGKILL timeout, and context-manager ZMQ server
5. Wired GitHub Actions CI: ruff ≥0.9.0 lint + pytest (unit tests only) on every push
6. Closed all integration wiring gaps and completed user documentation: quickstart, worked water example, thesis methods section

**Tech debt carried forward:**
- `compare`/`export` CLI commands are honest stubs (deferred to v2)
- `docs/ARCHITECTURE.md` and `docs/DEVELOPMENT.md` reference pre-refactor module names
- Phase 5 E2E GPU+Gaussian pipeline confirmation requires hardware access

**Archive:** `.planning/milestones/v1.0-ROADMAP.md`, `.planning/milestones/v1.0-REQUIREMENTS.md`, `.planning/milestones/v1.0-MILESTONE-AUDIT.md`

---

## v1.1 — Batch Benchmarking & Calculator Expansion (Shipped: 2026-03-28)

**Phases completed:** 8 phases (13–15), 16 plans
**Requirements satisfied:** Core calculator expansion, batch runner, SLURM integration, batch report
**Timeline:** Mar 2026 (~4 weeks)

**Key accomplishments:**
1. Added mace_off, mace_anicc, mace_polar as energy calculator options (7+ combos total)
2. Autograd dipole derivatives + polarizability passthrough in MACEDipoleCalculator
3. Scratch directory isolation for all Gaussian intermediate files
4. Hungarian optimal mode matching (bijective 1-to-1 ML-to-DFT pairing)
5. Frequency range coverage analysis with per-region error metrics
6. Batch runner with manifest-driven per-calculator restart and PubChem fetcher
7. SLURM DFT offloading via SSH/SCP with sacct polling

**Deferred to v1.3:**
- Phase 16: Benchmark campaign (~25 molecules × 7+ calculator combos)
- Phase 17: Architecture & development docs update

**Archive:** `.planning/milestones/v1.1-ROADMAP.md`, `.planning/milestones/v1.1-REQUIREMENTS.md`

---
