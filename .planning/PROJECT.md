# MACE-Gaussian Interface

## What This Is

A bridge between machine learning potentials (MACE) and Gaussian 16 for computing molecular IR spectra. It injects ML-calculated dipole derivatives into Gaussian's anharmonic frequency calculations in real-time via ZMQ, enabling fast spectral predictions without full DFT cost. Built as a research tool for a master thesis in theoretical chemistry. The codebase is now a proper `mace_gaussian/` Python package with CI, tests, and documentation, ready for other researchers to clone and reproduce.

## Core Value

Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data. A researcher can clone the repo, install the package (`pip install -e .`), run `mace-gaussian run water.xyz`, and get trustworthy spectral results.

## Requirements

### Validated (v1.0)

- ✓ ZMQ-based real-time ML dipole injection into Gaussian external interface — v1.0
- ✓ Geometry optimization with MACE energy calculators (mace_mp, mace_omol) — v1.0
- ✓ Harmonic and anharmonic frequency calculations via Gaussian + ML dipoles — v1.0
- ✓ Four ML model combinations: mace_mp/mace_omol × espaloma/mace_ml — v1.0
- ✓ DFT baseline calculations (B3LYP/6-31G(d,p)) for reference spectra — v1.0
- ✓ Mode matching via eigenvector dot product overlap — v1.0
- ✓ Statistical comparison metrics (MAE, RMSE, R², regression) — v1.0
- ✓ HTML report generation with embedded plots — v1.0
- ✓ Harmonic-only analysis mode with mode overlap heatmaps — v1.0
- ✓ CLI interface (`mace-gaussian`) for running workflows and diagnostics — v1.0
- ✓ 131-test pytest suite with GPU/Gaussian markers and regression baselines — v1.0
- ✓ mace_gaussian/ package with proper imports and `mace-gaussian` entry point — v1.0
- ✓ Safe MACE model loading via mace_loader.py (no sys.modules monkey-patching) — v1.0
- ✓ gaussian/ subpackage: LINGER=0 ZMQ fix, SIGKILL timeout, context-manager server — v1.0
- ✓ GitHub Actions CI: ruff ≥0.9.0 + pytest on push — v1.0
- ✓ User documentation: quickstart, worked water example, thesis methods section — v1.0
- ✓ calculation_parameters captured in results JSON — v1.0
- ✓ CUDANotAvailableWarning emitted via warnings.warn — v1.0
- ✓ Env vars resolved at CLI startup before prerequisite validation — v1.0

### Active

(No active requirements — define in next milestone with `/gsd:new-milestone`)

### Out of Scope

- Autograd-based dipole derivatives for MACE-ML — currently numerical finite differences (18 calls/water, 126/aspirin). Would be ~10× faster. Revisit if scaling to large molecules.
- Fixing ML intensity accuracy — research problem, not a code problem
- Adding new ML model backends — focus on existing 4 combinations
- Web interface or GUI — CLI is appropriate for research use
- Support for non-Gaussian QM packages — single-target for now
- Performance optimization for very large molecules (>100 atoms) — thesis molecules are smaller
- PyPI distribution — clone-and-reproduce is the target, not pip install from index
- `compare`/`export` CLI commands — stubs, deferred to v2

## Context

- **Current state (v1.0):** Distribution-ready Python package. 7,648 LOC source + 2,147 LOC tests. 131 tests pass. CI runs on every push. `mace-gaussian run water.xyz` works end-to-end.
- **Known limitations:** `compare`/`export` CLI stubs. docs/ARCHITECTURE.md references pre-refactor module names. E2E GPU+Gaussian test requires hardware.
- **Known bug:** Acetic acid (acoh) DFT frequency parsing fails for regression plots (commit a4384c4, xfail test documents it).
- **Thesis stage:** Method validated on water/CH4; expanding to larger molecules (BH3·NH3, acoh). mace_omol/mace_ml is the best-performing combination.
- **Target audience:** Computational chemistry researchers with Gaussian 16 who want ML-augmented spectral predictions.

## Constraints

- **External dependency**: Gaussian 16 must be installed and licensed — cannot be bundled or replaced
- **GPU**: CUDA-capable NVIDIA GPU required for practical performance (CPU fallback ~10× slower)
- **Custom MACE packages**: mace_torch and mace_dipole_core are local packages, not on PyPI
- **Platform**: Linux only (HPC environment assumed)
- **Python**: 3.9+ with uv as package manager

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Keep all 4 ML model combinations | Comparison across combinations is part of thesis methodology | ✓ Good — clean comparison data |
| Focus on refactoring before new features | Code quality and reproducibility are prerequisites for distribution | ✓ Good — solid foundation for v1.1 |
| Use quality model profile for GSD agents | Important project, want thorough analysis and planning | ✓ Good — thorough plans caught wiring gaps early |
| Reproducibility over pip-installability | Thesis needs clone-and-reproduce, not PyPI distribution | ✓ Good — uv.lock + documented install procedure |
| Safe MACE loading via pickle_module remapping | Eliminates sys.modules cleanup fragility | ✓ Good — clean, tested, no side effects |
| SIGKILL (not SIGTERM) for Gaussian timeout | Gaussian ignores SIGTERM | ✓ Good — reliable timeout behavior |
| LINGER=0 before bind() in ZMQ server | Prevents close() blocking on Gaussian crash | ✓ Good — race condition eliminated |
| mace_gaussian/ package with relative imports | Enables pip install -e . and entry-point CLI | ✓ Good — works end-to-end |
| No coverage threshold in CI | Informational only; focus on correctness | ✓ Good — 57% coverage without noise |
| ruff>=0.9.0 floor (not exact pin) in CI | Removes maintenance burden vs. exact pin | ✓ Good — aligned with dev dep declaration |

---
*Last updated: 2026-02-28 after v1.0 milestone*
