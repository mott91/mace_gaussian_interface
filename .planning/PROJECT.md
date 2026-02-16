# MACE-Gaussian Interface

## What This Is

A bridge between machine learning potentials (MACE) and Gaussian 16 for computing molecular IR spectra. It injects ML-calculated dipole derivatives into Gaussian's anharmonic frequency calculations in real-time via ZMQ, enabling fast spectral predictions without full DFT cost. Built as a research tool for a master thesis in theoretical chemistry, with the goal of being reproducible and eventually distributable to other researchers.

## Core Value

Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data. If nothing else works, a researcher must be able to clone the repo, run the pipeline on a molecule, and get trustworthy harmonic frequency results.

## Requirements

### Validated

- ZMQ-based real-time ML dipole injection into Gaussian external interface — existing
- Geometry optimization with MACE energy calculators (mace_mp, mace_omol) — existing
- Harmonic and anharmonic frequency calculations via Gaussian + ML dipoles — existing
- Four ML model combinations supported: mace_mp/mace_omol x espaloma/mace_ml — existing
- DFT baseline calculations (B3LYP/6-31G(d,p)) for reference spectra — existing
- Mode matching via eigenvector dot product overlap — existing
- Statistical comparison metrics (MAE, RMSE, R², regression) — existing
- HTML report generation with embedded plots — existing
- Harmonic-only analysis mode with mode overlap heatmaps — existing
- CLI interface for running workflows and diagnostics — existing

### Active

- [ ] Replace module monkey-patching with safe factory/isolation pattern
- [ ] Add subprocess timeout management for Gaussian processes
- [ ] Fix ZMQ socket cleanup race conditions
- [ ] Add unit tests for Gaussian parser (frequencies, overtones, combination bands)
- [ ] Add unit tests for checkpoint file handling and mode matching
- [ ] Validate environment variables and input files at startup
- [ ] Improve error messages (explicit errors instead of silent empty returns)
- [ ] Track optimization step counts properly
- [ ] Implement atomic file operations for checkpoint handling
- [ ] Log CUDA device selection (warn on CPU fallback)
- [ ] Clean up code to meet distribution-ready standards (consistent patterns, documentation)
- [ ] Ensure end-to-end reproducibility (someone can clone and reproduce results)

### Out of Scope

- Fixing ML intensity accuracy — research problem, not a code problem
- Adding new ML model backends — focus on cleaning existing 4 combinations first
- Web interface or GUI — CLI tool is appropriate for research use
- Supporting quantum chemistry packages other than Gaussian 16 — single-target for now
- Performance optimization for very large molecules (>100 atoms) — thesis molecules are smaller
- Automatic retry/restart logic — nice to have but not essential for thesis scope

## Context

- **Author background**: Not a programmer; codebase built with heavy AI assistance. Code quality and conventions may be inconsistent.
- **Current state**: Harmonic frequencies work well across molecule sizes. Anharmonic works but intensities are less reliable. mace_omol/mace_ml is the best-performing combination.
- **Known bug**: Acetic acid (acoh) DFT frequency parsing fails for regression plots (commit a4384c4).
- **Critical fragility**: Module monkey-patching for MACE is the biggest reliability risk. ZMQ socket cleanup and subprocess management are secondary concerns.
- **Thesis stage**: Early — method still being developed and validated. Plenty of time to get foundations right.
- **Target audience**: Eventually other computational chemistry researchers who have Gaussian 16 and want to try ML-augmented spectral predictions.

## Constraints

- **External dependency**: Gaussian 16 must be installed and licensed — cannot be bundled or replaced
- **GPU**: CUDA-capable NVIDIA GPU required for practical performance (CPU fallback exists but ~10x slower)
- **Custom MACE packages**: mace_torch and mace_dipole_core are local packages, not pip-installable from PyPI
- **Platform**: Linux only (HPC environment assumed)
- **Python**: 3.9+ with uv as package manager

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Keep all 4 ML model combinations | Comparison across combinations is part of thesis methodology | — Pending |
| Focus on refactoring before new features | Code quality and reproducibility are prerequisites for distribution | — Pending |
| Use quality model profile for GSD agents | Important project, want thorough analysis and planning | — Pending |
| Reproducibility over pip-installability | Thesis needs clone-and-reproduce, not package distribution yet | — Pending |

---
*Last updated: 2026-02-16 after initialization*
