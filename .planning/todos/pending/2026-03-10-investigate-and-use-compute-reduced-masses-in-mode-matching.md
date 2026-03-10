---
created: 2026-03-10T15:17:30.448Z
title: Investigate and use compute_reduced_masses in mode matching
area: general
files:
  - mace_gaussian/analysis/mode_matching.py:75-105
  - mace_gaussian/gaussian/fchk.py
---

## Problem

`compute_reduced_masses()` is defined in `mode_matching.py` but never called anywhere in the codebase. It computes μ_i = Σ m_atom * |L_{i,atom}|² for each mode.

The question is: should mode overlap be computed in mass-weighted coordinates or raw Cartesian? Gaussian's .fchk stores "Vib-Modes" as mass-weighted Cartesian displacement vectors (L matrix, already in √amu·Å units). If both DFT and ML modes come from Gaussian .fchk files in the same convention, raw dot product is already correct (they're both mass-weighted the same way). But if one side comes from a different source (e.g., direct MACE output), mass-weighting may need to be applied explicitly.

## Solution

1. Confirm what convention Gaussian uses in the .fchk "Vib-Modes" section — are vectors mass-weighted or not? Check against the Gaussian documentation or a known test case.
2. If mass-weighting is already baked into .fchk mode vectors: `compute_reduced_masses` is redundant for DFT-vs-ML comparison and can be documented as such (or removed).
3. If mass-weighting is NOT applied: apply it in `compute_mode_overlap()` before the dot product, using masses extracted from the .fchk.
4. Either way, document the conclusion clearly in the code.
