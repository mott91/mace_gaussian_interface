---
created: 2026-03-15T00:00:00.000Z
title: Frequency range coverage analysis for ML training set gaps
area: analysis
files:
  - mace_gaussian/analysis/analyze_spectra.py
---

## Problem

When ML models predict IR spectra, certain frequency ranges may be poorly represented in the training data, leading to systematically worse predictions in those regions. Currently there's no diagnostic to identify which frequency ranges or specific vibrational modes are underrepresented. If the same gaps appear across multiple different molecular species, that points to a systematic training set deficiency rather than molecule-specific issues.

## Solution

Build a frequency range coverage analyzer that:

1. Splits the IR spectrum into meaningful regions (e.g., fingerprint 400–1500 cm⁻¹, double bond 1500–2000 cm⁻¹, triple bond/X-H 2000–3500 cm⁻¹, O-H/N-H stretch 3500–4000 cm⁻¹) or uses finer bins.
2. For each region, computes per-region error metrics (MAE, R², RMSE) comparing ML vs DFT.
3. Flags regions where ML accuracy drops significantly below the molecule's average.
4. Cross-references across multiple molecules to identify systematic gaps: "Region X has poor accuracy across 8/10 molecules tested — likely undertrained."
5. Generates a heatmap or bar chart showing accuracy by frequency region × molecule.
6. Optionally suggests which functional groups / vibrational modes correspond to the weak regions.

This is valuable for the thesis: instead of reporting a single aggregate R², we can say "ML accuracy degrades in the X–Y cm⁻¹ range, consistent across molecule classes, suggesting training set gaps."
