---
created: 2026-03-10T14:14:33.440Z
title: Normalize intensity regression plots and add fit line
area: general
files:
  - run_analysis.py
  - run_analysis_harmonic.py
---

## Problem

Regression plots of IR intensities may be using raw (unnormalized) intensity values, making comparisons across molecules or methods harder to interpret. The plots should use normalized intensities. Additionally, a fit line (e.g., linear regression with displayed slope/R²/RMSE) should be added if not already present.

## Solution

1. Check `run_analysis.py` and `run_analysis_harmonic.py` to see how intensities are plotted.
2. If not already normalized: normalize intensities (e.g., to max = 1 or to sum = 1) before plotting regression.
3. Ensure a fit line is shown on the plot with R² and RMSE annotated (already a convention per CLAUDE.md).
4. Verify output looks correct across multiple molecules.
