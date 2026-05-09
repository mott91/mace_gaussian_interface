---
created: 2026-04-13T00:00:00.000Z
title: Add Pareto cost-accuracy plot back to HTML report
area: analysis
files:
  - mace_gaussian/analysis/html_report_generator.py
  - mace_gaussian/analysis/plotly_builders.py
---

## Problem

The Pareto (speedup vs MAE) scatter plot was removed from the HTML report because of persistent Plotly layout/sizing issues — the plot was displaced from its container, labels overlapped other sections, and the white box didn't match the plot dimensions.

## Solution

1. Fix the Plotly rendering issue (likely `responsive: True` fighting with fixed dimensions).
2. Re-integrate `build_pareto_figure` into the report. Filter to mace_ml dipole methods only (no espaloma).
3. Test that it renders correctly in the report container without overflowing or displacing.

Note: `build_pareto_figure` still exists in `plotly_builders.py` — just needs to be wired back in.
