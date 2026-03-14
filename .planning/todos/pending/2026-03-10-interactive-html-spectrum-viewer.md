---
created: 2026-03-10T16:05:00.000Z
title: Interactive HTML spectrum viewer with Plotly
area: analysis
files:
  - run_analysis.py
  - run_analysis_harmonic.py
---

## Problem

Static PNG plots can't be zoomed, peaks can't be hovered for assignments, and toggling between DFT/ML overlays requires regenerating plots. Crowded spectral regions are hard to inspect.

## Solution

1. Add Plotly-based HTML spectrum output alongside existing matplotlib PNGs.
2. Features: hover tooltips with peak frequency/intensity/assignment, zoom/pan, toggle DFT vs ML traces, toggle harmonic vs anharmonic.
3. Include mode matching info in tooltips (which DFT mode matches which ML mode, overlap score).
4. Generate as part of the existing analysis pipeline (e.g., `spectrum_interactive.html` next to the PNGs).
