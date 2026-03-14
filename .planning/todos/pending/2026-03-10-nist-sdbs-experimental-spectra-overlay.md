---
created: 2026-03-10T16:05:00.000Z
title: NIST/SDBS experimental spectra overlay
area: analysis
files:
  - run_analysis.py
  - run_analysis_harmonic.py
---

## Problem

Currently no way to visually compare computed spectra against real experimental data without manually downloading and overlaying reference spectra. This makes validation tedious.

## Solution

1. Build a fetcher that pulls experimental IR spectra from public databases (NIST WebBook, SDBS).
2. Parse the downloaded data (JCAMP-DX or digitized points) into frequency/intensity arrays.
3. Overlay experimental spectrum on the computed DFT and ML spectra plots.
4. Handle cases where experimental data isn't available (graceful skip with a warning).
5. Cache downloaded spectra locally to avoid repeated requests.
