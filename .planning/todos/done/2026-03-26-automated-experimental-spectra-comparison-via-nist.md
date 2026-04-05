---
created: 2026-03-26T18:35:00.000Z
title: Automated experimental spectra comparison via NIST API
area: analysis
files:
  - run_analysis.py
  - run_analysis_harmonic.py
---

## Problem

Currently no systematic comparison of computed spectra against experimental data. Comparing ML vs DFT only shows how well ML reproduces DFT, not how well either matches reality. For the thesis, experimental validation is critical. Doing this manually for each molecule is tedious.

## Solution

1. Use the NIST WebBook API (or web scraping if no public API) to automatically fetch gas-phase IR spectra for benchmark molecules.
2. Parse JCAMP-DX format (NIST's standard export format) into frequency/transmittance arrays.
3. Overlay experimental spectrum on computed DFT and ML spectra in analysis plots.
4. Compute quantitative metrics: peak position error (cm⁻¹) for matched peaks, overall spectral similarity score.
5. Cache fetched spectra locally to avoid repeated requests.
6. Include in batch report: per-molecule experimental comparison, aggregate accuracy stats.
7. Prioritize gas-phase data (closest to vacuum calculations). Flag molecules where only condensed-phase data exists and note expected systematic shifts.
8. Consider also SDBS (Spectral Database for Organic Compounds) as secondary source.
