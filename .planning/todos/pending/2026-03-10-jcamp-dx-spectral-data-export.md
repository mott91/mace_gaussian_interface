---
created: 2026-03-10T16:05:00.000Z
title: JCAMP-DX spectral data export
area: analysis
files: []
---

## Problem

Computed spectra are only available as PNG plots and internal JSON. There's no way to open results in standard spectroscopy software (SpectraGryph, MestReNova, etc.) or share data in a community-standard format.

## Solution

1. Implement JCAMP-DX (.jdx) export for computed IR spectra (both harmonic and anharmonic).
2. Include metadata (molecule name, method, basis set, calculator combination) in the JCAMP header.
3. Support both stick spectra and broadened spectra export.
4. Add as a step in the analysis pipeline or as a `mace-gaussian export` CLI command.
