---
created: 2026-03-10T16:05:00.000Z
title: Conformer-aware IR spectra via Boltzmann weighting
area: workflow
files:
  - mace_gaussian/workflow.py
---

## Problem

The pipeline currently optimizes to a single geometry and computes frequencies for that one conformer. For flexible molecules, different conformers can have significantly different IR spectra. The single-conformer result may not match experiment well.

## Solution

1. Add an optional conformer search step before frequency calculation (using RDKit ETKDG or CREST).
2. Optimize each conformer with MACE, compute relative energies.
3. Run frequency calculations on the low-energy conformers (within some kT window).
4. Boltzmann-weight the individual spectra to produce a population-averaged spectrum.
5. Report both individual conformer spectra and the weighted composite.
