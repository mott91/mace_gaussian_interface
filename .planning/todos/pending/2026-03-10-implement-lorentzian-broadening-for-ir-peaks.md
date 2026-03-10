---
created: 2026-03-10T14:14:33.443Z
title: Implement Lorentzian broadening for IR peaks
area: general
files:
  - run_analysis.py
  - run_analysis_harmonic.py
---

## Problem

Simulated IR spectra currently show stick spectra (discrete peaks). Real experimental spectra have Lorentzian (or Voigt) peak shapes due to lifetime broadening. Adding Lorentzian broadening would make simulated spectra directly visually comparable to experiment and is standard in computational IR spectroscopy.

## Solution

1. Implement a Lorentzian broadening function: L(ν) = (γ/π) / ((ν - ν₀)² + γ²), where γ is the half-width at half-maximum (HWHM).
2. Add a configurable HWHM parameter (default ~10 cm⁻¹ is typical).
3. Apply to both harmonic and anharmonic frequency outputs.
4. Plot the broadened spectrum overlaid with or replacing the stick spectrum.

Note: the user mentioned "Lorentz MoE" — confirm whether this means Lorentzian broadening or something else (mixture of experts?). Most likely it's Lorentzian peak broadening in the IR context.
