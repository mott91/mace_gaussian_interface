---
created: 2026-03-26T19:15:00.000Z
title: Reevaluate MACE-POLAR-1 as dipole calculator
area: research
files:
  - mace_gaussian/calculators/mace_loader.py
  - mace_gaussian/workflow.py
---

## Problem

MACE-POLAR-1 is currently wired as energy calculator only. It has learnable charge distributions and electric field response, so it might be able to predict dipole moments and dipole derivatives (needed for IR intensities). This hasn't been tested.

## Solution

1. Investigate MACE-POLAR-1 API: can it return dipole moments? (charges × positions)
2. Check if get_dielectric_derivatives() works with mace_polar calculator.
3. If dipoles are accessible, benchmark against mace_ml (MACE4IR) dipole predictions on water/formaldehyde.
4. If quality is comparable, wire as third dipole calculator option alongside espaloma and mace_ml.
5. Key question: does the polarizable charge equilibration give better dipole derivatives than the dedicated MACE4IR dipole model?
