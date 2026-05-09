---
created: 2026-03-10T14:14:33.437Z
title: Implement MACE polar model as calculator
area: general
files: []
---

## Problem

A new MACE polar model exists that may be capable of computing molecular polarizabilities and/or dipole derivatives. We don't yet know its full capabilities or how it fits into the existing calculator architecture. Need to research first, then integrate.

## Solution

1. Research the MACE polar model: what it predicts (polarizability tensor, dipole, dipole derivatives?), what inputs it requires, and what accuracy it achieves.
2. Determine whether it can replace or complement `mace_loader.py` / the MACE4IR dipole model.
3. Implement as a new calculator following the existing pattern in `calculators/`.
4. Wire it into the `--dipole-calculators` option in `mace-gaussian run`.
