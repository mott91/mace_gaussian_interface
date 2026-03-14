---
created: 2026-03-10T14:14:33.439Z
title: Research MACE4IR dipole derivative acceleration via autograd
area: general
files:
  - gm_main.py
  - calculators/mace_loader.py
---

## Problem

It is currently unclear how MACE4IR computes dipole derivatives and IR intensities. Dipole derivatives are typically computed via finite differences (expensive) but could potentially be accelerated using PyTorch autograd/backpropagation through the MACE model. This could significantly speed up the ZMQ-based injection into Gaussian.

## Solution

1. Read `gm_main.py` and `calculators/mace_loader.py` carefully to understand the current dipole derivative computation method.
2. Check whether the MACE4IR dipole model supports autograd (i.e., is the graph differentiable w.r.t. atomic positions?).
3. If yes, prototype a backprop-based dipole derivative calculator and benchmark against the current approach.
4. If feasible, replace or offer as an alternative computation path.
