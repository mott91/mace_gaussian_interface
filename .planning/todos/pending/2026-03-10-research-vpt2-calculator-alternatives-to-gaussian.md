---
created: 2026-03-10T14:14:33.441Z
title: Research VPT2 calculator alternatives to Gaussian
area: general
files: []
---

## Problem

The anharmonic frequency calculations currently depend entirely on Gaussian 16 for VPT2 (second-order vibrational perturbation theory). This creates a hard dependency on a commercial, licensed piece of software. Alternative open-source or more accessible VPT2 implementations may exist (e.g., in CFOUR, ORCA, PyVCI, or custom Python implementations).

## Solution

1. Survey existing VPT2 implementations outside of Gaussian (CFOUR, ORCA 5+, Psi4, PyVCI, etc.).
2. Assess feasibility of integration: API/interface compatibility, force field input format, license.
3. Identify the best candidate(s) and design an integration path that mirrors the existing Gaussian workflow.
4. Prototype if a strong candidate is found.
