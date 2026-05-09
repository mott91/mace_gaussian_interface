---
created: 2026-03-26T19:00:00.000Z
title: Integrate PyVPT2 and alternative anharmonic methods (eliminate Gaussian dependency)
area: workflow
files:
  - mace_gaussian/workflow.py
---

## Problem

The pipeline requires Gaussian 16 (proprietary, expensive license) for VPT2 anharmonic frequency calculations. This limits reproducibility, adoption, and who can use the tool. The ZMQ bridge adds complexity.

## Alternatives Found

### PyVPT2 (TOP PRIORITY)
- Pure Python VPT2 via QCEngine backend (already has ASE calculator support)
- Paper: Nelson, J. Chem. Phys. 2025
- GitHub: philipmnel/pyvpt2
- Integration: MACE ASE calculator → QCEngine ASE bridge → PyVPT2
- Would eliminate Gaussian, ZMQ, .gjf/.log parsing entirely

### iGVPT2
- C/C++ GVPT2 with plugin system (write bash/Python script that calls MACE)
- Better Fermi resonance handling (GVPT2, DCPT2 variants)
- GitHub: allouchear/iGVPT2

### PyVCI
- Python/C variational VCI (more accurate than VPT2 for strongly anharmonic modes)
- SourceForge: pyvci
- Needs pre-computed force constant tensors

### Psience/PyVibPTn
- Python VPT2 from McCoy group (U. Washington)
- GitHub: McCoyGroup/PyVibPTn

### Prior art: ML + VPT2
- PhysNet + VPT2 (Meuwly 2021, JCTC)
- ANI-1ccx-gelu + VPT2 (2025, JPCL)
- MACE-OFF23 composite IR (2024, JCTC)

## Solution

1. Prototype PyVPT2 integration: write QCEngine harness for MACE energy+dipole calculators.
2. Run on water/formaldehyde and compare against Gaussian VPT2 results (should be identical physics, different implementation).
3. If successful, add as alternative backend in workflow.py (--vpt2-engine=gaussian|pyvpt2).
4. Evaluate iGVPT2 for cases where Fermi resonance handling matters.
5. Long-term: PyVCI for strongly anharmonic systems where VPT2 breaks down.
