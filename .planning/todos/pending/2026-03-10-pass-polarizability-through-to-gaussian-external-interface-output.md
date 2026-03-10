---
created: 2026-03-10T15:17:30.443Z
title: Pass polarizability through to Gaussian external interface output
area: general
files:
  - mace_gaussian/gaussian/io.py:95-113
  - mace_gaussian/workflow.py:105-151
  - mace_gaussian/calculators/mace_loader.py
  - mace_dipole_pkg/build/lib/mace_dipole_core/calculators/mace.py:386-419
---

## Problem

`write_gaussian_output` in `io.py:95` hardcodes polarizability as zeros:
```python
polarizability = np.zeros(6)
```

The `AtomicDielectricMACE` model (type `DipolePolarizabilityMACE`) already computes the full polarizability tensor during its forward pass. `get_dielectric_derivatives()` returns both `dmu_dr` AND `dalpha_dr` when `use_polarizability=True`. We're discarding this output entirely.

Gaussian's external interface accepts polarizability (6 independent components of the symmetric 3×3 tensor, upper triangle: xx, xy, xz, yy, yz, zz) in the output file. If populated, Gaussian can use it for Raman activity calculations.

## Solution

1. Add a `calculate_polarizability(atoms)` method to `MACEDipoleCalculator` that calls `self.calc.get_dielectric_derivatives(atoms)` and extracts `dalpha_dr` (or separately calls the forward pass for just the polarizability tensor).
2. Thread the polarizability value through `calculate_dipole_properties()` in workflow.py (currently returns only dipole and dipole_derivatives).
3. Update `write_gaussian_output` signature and body to accept and write the real polarizability instead of zeros.
4. Verify units: Gaussian expects polarizability in Bohr³ (atomic units). Check what unit the model returns and convert if needed.

Note: this requires checking whether the loaded model file has `use_polarizability=True` (it's type `DipolePolarizabilityMACE` so it should). Verify at runtime.
