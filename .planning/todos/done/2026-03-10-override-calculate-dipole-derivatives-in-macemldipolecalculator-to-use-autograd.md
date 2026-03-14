---
created: 2026-03-10T15:17:30.442Z
title: Override calculate_dipole_derivatives in MACEMLDipoleCalculator to use autograd
area: general
files:
  - mace_gaussian/calculators/mace_ml.py
  - mace_gaussian/calculators/mace_loader.py
  - mace_gaussian/calculators/base.py
  - mace_dipole_pkg/build/lib/mace_dipole_core/calculators/mace.py:386-419
  - mace_dipole_pkg/build/lib/mace_dipole_core/modules/utils.py:422-464
---

## Problem

`MACEMLDipoleCalculator` inherits `calculate_dipole_derivatives()` from `DipoleCalculatorBase` which does finite differences: 2×3N forward passes of the MACE model. For a 10-atom molecule that's 60 forward passes.

The `AtomicDielectricMACE` model already computes the full Jacobian via autograd internally when called with `compute_dielectric_derivatives=True`. The ASE calculator exposes this as `get_dielectric_derivatives(atoms)` (mace.py:386), which calls the model with that flag and returns `dmu_dr` of shape `(3, N_atoms, 3)`. This gives the same result in ~3 backward passes instead of 6N forward passes.

We're not calling it. Instead we fall through to the base class finite differences.

## Solution

1. In `MACEDipoleCalculator` (mace_loader.py), add a `calculate_dipole_derivatives(atoms)` method that calls `self.calc.get_dielectric_derivatives(atoms)` and reshapes the output.
2. Output reshape: `dmu_dr` comes back as `(3, N_atoms, 3)` — (dipole_component, atom, spatial_dir). Required format for `write_gaussian_output` is `(3*N_atoms, 3)` — (atom*3 + dir, dipole_component). Conversion: `dmu_dr.transpose(1, 2, 0).reshape(3*N_atoms, 3)`.
3. Verify units are consistent with the finite-difference path (both should be e*Bohr / Angstrom).
4. Sanity-check on water: run both paths and assert the output matrices are close (within numerical tolerance).
5. `MACEMLDipoleCalculator` inherits from `MACEDipoleCalculator` indirectly, so the override propagates automatically.

Expected speedup: ~2N (where N = number of atoms). For a 10-atom molecule: ~20x. For a 30-atom molecule: ~60x.
