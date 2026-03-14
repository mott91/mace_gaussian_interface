---
created: 2026-03-10T15:17:30.447Z
title: Pass dalpha_dr Raman derivatives through to Gaussian output
area: general
files:
  - mace_gaussian/calculators/mace_loader.py
  - mace_gaussian/gaussian/io.py
  - mace_gaussian/workflow.py
  - mace_dipole_pkg/build/lib/mace_dipole_core/calculators/mace.py:386-419
---

## Problem

`get_dielectric_derivatives()` on the ASE calculator returns both `dmu_dr` (dipole Jacobian, used for IR intensities) AND `dalpha_dr` (polarizability Jacobian) when the model has `use_polarizability=True`. Our model is type `DipolePolarizabilityMACE` so it computes both. We discard `dalpha_dr` entirely.

`dalpha_dr` is the derivative of the polarizability tensor w.r.t. atomic displacements — this is exactly what's needed for Raman activity. Gaussian's external interface format has a slot for it (after the dipole derivatives section). If we pass it through, Gaussian would compute Raman spectra.

## Solution

1. Expose `dalpha_dr` alongside `dmu_dr` in `MACEDipoleCalculator.calculate_dipole_derivatives()` — return both or add a separate `calculate_polarizability_derivatives()` method.
2. Thread `dalpha_dr` through `calculate_dipole_properties()` in `workflow.py`.
3. Add it to `write_gaussian_output()` in `io.py` — check Gaussian's external interface spec for exact format (it follows the dipole derivative block).
4. Verify units: Gaussian expects polarizability derivatives in Bohr²/e or similar — confirm unit convention from the model.

Note: do this after the dipole derivative autograd fix (Phase 13.1) since the plumbing will be freshly laid.
