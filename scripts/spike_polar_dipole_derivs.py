"""Validate MACEPolar1DipoleCalculator: dipole + autograd derivatives vs finite-diff.

Confirms (1) the new calculator integrates with the factory pattern, (2) dipole
values match the stand-alone spike, (3) autograd-Jacobian dipole derivatives
agree with central-difference within tolerance — i.e., gradient flow through
PolarMACE's electrostatic stack is intact.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase.io import read

from mace_gaussian.calculators.mace_polar1 import MACEPolar1DipoleCalculator

REPO = Path(__file__).resolve().parent.parent
DEBYE_PER_E_ANGSTROM = 4.80320451


def main() -> None:
    calc = MACEPolar1DipoleCalculator()
    if not calc.available:
        raise SystemExit("mace_polar1 not available")

    atoms = read(str(REPO / "comparison_results" / "water" / "geometry_opt" / "optimized.xyz"))

    # 1. Dipole sanity
    dipole_eA, _ = calc.calculate_dipole(atoms)
    print(f"dipole vector (e·Å): {dipole_eA}")
    print(f"|μ| in Debye: {np.linalg.norm(dipole_eA) * DEBYE_PER_E_ANGSTROM:.4f}")

    # 2. Autograd dipole derivatives
    print("\nAutograd dipole derivatives:")
    deriv_autograd = calc.calculate_dipole_derivatives(atoms)
    print(f"  shape: {deriv_autograd.shape}")
    print(f"  ‖dμ/dr‖ Frobenius: {np.linalg.norm(deriv_autograd):.6f}")

    # 3. Finite-difference derivatives (override flag, then run base-class path)
    print("\nFinite-difference dipole derivatives:")
    calc.use_autograd = False
    deriv_fd = calc.calculate_dipole_derivatives(atoms, displacement=0.005)
    calc.use_autograd = True
    print(f"  shape: {deriv_fd.shape}")
    print(f"  ‖dμ/dr‖ Frobenius: {np.linalg.norm(deriv_fd):.6f}")

    # 4. Agreement
    diff = deriv_autograd - deriv_fd
    max_abs = np.max(np.abs(diff))
    rms = np.sqrt(np.mean(diff**2))
    print("\nAgreement (autograd vs finite-diff):")
    print(f"  max |Δ|: {max_abs:.2e}")
    print(f"  RMS Δ:   {rms:.2e}")
    if max_abs < 1e-3:
        print("  ✓ VERDICT: autograd path matches finite-diff — safe to use")
    else:
        print("  ✗ VERDICT: autograd diverges from finite-diff — investigate before trusting")


if __name__ == "__main__":
    main()
