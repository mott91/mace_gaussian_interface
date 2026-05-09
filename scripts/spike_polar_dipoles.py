"""Spike: can MACE-POLAR-1 replace the separate MACE4IR dipole model?

For each test molecule we compute the molecular dipole three ways:
  1. POLAR-1-L (energy calculator, dipole as byproduct of multipole expansion)
  2. MACE4IR (current dedicated dipole model — what we ship today)
  3. B3LYP/6-31G(d,p) reference parsed from gaussian_freq.log

If POLAR-1 agrees with B3LYP within or beyond MACE4IR's accuracy, we can retire
the separate dipole pipeline (and the pickle-class-remapping hack in mace_loader.py).

Run:
    micromamba run -n mace4ir_v2 python scripts/spike_polar_dipoles.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import torch
from ase.io import read

from mace.calculators import mace_polar
from mace_gaussian.calculators.mace_ml import MACEMLDipoleCalculator

REPO = Path(__file__).resolve().parent.parent
DEBYE_PER_E_ANGSTROM = 4.80320451  # ASE stores dipoles in e·Å

MOLECULES = [
    "water",              # polar, bent
    "methanol",           # polar, asymmetric, hydrogens of two types
    "hydrogen_fluoride",  # very polar diatomic
    "hydrogen_cyanide",   # polar, linear, C≡N
    "methane",            # zero dipole by Td symmetry — sanity check
    "ethane",             # zero dipole — larger sanity check
]


def parse_dft_dipole_debye(freq_log: Path) -> float:
    """Return the last total dipole magnitude (Debye) from a Gaussian freq log."""
    text = freq_log.read_text()
    # "Tot=              2.0408" — take the last occurrence
    matches = re.findall(r"Tot=\s+([\-\d.]+)", text)
    if not matches:
        raise ValueError(f"no dipole found in {freq_log}")
    return float(matches[-1])


def to_debye(dipole_e_angstrom: np.ndarray) -> tuple[np.ndarray, float]:
    """Convert (3,) dipole vector in e·Å to Debye and return (vector, magnitude)."""
    vec = np.asarray(dipole_e_angstrom).reshape(-1)[:3] * DEBYE_PER_E_ANGSTROM
    return vec, float(np.linalg.norm(vec))


def polar_dipole(calc, atoms) -> np.ndarray:
    """Run a forward pass through PolarMACE and pull dipole out of calc.results."""
    atoms = atoms.copy()
    atoms.calc = calc
    atoms.get_potential_energy()  # triggers calculate()
    dip = calc.results.get("dipole")
    if dip is None:
        raise RuntimeError(
            "POLAR-1 calculator did not populate results['dipole']. "
            f"Available keys: {sorted(calc.results.keys())}"
        )
    if isinstance(dip, torch.Tensor):
        dip = dip.detach().cpu().numpy()
    return np.asarray(dip).reshape(-1)[:3]


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")
    print("Loading POLAR-1-L (one-time)...")
    polar_calc = mace_polar(model="polar-1-l", device=device, default_dtype="float64")

    print("Loading MACE4IR dipole calculator (one-time)...")
    mace4ir = MACEMLDipoleCalculator()
    mace4ir._check_availability()
    if not mace4ir.available:
        raise SystemExit("MACE4IR dipole calculator unavailable — cannot compare")

    print()
    header = f"{'molecule':<18} {'B3LYP':>8} {'POLAR-1':>10} {'MACE4IR':>10} {'Δpolar':>8} {'Δm4ir':>8}"
    print(header)
    print("-" * len(header))

    for mol in MOLECULES:
        opt_xyz = REPO / "comparison_results" / mol / "geometry_opt" / "optimized.xyz"
        raw_xyz = REPO / "molecules" / f"{mol}.xyz"
        b3lyp_dir = REPO / "comparison_results" / mol / "b3lyp_6-31Gdp"
        candidate_logs = [
            b3lyp_dir / "gaussian_freq.log",
            b3lyp_dir / f"{mol}_freq_anharm.log",
        ]
        freq_log = next((p for p in candidate_logs if p.exists()), None)
        xyz = opt_xyz if opt_xyz.exists() else raw_xyz
        if not xyz.exists() or freq_log is None:
            print(f"{mol:<18}  (missing artifacts — skipping)")
            continue

        atoms = read(str(xyz))
        dft_mag = parse_dft_dipole_debye(freq_log)

        polar_vec_eA = polar_dipole(polar_calc, atoms)
        polar_vec_d, polar_mag = to_debye(polar_vec_eA)

        m4ir_vec_eA, _ = mace4ir.calculate_dipole(atoms)
        m4ir_vec_d, m4ir_mag = to_debye(m4ir_vec_eA)

        print(
            f"{mol:<18} {dft_mag:>8.4f} {polar_mag:>10.4f} {m4ir_mag:>10.4f} "
            f"{polar_mag - dft_mag:>+8.4f} {m4ir_mag - dft_mag:>+8.4f}"
        )

    print()
    print("Magnitudes in Debye. Δ columns are signed errors vs B3LYP/6-31G(d,p).")
    print("If |Δpolar| ≲ |Δm4ir|, POLAR-1 is a viable drop-in for the dipole pipeline.")


if __name__ == "__main__":
    main()
