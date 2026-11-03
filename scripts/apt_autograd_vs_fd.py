#!/usr/bin/env python3
"""Autograd dipole derivatives (atomic polar tensors) against central finite differences.

For each dipole model of the campaign, the derivatives the pipeline hands to Gaussian
(autograd) are compared with central differences of the same model's dipole, computed by
the base-class method (``DipoleCalculatorBase.calculate_dipole_derivatives``), at several
step sizes. Units: e (atomic units), as passed to Gaussian.

Usage::

    python scripts/apt_autograd_vs_fd.py OUT_DIR [molecule ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mace_gaussian.calculators import dipole_factory
from mace_gaussian.calculators.base import DipoleCalculatorBase

CAMPAIGN = Path(__file__).resolve().parent.parent / "campaigns" / "2026" / "comparison_results"
MODELS = ["mace_ml", "mace_mdp", "mace_polar1"]
STEPS_A = [1e-2, 5e-3, 1e-3]
REL_FLOOR = 0.05  # e; relative deviation only over entries larger than this


def compare(calc, atoms) -> dict:
    auto = np.asarray(calc.calculate_dipole_derivatives(atoms.copy()))
    out = {"max_abs_entry": float(np.abs(auto).max())}
    big = np.abs(auto) > REL_FLOOR
    for h in STEPS_A:
        fd = DipoleCalculatorBase.calculate_dipole_derivatives(calc, atoms.copy(), displacement=h)
        d = np.abs(np.asarray(fd) - auto)
        out[f"{h:g}"] = {
            "max_abs_dev_e": float(d.max()),
            "max_rel_dev": float((d[big] / np.abs(auto[big])).max()) if big.any() else None,
            "sum_rule_fd_e": float(np.abs(np.asarray(fd).reshape(-1, 3, 3).sum(axis=0)).max()),
        }
    # Charge-neutral molecule: the APTs of all atoms sum to zero (translation invariance).
    out["sum_rule_autograd_e"] = float(np.abs(auto.reshape(-1, 3, 3).sum(axis=0)).max())
    return out


def main(out_dir: str, molecules: list[str]) -> None:
    results = {}
    for mol in molecules:
        atoms = read(str(CAMPAIGN / mol / "geometry_opt" / "optimized.xyz"))
        atoms.info["charge"], atoms.info["spin"] = 0.0, 1.0
        for m in MODELS:
            res = compare(dipole_factory.get_calculator(m), atoms)
            results[f"{mol}/{m}"] = res
            parts = []
            for h in STEPS_A:
                r = res[f"{h:g}"]
                parts.append(f"h={h:g}: abs {r['max_abs_dev_e']:.1e} rel {r['max_rel_dev']:.1e}")
            steps = "  ".join(parts)
            print(f"{mol:10s} {m:12s} max|APT| {res['max_abs_entry']:.2f} e  {steps}  "
                  f"sum rule {res['sum_rule_autograd_e']:.1e}")  # fmt: skip
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "apt_autograd_vs_fd.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or ["carbon_dioxide", "methanol", "glycine"])
