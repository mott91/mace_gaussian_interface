"""Follow-up 2, B6: smoothness of the dipole and energy surfaces on the finite-difference scale.

At the B3LYP/6-31G(d,p) minimum, along the B3LYP normal modes (rebuilt as in b1_b2_analyse.py):
each mode k gives a Cartesian unit vector e_k; displacements are x0 +/- h e_k with
h = 0.0025, 0.005, 0.01, 0.02, 0.04 Angstrom.

- dipole models (APT from the pipeline routine at every displaced point):
    M_kj(h) = e_j . [APT(+h e_k) - APT(-h e_k)] / 2h      (3-vector, d2mu / dq_k dq_j)
  This is how VPT2 obtains the second dipole derivatives (finite differences of APTs).
- energy models (analytic Hessian at every displaced point):
    C_kjj(h) = e_j . [H(+h e_k) - H(-h e_k)] / 2h . e_j   (cubic constant)
  and C_jkj from displacing j instead.
- Noise measure: the same mixed derivative from two different displacements must agree.
  asym(h) = sum |M_kj - M_jk| / sum |M_kj + M_jk|  (dipole)
  asym(h) = sum |T_kab - T_akb| / sum |T_kab + T_akb| (energy, T_kab = e_a.dH/dq_k.e_b), k != a.
  A smooth surface gives asym -> O(h^2) as h -> 0; noise gives asym growing as 1/h.

Usage: python b6_smoothness.py dipole|energy <model> <mol> [<mol> ...]
Appends one JSON line per (molecule, model) to b6_results.jsonl.
"""

from __future__ import annotations

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")
from ase import Atoms  # noqa: E402
from b1_b2_analyse import fchk, normal_modes  # noqa: E402

from mace_gaussian import workflow as wf  # noqa: E402
from mace_gaussian.calculators import dipole_factory  # noqa: E402

BOHR = 0.529177210903
STEPS = (0.0025, 0.005, 0.01, 0.02, 0.04)
Z = {"methanol": None}


def unit_modes(mol):
    f = fchk(mol)
    _freqs, disp = normal_modes(f["x"], f["m"], f["hess"])
    e = disp / np.linalg.norm(disp, axis=0)  # Cartesian unit vectors (3N, nm)
    t = (ROOT / "campaigns/2026" / mol / f"{mol}_freq_anharm.fchk").read_text()
    i = t.index("Atomic numbers")
    n = int(t[i : t.index("\n", i)].split("N=")[1])
    z = np.array([int(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])
    return z, f["x"] * BOHR, e


def main():
    kind, model, *mols = sys.argv[1:]
    if kind == "dipole":
        calc = dipole_factory.get_calculator(model)

        def evaluate(atoms):
            return np.asarray(wf.calculate_dipole_properties(atoms, calc, 2, True)[1])  # (3N, 3)
    else:
        calc = wf.calculator(model)

        def evaluate(atoms):
            atoms.calc = calc
            return np.asarray(wf.calculate_hessian(atoms, calc, len(atoms)))  # (3N, 3N)

    for mol in mols:
        z, x0, e = unit_modes(mol)
        nm = e.shape[1]
        res = dict(mol=mol, kind=kind, model=model, steps=list(STEPS), asym=[], diag_rel_change=[])
        diag_ref = None
        per_h = {}
        for h in STEPS:
            d = []
            for k in range(nm):
                vals = []
                for s in (+1, -1):
                    a = Atoms(numbers=z, positions=x0 + s * h * e[:, k].reshape(-1, 3))
                    a.info.update(charge=0.0, spin=1.0)
                    vals.append(evaluate(a))
                d.append((vals[0] - vals[1]) / (2 * h))
            d = np.array(d)  # dipole: (nm, 3N, 3); energy: (nm, 3N, 3N)
            if kind == "dipole":
                m = np.einsum("kic,ij->kjc", d, e)  # M[k, j, :]
                num = np.abs(m - m.transpose(1, 0, 2)).sum()
                den = np.abs(m + m.transpose(1, 0, 2)).sum()
                diag = np.array([m[k, k] for k in range(nm)])
            else:
                # T[k, a, b] = e_a . dH/dq_k . e_b is the cubic constant Phi_kab obtained by
                # displacing along k. The same constant from displacing along a is T[a, k, b];
                # the two must agree (fully symmetric tensor). Compare k <-> a, off-diagonal.
                t = np.einsum("kxy,xa,yb->kab", d, e, e)
                ts = t.transpose(1, 0, 2)
                off = ~np.eye(nm, dtype=bool)[:, :, None].repeat(nm, axis=2)
                num = np.abs(t - ts)[off].sum()
                den = np.abs(t + ts)[off].sum()
                diag = np.array([t[k, k, k] for k in range(nm)])
            per_h[h] = diag
            res["asym"].append(float(num / den))
        diag_ref = per_h[0.02]
        for h in STEPS:
            res["diag_rel_change"].append(float(np.abs(per_h[h] - diag_ref).sum() / np.abs(diag_ref).sum()))
        print(json.dumps(res), flush=True)
        with (HERE / "b6_results.jsonl").open("a") as fh:
            fh.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
