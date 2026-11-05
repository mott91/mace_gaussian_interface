"""Follow-up 2, B6 (addendum): step dependence of the X-H cubic and quartic constants.

The consistency test in b6_smoothness.py detects noise in first differences of the Hessian
(cubic constants). The x_ii of VPT2 also contain the quartic constants, which come from second
differences, H(+h) + H(-h) - 2 H(0), divided by h^2: an irregularity eps of the surface enters
them as eps / h^2. Here, along each X-H stretch of the B3LYP modes (harmonic > 2700 cm-1) at the
B3LYP minimum, for h = 0.0025 ... 0.04 Angstrom:
    cubic   Phi_kkk(h)  = e_k . [H(+h) - H(-h)] . e_k / 2h
    quartic Phi_kkkk(h) = e_k . [H(+h) + H(-h) - 2 H(0)] . e_k / h^2
(units: Hessian as returned by the pipeline, per Angstrom of displacement along the unit vector).
A smooth surface gives values that converge as h -> 0; irregularities make the quartic scatter
or blow up at small h.

Usage: python b6_quartic.py <energy model> <mol> [<mol> ...]; appends to b6_quartic.jsonl.
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
from b6_smoothness import STEPS, unit_modes  # noqa: E402

from mace_gaussian import workflow as wf  # noqa: E402


def main():
    model, *mols = sys.argv[1:]
    calc = wf.calculator(model)

    def hess(z, x):
        a = Atoms(numbers=z, positions=x)
        a.info.update(charge=0.0, spin=1.0)
        a.calc = calc
        return np.asarray(wf.calculate_hessian(a, calc, len(a)))

    for mol in mols:
        z, x0, e = unit_modes(mol)
        f = fchk(mol)
        freqs, _ = normal_modes(f["x"], f["m"], f["hess"])
        xh = [k for k in range(len(freqs)) if freqs[k] > 2700]
        h0 = hess(z, x0)
        res = dict(mol=mol, model=model, steps=list(STEPS), modes=[float(freqs[k]) for k in xh],
                   cubic=[], quartic=[])  # fmt: skip
        for k in xh:
            ek = e[:, k]
            cub, qua = [], []
            for h in STEPS:
                hp = hess(z, x0 + h * ek.reshape(-1, 3))
                hm = hess(z, x0 - h * ek.reshape(-1, 3))
                cub.append(float(ek @ (hp - hm) @ ek / (2 * h)))
                qua.append(float(ek @ (hp + hm - 2 * h0) @ ek / h**2))
            res["cubic"].append(cub)
            res["quartic"].append(qua)
        print(json.dumps(res), flush=True)
        with (HERE / "b6_quartic.jsonl").open("a") as fh:
            fh.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
