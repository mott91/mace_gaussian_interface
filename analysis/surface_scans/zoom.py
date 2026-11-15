"""Very fine curvature scan near the minimum, to resolve structure below the 0.002 A grid.

For one target and MACE model: k(Q) = d^T H d on dr = -0.012 ... +0.012 A in steps of 0.0002 A
(121 points), same machinery as scan_mace.py.   python zoom.py <target> <model> [...]
"""

from __future__ import annotations

import sys
import time

import numpy as np

import scan_mace as sm
from scan_mace import Atoms, c, wf

DR = np.round(np.arange(-0.012, 0.01201, 0.0002), 5)


def main(pairs):
    for key, model in pairs:
        t0 = time.time()
        calc = wf.calculator(model)
        mode = c.load_mode(key, model)
        q = mode.q_for_bond(DR)
        d = mode.d.ravel()
        k = []
        for qq in q:
            a = Atoms(numbers=mode.z, positions=mode.geometry(qq))
            a.info.update(charge=0.0, spin=1.0)
            a.calc = calc
            k.append(float(d @ (np.asarray(wf.calculate_hessian(a, calc, len(a))) * sm.TO_AJ_A2) @ d))
        c.save(f"zoom_{key}_{model}.json", {"key": key, "method": model, "q": q.tolist(), "dr": mode.bond(q).tolist(),
                                            "k": k, "seconds": time.time() - t0})  # fmt: skip
        print(f"zoom {key} {model} {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    main(list(zip(args[::2], args[1::2])))
