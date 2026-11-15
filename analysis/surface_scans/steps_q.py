"""Curvature at Gaussian's VPT2 displacement points, with the step in the normal coordinate.

The wobbly MACE-OMOL-0 curvature (zoom.py) shows that Gaussian's printed quartic constants are
reproduced with a displacement of 0.010 in the mass-weighted normal coordinate Q (amu^1/2 A), not
with 0.010 A of Cartesian displacement. This script evaluates k = d^T H d at Q = +/- s for
s = 0.005, 0.010, 0.020, 0.040 amu^1/2 A for every target and method, so that the finite-difference
constants can be compared with the VPT2 logs at exactly Gaussian's points.

python steps_q.py b3lyp [workers] | mace     -> data/stepsq_<method>_<target>.json
"""

from __future__ import annotations

import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

SQ = (0.005, 0.010, 0.020, 0.040)  # amu^1/2 Angstrom
QS = np.array([s * x for x in SQ for s in (-1, 1)])


def save(method, key, mode, k, k0, seconds):
    c.save(f"stepsq_{method}_{key}.json", {"key": key, "method": method, "steps_q": list(SQ), "q": QS.tolist(),
                                           "dr": mode.bond(QS).tolist(), "k": k, "k0": k0, "seconds": seconds})  # fmt: skip


def run_b3lyp(workers=3):
    import scan_b3lyp as sb

    for key, _mol, _i in c.TARGETS:
        t0 = time.time()
        mode = c.load_mode(key, "b3lyp")
        wd = sb.WORK / key
        jobs = []
        for j, q in enumerate(QS):
            p = wd / f"stepsq_{j:03d}.gjf"
            sb.write_input(p, mode, q, True)
            jobs.append(p)
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(lambda p: sb.run(p, True), jobs))
        k = [sb.parse(p, mode, True)[1] for p in jobs]
        k0 = sb.parse(wd / "zero_000.gjf", mode, True)[1]
        save("b3lyp", key, mode, k, k0, time.time() - t0)
        print(f"stepsq b3lyp {key} {time.time() - t0:.0f} s", flush=True)


def run_mace():
    import scan_mace as sm

    for model in c.MODELS:
        calc = sm.wf.calculator(model)
        for key, _mol, _i in c.TARGETS:
            t0 = time.time()
            mode = c.load_mode(key, model)
            d = mode.d.ravel()
            k = []
            for q in [0.0, *QS]:
                a = sm.Atoms(numbers=mode.z, positions=mode.geometry(q))
                a.info.update(charge=0.0, spin=1.0)
                a.calc = calc
                k.append(float(d @ (np.asarray(sm.wf.calculate_hessian(a, calc, len(a))) * sm.TO_AJ_A2) @ d))
            save(model, key, mode, k[1:], k[0], time.time() - t0)
            print(f"stepsq {model} {key} {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "b3lyp":
        run_b3lyp(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
    else:
        run_mace()
