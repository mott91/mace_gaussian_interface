"""2D map of one coupling: methanol O-H stretch against the C-O-H bend (B3LYP mode 5, the mode
with the largest |chi_ij| to the O-H stretch in B3LYP, -8.9 cm-1).

Each method on its own minimum and its own two normal modes (as in the 1D scans):
    x(Q1, Q2) = x0 + Q1 d1 + Q2 d2
Q1 on 21 points from dr(O-H) = -0.20 to +0.30 A; Q2 on 21 points from -0.45 to +0.45 amu^1/2 A
(about three times the zero-point turning point of the bend), with +Q2 opening the C-O-H angle.
Energies only: MACE in the campaign environment, B3LYP as single points with scf=tight.

python scan_2d.py mace | b3lyp [workers]   -> data/map2d_<method>.json
"""

from __future__ import annotations

import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

N = 21
DR1 = np.linspace(-0.20, 0.30, N)
Q2 = np.linspace(-0.45, 0.45, N)
BEND = ("methanol_bend", "methanol", 5)


def modes(method):
    m1 = c.load_mode("methanol_oh", method)
    m2 = c.load_mode(*BEND[:1], method, *BEND[1:])
    # +Q2 opens the C-O-H angle (the sign of a normal mode is arbitrary)
    o, h = m1.x_atom, m1.h_atom
    cc = int(np.argmin(np.where(m1.z == 6, np.linalg.norm(m1.x0 - m1.x0[o], axis=1), np.inf)))

    def angle(x):
        a, b = x[cc] - x[o], x[h] - x[o]
        return np.degrees(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b)))

    if angle(m2.x0 + 1e-3 * m2.d) < angle(m2.x0):
        m2.d = -m2.d
    assert np.allclose(m1.x0, m2.x0)
    return m1, m2


def grid(method):
    m1, m2 = modes(method)
    q1 = m1.q_for_bond(DR1)
    pts = [(i, j, m1.x0 + a * m1.d + b * m2.d) for i, a in enumerate(q1) for j, b in enumerate(Q2)]
    return m1, m2, q1, pts


def run_mace():
    sys.path.insert(0, str(c.ROOT))
    import logging
    import warnings

    logging.disable(logging.CRITICAL)
    warnings.filterwarnings("ignore")
    from ase import Atoms

    from mace_gaussian import workflow as wf

    for model in c.MODELS:
        t0 = time.time()
        calc = wf.calculator(model)
        m1, m2, q1, pts = grid(model)
        e = np.zeros((N, N))
        for i, j, x in pts:
            a = Atoms(numbers=m1.z, positions=x)
            a.info.update(charge=0.0, spin=1.0)
            a.calc = calc
            e[i, j] = a.get_potential_energy() * c.EV_AJ
        save(model, m1, m2, q1, e, time.time() - t0)


def run_b3lyp(workers=2):
    import scan_b3lyp as sb

    t0 = time.time()
    m1, m2, q1, pts = grid("b3lyp")
    wd = sb.WORK / "map2d"
    wd.mkdir(parents=True, exist_ok=True)
    jobs = []
    for i, j, x in pts:
        p = wd / f"g_{i:02d}_{j:02d}.gjf"
        lines = [f"%chk={p.stem}.chk", "%mem=2GB", f"%NProcShared={sb.NPROC}", sb.ROUTE_SP, "", f"map2d {i} {j}", "", "0 1"]
        lines += [f"{c.SYMBOL[z]} {a:16.10f} {b:16.10f} {cc:16.10f}" for z, (a, b, cc) in zip(m1.z, x)]
        p.write_text("\n".join(lines) + "\n\n")
        jobs.append(p)
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda p: sb.run(p, False), jobs))
    e = np.zeros((N, N))
    for i, j, _x in pts:
        e[i, j] = sb.parse(wd / f"g_{i:02d}_{j:02d}.gjf", m1, False)[0]
    save("b3lyp", m1, m2, q1, e, time.time() - t0)


def save(method, m1, m2, q1, e, seconds):
    c.save(f"map2d_{method}.json", {"method": method, "q1": q1.tolist(), "dr1": DR1.tolist(), "q2": Q2.tolist(),
                                    "E": e.tolist(), "freq1": m1.freq, "freq2": m2.freq, "own2": m2.own_mode,
                                    "seconds": seconds})  # fmt: skip
    print(f"map2d {method:11s} {seconds:.0f} s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "mace":
        run_mace()
    else:
        run_b3lyp(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
