"""Which unit Gaussian's VPT2 "Nuclear step= 0.010000 Angstroms" really has, from the log itself.

Gaussian prints the nuclear repulsion energy V_nn (10 decimals, Hartree) at every displaced point of
the anharmonic run ("Re-enter D2Numr: normal mode n step-up/step-down"). V_nn depends only on the
geometry, so the printed values can be compared with V_nn at x0 +/- s d for the two readings of the
printed step s (0.010 by default, 1e-4 N with freq=(anharm,step=N)):
  - Q: s amu^1/2 A in the mass-weighted normal coordinate, d = L / sqrt(m) (L normalized);
  - cart: s A of Cartesian displacement along the unit vector d / |d|.
x0 and the mode vectors are the run's fchk (Current cartesian coordinates, Vib-Modes).

python step_unit_check.py <log> <fchk>   (B3LYP runs only: ML runs via `external` print no V_nn)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402


def vnn(z, x_ang):
    x = x_ang / c.BOHR_A
    e = 0.0
    for a in range(len(z)):
        for b in range(a):
            e += z[a] * z[b] / np.linalg.norm(x[a] - x[b])
    return e


def printed(log: Path) -> tuple[float, float, dict]:
    """The printed step, V_nn of the reference point (first one in the freq step) and
    {(mode, 'up'/'down'): V_nn}."""
    t = log.read_text()
    t = t[t.index("Numerical evaluation of third derivatives") :]
    step = float(re.search(r"Nuclear step=\s*(\S+) Angstroms", t).group(1))
    ref = float(re.search(r"nuclear repulsion energy\s+(\S+) Hartrees", t).group(1))
    # "Re-enter D2Numr: normal mode n step-up" is printed when that point is done: its V_nn is the last
    # one printed before the line
    out = {}
    for mt in re.finditer(r"Re-enter D2Numr:\s+normal mode\s+(\d+) step-(up|down)\s+number 1", t):
        e = float(re.findall(r"nuclear repulsion energy\s+(\S+) Hartrees", t[: mt.start()])[-1])
        out[(int(mt.group(1)), mt.group(2))] = e
    return step, ref, out


def main(log: Path, fchk: Path) -> None:
    a = c.fchk_arrays(fchk, ["Atomic numbers", "Current cartesian coordinates", "Real atomic weights", "Vib-E2", "Vib-Modes"])
    z = a["Atomic numbers"].astype(int)
    n = len(z)
    x0 = a["Current cartesian coordinates"].reshape(n, 3) * c.BOHR_A
    m = a["Real atomic weights"]
    nm = 3 * n - 6
    freqs = a["Vib-E2"][:nm]
    modes = a["Vib-Modes"].reshape(nm, n, 3)
    step, ref, pr = printed(log)
    v0 = vnn(z, x0)
    print(f"{log.name}: printed Nuclear step= {step:.6f} Angstroms; V_nn reference printed {ref:.10f}, from fchk geometry {v0:.10f}")
    print(f"{'mode':>4} {'cm-1':>7} {'|d|':>6}  {'printed dV (up)':>15}  {'Q':>12}  {'cart':>12}   "
          f"{'printed dV (down)':>17}  {'Q':>12}  {'cart':>12}")
    dev = {"Q": [], "cart": []}
    for k in range(nm):
        L = modes[k] * np.sqrt(m)[:, None]
        L /= np.linalg.norm(L)
        d = L / np.sqrt(m)[:, None]
        e = d / np.linalg.norm(d)
        # Gaussian's mode numbering in D2Numr: ascending (checked by the match below); sign of its vector
        # unknown, so compare up/down as a pair against +/- in both orders
        cand = {"Q": [vnn(z, x0 + s * step * d) - v0 for s in (1, -1)],
                "cart": [vnn(z, x0 + s * step * e) - v0 for s in (1, -1)]}
        up, dn = pr[(k + 1, "up")] - ref, pr[(k + 1, "down")] - ref
        for key, (p, q) in cand.items():
            dev[key].append(min(abs(up - p) + abs(dn - q), abs(up - q) + abs(dn - p)))
        print(f"{k + 1:4d} {freqs[k]:7.1f} {np.linalg.norm(d):6.3f}  {up:15.10f}  {cand['Q'][0]:12.10f}  "
              f"{cand['cart'][0]:12.10f}   {dn:17.10f}  {cand['Q'][1]:12.10f}  {cand['cart'][1]:12.10f}")
    for key, v in dev.items():
        print(f"{key:5s}: max |printed - computed| over up+down = {max(v):.2e} Hartree")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
