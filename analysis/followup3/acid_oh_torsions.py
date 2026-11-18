"""The O-H torsion of each carboxylic acid (and glycine's COOH): the B3LYP mode that turns the acid
hydrogen about the C-O bond most, with its harmonic wavenumber and region.

For every B3LYP mode, the mass-weighted mode is turned into a Cartesian displacement of 0.02 A norm
(as rotor_check.py); the torsion change is the mean change of the dihedrals X-C-O-H (X = the other
neighbours of the carboxyl C), and the kinetic-energy share on the acid hydrogen is |L_H|^2.
Output: out/acid_oh_torsions.csv and a printout.   python analysis/followup3/acid_oh_torsions.py
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
logging.disable(logging.WARNING)
import thesis_stats as ts  # noqa: E402
from mace_gaussian.analysis.mode_matching import extract_mode_data_from_checkpoint  # noqa: E402

ACIDS = ts.GROUPS["acids"] + ["glycine"]


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return np.arctan2(np.cross(b1, v) @ w, v @ w)


def main():
    rows = []
    for mol in ACIDS:
        fchk = next(iter(sorted((ts.CAMP / "comparison_results" / mol).glob("b3lyp*/*.fchk"))))
        modes, freqs, x, m, n = extract_mode_data_from_checkpoint(str(fchk), force_harmonic=True)
        x, z = np.asarray(x).reshape(n, 3), np.round(m).astype(int)
        d = np.linalg.norm(x[:, None] - x[None], axis=-1)
        # the acid hydrogen: H on an O whose other neighbour is a C bonded to two O (the carboxyl C)
        found = None
        for h in np.where(z == 1)[0]:
            o = int(np.argmin(np.where(np.arange(n) == h, 9, d[h])))
            if z[o] != 16:
                continue
            for c in np.where((z == 12) & (d[o] < 1.6))[0]:
                if ((z == 16) & (d[c] < 1.6)).sum() == 2:
                    found = (h, o, int(c))
        h, o, c = found
        xs = [a for a in range(n) if a not in (o,) and d[c, a] < 1.6 and z[a] != 1]
        out = []
        for k, f in enumerate(freqs):
            L = np.asarray(modes[k]).reshape(n, 3)
            v = L / np.sqrt(m)[:, None]
            v = v / np.linalg.norm(v) * 0.02
            t = np.mean([abs((dih(*(x + v)[[a, c, o, h]]) - dih(*x[[a, c, o, h]]) + np.pi) % (2 * np.pi) - np.pi)
                         for a in xs])  # fmt: skip
            ke = (L[h] ** 2).sum() / (L**2).sum()
            out.append((k + 1, f, np.degrees(t), ke))
        o_ = sorted(out, key=lambda r: -r[2])
        best, second = o_[0], o_[1]
        rows.append(dict(molecule=mol, mode=best[0], freq=round(best[1], 1), region=ts.region(best[1]),
                         torsion_deg=round(best[2], 2), ke_on_H=round(best[3], 2),
                         next_mode=second[0], next_freq=round(second[1], 1), next_torsion_deg=round(second[2], 2)))  # fmt: skip
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "out" / "acid_oh_torsions.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
