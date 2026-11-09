"""Hand check of the X-H internal-rotation labels below 500 cm-1 (B3LYP harmonic).

For every B3LYP mode below 500 cm-1 of the 20 pooled molecules: the bond whose torsion changes most
(summed over all dihedrals around it), the groups at its two ends, and where the kinetic energy sits
(H on O/N, H on C, heavy atoms). Printed for reading; the decision per band is written by hand into
rotor_labels.csv (column 'rotor', plus a short reason). Output: rotor_check.txt.
"""
import itertools, logging, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
logging.disable(logging.WARNING)
from mace_gaussian.analysis.mode_matching import extract_mode_data_from_checkpoint
R = Path(__file__).resolve().parents[2] / "campaigns/2026/comparison_results"
EL = {1: "H", 12: "C", 14: "N", 16: "O", 31: "P", 32: "S"}
COV = {"H": .31, "C": .76, "N": .71, "O": .66, "P": 1.07, "S": 1.05}
def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2; b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1; w = b2 - (b2 @ b1) * b1
    return np.arctan2(np.cross(b1, v) @ w, v @ w)
lines = []
for mol in sorted(p.name for p in R.iterdir() if p.is_dir() and p.name not in ("ammonia", "phosphine")):
    f = next(iter(sorted((R / mol).glob("b3lyp*/*_freq_anharm.fchk"))), None)
    if f is None: continue
    modes, freqs, x, m, n = extract_mode_data_from_checkpoint(str(f), force_harmonic=True)
    el = [EL[round(v)] for v in m]
    bonds = [(i, j) for i, j in itertools.combinations(range(n), 2) if np.linalg.norm(x[i]-x[j]) < 1.25*(COV[el[i]]+COV[el[j]])]
    nb = {i: {j for a, b in bonds for i_, j in ((a, b), (b, a)) if i_ == i} for i in range(n)}
    def group(a, b):  # what hangs on atom b, seen from a
        hs = sum(el[j] == "H" for j in nb[b] - {a}); heavy = [el[j] for j in nb[b] - {a} if el[j] != "H"]
        return f"{el[b]}H{hs if hs > 1 else ''}" if hs and not heavy else f"{el[b]}({''.join(sorted(heavy))}{'H'*hs})"
    hO = [i for i in range(n) if el[i] == "H" and any(el[j] in "ON" for j in nb[i])]
    hC = [i for i in range(n) if el[i] == "H" and i not in hO]
    for k, fr in enumerate(freqs):
        if fr >= 500: continue
        v = np.asarray(modes[k]).reshape(n, 3) / np.sqrt(m)[:, None]; v = v / np.linalg.norm(v) * 0.02
        tor = {}
        for b, c in bonds:
            if el[b] == "H" or el[c] == "H": continue
            s = 0.0
            for a in nb[b] - {c}:
                for e in nb[c] - {b}:
                    if a == e: continue
                    s += abs((dih(*(x + v)[[a, b, c, e]]) - dih(*x[[a, b, c, e]]) + np.pi) % (2*np.pi) - np.pi)
            tor[(b, c)] = s / max(1, (len(nb[b]) - 1) * (len(nb[c]) - 1))
        (b, c), t = max(tor.items(), key=lambda kv: kv[1])
        second = sorted(tor.values())[-2] if len(tor) > 1 else 0
        w = (v ** 2).sum(1); w /= w.sum()
        lines.append(f"{mol:17s} {k+1:3d} {fr:6.1f}  torsion {el[b]}{b+1}-{el[c]}{c+1} [{group(c, b)} | {group(b, c)}] "
                     f"{np.degrees(t):5.1f}° (2nd {np.degrees(second):4.1f}°)  KE: H(O/N) {w[hO].sum():.0%} H(C) {w[hC].sum():.0%} heavy {1-w[hO].sum()-w[hC].sum():.0%}")
Path(__file__).with_name("out").joinpath("rotor_check.txt").write_text("\n".join(lines))
print("\n".join(lines))
