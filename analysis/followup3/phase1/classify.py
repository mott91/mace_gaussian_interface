"""Rough mode-type label for B3LYP modes: dominant internal-coordinate change (torsion vs bend vs stretch)."""
import sys, itertools, numpy as np, pandas as pd
sys.path.insert(0, "/home/mot/mace_gaussian")
from mace_gaussian.analysis.mode_matching import extract_mode_data_from_checkpoint
from pathlib import Path
R = Path("/home/mot/mace_gaussian/campaigns/2026/comparison_results")
COV = {1: .31, 6: .76, 7: .71, 8: .66, 15: 1.07, 16: 1.05}
from mace_gaussian.gaussian.fchk import parse_fchk_section
def ang(a, b, c):
    u, v = a - b, c - b; return np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2; b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1; w = b2 - (b2 @ b1) * b1
    return np.arctan2(np.cross(b1, v) @ w, v @ w)
rows = []
for mol in sorted(p.name for p in R.iterdir() if p.is_dir()):
    d = next(iter(sorted((R / mol).glob("b3lyp*/*.fchk"))), None)
    if d is None: continue
    modes, freqs, xyz, masses, n = extract_mode_data_from_checkpoint(str(d), force_harmonic=True)
    z = np.array([min({1:1.008,6:12.011,7:14.007,8:15.999,15:30.974,16:32.06}.items(), key=lambda kv: abs(kv[1]-mm))[0] for mm in masses])
    xyz = np.asarray(xyz).reshape(n, 3)
    if xyz.max() > 20: xyz = xyz * 0.529177
    bonds = [(i, j) for i, j in itertools.combinations(range(n), 2)
             if np.linalg.norm(xyz[i] - xyz[j]) < 1.25 * (COV[z[i]] + COV[z[j]])]
    nb = {i: set() for i in range(n)}
    for i, j in bonds: nb[i].add(j); nb[j].add(i)
    angles = [(a, b, c) for b in range(n) for a, c in itertools.combinations(sorted(nb[b]), 2)]
    tors = [(a, b, c, e) for b, c in bonds for a in nb[b] - {c} for e in nb[c] - {b} if a != e]
    m = np.sqrt(masses)
    for k, f in enumerate(freqs):
        if f >= 1000: continue
        v = np.asarray(modes[k]).reshape(n, 3) / m[:, None]   # cartesian displacement
        v = v / np.linalg.norm(v) * 0.05
        def d_of(fn, idx, wrap=False):
            out = []
            for t in idx:
                x0 = fn(*[xyz[i] for i in t]); x1 = fn(*[(xyz + v)[i] for i in t])
                dd = x1 - x0
                if wrap: dd = (dd + np.pi) % (2 * np.pi) - np.pi
                out.append(abs(dd))
            return np.array(out) if out else np.zeros(1)
        dt = d_of(dih, tors, True); da = d_of(ang, angles); 
        dr = np.array([abs(np.linalg.norm((xyz+v)[i]-(xyz+v)[j]) - np.linalg.norm(xyz[i]-xyz[j])) for i, j in bonds])
        # torsion around bond b-c: max over that bond; H rotor if the moving end is terminal H
        best = int(np.argmax(dt)) if len(tors) else None
        kind = "bend"
        if len(tors) and dt.max() > 1.5 * da.max():
            a, b, c, e = tors[best]
            hrot = z[a] == 1 or z[e] == 1
            kind = "torsion (X-H rotor)" if hrot else "torsion (skeletal)"
        if dr.max() > 0.02: kind += "/stretch"
        rows.append(dict(molecule=mol, i=k + 1, freq=f, kind=kind, dtors=dt.max(), dang=da.max()))
pd.DataFrame(rows).to_csv("modes1000.csv", index=False)
print(pd.DataFrame(rows).groupby("kind").size())
