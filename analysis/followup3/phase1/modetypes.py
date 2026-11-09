"""Rough type of every B3LYP mode (all molecules): dominant bond stretch (with atom types) or bend/torsion.
Stretch if the largest bond-length change is > 0.5 of the largest atom displacement. Output modetypes.csv."""
import sys, itertools, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, "/home/mot/mace_gaussian")
from mace_gaussian.analysis.mode_matching import extract_mode_data_from_checkpoint
R = Path("/home/mot/mace_gaussian/campaigns/2026/comparison_results")
EL = {1: "H", 12: "C", 14: "N", 16: "O", 31: "P", 32: "S"}
COV = {"H": .31, "C": .76, "N": .71, "O": .66, "P": 1.07, "S": 1.05}
rows = []
for mol in sorted(p.name for p in R.iterdir() if p.is_dir()):
    f = next(iter(sorted((R / mol).glob("b3lyp*/*_freq_anharm.fchk"))), None)
    if f is None: continue
    modes, freqs, x, m, n = extract_mode_data_from_checkpoint(str(f), force_harmonic=True)
    el = [EL[round(v)] for v in m]
    bonds = [(i, j) for i, j in itertools.combinations(range(n), 2)
             if np.linalg.norm(x[i] - x[j]) < 1.25 * (COV[el[i]] + COV[el[j]])]
    for k, fr in enumerate(freqs):
        v = np.asarray(modes[k]).reshape(n, 3) / np.sqrt(m)[:, None]; v = v / np.abs(v).max() * 0.05
        dr = [abs(np.linalg.norm((x + v)[i] - (x + v)[j]) - np.linalg.norm(x[i] - x[j])) for i, j in bonds]
        b = int(np.argmax(dr)); i, j = bonds[b]
        stretch = dr[b] > 0.5 * np.linalg.norm(v, axis=1).max()
        kind = "-".join(sorted((el[i], el[j]), key=lambda e: "CNOPSH".index(e))) + " str" if stretch else "bend/tors"
        rows.append(dict(molecule=mol, i=k + 1, freq=fr, kind=kind))
pd.DataFrame(rows).to_csv("modetypes.csv", index=False)
