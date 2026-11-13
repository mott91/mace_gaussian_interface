import sys
from pathlib import Path
import numpy as np, torch
sys.path.insert(0, "/home/mot/mace_gaussian"); sys.path.insert(0, "/home/mot/mace_gaussian/analysis/followup2")
from ase import Atoms
from b5_decomposition import charges_dipoles
from c60_dipoles import fchk
from mace_gaussian.calculators.mace_loader import MACEDipoleCalculator
from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL
d = MACEDipoleCalculator(DEFAULT_MACE_DIPOLE_MODEL, device="cpu"); d._ensure_calculator(); calc = d.calc
print("cutoff", float(calc.models[0].r_max), "interactions", int(calc.models[0].num_interactions))
for name, p in {"C60": "/home/mot/mace_gaussian/campaigns/c60/b3lyp_harmonic/c60_b3lyp_harm.fchk",
                "glycine": "/home/mot/mace_gaussian/campaigns/2026/comparison_results/glycine/b3lyp_6-31Gdp/glycine_freq_anharm.fchk"}.items():
    f = fchk(Path(p)); x = f["x"] * 0.529177210903; x -= x.mean(0)
    a = Atoms(numbers=f["z"], positions=x); a.info["charge"], a.info["spin"] = 0.0, 1.0
    pos, q, _, _ = charges_dipoles(calc, "mace_ml", a); n = len(a)
    dq = np.stack([torch.autograd.grad(q[i], pos, retain_graph=True)[0].numpy() for i in range(n)])  # [i, j, a]
    g = np.linalg.norm(dq, axis=2)  # |dq_i/dr_j|
    dist = np.linalg.norm(x[:, None] - x[None], axis=2)
    j = 0
    print(f"== {name}: max |dq_i/dr_j| {g.max():.3f} e/A; sum_i |dq_i/dr_j| per moved atom: median {np.median(g.sum(0)):.2f}, max {g.sum(0).max():.2f} e/A")
    for lo, hi in [(0, 1.6), (1.6, 3.0), (3.0, 5.0), (5.0, 8.0)]:
        m = (dist >= lo) & (dist < hi)
        print(f"   pairs at {lo}-{hi} A: n {m.sum():5d}, mean |dq_i/dr_j| {g[m].mean() if m.any() else 0:.3f}, max {g[m].max() if m.any() else 0:.3f}")
