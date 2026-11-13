"""Where does MACE4IR's huge C60 APT come from? B5 split (charge / charge flux / dipole flux)."""

import sys
from pathlib import Path

import numpy as np

F2 = Path("/home/mot/mace_gaussian/analysis/followup2")
sys.path.insert(0, "/home/mot/mace_gaussian")
sys.path.insert(0, str(F2))
from ase import Atoms  # noqa: E402
from b5_decomposition import decompose  # noqa: E402
from c60_dipoles import fchk  # noqa: E402

from mace_gaussian.calculators.mace_loader import MACEDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL  # noqa: E402

BOHR = 0.529177210903
d = MACEDipoleCalculator(DEFAULT_MACE_DIPOLE_MODEL, device="cpu")
d._ensure_calculator()
calc = d.calc

cases = {
    "C60 (B3LYP geometry)": Path("/home/mot/mace_gaussian/campaigns/c60/b3lyp_harmonic/c60_b3lyp_harm.fchk"),
    "glycine (B3LYP geometry)": Path(
        "/home/mot/mace_gaussian/campaigns/2026/comparison_results/glycine/b3lyp_6-31Gdp/glycine_freq_anharm.fchk"),
}
for name, p in cases.items():
    f = fchk(p)
    x = f["x"] * BOHR
    x = x - x.mean(0)
    atoms = Atoms(numbers=f["z"], positions=x)
    atoms.info["charge"], atoms.info["spin"] = 0.0, 1.0
    r = decompose(calc, "mace_ml", atoms)
    full_model = np.asarray(d.calculate_dipole_derivatives(atoms.copy()))
    tot = r["charge"] + r["cflux"] + r["dflux"]
    ref = f["apt"]
    print(f"== {name}: n={len(atoms)}, radius max {np.linalg.norm(x, axis=1).max():.2f} A")
    print(f"   q: min {r['q'].min():+.4f} max {r['q'].max():+.4f} sum {r['q'].sum():+.2e} e")
    print(f"   check: |split - model APT| max {np.abs(tot - full_model).max():.1e}, mu_err {r['mu_err']:.1e}, sum dq {r['sum_dq']:.1e}")
    for t in ("charge", "cflux", "dflux"):
        print(f"   {t:6s}: max |elem| {np.abs(r[t]).max():8.3f} e   Frobenius {np.linalg.norm(r[t]):8.2f}")
    print(f"   total : max |elem| {np.abs(tot).max():8.3f} e   Frobenius {np.linalg.norm(tot):8.2f}")
    if ref is not None:
        print(f"   B3LYP : max |elem| {np.abs(ref).max():8.3f} e   Frobenius {np.linalg.norm(ref):8.2f}")
    # charge flux per atom pair: how far the charge moves. max over i,j of |dq_i/dr_j| and of r_i * dq
    # Also: sum_i dq_i/dr_j = 0 (neutrality) but the lever arm is r_i: measure spread of dq over the molecule
    # dq is not returned; recompute magnitude via cflux / radius
