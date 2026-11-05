"""Follow-up 2, B3 (analysis): how good is the reference, and is the ranking a training-level effect?

Every level's APT at the B3LYP/6-31G(d,p) minimum (same geometry, same frame), projected on the
same B3LYP/6-31G(d,p) normal modes as in B1. Intensities then differ only by the level of theory
(reference rows) or by the dipole model (ML rows, from b1_apt_<model>.npz).

Levels (b3_reference/<level>/<mol>/):
  b3lyp_631gdp   campaign fchk (the thesis reference)
  b3lyp_avtz     g_b3lyp_avtz_fix  (fchk "Dipole Derivatives")
  pbe0_tzvp      g_pbe0_tzvp_fix   (MACE4IR stand-in)
  wb97mv_tzvpd   o_wb97mv_fix      (ORCA .hess $dipole_derivatives; POLAR-1 / OMol25 level)
  wb97md3_tzvppd o_wb97md3bj_fix   (MDP / SPICE level)
Also g_b3lyp_avtz_opt: B3LYP/aug-cc-pVTZ at its own minimum, harmonic intensities from its log,
paired with the 6-31G(d,p) modes by Kabsch-aligned mode overlap (Hungarian), as in the campaign.

Writes b3_tables.md.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
logging.disable(logging.CRITICAL)
from b1_b2_analyse import KMMOL, fchk, merge_degenerate, normal_modes, region  # noqa: E402

from mace_gaussian.analysis.mode_matching import align_modes_to_reference, match_modes  # noqa: E402

B3 = HERE / "b3_reference"
MOLS = ["methanol", "formic_acid_trans", "acetic_acid", "carbon_dioxide", "nitric_acid",
        "carbonic_acid", "hydrogen_peroxide", "glycine"]  # fmt: skip
LEVELS = {"b3lyp_631gdp": "B3LYP/6-31G(d,p)", "b3lyp_avtz": "B3LYP/aug-cc-pVTZ",
          "pbe0_tzvp": "PBE0/def2-TZVP", "wb97mv_tzvpd": "ωB97M-V/def2-TZVPD",
          "wb97md3_tzvppd": "ωB97M-D3(BJ)/def2-TZVPPD"}  # fmt: skip
ML = {"mace_polar1": "POLAR-1", "polar1_field": "POLAR-1 (field)", "mace_mdp": "MDP",
      "mace_ml": "MACE4IR"}  # fmt: skip


def fchk_any(path, key):
    t = path.read_text()
    i = t.index(key + " ")
    n = int(t[i : t.index("\n", i)].split("N=")[1])
    return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])


def orca_hess(path):
    """$dipole_derivatives (3N x 3) and $atoms coordinates (Bohr) from an ORCA .hess file."""
    lines = path.read_text().split("\n")
    i = lines.index("$dipole_derivatives")
    n3 = int(lines[i + 1])
    apt = np.array([[float(v) for v in lines[i + 2 + k].split()] for k in range(n3)])
    j = lines.index("$atoms")
    na = int(lines[j + 1])
    x = np.array([[float(v) for v in lines[j + 2 + k].split()[2:5]] for k in range(na)])
    return apt, x


def level_apt(level, mol, f0):
    """APT (3N, 3) in e at the B3LYP/6-31G(d,p) geometry, plus max coordinate mismatch (Bohr)."""
    if level == "b3lyp_631gdp":
        return f0["apt"], 0.0
    if level.startswith("wb97"):
        d = {"wb97mv_tzvpd": "o_wb97mv_fix", "wb97md3_tzvppd": "o_wb97md3bj_fix"}[level]
        p = B3 / d / mol / f"{mol}.hess"
        if not p.exists():
            return None, None
        apt, x = orca_hess(p)
        # ORCA's analytic Freq translates the molecule (centre of mass); it does not rotate it.
        # The APT of a neutral molecule is translation invariant, so compare centred coordinates.
        return apt, float(np.abs((x - x.mean(0)) - (f0["x"] - f0["x"].mean(0))).max())
    d = {"b3lyp_avtz": "g_b3lyp_avtz_fix", "pbe0_tzvp": "g_pbe0_tzvp_fix"}[level]
    p = B3 / d / mol / f"{mol}.fchk"
    if not p.exists():
        return None, None
    n = len(f0["m"])
    apt = fchk_any(p, "Dipole Derivatives").reshape(3 * n, 3)
    x = fchk_any(p, "Current cartesian coordinates").reshape(-1, 3)
    return apt, float(np.abs(x - f0["x"]).max())


def own_geometry_avtz(mol, f0, freqs0, disp0):
    """B3LYP/aug-cc-pVTZ own-minimum harmonic intensities paired to the 6-31G(d,p) modes."""
    p = B3 / "g_b3lyp_avtz_opt" / mol / f"{mol}.fchk"
    if not p.exists():
        return None
    n = len(f0["m"])
    x = fchk_any(p, "Current cartesian coordinates").reshape(-1, 3)
    h = np.zeros((3 * n, 3 * n))
    h[np.tril_indices(3 * n)] = fchk_any(p, "Cartesian Force Constants")
    h = h + h.T - np.diag(np.diag(h))
    apt = fchk_any(p, "Dipole Derivatives").reshape(3 * n, 3)
    fr, disp = normal_modes(x, f0["m"], h)
    inten = KMMOL * np.sum((apt.T @ disp) ** 2, 0)
    # pair modes: Cartesian mode vectors (n_modes, n, 3), aligned onto the 6-31G(d,p) frame
    mc = (disp / np.linalg.norm(disp, axis=0)).T.reshape(-1, n, 3)
    mr = (disp0 / np.linalg.norm(disp0, axis=0)).T.reshape(-1, n, 3)
    mc, _rmsd = align_modes_to_reference(mc, x, f0["x"])
    m = match_modes(mc, mr)
    out = np.full(len(freqs0), np.nan)
    ov = np.full(len(freqs0), np.nan)
    for i_c, (i_r, o) in m.items():
        if i_r is not None:
            out[i_r] = inten[i_c]
            ov[i_r] = o
    return out, ov


def main():
    ml = {m: np.load(HERE / f"b1_apt_{'mace_polar1' if m == 'polar1_field' else m}.npz") for m in ML}
    rows, apt_rows, mism = [], [], []
    for mol in MOLS:
        f0 = fchk(mol)
        freqs, disp = normal_modes(f0["x"], f0["m"], f0["hess"])
        ints, apts = {}, {}
        for lev in LEVELS:
            a, dx = level_apt(lev, mol, f0)
            if a is None:
                continue
            mism.append((lev, mol, dx))
            apts[lev] = a
            ints[lev] = KMMOL * np.sum((a.T @ disp) ** 2, 0)
        for m in ML:
            if m == "polar1_field":
                a = -ml["mace_polar1"][f"{mol}/apt_field_0.001"]
            else:
                a = ml[m][f"{mol}/apt"]
            apts[m] = a
            ints[m] = KMMOL * np.sum((a.T @ disp) ** 2, 0)
        own = own_geometry_avtz(mol, f0, freqs, disp)
        if own is not None:
            ints["b3lyp_avtz_own"] = own[0]
        keys = list(ints)
        fm, merged = merge_degenerate(freqs, *[np.nan_to_num(ints[k], nan=-1.0) for k in keys])
        for k in range(len(fm)):
            rows.append(dict(mol=mol, freq=fm[k], **{key: merged[j][k] for j, key in enumerate(keys)}))
        for a in apts:
            for b in apts:
                apt_rows.append(dict(mol=mol, a=a, b=b, rel=float(np.linalg.norm(apts[a] - apts[b]) / np.linalg.norm(apts[b]))))

    out = ["# B3 tables (generated by b3_analyse.py)\n"]
    w = out.append
    w("Frame check, max |x_level - x_B3LYP/6-31G(d,p)| after centring (Bohr): " + ", ".join(
        f"{LEVELS[lev]} {max(d for lv, _m, d in mism if lv == lev):.1e}" for lev in LEVELS
        if any(lv == lev for lv, _m, _d in mism)) + "\n")  # fmt: skip
    have = [lev for lev in LEVELS if all(lev in r for r in rows)]
    missing = [lev for lev in LEVELS if lev not in have]
    if missing:
        w(f"**Missing (not all molecules finished): {', '.join(missing)}**\n")

    sel = [r for r in rows if r["b3lyp_631gdp"] > 1]
    # 1. how far the reference moves
    w("\n## 1. How far the reference moves (bands with B3LYP/6-31G(d,p) > 1 km/mol)\n")
    w("Ratio I(level) / I(B3LYP/6-31G(d,p)), median over the bands of each region, same "
      "geometry and modes; last column: B3LYP/aug-cc-pVTZ at its own minimum.\n")  # fmt: skip
    cols = [c for c in have if c != "b3lyp_631gdp"] + (["b3lyp_avtz_own"] if all("b3lyp_avtz_own" in r for r in rows) else [])
    w("| region | bands | " + " | ".join(LEVELS.get(c, "B3LYP/aug-cc-pVTZ own geometry") for c in cols) + " |")
    w("|---" * (len(cols) + 2) + "|")
    for g in ("X-H >2700", "1700-2700", "500-1700", "<500", "all"):
        b = [r for r in sel if g == "all" or region(r["freq"]) == g]
        if not b:
            continue
        cells = []
        for c in cols:
            v = np.array([r[c] / r["b3lyp_631gdp"] for r in b if r[c] >= 0])
            cells.append(f"{np.median(v):.2f} ({v.min():.2f}–{v.max():.2f})")
        w(f"| {g} | {len(b)} | " + " | ".join(cells) + " |")
    w("\nPer band, 1700-2700 cm-1 (km/mol):\n")
    w("| molecule | freq | " + " | ".join(LEVELS[c] for c in have) + " | " + " | ".join(ML.values()) + " |")
    w("|---" * (2 + len(have) + len(ML)) + "|")
    for r in sel:
        if 1700 < r["freq"] <= 2700:
            w(f"| {r['mol']} | {r['freq']:.0f} | " + " | ".join(f"{r[c]:.0f}" for c in have)
              + " | " + " | ".join(f"{r[m]:.0f}" for m in ML) + " |")  # fmt: skip

    # 2. ML models against each reference level
    def mae_table(bs, title):
        w(f"\n### {title}\n")
        w("| dipole model | " + " | ".join(LEVELS[c] for c in have) + " |")
        w("|---" * (len(have) + 1) + "|")
        for m in ML:
            cells = []
            vals = {}
            for c in have:
                d = np.array([abs(r[m] - r[c]) for r in bs])
                lr = np.array([abs(np.log10(max(r[m], 1e-6) / r[c])) for r in bs])
                vals[c] = (d.mean(), np.median(lr))
            for c in have:
                cells.append(f"{vals[c][0]:.1f} / {vals[c][1]:.3f}")
            w(f"| {ML[m]} | " + " | ".join(cells) + " |")
        w("\nBest charge-dipole model (POLAR-1 / MDP / MACE4IR) per reference, MAE | median log: " + "; ".join(
            f"{LEVELS[c]}: {ML[min(['mace_polar1', 'mace_mdp', 'mace_ml'], key=lambda m: np.mean([abs(r[m] - r[c]) for r in bs]))]} | "
            f"{ML[min(['mace_polar1', 'mace_mdp', 'mace_ml'], key=lambda m: np.median([abs(np.log10(max(r[m], 1e-6) / r[c])) for r in bs]))]}"
            for c in have))  # fmt: skip

    w("\n## 2. Each dipole model against each reference level: MAE (km/mol) / median abs(log10 ratio)\n")
    w("Bands: B3LYP/6-31G(d,p) harmonic > 1 km/mol, on the 6-31G(d,p) modes at its geometry.\n")
    mae_table(sel, f"all bands ({len(sel)}, {len(MOLS)} molecules)")
    mae_table([r for r in sel if not 1700 < r["freq"] <= 2700], "without the 1700-2700 cm-1 bands")
    mae_table([r for r in sel if region(r["freq"]) == "X-H >2700"], "X-H stretches only")

    # 3. element-wise APT
    w("\n## 3. Element-wise APT: relative Frobenius difference, mean over molecules\n")
    w("| | " + " | ".join(LEVELS[c] for c in have) + " |")
    w("|---" * (len(have) + 1) + "|")
    for a in list(ML) + have:
        cells = []
        for b in have:
            v = [r["rel"] for r in apt_rows if r["a"] == a and r["b"] == b]
            cells.append(f"{np.mean(v):.3f}")
        w(f"| {ML.get(a, LEVELS.get(a))} | " + " | ".join(cells) + " |")
    (HERE / "b3_tables.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
