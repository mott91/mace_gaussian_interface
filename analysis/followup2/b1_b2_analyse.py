"""Follow-up 2, B1 + B2 (analysis): dipole models alone, at the B3LYP geometry and modes.

Inputs: b1_apt_<model>.npz (b1_b2_compute.py) and the campaign B3LYP fchk files.
- B3LYP normal modes are rebuilt from the fchk "Cartesian Force Constants" and "Real atomic
  weights" (mass-weighted, translations/rotations projected out). Unit test: projecting the
  B3LYP APT ("Dipole Derivatives") on these modes must reproduce Gaussian's own harmonic
  intensities from the log.
- Every model's APT is projected on the same B3LYP modes: intensities then differ only by the
  dipole model. Band selection as in A1 (B3LYP harmonic > cutoff, degenerate levels within
  0.5 cm-1 summed, ML not filtered).
- Element-wise APT comparison (no modes): RMS of APT_model - APT_B3LYP over all 9N elements.
- mu0 at the same geometry: no frame or geometry effect.
- B2: POLAR-1 field response. The model adds +E.mu to its energy, so the physical dipole is
  +dE/dE_field: mu_field = -(stored), APT_field = -(stored) (sign fixed here, see report).

Writes b1_tables.md next to this script.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import logging  # noqa: E402

logging.disable(logging.CRITICAL)
from mace_gaussian.gaussian.parser import parse_gaussian_log  # noqa: E402

CAMP = ROOT / "campaigns/2026"
AU2CM = 5140.48706  # sqrt(Eh / (a0^2 amu)) in cm-1
KMMOL = 974.8801  # km/mol per (e^2 / amu)
AU2D = 2.541746
SYMTOP = {"ammonia", "phosphine"}
MODELS = ["mace_polar1", "polar1_field", "mace_mdp", "mace_ml"]
LAB = {"mace_polar1": "POLAR-1", "polar1_field": "POLAR-1 (field)", "mace_mdp": "MDP",
       "mace_ml": "MACE4IR"}  # fmt: skip
FAMILY = {
    **dict.fromkeys(["formic_acid_trans", "formic_acid_cis", "acetic_acid", "propionic_acid",
                     "butyric_acid", "valeric_acid", "caproic_acid"], "acids"),
    **dict.fromkeys(["methanol", "ethanol", "1-propanol", "isopropanol", "1-butanol",
                     "tert-butanol"], "alcohols"),
    **dict.fromkeys(["carbon_dioxide", "sulfur_dioxide", "ammonia", "phosphine",
                     "hydrogen_peroxide", "nitric_acid", "carbonic_acid"], "inorganic"),
    **dict.fromkeys(["glycine", "beta_glucose"], "glycine/glucose"),
}  # fmt: skip
RNG = np.random.default_rng(20261105)


def fchk(mol):
    t = (CAMP / mol / f"{mol}_freq_anharm.fchk").read_text()

    def arr(key):
        i = t.index(key + " ")
        n = int(t[i : t.index("\n", i)].split("N=")[1])
        return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])

    n = len(arr("Atomic numbers"))
    fc = arr("Cartesian Force Constants")
    h = np.zeros((3 * n, 3 * n))
    h[np.tril_indices(3 * n)] = fc
    h = h + h.T - np.diag(np.diag(h))
    return dict(x=arr("Current cartesian coordinates").reshape(-1, 3), m=arr("Real atomic weights"),
                hess=h, apt=arr("Dipole Derivatives").reshape(3 * n, 3), mu=arr("Dipole Moment"))  # fmt: skip


def normal_modes(x, m, h):
    """Mass-weighted normal modes with translations/rotations projected out.
    Returns freqs (cm-1) and the Cartesian displacement per unit normal coordinate
    (3N, nmodes), i.e. L / sqrt(m)."""
    n = len(m)
    sm = np.repeat(np.sqrt(m), 3)
    hw = h / np.outer(sm, sm)
    com = (m[:, None] * x).sum(0) / m.sum()
    r = x - com
    vecs = []
    for k in range(3):
        v = np.zeros((n, 3))
        v[:, k] = 1
        vecs.append((v * np.sqrt(m)[:, None]).ravel())
    for k in range(3):
        e = np.zeros(3)
        e[k] = 1
        vecs.append((np.cross(e, r) * np.sqrt(m)[:, None]).ravel())
    v = np.array(vecs).T
    u, s, _ = np.linalg.svd(v, full_matrices=False)
    u = u[:, s > 1e-6 * s.max()]
    p = np.eye(3 * n) - u @ u.T
    w, lv = np.linalg.eigh(p @ hw @ p)
    keep = np.argsort(np.abs(w))[u.shape[1]:]  # drop the TR null space
    keep = keep[np.argsort(w[keep])]
    freqs = np.sign(w[keep]) * np.sqrt(np.abs(w[keep])) * AU2CM
    return freqs, lv[:, keep] / sm[:, None]


def intensities(apt, disp):
    dmu = apt.T @ disp  # (3, nmodes), e / sqrt(amu)
    return KMMOL * np.sum(dmu**2, axis=0)


def merge_degenerate(freqs, *ints):
    order = np.argsort(freqs)
    groups, cur = [], [order[0]]
    for i in order[1:]:
        if freqs[i] - freqs[cur[-1]] < 0.5:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    f = np.array([freqs[g].mean() for g in groups])
    return f, [np.array([iv[g].sum() for g in groups]) for iv in ints]


def region(f):
    return ("X-H >2700" if f > 2700 else "1700-2700" if f > 1700 else "500-1700" if f >= 500
            else "<500")  # fmt: skip


def main():
    data = {m: np.load(HERE / f"b1_apt_{'mace_polar1' if m == 'polar1_field' else m}.npz")
            for m in MODELS}  # fmt: skip
    mols = sorted({k.split("/")[0] for k in data["mace_ml"].files})
    bands, apt_err, mu_rows, unit_test, fr = [], [], [], [], []
    for mol in mols:
        f = fchk(mol)
        freqs, disp = normal_modes(f["x"], f["m"], f["hess"])
        ib = intensities(f["apt"], disp)
        # unit test against Gaussian's log
        g = parse_gaussian_log(str(CAMP / mol / f"{mol}_freq_anharm.log"))["harmonic"]
        gf = np.array([r["freq_cm"] for r in g])
        gi = np.array([r["ir_intensity"] for r in g])
        o = np.argsort(freqs)
        unit_test.append((mol, np.max(np.abs(np.sort(gf) - freqs[o])),
                          np.max(np.abs(gi[np.argsort(gf)] - ib[o]) / np.maximum(gi[np.argsort(gf)], 1))))  # fmt: skip
        apts = {}
        for m in MODELS:
            if m == "polar1_field":
                apts[m] = -data["mace_polar1"][f"{mol}/apt_field_0.001"]
                mu = -data["mace_polar1"][f"{mol}/mu_field_0.001"] / 0.529177210903  # e*A -> e*Bohr
            else:
                apts[m] = data[m][f"{mol}/apt"]
                mu = data[m][f"{mol}/mu"]
            nb = np.linalg.norm(f["mu"])
            ang = (np.degrees(np.arccos(np.clip(mu @ f["mu"] / (np.linalg.norm(mu) * nb), -1, 1)))
                   if nb * AU2D > 0.1 else np.nan)  # fmt: skip
            mu_rows.append(dict(mol=mol, model=m, err=(np.linalg.norm(mu) - nb) * AU2D, angle=ang,
                                vec=np.linalg.norm(mu - f["mu"]) * AU2D))  # fmt: skip
            apt_err.append(dict(mol=mol, model=m, rms=float(np.sqrt(np.mean((apts[m] - f["apt"]) ** 2))),
                                rel=float(np.linalg.norm(apts[m] - f["apt"]) / np.linalg.norm(f["apt"]))))  # fmt: skip
        ims = [intensities(apts[m], disp) for m in MODELS]
        fm, merged = merge_degenerate(freqs, ib, *ims)
        for k in range(len(fm)):
            bands.append(dict(mol=mol, freq=fm[k], dft=merged[0][k],
                              **{m: merged[1 + j][k] for j, m in enumerate(MODELS)}))  # fmt: skip
        # B2 detail: field vs charge, all field steps
        for d in (1e-4, 1e-3, 1e-2):
            a_f = -data["mace_polar1"][f"{mol}/apt_field_{d:g}"]
            a_c = data["mace_polar1"][f"{mol}/apt"]
            mu_f = -data["mace_polar1"][f"{mol}/mu_field_{d:g}"]
            mu_c = data["mace_polar1"][f"{mol}/mu_charge_eA"]
            fr.append(dict(mol=mol, d=d, apt_rel=np.linalg.norm(a_f - a_c) / np.linalg.norm(a_c),
                           mu_ratio=np.linalg.norm(mu_f) / max(np.linalg.norm(mu_c), 1e-9),
                           mu_ang=np.degrees(np.arccos(np.clip(mu_f @ mu_c / (np.linalg.norm(mu_f) * np.linalg.norm(mu_c) + 1e-12), -1, 1)))))  # fmt: skip

    out = ["# B1 + B2 tables (generated by b1_b2_analyse.py)\n"]
    w = out.append
    ut = np.array([(u[1], u[2]) for u in unit_test])
    w(f"Unit test (22 molecules): B3LYP frequencies from the fchk Hessian agree with Gaussian's "
      f"log to {ut[:, 0].max():.3f} cm-1; intensities from APT x modes to "
      f"{100 * ut[:, 1].max():.3f} % (relative to max(I, 1 km/mol)).\n")  # fmt: skip

    def sel(cut, excl, extra=lambda b: True):
        return [b for b in bands if b["dft"] > cut and b["mol"] not in excl and extra(b)]

    def met(bs, m):
        d = np.array([b["dft"] for b in bs])
        x = np.array([b[m] for b in bs])
        lr = np.abs(np.log10(np.maximum(x, 1e-6) / d))
        return np.mean(np.abs(x - d)), np.median(lr), 100 * np.mean(lr < np.log10(2)), len(bs)

    def row(lab, bs):
        v = {m: met(bs, m) for m in MODELS}
        best_mae = min(MODELS[0:1] + MODELS[2:], key=lambda m: v[m][0])
        best_lr = min(MODELS[0:1] + MODELS[2:], key=lambda m: v[m][1])
        cells = []
        for m in MODELS:
            a = f"{v[m][0]:.1f}"
            b = f"{v[m][1]:.3f}"
            cells.append(f"{'**' + a + '**' if m == best_mae else a} / {'**' + b + '**' if m == best_lr else b} / {v[m][2]:.0f} %")
        return f"| {lab} | " + " | ".join(cells) + f" | {v[MODELS[0]][3]} |"

    hdr = "| set | " + " | ".join(LAB[m] for m in MODELS) + " | bands |"
    w("\n## Intensities on B3LYP modes: MAE (km/mol) / median abs(log10 ratio) / within x2\n")
    w("Bold: best of the three charge-dipole models (POLAR-1, MDP, MACE4IR).\n")
    w(hdr)
    w("|---" * 6 + "|")
    for cut in (0.1, 1.0, 10.0):
        w(row(f"headline, > {cut:g} km/mol", sel(cut, SYMTOP)))
    w(row("all 22 molecules, > 1", sel(1.0, set())))
    w(row("headline > 1, without 1700-2700", sel(1.0, SYMTOP, lambda b: not 1700 < b["freq"] <= 2700)))
    for g in ("X-H >2700", "1700-2700", "500-1700", "<500"):
        w(row(f"region {g}", sel(1.0, SYMTOP, lambda b, g=g: region(b["freq"]) == g)))
    for fam in ("acids", "alcohols", "glycine/glucose", "inorganic"):
        w(row(f"family {fam}", sel(1.0, SYMTOP, lambda b, f=fam: FAMILY[b["mol"]] == f)))

    # per molecule best + bootstrap (three charge-dipole models)
    three = ["mace_polar1", "mace_mdp", "mace_ml"]
    hb = sel(1.0, SYMTOP)
    hm = sorted({b["mol"] for b in hb})
    per = {mol: {m: np.mean([abs(b[m] - b["dft"]) for b in hb if b["mol"] == mol]) for m in three}
           for mol in hm}  # fmt: skip
    cnt = {m: sum(min(three, key=lambda k: per[mol][k]) == m for mol in hm) for m in three}
    w(f"\nBest per molecule (headline, > 1 km/mol, {len(hm)} molecules): "
      + ", ".join(f"{LAB[m]} {cnt[m]}" for m in three))  # fmt: skip
    acc = {mol: {m: (sum(abs(b[m] - b["dft"]) for b in hb if b["mol"] == mol),
                     sum(1 for b in hb if b["mol"] == mol)) for m in three} for mol in hm}  # fmt: skip
    diffs = {"mace_mdp": [], "mace_ml": []}
    first = 0
    for _ in range(20000):
        dr = RNG.choice(hm, size=len(hm))
        mae = {m: sum(acc[x][m][0] for x in dr) / sum(acc[x][m][1] for x in dr) for m in three}
        for o in diffs:
            diffs[o].append(mae["mace_polar1"] - mae[o])
        first += min(three, key=lambda m: mae[m]) == "mace_polar1"
    for o, v in diffs.items():
        v = np.array(v)
        w(f"- bootstrap POLAR-1 − {LAB[o]}: {v.mean():+.2f} [{np.percentile(v, 2.5):+.2f}, "
          f"{np.percentile(v, 97.5):+.2f}] km/mol")  # fmt: skip
    w(f"- P(POLAR-1 first): {100 * first / 20000:.1f} %")

    # 1700-2700 bands
    w("\n## The 1700-2700 cm-1 bands on B3LYP modes (km/mol)\n")
    w("| molecule | freq | B3LYP | " + " | ".join(LAB[m] for m in MODELS) + " |")
    w("|---" * 7 + "|")
    for b in sel(1.0, SYMTOP, lambda b: 1700 < b["freq"] <= 2700):
        w(f"| {b['mol']} | {b['freq']:.0f} | {b['dft']:.0f} | " + " | ".join(f"{b[m]:.0f}" for m in MODELS) + " |")

    # element-wise APT
    w("\n## Element-wise APT against B3LYP (no modes): RMS (e) / relative Frobenius error\n")
    w("| set | " + " | ".join(LAB[m] for m in MODELS) + " |")
    w("|---" * 5 + "|")
    for lab, ex in (("headline (20)", SYMTOP), ("all 22", set())):
        cells = []
        for m in MODELS:
            r = [a for a in apt_err if a["model"] == m and a["mol"] not in ex]
            cells.append(f"{np.mean([a['rms'] for a in r]):.4f} / {np.mean([a['rel'] for a in r]):.3f}")
        w(f"| {lab} | " + " | ".join(cells) + " |")
    bestapt = {m: 0 for m in three}
    for mol in mols:
        if mol in SYMTOP:
            continue
        r = {a["model"]: a["rms"] for a in apt_err if a["mol"] == mol}
        bestapt[min(three, key=lambda m: r[m])] += 1
    w("\nBest APT per molecule (headline): " + ", ".join(f"{LAB[m]} {bestapt[m]}" for m in three))

    # mu0 at the B3LYP geometry
    w("\n## mu0 at the B3LYP geometry (no geometry or frame effect), headline\n")
    w("| metric | " + " | ".join(LAB[m] for m in MODELS) + " |")
    w("|---" * 5 + "|")
    for met_, fn, lab in (("err", lambda v: np.mean(np.abs(v)), "MAE of abs(mu) (D)"),
                          ("err", np.mean, "mean signed error (D)"),
                          ("vec", np.mean, "mean vector error (D)"),
                          ("angle", np.nanmedian, "median angle (deg)")):  # fmt: skip
        w(f"| {lab} | " + " | ".join(
            f"{fn(np.array([r[met_] for r in mu_rows if r['model'] == m and r['mol'] not in SYMTOP])):.3f}"
            for m in MODELS) + " |")  # fmt: skip

    # B2 detail
    w("\n## B2: POLAR-1 field response vs charge distribution (all 22 molecules)\n")
    w("| field step (V/A) | median abs(mu_field)/abs(mu_charge) | range | median angle (deg) | "
      "median relative APT difference | max |")  # fmt: skip
    w("|---|---|---|---|---|---|")
    for d in (1e-4, 1e-3, 1e-2):
        r = [x for x in fr if x["d"] == d and x["mol"] != "carbon_dioxide"]
        q = np.array([x["mu_ratio"] for x in r])
        a = np.array([x["apt_rel"] for x in fr if x["d"] == d])
        w(f"| {d:g} | {np.median(q):.3f} | {q.min():.3f}–{q.max():.3f} | "
          f"{np.median([x['mu_ang'] for x in r]):.2f} | {np.median(a):.3f} | {a.max():.3f} |")  # fmt: skip
    (HERE / "b1_tables.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))
    print("\nunit test per molecule:", [(u[0], round(u[1], 3), round(100 * u[2], 3)) for u in unit_test])


if __name__ == "__main__":
    main()
