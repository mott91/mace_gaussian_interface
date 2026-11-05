"""Follow-up 2, A1: robustness of the dipole-model ranking on harmonic IR intensities.

Data (read only): campaigns/2026/analysis_results_harmonic/<mol>/data/comparison_<E>_<D>.csv,
written by run_analysis_harmonic.py. Each row pairs one B3LYP harmonic mode with one ML harmonic
mode of run <E>_<D>: Hungarian assignment on the overlap of the Kabsch-aligned Cartesian
normal-mode vectors (mace_gaussian/analysis/mode_matching.py). Intensities are the harmonic
(double-harmonic) IR intensities in km/mol, B3LYP from its own log, ML from the run's log.

Band selection, everything below:
- one band = one B3LYP harmonic fundamental with its matched ML mode;
- degenerate B3LYP levels (harmonic frequencies within 0.5 cm-1 of each other) are merged into one
  band, intensities summed on both sides (only CO2's bend, NH3's and PH3's e modes);
- a band is kept if its B3LYP harmonic intensity exceeds the cutoff (0.1, 1 or 10 km/mol);
  the ML intensity is not filtered;
- molecule sets: "headline" = without NH3 and PH3; "all" = every molecule the run exists for;
  ANI-cc has no SO2 and PH3 runs (no S/P parameters); "common" = molecules all 12 runs share.

Outputs: a1_tables.md (all tables) and a1_bands.csv (the selected bands) next to this script.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, wilcoxon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HARM = ROOT / "campaigns/2026/analysis_results_harmonic"
ENERGY = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]
DIPOLE = ["mace_polar1", "mace_mdp", "mace_ml"]
ELAB = {"mace_anicc": "ANI-cc", "mace_omol": "OMOL", "mace_off": "OFF", "mace_polar": "POLAR"}
DLAB = {"mace_polar1": "POLAR-1", "mace_mdp": "MDP", "mace_ml": "MACE4IR"}
FAMILY = {
    **dict.fromkeys(
        ["formic_acid_trans", "formic_acid_cis", "acetic_acid", "propionic_acid",
         "butyric_acid", "valeric_acid", "caproic_acid"], "acids"),  # fmt: skip
    **dict.fromkeys(
        ["methanol", "ethanol", "1-propanol", "isopropanol", "1-butanol", "tert-butanol"],
        "alcohols"),  # fmt: skip
    **dict.fromkeys(
        ["carbon_dioxide", "sulfur_dioxide", "ammonia", "phosphine", "hydrogen_peroxide",
         "nitric_acid", "carbonic_acid"], "inorganic"),  # fmt: skip
    **dict.fromkeys(["glycine", "beta_glucose"], "glycine/glucose"),
}
SYMTOP = {"ammonia", "phosphine"}
RNG = np.random.default_rng(20261105)
NBOOT = 20000


def region(f: float) -> str:
    if f > 2700:
        return "X-H stretch >2700"
    if f > 1700:
        return "1700-2700 (C=O, NO2, CO2 asym)"
    if f >= 500:
        return "500-1700"
    return "<500"


def load():
    """bands[(mol, E, D)] = list of dicts (one per band, degenerate levels merged)."""
    out = {}
    merged = 0
    for mol_dir in sorted(HARM.iterdir()):
        mol = mol_dir.name
        for E in ENERGY:
            for D in DIPOLE:
                f = mol_dir / "data" / f"comparison_{E}_{D}.csv"
                if not f.exists():
                    continue
                rows = sorted(
                    (
                        (float(r["DFT_Frequency_cm"]), float(r["DFT_Intensity_km_mol"]),
                         float(r["ML_Intensity_km_mol"]), float(r["Mode_Overlap"]),
                         float(r["ML_Frequency_cm"]))
                        for r in csv.DictReader(f.open())
                    )
                )  # fmt: skip
                groups, cur = [], [rows[0]]
                for r in rows[1:]:
                    if r[0] - cur[-1][0] < 0.5:
                        cur.append(r)
                    else:
                        groups.append(cur)
                        cur = [r]
                groups.append(cur)
                merged += sum(len(g) > 1 for g in groups)
                out[(mol, E, D)] = [
                    dict(mol=mol, E=E, D=D, freq=np.mean([g[0] for g in grp]),
                         dft=sum(g[1] for g in grp), ml=sum(g[2] for g in grp),
                         overlap=min(g[3] for g in grp), ml_freq=np.mean([g[4] for g in grp]),
                         ndeg=len(grp))
                    for grp in groups
                ]  # fmt: skip
    return out, merged


def metrics(b: list[dict]) -> dict:
    if not b:
        return dict(n=0, mae=np.nan, medlog=np.nan, f2=np.nan)
    d = np.array([x["dft"] for x in b])
    m = np.array([x["ml"] for x in b])
    lr = np.abs(np.log10(np.maximum(m, 1e-6) / d))
    return dict(n=len(b), mae=float(np.mean(np.abs(m - d))), medlog=float(np.median(lr)),
                f2=float(np.mean(lr < np.log10(2))))  # fmt: skip


def select(bands, mols, E, D, cut):
    return [x for mol in mols for x in bands.get((mol, E, D), []) if x["dft"] > cut]


def mols_for(bands, E, excl):
    have = {m for (m, e, d) in bands if e == E}
    have = [m for m in sorted(have) if all((m, E, D) in bands for D in DIPOLE)]
    return [m for m in have if m not in excl]


def table(bands, cut, excl, common=None, metric="mae", fmt="{:.1f}"):
    lines = ["| energy model | " + " | ".join(DLAB[D] for D in DIPOLE) + " | bands / molecules |",
             "|---" * 5 + "|"]  # fmt: skip
    for E in ENERGY:
        mols = common if common is not None else mols_for(bands, E, excl)
        vals = {D: metrics(select(bands, mols, E, D, cut)) for D in DIPOLE}
        best = (min if metric in ("mae", "medlog") else max)(DIPOLE, key=lambda D: vals[D][metric])
        cells = []
        for D in DIPOLE:
            s = fmt.format(vals[D][metric] * (100 if metric == "f2" else 1))
            cells.append(f"**{s}**" if D == best else s)
        lines.append(f"| {ELAB[E]} | " + " | ".join(cells)
                     + f" | {vals[DIPOLE[0]]['n']} / {len(mols)} |")  # fmt: skip
    return "\n".join(lines)


def per_mol_mae(bands, mols, E, cut):
    """{mol: {D: MAE}} for molecules with at least one band above cut."""
    out = {}
    for mol in mols:
        r = {D: metrics(select(bands, [mol], E, D, cut)) for D in DIPOLE}
        if r[DIPOLE[0]]["n"]:
            out[mol] = {D: r[D]["mae"] for D in DIPOLE}
    return out


def pooled_mae_from(bandlist_by_mol, mols_draw, D):
    num = sum(bandlist_by_mol[m][D][0] for m in mols_draw)
    den = sum(bandlist_by_mol[m][D][1] for m in mols_draw)
    return num / den


def bootstrap(bands, Es, mols, cut):
    """Resample molecules (with replacement); statistic = band-pooled MAE per dipole model,
    averaged over the energy models in Es (all on the same molecule draw)."""
    acc = {(m, E, D): (0.0, 0) for m in mols for E in Es for D in DIPOLE}
    for m in mols:
        for E in Es:
            for D in DIPOLE:
                b = select(bands, [m], E, D, cut)
                acc[(m, E, D)] = (float(sum(abs(x["ml"] - x["dft"]) for x in b)), len(b))
    mols = np.array(mols)
    diffs = {("mace_polar1", "mace_mdp"): [], ("mace_polar1", "mace_ml"): []}
    first = 0
    for _ in range(NBOOT):
        draw = RNG.choice(mols, size=len(mols), replace=True)
        mae = {}
        for D in DIPOLE:
            v = []
            for E in Es:
                num = sum(acc[(m, E, D)][0] for m in draw)
                den = sum(acc[(m, E, D)][1] for m in draw)
                v.append(num / den)
            mae[D] = float(np.mean(v))
        for a, b in diffs:
            diffs[(a, b)].append(mae[a] - mae[b])
        first += min(DIPOLE, key=lambda D: mae[D]) == "mace_polar1"
    res = {}
    for k, v in diffs.items():
        v = np.array(v)
        res[k] = (float(np.mean(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    return res, first / NBOOT


def main():
    bands, merged = load()
    out = []
    w = out.append
    w("# A1 tables (generated by a1_harmonic_intensities.py)\n")
    w(f"Degenerate B3LYP levels merged: {merged} (summed over all 255 runs).\n")
    common = sorted(set.intersection(*[set(mols_for(bands, E, SYMTOP)) for E in ENERGY]))
    w(f"Common molecule set (all 12 runs, without NH3/PH3): {len(common)} molecules.\n")

    for cut in (1.0, 0.1, 10.0):
        w(f"\n## Cutoff: B3LYP harmonic intensity > {cut:g} km/mol\n")
        for name, excl in (("headline, without NH3/PH3", SYMTOP), ("all molecules", set())):
            w(f"\n### MAE (km/mol), {name}\n")
            w(table(bands, cut, excl))
        w("\n### MAE (km/mol), common set (same molecules in every row)\n")
        w(table(bands, cut, SYMTOP, common=common))
        w("\n### median |log10(I_ML / I_B3LYP)|, headline\n")
        w(table(bands, cut, SYMTOP, metric="medlog", fmt="{:.3f}"))
        w("\n### % of bands within a factor 2, headline\n")
        w(table(bands, cut, SYMTOP, metric="f2", fmt="{:.0f}"))

    cut = 1.0
    # best per molecule (headline set, per energy model)
    w("\n## Per molecule, cutoff 1 km/mol, headline\n")
    w("\nBest dipole model per (molecule, energy model), counts:\n")
    w("| energy model | " + " | ".join(DLAB[D] for D in DIPOLE) + " | molecules |")
    w("|---" * 5 + "|")
    allcount = defaultdict(int)
    pm_all = {}
    for E in ENERGY:
        pm = per_mol_mae(bands, mols_for(bands, E, SYMTOP), E, cut)
        pm_all[E] = pm
        c = defaultdict(int)
        for mol, v in pm.items():
            c[min(DIPOLE, key=lambda D: v[D])] += 1
        for D in DIPOLE:
            allcount[D] += c[D]
        w(f"| {ELAB[E]} | " + " | ".join(str(c[D]) for D in DIPOLE) + f" | {len(pm)} |")
    w("| all four | " + " | ".join(str(allcount[D]) for D in DIPOLE) + " | |")

    # paired tests over molecules
    w("\n## Paired tests over molecules (per-molecule MAE, cutoff 1, headline)\n")
    w("| energy model | pair | median diff | POLAR-1 better in | Wilcoxon p | sign-test p |")
    w("|---|---|---|---|---|---|")
    for E in ENERGY:
        pm = pm_all[E]
        for other in ("mace_mdp", "mace_ml"):
            d = np.array([v["mace_polar1"] - v[other] for v in pm.values()])
            k = int(np.sum(d < 0))
            nz = int(np.sum(d != 0))
            pw = wilcoxon(d).pvalue if nz else np.nan
            ps = binomtest(k, nz, 0.5).pvalue if nz else np.nan
            w(f"| {ELAB[E]} | POLAR-1 − {DLAB[other]} | {np.median(d):+.2f} | {k} of {len(d)} | "
              f"{pw:.3f} | {ps:.3f} |")  # fmt: skip

    # bootstrap
    w(f"\n## Bootstrap over molecules ({NBOOT} draws, cutoff 1, headline)\n")
    w("Statistic: band-pooled MAE of each dipole model; for 'all four' the mean over the four "
      "energy models on the same draw of the common molecule set.\n")  # fmt: skip
    w("| energy model(s) | molecules | POLAR-1 − MDP [95 %] | POLAR-1 − MACE4IR [95 %] | "
      "P(POLAR-1 first) |")  # fmt: skip
    w("|---|---|---|---|---|")
    for Es, mols, lab in (
        *[([E], [m for m in pm_all[E]], ELAB[E]) for E in ENERGY],
        (ENERGY, [m for m in common if all(m in pm_all[E] for E in ENERGY)], "all four"),
    ):
        res, pf = bootstrap(bands, Es, mols, cut)
        a, b = res[("mace_polar1", "mace_mdp")], res[("mace_polar1", "mace_ml")]
        w(f"| {lab} | {len(mols)} | {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}] | "
          f"{b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] | {100 * pf:.1f} % |")  # fmt: skip

    # splits
    for key, fn, lab in (("region", lambda x: region(x["freq"]), "band region (B3LYP harmonic)"),
                         ("family", lambda x: FAMILY[x["mol"]], "molecule family")):  # fmt: skip
        w(f"\n## Split by {lab}, cutoff 1, headline: MAE (km/mol) / bands\n")
        groups = sorted({fn(x) for v in bands.values() for x in v})
        w("| energy model | group | " + " | ".join(DLAB[D] for D in DIPOLE) + " | bands |")
        w("|---" * 6 + "|")
        for E in ENERGY:
            mols = mols_for(bands, E, SYMTOP)
            for g in groups:
                vals = {D: metrics([x for x in select(bands, mols, E, D, cut) if fn(x) == g])
                        for D in DIPOLE}  # fmt: skip
                if not vals[DIPOLE[0]]["n"]:
                    continue
                best = min(DIPOLE, key=lambda D: vals[D]["mae"])
                cells = [(f"**{vals[D]['mae']:.1f}**" if D == best else f"{vals[D]['mae']:.1f}")
                         for D in DIPOLE]  # fmt: skip
                w(f"| {ELAB[E]} | {g} | " + " | ".join(cells) + f" | {vals[DIPOLE[0]]['n']} |")

    # concentration: how much of the POLAR-1 lead comes from the largest band differences
    w("\n## Concentration of the lead (cutoff 1, headline, all energy models pooled)\n")
    for other in ("mace_mdp", "mace_ml"):
        gains = []
        for E in ENERGY:
            mols = mols_for(bands, E, SYMTOP)
            p = select(bands, mols, E, "mace_polar1", cut)
            o = select(bands, mols, E, other, cut)
            for x, y in zip(p, o):
                assert x["mol"] == y["mol"] and abs(x["freq"] - y["freq"]) < 1e-6
                gains.append((abs(y["ml"] - y["dft"]) - abs(x["ml"] - x["dft"]), x["mol"],
                              x["freq"], E))  # fmt: skip
        g = np.array([t[0] for t in gains])
        tot = g.sum()
        order = np.argsort(-g)
        top = g[order[:10]].sum()
        w(f"- vs {DLAB[other]}: total lead {tot:.0f} km/mol over {len(g)} band×energy-model "
          f"pairs; POLAR-1 better in {np.mean(g > 0) * 100:.0f} %; the 10 largest gains sum to "
          f"{top:.0f} km/mol ({100 * top / tot:.0f} % of the total). Largest: "
          + "; ".join(f"{gains[i][1]} {gains[i][2]:.0f} cm-1 ({ELAB[gains[i][3]]}) {g[i]:+.0f}"
                      for i in order[:5]))  # fmt: skip

    # the 1700-2700 cm-1 bands: what they are, and the ranking without them
    w("\n## The 1700-2700 cm-1 bands (OMOL modes, cutoff 1, headline), km/mol\n")
    w("| molecule | B3LYP freq | B3LYP | POLAR-1 | MDP | MACE4IR |")
    w("|---|---|---|---|---|---|")
    E = "mace_omol"
    mols = mols_for(bands, E, SYMTOP)
    sel = {D: [x for x in select(bands, mols, E, D, cut) if 1700 < x["freq"] <= 2700]
           for D in DIPOLE}  # fmt: skip
    for i, x in enumerate(sel["mace_polar1"]):
        w(f"| {x['mol']} | {x['freq']:.0f} | {x['dft']:.0f} | "
          + " | ".join(f"{sel[D][i]['ml']:.0f}" for D in DIPOLE) + " |")  # fmt: skip
    for lab, keep in (("without the 1700-2700 cm-1 bands", lambda x: not 1700 < x["freq"] <= 2700),
                      ("without CO2", lambda x: x["mol"] != "carbon_dioxide")):  # fmt: skip
        w(f"\n### MAE (km/mol) {lab} (cutoff 1, headline)\n")
        w(table({k: [x for x in v if keep(x)] for k, v in bands.items()}, cut, SYMTOP))
        w(f"\n### median |log10 ratio| {lab}\n")
        w(table({k: [x for x in v if keep(x)] for k, v in bands.items()}, cut, SYMTOP,
                metric="medlog", fmt="{:.3f}"))  # fmt: skip

    # overlap quality of the selected bands
    sel = [x for E in ENERGY for D in DIPOLE
           for x in select(bands, mols_for(bands, E, SYMTOP), E, D, cut)]  # fmt: skip
    ov = np.array([x["overlap"] for x in sel])
    w(f"\n## Mode-matching quality of the selected bands\n\nOverlap < 0.5: {np.mean(ov < 0.5)*100:.1f} %, "
      f"< 0.8: {np.mean(ov < 0.8)*100:.1f} %, median {np.median(ov):.3f} (n = {len(ov)}).\n")  # fmt: skip
    w("\n### MAE with overlap >= 0.8 only (cutoff 1, headline)\n")
    filt = {k: [x for x in v if x["overlap"] >= 0.8] for k, v in bands.items()}
    w(table(filt, cut, SYMTOP))

    (HERE / "a1_tables.md").write_text("\n".join(out) + "\n")
    with (HERE / "a1_bands.csv").open("w", newline="") as fh:
        fields = ["mol", "E", "D", "freq", "ml_freq", "dft", "ml", "overlap", "ndeg"]
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for v in bands.values():
            for x in v:
                wr.writerow({k: x[k] for k in fields})
    print("\n".join(out))


if __name__ == "__main__":
    main()
