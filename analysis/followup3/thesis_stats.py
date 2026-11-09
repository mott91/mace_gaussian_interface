"""Follow-up 3, phase 2: every statistic of the results chapter, as defined in Methods 3.6.

One script, one data layer: the band table of the campaign (campaigns/2026/figures/supervisor/data/
bands.csv and chi.csv, written by mace_gaussian/analysis/band_table.py), extended here with the
harmonic IR intensities read from each run's results.json through the same eigenvector assignment.

Rules (Methods 3.6, thesis repo latex/chapters/03-methods.tex, as of 2026-11-09):
- errors are MACE minus B3LYP; frequencies: median |error|, MAE, bias, RMSE, max (worst band named),
  share > 100 cm-1, R2 against y = x (reported, not used to rank);
- degenerate B3LYP levels (harmonic frequencies within 0.5 cm-1; NH3/PH3 E modes and the CO2 bend)
  count as one band: mean frequency, summed intensity;
- NH3 and PH3 are left out of every anharmonic quantity (shown alone, and as [with] values for the
  fundamentals); harmonic quantities keep them;
- CO2 bend overtones and combination bands are not compared;
- overtones/combinations over all paired bands and over bands with B3LYP VPT2 intensity >= 1 km/mol;
- five regions by B3LYP harmonic frequency (<500, 500-1600, 1600-2700, 2700-3200, >3200); below 500
  split into X-H internal rotations and other modes (hand-checked labels, rotor_labels.csv);
- intensities: B3LYP >= 1 km/mol in the respective comparison (harmonic / VPT2); MAE, share within a
  factor of two, R2 of log10 against y = x; per energy model and pooled over the four energy models;
- pooled over each model's own molecules, each band once; groups of Table 3.3; per molecule; and the
  common set of all models as a check;
- per-molecule win counts (molecules in which the row model has the lower MAE than the column model);
- scaled harmonic: lambda = sum(w nu) / sum(w^2) against B3LYP VPT2 fundamentals, leave-one-molecule-out
  (the text) and fitted on all molecules (for comparison); B3LYP harmonic the same way;
- experiment: band_origins.json (explicit modes where given, else nearest B3LYP VPT2 within 200 cm-1,
  excluded bands dropped), B3LYP harmonic scaled by 0.9614, headline without NH3/PH3 and CO2 nu1;
- VCI: the character-paired table of scripts/vci_comparison.py, clean rows only.

Outputs (next to this script, in out/): stats.md with every table, and one CSV per table.
Run:  python analysis/followup3/thesis_stats.py   (environment with pandas; reads campaign data only)
"""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from mace_gaussian.analysis.master_table import assign_band_origins, load_band_origins  # noqa: E402

CAMP = ROOT / "campaigns" / "2026"
DATA = CAMP / "figures" / "supervisor" / "data"
OUT = HERE / "out"

E_MODELS = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]
D_MODELS = ["mace_polar1", "mace_mdp", "mace_ml"]
LAB = {"mace_anicc": "ANI-cc", "mace_omol": "OMOL", "mace_off": "OFF", "mace_polar": "POLAR",
       "mace_polar1": "POLAR-1", "mace_mdp": "MDP", "mace_ml": "MACE4IR"}  # fmt: skip
SYMTOP = {"ammonia", "phosphine"}
GROUPS = {
    "acids": ["formic_acid_trans", "formic_acid_cis", "acetic_acid", "propionic_acid",
              "butyric_acid", "valeric_acid", "caproic_acid"],
    "amino acid": ["glycine"],
    "alcohols": ["methanol", "ethanol", "1-propanol", "isopropanol", "1-butanol", "tert-butanol"],
    "small molecules": ["carbon_dioxide", "sulfur_dioxide", "ammonia", "phosphine",
                        "hydrogen_peroxide", "nitric_acid", "carbonic_acid"],
    "sugar": ["beta_glucose"],
}  # fmt: skip
GROUP_OF = {m: g for g, ms in GROUPS.items() for m in ms}
REGIONS = [(0, 500, "< 500"), (500, 1600, "500-1600"), (1600, 2700, "1600-2700"),
           (2700, 3200, "2700-3200"), (3200, 1e9, "> 3200")]  # fmt: skip
DEGENERATE_CM = 0.5
INT_MIN = 1.0  # km/mol
LARGE = 100.0  # cm-1
SCALE_B3LYP_EXP = 0.9614
CO2_BEND = {1, 2}  # B3LYP modes of the CO2 bend (checkpoint order)

TABLES: list[tuple[str, pd.DataFrame | str]] = []


def table(name: str, df: pd.DataFrame | str) -> None:
    TABLES.append((name, df))
    if isinstance(df, pd.DataFrame):
        df.to_csv(OUT / f"{name}.csv", index=True)


def region(f: float) -> str:
    return next(lab for lo, hi, lab in REGIONS if lo <= f < hi)


# ------------------------------------------------------------------------------
# Metrics
# ------------------------------------------------------------------------------


def r2(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ss = ((x - x.mean()) ** 2).sum()
    return float(1 - ((y - x) ** 2).sum() / ss) if len(x) > 1 and ss > 0 else np.nan


def freq_metrics(df: pd.DataFrame, ref: str, ml: str) -> dict:
    if df.empty:
        return {"n": 0}
    e = (df[ml] - df[ref]).to_numpy()
    a = np.abs(e)
    w = df.iloc[int(np.argmax(a))]
    worst = f"{w.molecule} {w.get('band', '')}".strip()
    return {"n": len(e), "median": np.median(a), "MAE": a.mean(), "bias": e.mean(),
            "RMSE": np.sqrt((e ** 2).mean()), "max": a.max(), "worst": worst,
            ">100 %": 100 * (a > LARGE).mean(), "R2": r2(df[ref], df[ml])}  # fmt: skip


def int_metrics(df: pd.DataFrame, ref: str, ml: str) -> dict:
    df = df[df[ref] >= INT_MIN]
    if df.empty:
        return {"n": 0}
    x, y = df[ref].to_numpy(), df[ml].to_numpy()
    ratio = np.where(y > 0, y / x, 0)
    yl = np.log10(np.clip(y, 1e-4, None))  # a zero ML intensity is set to 1e-4 km/mol for the log
    return {"n": len(x), "MAE": np.abs(y - x).mean(), "within 2x %": 100 * ((ratio >= 0.5) & (ratio <= 2)).mean(),
            "R2 log10": r2(np.log10(x), yl), "ML zero": int((y <= 0).sum())}  # fmt: skip


def fmt(v, nd=1):
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    return f"{v:.{nd}f}" if isinstance(v, float | np.floating) else str(v)


def md(df: pd.DataFrame, nd=1) -> str:
    def dec(col):
        c = str(col)
        if c in ("n", "molecules", "ML zero"):
            return 0
        if "lambda" in c:
            return 4
        return 3 if ("R2" in c or "r2" in c) else nd
    cols = [str(c) for c in df.columns]
    name = df.index.name if isinstance(df.index.name, str) else " / ".join(df.index.names or [])
    out = ["| " + " | ".join([name or ""] + cols) + " |", "|" + "---|" * (len(cols) + 1)]
    for idx, row in df.iterrows():
        idx = " / ".join(idx) if isinstance(idx, tuple) else str(idx)
        cells = []
        for c, v in zip(df.columns, row):
            d = dec(c)
            cells.append(str(int(v)) if d == 0 and isinstance(v, float | int | np.floating) and v == v else fmt(v, d))
        out.append("| " + " | ".join([idx] + cells) + " |")
    return "\n".join(out)


# ------------------------------------------------------------------------------
# Data: band table + harmonic intensities, degenerate merge, regions, rotors
# ------------------------------------------------------------------------------


def harmonic_table() -> pd.DataFrame:
    """Harmonic frequency + intensity of every B3LYP mode and its assigned MACE mode, per run."""
    rows = []
    for rep_path in sorted((CAMP / "analysis_results").glob("*/report_data.json")):
        mol = rep_path.parent.name
        rep = json.loads(rep_path.read_text())
        dft_dir = next(p for p in sorted((CAMP / "comparison_results" / mol).glob("b3lyp*")) if p.is_dir())
        dft = sorted(json.loads((dft_dir / "results.json").read_text())["frequencies"]["harmonic"],
                     key=lambda e: e["freq_cm"])  # fmt: skip
        for comp in rep["comparisons"]:
            run = CAMP / "comparison_results" / mol / comp["name"] / "results.json"
            if not run.exists():
                continue
            energy, dipole = next((comp["name"][: -len(d) - 1], d) for d in D_MODELS if comp["name"].endswith("_" + d))
            ml = sorted(json.loads(run.read_text())["frequencies"]["harmonic"], key=lambda e: e["freq_cm"])
            for m_idx, d_idx in comp["mode_mapping"].items():
                rows.append(dict(molecule=mol, energy_model=energy, dipole_model=dipole, i=int(d_idx) + 1,
                                 dft_harmonic=dft[int(d_idx)]["freq_cm"], ml_harmonic=ml[int(m_idx)]["freq_cm"],
                                 dft_int=dft[int(d_idx)]["ir_intensity"], ml_int=ml[int(m_idx)]["ir_intensity"]))  # fmt: skip
    return pd.DataFrame(rows)


def degenerate_groups(freqs: pd.Series) -> dict[int, int]:
    """{mode i: group id} for one molecule; modes within 0.5 cm-1 share the id of the lowest one."""
    s = freqs.sort_values()
    gid, prev, out = None, None, {}
    for i, f in s.items():
        if prev is None or f - prev > DEGENERATE_CM:
            gid = i
        out[i] = gid
        prev = f
    return out


def merge_degenerate(df: pd.DataFrame, mean_cols: list[str], sum_cols: list[str], keys: list[str]) -> pd.DataFrame:
    agg = {c: "mean" for c in mean_cols} | {c: "sum" for c in sum_cols} | {"i": "min", "ndeg": "size"}
    return df.assign(ndeg=1).groupby(keys + ["grp"], as_index=False).agg(agg)


def load():
    bands = pd.read_csv(DATA / "bands.csv")
    chi = pd.read_csv(DATA / "chi.csv")
    harm = harmonic_table()
    # degenerate groups from the B3LYP harmonic frequencies of each molecule
    grp = {}
    for mol, g in harm.drop_duplicates(["molecule", "i"]).groupby("molecule"):
        for i, gid in degenerate_groups(g.set_index("i").dft_harmonic).items():
            grp[(mol, i)] = gid
    harm["grp"] = [grp[(m, i)] for m, i in zip(harm.molecule, harm.i)]
    fund = bands[bands.band_type == "fundamental"].copy()
    fund["grp"] = [grp[(m, int(i))] for m, i in zip(fund.molecule, fund.i)]

    # checks: the joined harmonic table and the band table pair the same modes
    chk = fund.merge(harm, on=["molecule", "energy_model", "dipole_model", "i"], suffixes=("", "_h"))
    assert len(chk) == len(fund), (len(chk), len(fund))
    assert (chk.dft_harmonic - chk.dft_harmonic_h).abs().max() < 0.05
    assert (chk.ml_harmonic - chk.ml_harmonic_h).abs().max() < 0.05

    keys = ["molecule", "energy_model", "dipole_model"]
    fund_m = merge_degenerate(fund, ["dft_harmonic", "dft_anharmonic", "ml_harmonic", "ml_anharmonic", "overlap"],
                              ["dft_intensity", "ml_intensity"], keys)  # fmt: skip
    harm_m = merge_degenerate(harm, ["dft_harmonic", "ml_harmonic"], ["dft_int", "ml_int"], keys)
    rot = pd.read_csv(HERE / "rotor_labels.csv")
    rotors = {(r.molecule, r.i) for r in rot.itertuples() if r.rotor}
    for df in (fund_m, harm_m):
        df["region"] = df.dft_harmonic.map(region)
        df["rotor"] = [(m, i) in rotors for m, i in zip(df.molecule, df.i)]
        df["group"] = df.molecule.map(GROUP_OF)
        df["band"] = [f"mode {i}" + (f" (+{n - 1})" if n > 1 else "") + f", {f:.0f} cm-1"
                      for i, n, f in zip(df.i, df.ndeg, df.dft_harmonic)]  # fmt: skip
    return bands, chi, fund_m, harm_m, rotors


# ------------------------------------------------------------------------------
# Sections
# ------------------------------------------------------------------------------


def pooled_freq(df, ref, ml, label, models=E_MODELS):
    rows = {LAB[m]: freq_metrics(df[df.energy_model == m], ref, ml) for m in models}
    return pd.DataFrame(rows).T.rename_axis(label)


def by(df, ref, ml, key, models=E_MODELS, what="MAE"):
    """MAE (or another metric) per key value and model, with n, as 'value (median) n'."""
    out = {}
    for k, g in df.groupby(key, sort=False):
        out[k] = {}
        for m in models:
            s = freq_metrics(g[g.energy_model == m], ref, ml)
            out[k][LAB[m]] = f"{s['MAE']:.1f} ({s['median']:.1f}) n{s['n']}" if s["n"] else "–"
    return pd.DataFrame(out).T


def wins(per_mol: pd.DataFrame) -> pd.DataFrame:
    """per_mol: index molecule, columns models, values MAE. Entry [a, b] = 'w/n': a lower in w of n."""
    cols = list(per_mol.columns)
    out = pd.DataFrame("", index=cols, columns=cols)
    for a in cols:
        for b in cols:
            if a == b:
                out.loc[a, b] = "–"
                continue
            both = per_mol[[a, b]].dropna()
            out.loc[a, b] = f"{int((both[a] < both[b]).sum())}/{len(both)}"
    return out


def section_frequencies(fund, harm):
    f = fund[fund.dipole_model == "mace_ml"]
    h = harm[harm.dipole_model == "mace_ml"]
    anh, anh_all = f[~f.molecule.isin(SYMTOP)], f
    table("freq_harmonic_pooled", pooled_freq(h, "dft_harmonic", "ml_harmonic", "harmonic (all molecules)"))
    table("freq_vpt2_pooled", pooled_freq(anh, "dft_anharmonic", "ml_anharmonic", "VPT2, NH3/PH3 out"))
    table("freq_vpt2_pooled_with_nh3_ph3", pooled_freq(anh_all, "dft_anharmonic", "ml_anharmonic", "VPT2, with NH3/PH3"))
    common = set.intersection(*[set(anh[anh.energy_model == m].molecule) for m in E_MODELS])
    common_h = set.intersection(*[set(h[h.energy_model == m].molecule) for m in E_MODELS])
    table("freq_vpt2_common", pooled_freq(anh[anh.molecule.isin(common)], "dft_anharmonic", "ml_anharmonic",
                                          f"VPT2, common set ({len(common)} molecules)"))  # fmt: skip
    table("freq_harmonic_common", pooled_freq(h[h.molecule.isin(common_h)], "dft_harmonic", "ml_harmonic",
                                              f"harmonic, common set ({len(common_h)} molecules)"))  # fmt: skip
    # NH3 / PH3 alone (VPT2)
    sym = {}
    for mol in sorted(SYMTOP):
        for m in E_MODELS:
            s = freq_metrics(f[(f.molecule == mol) & (f.energy_model == m)], "dft_anharmonic", "ml_anharmonic")
            sym.setdefault(mol, {})[LAB[m]] = s.get("MAE", np.nan)
    table("freq_vpt2_nh3_ph3_alone", pd.DataFrame(sym).T.rename_axis("VPT2 MAE"))
    # regions, with the < 500 split
    for df, ref, ml, name in ((anh, "dft_anharmonic", "ml_anharmonic", "vpt2"), (h, "dft_harmonic", "ml_harmonic", "harmonic")):
        d = df.copy()
        d["region5"] = d.region
        d.loc[d.region == "< 500", "region5"] = np.where(d.loc[d.region == "< 500", "rotor"], "< 500: X-H rotations", "< 500: other")
        order = ["< 500", "< 500: X-H rotations", "< 500: other", "500-1600", "1600-2700", "2700-3200", "> 3200"]
        t = pd.concat([by(d, ref, ml, "region"), by(d[d.region == "< 500"], ref, ml, "region5")])
        table(f"freq_{name}_regions", t.loc[[o for o in order if o in t.index]].rename_axis(f"{name}: MAE (median) n"))
        # full metric set per region (csv only)
        full = {(r, LAB[m]): freq_metrics(g[g.energy_model == m], ref, ml) for r, g in d.groupby("region5") for m in E_MODELS}
        pd.DataFrame(full).T.to_csv(OUT / f"freq_{name}_regions_full.csv")
        t = by(d, ref, ml, "group").loc[list(GROUPS)] if name == "harmonic" else by(d, ref, ml, "group").reindex(list(GROUPS)).dropna(how="all")
        table(f"freq_{name}_groups", t.rename_axis(f"{name}: MAE (median) n"))
        per_mol = d.groupby(["molecule", "energy_model"]).apply(lambda g: (g[ml] - g[ref]).abs().mean()).unstack()[E_MODELS]
        table(f"freq_{name}_per_molecule", per_mol.rename(columns=LAB).rename_axis(f"{name} MAE"))
        table(f"freq_{name}_wins", wins(per_mol.rename(columns=LAB)).rename_axis(f"{name}: row lower in"))
    return anh


def section_anharmonicity(fund, chi):
    f = fund[(fund.dipole_model == "mace_ml") & ~fund.molecule.isin(SYMTOP)].copy()
    f["dft_shift"], f["ml_shift"] = f.dft_anharmonic - f.dft_harmonic, f.ml_anharmonic - f.ml_harmonic
    out = {}
    for m in E_MODELS:
        g = f[f.energy_model == m]
        out[LAB[m]] = {"n": len(g), "MAE": (g.ml_shift - g.dft_shift).abs().mean(), "R2": r2(g.dft_shift, g.ml_shift)}
    table("shift_pooled", pd.DataFrame(out).T.rename_axis("anharmonic shift"))
    reg = {}
    for r, g in f.groupby("region"):
        reg[r] = {LAB[m]: f"{(x.ml_shift - x.dft_shift).abs().mean():.1f} / R2 {r2(x.dft_shift, x.ml_shift):.2f}"
                  for m in E_MODELS for x in [g[g.energy_model == m]] if len(x)}  # fmt: skip
    table("shift_regions", pd.DataFrame(reg).T.reindex([lab for *_, lab in REGIONS]).rename_axis("shift MAE / R2"))
    grp = {}
    for k, g in f.groupby("group"):
        grp[k] = {LAB[m]: f"{(x.ml_shift - x.dft_shift).abs().mean():.1f} / R2 {r2(x.dft_shift, x.ml_shift):.2f}"
                  for m in E_MODELS for x in [g[g.energy_model == m]] if len(x)}  # fmt: skip
    table("shift_groups", pd.DataFrame(grp).T.reindex([g for g in GROUPS if g in grp]).rename_axis("shift MAE / R2"))
    # x_ii from the X matrix (degenerate modes are already left out of chi.csv)
    c = chi[(chi.kind == "diagonal") & ~chi.molecule.isin(SYMTOP)].dropna(subset=["dft_xmatrix", "ml_xmatrix"]).copy()
    harmf = f.drop_duplicates(["molecule", "i"]).set_index(["molecule", "i"]).dft_harmonic
    c["region"] = [region(harmf.get((m, i), np.nan)) if (m, i) in harmf.index else "?" for m, i in zip(c.molecule, c.i)]
    c["group"] = c.molecule.map(GROUP_OF)
    out = {}
    for m in E_MODELS:
        g = c[c.energy_model == m]
        out[LAB[m]] = {"n": len(g), "MAE": (g.ml_xmatrix - g.dft_xmatrix).abs().mean(), "R2": r2(g.dft_xmatrix, g.ml_xmatrix),
                       "Pearson r2": np.corrcoef(g.dft_xmatrix, g.ml_xmatrix)[0, 1] ** 2}  # fmt: skip
    table("xii_pooled", pd.DataFrame(out).T.rename_axis("x_ii (X matrix)"))
    for key, name in (("region", "xii_regions"), ("group", "xii_groups")):
        t = {}
        for k, g in c.groupby(key):
            t[k] = {LAB[m]: f"{(x.ml_xmatrix - x.dft_xmatrix).abs().mean():.1f} / R2 {r2(x.dft_xmatrix, x.ml_xmatrix):.2f} n{len(x)}"
                    for m in E_MODELS for x in [g[g.energy_model == m]] if len(x)}  # fmt: skip
        order = [lab for *_, lab in REGIONS] if key == "region" else list(GROUPS)
        table(name, pd.DataFrame(t).T.reindex([o for o in order if o in t]).rename_axis("x_ii MAE / R2"))


def section_overtones(bands):
    b = bands[bands.band_type.isin(["overtone", "combination"]) & (bands.dipole_model == "mace_ml")
              & ~bands.molecule.isin(SYMTOP)].copy()  # fmt: skip
    co2_bend = (b.molecule == "carbon_dioxide") & (b.i.isin(CO2_BEND) | b.j.isin(CO2_BEND))
    b = b[~co2_bend]
    b["region"] = b.dft_harmonic.map(region)
    b["group"] = b.molecule.map(GROUP_OF)
    b["band"] = [f"{t[:4]} {int(i)}" + (f"+{int(j)}" if j == j else "") + f", {f:.0f} cm-1"
                 for t, i, j, f in zip(b.band_type, b.i, b.j, b.dft_anharmonic)]  # fmt: skip
    for vis, name in ((False, "all"), (True, "visible")):
        d = b[b.dft_intensity >= INT_MIN] if vis else b
        for bt in ("overtone", "combination"):
            x = d[d.band_type == bt]
            table(f"{bt}s_{name}_pooled", pooled_freq(x, "dft_anharmonic", "ml_anharmonic", f"{bt}s, {name}"))
            table(f"{bt}s_{name}_regions", by(x, "dft_anharmonic", "ml_anharmonic", "region")
                  .reindex([lab for *_, lab in REGIONS]).dropna(how="all").rename_axis(f"{bt}s {name}: MAE (median) n"))  # fmt: skip
            table(f"{bt}s_{name}_groups", by(x, "dft_anharmonic", "ml_anharmonic", "group")
                  .reindex(list(GROUPS)).dropna(how="all").rename_axis(f"{bt}s {name}: MAE (median) n"))  # fmt: skip
    # regions by parent mode: overtone -> region of its fundamental, combination -> region of the higher mode
    fh = bands[(bands.band_type == "fundamental") & (bands.dipole_model == "mace_ml")].drop_duplicates(["molecule", "i"])
    fh = fh.set_index(["molecule", "i"]).dft_harmonic
    b["parent_region"] = [region(max(fh[(m, int(i))], fh[(m, int(j))] if j == j else 0)) for m, i, j in zip(b.molecule, b.i, b.j)]
    for vis, name in ((False, "all"), (True, "visible")):
        d = b[b.dft_intensity >= INT_MIN] if vis else b
        for bt in ("overtone", "combination"):
            x = d[d.band_type == bt]
            table(f"{bt}s_{name}_regions_by_parent", by(x, "dft_anharmonic", "ml_anharmonic", "parent_region")
                  .reindex([lab for *_, lab in REGIONS]).dropna(how="all")
                  .rename_axis(f"{bt}s {name}, region of the parent mode: MAE (median) n"))  # fmt: skip
    n_removed = int(co2_bend.sum()) // 4 if co2_bend.any() else 0
    table("overtones_note", f"CO2 bend overtones/combinations removed: {n_removed} per energy model.")


def section_intensities(fund, harm):
    for df, ref, ml, name, keep_sym in ((harm, "dft_int", "ml_int", "harmonic", True),
                                        (fund, "dft_intensity", "ml_intensity", "vpt2", False)):  # fmt: skip
        d = df if keep_sym else df[~df.molecule.isin(SYMTOP)]
        rows = {}
        for e in E_MODELS + ["pooled"]:
            for dp in D_MODELS:
                x = d[d.dipole_model == dp] if e == "pooled" else d[(d.dipole_model == dp) & (d.energy_model == e)]
                rows[(LAB.get(e, e), LAB[dp])] = int_metrics(x, ref, ml)
        table(f"int_{name}", pd.DataFrame(rows).T.rename_axis(("energy model", "dipole model")))
        for key, order in (("region", [lab for *_, lab in REGIONS]), ("group", list(GROUPS))):
            t = {}
            for k, g in d.groupby(key):
                t[k] = {LAB[dp]: (lambda s: f"{s['MAE']:.1f} / {s['within 2x %']:.0f} % n{s['n']}" if s["n"] else "–")(
                    int_metrics(g[g.dipole_model == dp], ref, ml)) for dp in D_MODELS}  # fmt: skip
            table(f"int_{name}_{key}s", pd.DataFrame(t).T.reindex([o for o in order if o in t])
                  .rename_axis(f"{name} intensity, pooled over energy models: MAE / within 2x"))  # fmt: skip
        vis = d[d[ref] >= INT_MIN].assign(ae=lambda x: (x[ml] - x[ref]).abs())
        pm = vis.groupby(["molecule", "dipole_model"]).ae.mean().unstack()[D_MODELS].rename(columns=LAB)
        pm["n (per energy model)"] = vis[vis.dipole_model == "mace_ml"].groupby("molecule").size() / vis[vis.dipole_model == "mace_ml"].groupby("molecule").energy_model.nunique()
        table(f"int_{name}_per_molecule", pm.rename_axis(f"{name} intensity MAE, pooled over energy models"))
        for e in E_MODELS + ["pooled"]:
            x = vis if e == "pooled" else vis[vis.energy_model == e]
            per_mol = x.groupby(["molecule", "dipole_model"]).ae.mean().unstack()[D_MODELS].rename(columns=LAB)
            table(f"int_{name}_wins_{LAB.get(e, e)}", wins(per_mol).rename_axis(f"{name}, {LAB.get(e, e)}: row lower in"))


def lam(w, nu):
    return float((w * nu).sum() / (w ** 2).sum())


def section_scaling(fund):
    f = fund[(fund.dipole_model == "mace_ml") & ~fund.molecule.isin(SYMTOP)]
    rows, lams = {}, {}
    sources = [(LAB[m], f[f.energy_model == m], "ml_harmonic") for m in E_MODELS]
    sources.append(("B3LYP", f[f.energy_model == "mace_omol"], "dft_harmonic"))  # B3LYP: all 20 molecules
    for name, g, col in sources:
        lam_all = lam(g[col], g.dft_anharmonic)
        loo = []
        for mol in g.molecule.unique():
            tr, te = g[g.molecule != mol], g[g.molecule == mol]
            loo.append(lam(tr[col], tr.dft_anharmonic) * te[col] - te.dft_anharmonic)
        e_loo, e_all = pd.concat(loo).abs(), (lam_all * g[col] - g.dft_anharmonic).abs()
        vpt2 = (g.ml_anharmonic - g.dft_anharmonic).abs() if name != "B3LYP" else pd.Series(dtype=float)
        rows[name] = {"n": len(g), "lambda (all)": lam_all, "scaled MAE, LOO": e_loo.mean(), "median, LOO": e_loo.median(),
                      "scaled MAE, all-fit": e_all.mean(), "median, all-fit": e_all.median(),
                      "VPT2 MAE": vpt2.mean() if len(vpt2) else np.nan, "VPT2 median": vpt2.median() if len(vpt2) else np.nan}  # fmt: skip
        lams[name] = lam_all
    table("scaled_harmonic", pd.DataFrame(rows).T.rename_axis("against B3LYP VPT2 fundamentals"))
    # by region (LOO)
    reg = {}
    for name, g, col in sources:
        parts = []
        for mol in g.molecule.unique():
            tr, te = g[g.molecule != mol], g[g.molecule == mol]
            parts.append(te.assign(err=lam(tr[col], tr.dft_anharmonic) * te[col] - te.dft_anharmonic))
        e = pd.concat(parts)
        for r, x in e.groupby("region"):
            reg.setdefault(r, {})[name] = f"{x.err.abs().mean():.1f}"
    table("scaled_harmonic_regions", pd.DataFrame(reg).T.reindex([lab for *_, lab in REGIONS]).rename_axis("scaled (LOO) MAE"))


def section_experiment(fund):
    f = fund[fund.dipole_model == "mace_ml"]
    rows = []
    for mol in sorted(f.molecule.unique()):
        entry = load_band_origins(mol)
        if entry is None:
            continue
        d = f[f.molecule == mol]
        # the experiment is matched mode by mode (a degenerate band takes both of its modes), so use
        # the unmerged band table here: B3LYP VPT2 per mode
        raw = pd.read_csv(DATA / "bands.csv")
        raw = raw[(raw.molecule == mol) & (raw.band_type == "fundamental") & (raw.dipole_model == "mace_ml")]
        ref = raw.drop_duplicates("i").set_index("i").dft_anharmonic.to_dict()
        taken = assign_band_origins(entry["bands"], ref)
        groups: dict[str, list[int]] = {}
        info = {}
        for k, band in taken.items():
            groups.setdefault(band["label"], []).append(k)
            info[band["label"]] = band
        for label, ks in groups.items():
            band = info[label]
            dd = raw[raw.energy_model == "mace_omol"].set_index("i")
            row = dict(molecule=mol, label=label, description=band.get("description", ""), exp=float(band["freq_cm"]),
                       modes="+".join(map(str, sorted(ks))), assignment=band["assignment"],
                       B3LYP_vpt2=dd.loc[ks, "dft_anharmonic"].mean(),
                       B3LYP_harm_scaled=SCALE_B3LYP_EXP * dd.loc[ks, "dft_harmonic"].mean())  # fmt: skip
            for m in E_MODELS:
                mm = raw[raw.energy_model == m].set_index("i")
                row[LAB[m]] = mm.loc[ks, "ml_anharmonic"].mean() if set(ks) <= set(mm.index) else np.nan
            rows.append(row)
        for band in entry["bands"]:
            if band.get("exclude"):
                rows.append(dict(molecule=mol, label=band["label"], description=band.get("description", ""),
                                 exp=float(band["freq_cm"]), modes="", assignment="excluded: " + band["exclude"]))  # fmt: skip
    ex = pd.DataFrame(rows)
    ex.to_csv(OUT / "experiment_bands.csv", index=False)
    used = ex[~ex.assignment.str.startswith("excluded")]
    co2_nu1 = (used.molecule == "carbon_dioxide") & (used.label == "ν1")
    head = used[~used.molecule.isin(SYMTOP) & ~co2_nu1]
    with_sym = used[~co2_nu1]
    methods = ["B3LYP_vpt2", "B3LYP_harm_scaled"] + [LAB[m] for m in E_MODELS]
    for d, name in ((head, "headline"), (with_sym, "with_nh3_ph3")):
        t = {}
        for mth in methods:
            x = d.dropna(subset=[mth])
            e = x[mth] - x.exp
            w = x.iloc[int(np.argmax(e.abs()))]
            t[mth] = {"n": len(x), "molecules": x.molecule.nunique(), "MAE": e.abs().mean(), "median": e.abs().median(),
                      "bias": e.mean(), "max": e.abs().max(), "worst": f"{w.molecule} {w.label}"}  # fmt: skip
        table(f"experiment_{name}", pd.DataFrame(t).T.rename_axis(f"experiment, {name}"))
    per_mol = head.groupby("molecule").apply(lambda x: pd.Series({m: (x[m] - x.exp).abs().mean() for m in methods} | {"n": len(x)}))
    table("experiment_per_molecule", per_mol.rename_axis("experiment MAE"))
    counts = ex.assignment.str.split(":").str[0].value_counts().to_dict()
    table("experiment_note", f"bands: {len(ex)} listed, {len(used)} assigned ({counts}); headline {len(head)} bands in "
          f"{head.molecule.nunique()} molecules; with NH3/PH3 {len(with_sym)}.")  # fmt: skip


def section_vci():
    p = CAMP / "experiment_2026-11-03" / "vci_comparison_character_paired.csv"
    v = pd.read_csv(p)
    clean = v[v.flag.isna()]
    t = {}
    for name in ["B3LYP"] + E_MODELS:
        col = f"{name}_vpt2"
        x = clean.dropna(subset=[col])
        e = x[col] - x.vci
        t[LAB.get(name, name)] = {"n": len(x), "MAE": e.abs().mean(), "median": e.abs().median(), "bias": e.mean(), "max": e.abs().max()}
    table("vci_vpt2", pd.DataFrame(t).T.rename_axis(f"VPT2 vs VCI ({len(clean)} clean of {len(v)})"))
    hh = v.dropna(subset=["harm_cc"])
    hh = hh[~((hh.molecule == "methanol") & (hh.label == "nu12"))]
    t = {}
    for name in ["B3LYP"] + E_MODELS:
        col = f"{name}_h"
        x = hh.dropna(subset=[col])
        e = x[col] - x.harm_cc
        t[LAB.get(name, name)] = {"n": len(x), "MAE": e.abs().mean(), "median": e.abs().median(), "bias": e.mean(), "max": e.abs().max()}
    table("vci_harmonic", pd.DataFrame(t).T.rename_axis("harmonic vs CCSD(T)-F12 harmonic"))


def main():
    OUT.mkdir(exist_ok=True)
    bands, chi, fund, harm, rotors = load()
    section_frequencies(fund, harm)
    section_anharmonicity(fund, chi)
    section_overtones(bands)
    section_intensities(fund, harm)
    section_scaling(fund)
    section_experiment(fund)
    section_vci()
    text = ["# Follow-up 3 statistics (generated by thesis_stats.py)", "",
            f"Rotor labels (hand-checked): {len(rotors)} X-H internal rotations below 500 cm-1.", ""]  # fmt: skip
    for name, df in TABLES:
        text += [f"## {name}", "", df if isinstance(df, str) else md(df), ""]
    (OUT / "stats.md").write_text("\n".join(text))
    print(f"wrote {len(TABLES)} tables to {OUT}")


if __name__ == "__main__":
    main()
