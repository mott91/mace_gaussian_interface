"""Cost figures of the results chapter, drawn from the raw timings.

  cost_scaling    carboxylic acid series: (a) wall time against the number of atoms N, log-log, power-law
                  fits, 6-31G(d,p) basis functions on the upper axis; (b) speedup over B3LYP
  cost_accuracy   mean wall time per molecule against MAE, each panel with the time taken over the same
                  molecules as its errors: (a) coupled-cluster set, (b) experiment set (the panel against
                  B3LYP on the common set was removed in figures-ch4-tweaks)

Times: B3LYP = total_elapsed_s of the baseline results.json (both Gaussian links, optimisation and
VPT2); MACE = mean wall_s of the three dipole-model runs (analysis/ch5_timing/timing_runs.csv, which is
results.json runtime_s). Basis functions: NBasis of the B3LYP Gaussian logs (Cartesian d functions,
Gaussian's default for 6-31G(d,p)). MAEs from the stats tables of analysis/followup3/out.

    python analysis/results_figures/fig_cost.py [output_dir]
"""

import importlib.util
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("thesis_style", ROOT / "mace_gaussian" / "analysis" / "thesis_style.py")
style = importlib.util.module_from_spec(spec)
sys.modules["thesis_style"] = style
spec.loader.exec_module(style)

CAMP = ROOT / "campaigns" / "2026"
STATS = ROOT / "analysis" / "followup3" / "out"
ELAB = {"mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23", "mace_polar": "MACE-POLAR-1"}
E = list(ELAB)
KEY = {"b3lyp": "B3LYP", "mace_anicc": "ANI-cc", "mace_omol": "OMOL", "mace_off": "OFF", "mace_polar": "POLAR"}
ACIDS = ["formic_acid_trans", "acetic_acid", "propionic_acid", "butyric_acid", "valeric_acid", "caproic_acid"]
CC_SET = ["methanol", "formic_acid_trans", "carbon_dioxide", "carbonic_acid"]  # vci_vpt2 (30 bands)
SYMTOP = ["ammonia", "phosphine"]  # out of the pooled VPT2 and experiment statistics


def color(k: str) -> str:
    return style.INK if k == "b3lyp" else style.model_color(k)


def label(k: str) -> str:
    return "B3LYP" if k == "b3lyp" else ELAB[k]


def find(d, key):
    if isinstance(d, dict):
        if d.get(key) is not None:
            return d[key]
        for v in d.values():
            r = find(v, key)
            if r is not None:
                return r
    return None


def timings() -> pd.DataFrame:
    """Wall time in minutes, one row per molecule, columns b3lyp and the four energy models."""
    t = pd.read_csv(ROOT / "analysis" / "ch5_timing" / "timing_runs.csv")
    g = t.groupby(["molecule", "energy_model"]).wall_s.agg(["count", "mean"])
    assert (g["count"] == 3).all(), g[g["count"] != 3]
    out = g["mean"].unstack() / 60
    b = {}
    for mol in out.index:
        res = next((CAMP / "comparison_results" / mol).glob("b3lyp*/results.json"))
        b[mol] = float(find(json.loads(res.read_text()), "total_elapsed_s")) / 60
    out.insert(0, "b3lyp", pd.Series(b))
    out["atoms"] = t.groupby("molecule").atoms.first()
    return out


def nbasis(mol: str) -> int:
    log = next((CAMP / "comparison_results" / mol).glob("b3lyp*/*.log")).read_text()
    nb, cart = re.search(r"(\d+) basis functions,\s+\d+ primitive gaussians,\s+(\d+) cartesian", log).groups()
    assert nb == cart, (mol, nb, cart)  # Cartesian d functions (6D)
    return int(nb)


def cost_scaling(t: pd.DataFrame, out: Path) -> None:
    a = t.loc[ACIDS]
    n = a.atoms.to_numpy(float)
    nb = [nbasis(m) for m in ACIDS]
    print("acid series (min):")
    print(a.assign(nbasis=nb).round(2).to_string())
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(style.figwidth(), 2.75))
    handles = []
    for k in ["b3lyp", *E]:
        y = a[k].to_numpy()
        slope, icpt = np.polyfit(np.log(n), np.log(y), 1)
        xs = np.array([4.6, 21.5])
        ax.plot(xs, np.exp(icpt) * xs**slope, color=color(k), lw=1.0, zorder=2)
        ax.plot(n, y, "o", color=color(k), ms=3.6, mec="none", zorder=3)
        handles.append(Line2D([], [], color=color(k), marker="o", ms=3.6, mec="none", lw=1.0,
                              label=f"{label(k)} ($N^{{{slope:.1f}}}$)"))  # fmt: skip
        print(f"exponent vs N  {k:11s} {slope:.3f}")
        if k != "b3lyp":
            bx.plot(n, a.b3lyp / a[k], "o-", color=color(k), ms=3.6, mec="none", lw=1.0)
    sb = np.polyfit(np.log(nb), np.log(a.b3lyp), 1)[0]
    print(f"exponent vs NBasis  b3lyp {sb:.3f}")
    print("speedup over B3LYP:")
    print((1 / a[E].div(a.b3lyp, axis=0)).round(2).to_string())
    for x in (ax, bx):
        style.clean_axes(x)
        x.set_xscale("log")
        x.set_xlim(4.6, 21.5)
        x.set_xticks(n, [f"{v:.0f}" for v in n])
        x.xaxis.set_minor_locator(plt.NullLocator())
        x.set_xlabel(r"number of atoms $N$")
    ax.set_yscale("log"), ax.set_ylim(0.3, 400)
    ax.set_ylabel("wall time / min")
    top = ax.secondary_xaxis("top")
    top.set_xticks(n, [str(v) for v in nb])
    top.xaxis.set_minor_locator(plt.NullLocator())
    top.set_xlabel("6-31G(d,p) basis functions", fontsize=style.ANNOT_FS)
    top.tick_params(labelsize=style.ANNOT_FS)
    bx.axhline(1, color=style.INK_MUTED, lw=0.7, ls=":", zorder=1)
    bx.set_ylim(0, 16), bx.set_yticks([0, 5, 10, 15])
    bx.set_ylabel("speedup over B3LYP")
    for x, s in ((ax, "(a)"), (bx, "(b)")):
        x.text(-0.2, 1.0, s, transform=x.transAxes, ha="right", va="bottom", fontsize=style.BASE_FS)
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.0), frameon=False,
               fontsize=style.ANNOT_FS, handlelength=1.6, columnspacing=1.4, borderpad=0)  # fmt: skip
    style.tight(fig, pad=0.4, w_pad=1.5)
    style.save(fig, out, "cost_scaling")


def cost_accuracy(t: pd.DataFrame, out: Path) -> None:
    exp_bands = pd.read_csv(STATS / "experiment_bands.csv")
    exp_set = sorted(set(exp_bands[~exp_bands.assignment.astype(str).str.startswith("excluded")].molecule) - set(SYMTOP))
    exp = pd.read_csv(STATS / "experiment_headline.csv", index_col=0)
    cc = pd.read_csv(STATS / "vci_vpt2.csv", index_col=0)
    assert exp.loc["ANI-cc", "molecules"] == len(exp_set) == 5 and exp.loc["ANI-cc", "n"] == 46
    assert cc.loc["ANI-cc", "n"] == 30
    panels = [
        ("(a) against coupled cluster", CC_SET, lambda k: cc.loc[KEY[k], "MAE"], [0.5, 1, 2]),
        ("(b) against experiment", exp_set, lambda k: exp.loc["B3LYP_vpt2" if k == "b3lyp" else KEY[k], "MAE"],
         [0.5, 1, 2, 3]),
    ]  # fmt: skip
    fig, axes = plt.subplots(1, 2, figsize=(style.figwidth(), 2.35), sharey=True)
    for ax, (title, mols, mae, ticks) in zip(axes, panels):
        tm = t.loc[mols, ["b3lyp", *E]].mean()
        print(f"{title}: {len(mols)} molecules {mols}")
        for k in ["b3lyp", *E]:
            y = mae(k)
            print(f"   {k:11s} mean wall time {tm[k]:7.2f} min  MAE {y:5.1f}")
            ax.plot([tm[k]], [y], "o", ms=5, mec="none", color=color(k), zorder=3)
        lo, hi = tm.min(), tm.max()
        ax.set_xscale("log"), ax.set_xlim(lo / 1.5, hi * 1.5)
        assert ticks[0] > lo / 1.5 and ticks[-1] < hi * 1.5, (title, lo, hi)
        ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        style.clean_axes(ax)
        ax.set_title(title, fontsize=style.BASE_FS, loc="left")
    axes[0].set_ylim(0, 62), axes[0].set_yticks([0, 20, 40, 60])
    axes[0].set_ylabel("MAE of the VPT2\nfundamentals / cm$^{-1}$")
    fig.supxlabel("mean wall time per molecule / min", fontsize=style.BASE_FS)
    handles = [Line2D([], [], color=color(k), marker="o", ms=5, mec="none", ls="none", label=label(k)) for k in ["b3lyp", *E]]
    fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 1.0), frameon=False,
               fontsize=style.ANNOT_FS, handletextpad=0.2, columnspacing=1.0, borderpad=0)  # fmt: skip
    style.tight(fig, pad=0.4, w_pad=0.8)
    style.save(fig, out, "cost_accuracy")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out"
    style.apply_style()
    t = timings()
    cost_scaling(t, out)
    cost_accuracy(t, out)
