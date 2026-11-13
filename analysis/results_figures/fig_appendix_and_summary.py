"""Appendix gallery (stacked spectra and mode-overlap matrices of all 22 molecules) and two summary figures.

  <out>/../results_appendix/spectra_stacked_<mol>, overlap_matrices_<mol>   via scripts/make_thesis_figures.py
  <out>/band_types      MAE of fundamentals, overtones and combination bands (stats CSVs of followup3)
  <out>/cost_accuracy   wall time for caproic acid against MAE vs experiment and vs coupled cluster

Run from a scratch directory:  python analysis/results_figures/fig_appendix_and_summary.py <out_dir>
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig_spectra_cost_c60 as base

mtf, style, ROOT, CAMP, ELAB, E = base.mtf, base.style, base.ROOT, base.CAMP, base.ELAB, base.E
STATS = ROOT / "analysis" / "followup3" / "out"
KEY = {"mace_anicc": "ANI-cc", "mace_omol": "OMOL", "mace_off": "OFF", "mace_polar": "POLAR"}


def gallery(out: Path) -> None:
    style.MODEL_LABEL.update(ELAB)
    mtf.ANALYSIS, mtf.COMPARISON, mtf.OUT = CAMP / "analysis_results", CAMP / "comparison_results", out
    mtf.DIPOLE_PREF[:] = ["mace_polar1", "mace_mdp", "mace_ml"]
    mtf.MODEL_ORDER[:] = ["mace_omol", "mace_off", "mace_anicc", "mace_polar"]
    mols = sorted(p.name for p in (CAMP / "analysis_results").iterdir() if (p / "report_data.json").exists())
    for mol in mols:
        mtf.fig_stacked_spectra(mol)
        mtf.fig_overlap_matrices(mol)
        plt.close("all")
    made = sorted(p.stem for p in out.glob("*.pdf"))
    print(len(mols), "molecules;", sum(m.startswith("spectra") for m in made), "spectra,", sum(m.startswith("overlap") for m in made), "overlap figures")


def band_types(out: Path) -> None:
    sets = [("fundamentals", "freq_vpt2_pooled"), ("overtones,\nall", "overtones_all_pooled"), ("overtones,\nvisible", "overtones_visible_pooled"),
            ("combinations,\nall", "combinations_all_pooled"), ("combinations,\nvisible", "combinations_visible_pooled")]
    fig, ax = plt.subplots(figsize=(style.figwidth(), 2.9))
    style.clean_axes(ax)
    w = 0.19
    for j, m in enumerate(E):
        vals = [pd.read_csv(STATS / f"{f}.csv", index_col=0).loc[KEY[m], "MAE"] for _, f in sets]
        ax.bar(np.arange(len(sets)) + (j - 1.5) * w, vals, width=w, color=style.model_color(m), label=ELAB[m])
        print(m, [round(float(v), 1) for v in vals])
    ax.set_xticks(range(len(sets)), [s for s, _ in sets])
    ax.set_ylabel("MAE / cm$^{-1}$")
    ax.legend(frameon=False, fontsize=style.ANNOT_FS, ncol=2, loc="upper left")
    fig.tight_layout()
    style.save(fig, out, "band_types")


def cost_accuracy(out: Path) -> None:
    info = json.loads((CAMP / "comparison_results" / "batch_manifest.json").read_text())["molecules"]["caproic_acid"]
    b = next((CAMP / "comparison_results" / "caproic_acid").glob("b3lyp*/results.json"))
    t = {"b3lyp": base.find(json.loads(b.read_text()), "total_elapsed_s") / 60}
    for e in E:
        t[e] = np.mean([c["runtime_s"] for k, c in info["combinations"].items() if k.startswith(e + "_")]) / 60
    exp = pd.read_csv(STATS / "experiment_headline.csv", index_col=0)["MAE"]
    vci = pd.read_csv(STATS / "vci_vpt2.csv", index_col=0)["MAE"]
    fig, ax = plt.subplots(figsize=(style.figwidth(0.8), 3.2))
    style.clean_axes(ax)
    for key in ["b3lyp", *E]:
        color = style.INK if key == "b3lyp" else style.model_color(key)
        ye, yc = exp["B3LYP_vpt2" if key == "b3lyp" else KEY[key]], vci["B3LYP" if key == "b3lyp" else KEY[key]]
        ax.plot([t[key]], [ye], "o", ms=6, color=color), ax.plot([t[key]], [yc], "o", ms=6, mfc="none", mec=color, mew=1.1)
        ax.annotate("B3LYP" if key == "b3lyp" else ELAB[key], (t[key], max(ye, yc)), xytext=(0, 6), textcoords="offset points",
                    ha="center", fontsize=style.ANNOT_FS, color=color)
        print(f"{key:11s} {t[key]:6.1f} min  exp {ye:5.1f}  cc {yc:5.1f}")
    ax.plot([], [], "o", ms=6, color=style.INK_MUTED, label="against experiment")
    ax.plot([], [], "o", ms=6, mfc="none", mec=style.INK_MUTED, mew=1.1, label="against coupled cluster")
    ax.set_xscale("log"), ax.set_xlim(8, 500), ax.set_ylim(0, 66)
    ax.set_xticks([10, 30, 100, 300]), ax.set_xticklabels(["10", "30", "100", "300"])
    ax.set_xlabel("wall time for caproic acid / min"), ax.set_ylabel("MAE of the VPT2 fundamentals / cm$^{-1}$")
    ax.legend(frameon=False, fontsize=style.ANNOT_FS, loc="upper right")
    fig.tight_layout()
    style.save(fig, out, "cost_accuracy")


if __name__ == "__main__":
    out = Path(sys.argv[1])
    style.apply_style()
    band_types(out), cost_accuracy(out)
    gallery(out.parent / "results_appendix")
