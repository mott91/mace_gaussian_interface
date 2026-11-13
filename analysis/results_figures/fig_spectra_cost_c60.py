"""Stacked spectra (methanol, glycine, acetic acid), cost scaling and the C60 spectrum for the results chapter.

The stacked spectra reuse fig_stacked_spectra of scripts/make_thesis_figures.py on the 2026 campaign,
with the model names of the thesis text and MACE-POLAR-1 as the dipole model.
Run from a scratch directory (importing mace_gaussian writes .model.pt into the cwd):
    python analysis/results_figures/fig_spectra_cost_c60.py [output_dir]
"""

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("make_thesis_figures", ROOT / "scripts" / "make_thesis_figures.py")
mtf = importlib.util.module_from_spec(spec)
sys.modules["make_thesis_figures"] = mtf
spec.loader.exec_module(mtf)
import mace_gaussian.analysis.thesis_style as style  # noqa: E402

CAMP = ROOT / "campaigns" / "2026"
ELAB = {"mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23", "mace_polar": "MACE-POLAR-1"}
E = list(ELAB)
ACIDS = ["formic_acid_trans", "acetic_acid", "propionic_acid", "butyric_acid", "valeric_acid", "caproic_acid"]


def find(d, key):
    if isinstance(d, dict):
        if d.get(key) is not None:
            return d[key]
        for v in d.values():
            r = find(v, key)
            if r is not None:
                return r
    return None


def spectra(out: Path) -> None:
    style.MODEL_LABEL.update(ELAB)
    mtf.ANALYSIS, mtf.COMPARISON, mtf.OUT = CAMP / "analysis_results", CAMP / "comparison_results", out
    mtf.DIPOLE_PREF[:] = ["mace_polar1", "mace_mdp", "mace_ml"]
    mtf.MODEL_ORDER[:] = ["mace_omol", "mace_off", "mace_anicc", "mace_polar"]
    for mol in ("methanol", "glycine", "acetic_acid"):
        reps = mtf.representative_runs(mtf.load_report(mol))
        print(mol, {k: v["name"] for k, v in reps.items()})
        mtf.fig_stacked_spectra(mol)


def cost(out: Path) -> None:
    man = json.loads((CAMP / "comparison_results" / "batch_manifest.json").read_text())["molecules"]
    pts = {m: [] for m in ["b3lyp", *E]}
    for mol, info in man.items():
        n = int(Path(info["xyz_path"]).read_text().splitlines()[0])
        b = next((CAMP / "comparison_results" / mol).glob("b3lyp*/results.json"))
        t = find(json.loads(b.read_text()), "total_elapsed_s")
        if t:
            pts["b3lyp"].append((mol, n, float(t)))
        for e in E:
            ts_ = [c["runtime_s"] for k, c in info["combinations"].items()
                   if k.startswith(e + "_") and c.get("status") == "complete" and c.get("runtime_s")]
            if len(ts_) == 3:
                pts[e].append((mol, n, float(np.mean(ts_))))
    fig, ax = plt.subplots(figsize=(style.figwidth(0.8), 3.4))
    style.clean_axes(ax)
    for key, rows in pts.items():
        color = style.INK if key == "b3lyp" else style.model_color(key)
        acid = sorted((n, t) for m, n, t in rows if m in ACIDS)
        rest = [(n, t) for m, n, t in rows if m not in ACIDS]
        slope = np.polyfit(np.log([a[0] for a in acid]), np.log([a[1] for a in acid]), 1)[0]
        ax.plot(*zip(*rest), "o", ms=3.6, mfc="none", mec=color, mew=0.8, alpha=0.8)
        ax.plot(*zip(*acid), "o-", ms=3.8, lw=1.0, color=color, label=f"{'B3LYP' if key == 'b3lyp' else ELAB[key]} ($N^{{{slope:.1f}}}$)")
        print(f"{key:11s} exponent {slope:.2f}  caproic {dict((m, t) for m, n, t in rows)['caproic_acid'] / 60:.0f} min  n {len(rows)}")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xticks([3, 5, 10, 20]), ax.set_xticklabels(["3", "5", "10", "20"])
    ax.set_xlabel(r"number of atoms $N$"), ax.set_ylabel("wall time / s")
    ax.legend(frameon=False, fontsize=style.ANNOT_FS, loc="upper left")
    fig.tight_layout()
    style.save(fig, out, "cost_scaling")


def c60(out: Path) -> None:
    r = json.loads((ROOT / "campaigns/c60/comparison_results/c60/mace_omol_mace_polar1/results.json").read_text())["frequencies"]
    grid = np.arange(400.0, 1700.0, 0.5)
    hf, hi = zip(*[(b["freq_cm"], b["ir_intensity"]) for b in r["harmonic"]])
    bands = ([(b["freq_cm"], b["ir_intensity"]) for b in r["anharmonic"]]
             + [(b["freq_anharmonic"], b["ir_intensity"]) for b in r["overtones"]]
             + [(b["freq_anharmonic"], b["ir_intensity"]) for b in r["combination_bands"]])
    bands = [(f, i) for f, i in bands if f is not None and 380 < f < 1720]
    vf, vi = zip(*bands)
    fund = sum(b["ir_intensity"] for b in r["anharmonic"] if 400 < b["freq_cm"] < 1700)
    print(f"C60 window 400-1700: VPT2 bands {len(bands)}, summed intensity {sum(vi):.0f}, of which fundamentals {fund:.0f}")
    yh, ph = mtf.broaden(list(hf), list(hi), grid)
    yv, _ = mtf.broaden(list(vf), list(vi), grid, scale=ph)
    fig, ax = plt.subplots(figsize=(style.figwidth(), 2.9))
    style.clean_axes(ax, grid=False)
    off = 1.15 * max(1.0, float(np.max(yh)))
    for x in (527, 576, 1183, 1429):
        ax.axvline(x, color=style.INK_MUTED, lw=0.7, ls=":")
    c = style.model_color("mace_omol")
    ax.plot(grid, yh, color=style.INK_MUTED, lw=1.2), ax.plot(grid, yv + off, color=c, lw=1.2)
    ax.text(1690, 0.35, "harmonic", ha="right", fontsize=style.ANNOT_FS, color=style.INK_MUTED)
    ax.text(1690, off + 0.35, "VPT2", ha="right", fontsize=style.ANNOT_FS, color=c)
    ax.set_xlim(400, 1700), ax.set_yticks([0, off]), ax.set_yticklabels(["0", ""])
    ax.set_xlabel(r"wavenumber $\tilde{\nu}$ / cm$^{-1}$"), ax.set_ylabel("absorbance (harmonic peak $=1$)")
    fig.tight_layout()
    style.save(fig, out, "c60_spectrum")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out"
    style.apply_style()
    spectra(out), cost(out), c60(out)
