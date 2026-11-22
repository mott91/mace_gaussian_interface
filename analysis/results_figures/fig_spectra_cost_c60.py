"""Stacked spectra (caproic acid = opening figure of Chapter 4, methanol, glycine, acetic acid) and the C60 spectrum for the results chapter (the cost figures: fig_cost.py).

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
# stacked spectra, from the bottom up: experiment, B3LYP (VPT2), then the fixed chapter order, so that
# MACE-ANICC sits directly above the reference (decided 2026-11-13); fig_stacked_spectra draws its
# MODEL_ORDER from the bottom up above B3LYP
SPECTRA_ORDER = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]
EXPERIMENT_BANDS = ROOT / "analysis" / "followup3" / "out" / "experiment_bands.csv"


def experiment_lines() -> dict[str, list[tuple[str, float]]]:
    """Dotted lines = exactly the bands of the experiment comparison (thesis_stats.section_experiment,
    'with NH3/PH3' set): excluded bands and CO2 nu1 (Fermi dyad, dropped there as co2_nu1) get none;
    a degenerate band is one row."""
    import pandas as pd

    e = pd.read_csv(EXPERIMENT_BANDS)
    e = e[~e.assignment.astype(str).str.startswith("excluded")]
    e = e[~((e.molecule == "carbon_dioxide") & (e.label == "ν1"))]
    return {m: list(zip(g.label, g.exp.astype(float))) for m, g in e.groupby("molecule")}


# dotted lines of the C60 figure: the four T1u values of Table 4.14 (room-temperature film values).
# The para-H2 matrix values of wakabayashi2019c60 are 529.77, 578.24, 1184.7, 1432 (followup-1 report);
# change here if the table switches (figures-ch4-final report, item 4)
C60_EXP = (527, 576, 1183, 1429)


def spectra(out: Path) -> None:
    style.MODEL_LABEL.update(ELAB)
    mtf.ANALYSIS, mtf.COMPARISON, mtf.OUT = CAMP / "analysis_results", CAMP / "comparison_results", out
    mtf.DIPOLE_PREF[:] = ["mace_polar1", "mace_mdp", "mace_ml"]
    mtf.MODEL_ORDER[:] = SPECTRA_ORDER
    mtf.BAND_LINES = experiment_lines()
    for mol in ("caproic_acid", "methanol", "glycine", "acetic_acid"):
        reps = mtf.representative_runs(mtf.load_report(mol))
        print(mol, {k: v["name"] for k, v in reps.items()})
        mtf.fig_stacked_spectra(mol)


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
    for x in C60_EXP:
        ax.axvline(x, color=style.INK_MUTED, lw=0.7, ls=":")
    c = style.model_color("mace_omol")
    ax.plot(grid, yh, color=style.INK_MUTED, lw=1.2), ax.plot(grid, yv + off, color=c, lw=1.2)
    ax.text(1690, 0.35, "harmonic", ha="right", fontsize=style.ANNOT_FS, color=style.INK_MUTED)
    ax.text(1690, off + 0.35, "VPT2", ha="right", fontsize=style.ANNOT_FS, color=c)
    ax.set_xlim(400, 1700), ax.set_yticks([0, off]), ax.set_yticklabels(["0", ""])
    ax.set_xlabel(r"wavenumber $\tilde{\nu}$ / cm$^{-1}$"), ax.set_ylabel("relative absorbance")
    style.tight(fig)
    style.save(fig, out, "c60_spectrum")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out"
    style.apply_style()
    spectra(out), c60(out)
