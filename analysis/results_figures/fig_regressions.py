"""Regression plots (MACE value against B3LYP value, parity line) for the results chapter.

  regression_harmonic        harmonic frequencies, one panel per energy model (MAE and R2)
  regression_vpt2            VPT2 fundamentals, one panel per energy model (MAE and R2)
  regression_anharmonicity   anharmonic shift (top) and chi_ii (bottom), one column per energy model
  regression_shift           the top row of regression_anharmonicity on its own
  regression_xii             the bottom row of regression_anharmonicity on its own
  regression_intensities     harmonic intensities, one panel per dipole model (pooled over energy models)

Same band sets as the tables of analysis/followup3/thesis_stats.py; the MAE in each panel must equal
the table value.   python analysis/results_figures/fig_regressions.py [output_dir]
"""

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ts = _load("thesis_stats", ROOT / "analysis" / "followup3" / "thesis_stats.py")
style = _load("thesis_style", ROOT / "mace_gaussian" / "analysis" / "thesis_style.py")

ELAB = {"mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23", "mace_polar": "MACE-POLAR-1"}
DLAB = {"mace_polar1": "MACE-POLAR-1", "mace_mdp": "MACE-MDP", "mace_ml": "MACE4IR"}
DOT = dict(ms=2.6, mec="none", alpha=0.6)
# harmonic MAE order of Table 4.1 (OMOL and POLAR tie at 12.5 and keep the chapter order)
H_ORDER = ["mace_anicc", "mace_omol", "mace_polar", "mace_off"]
# R2 decimals in the two frequency plots: two give 1.00 and three give 1.000 for three of the four
# harmonic panels (R2 0.9995-0.9998), so four, in both plots (report figures-4-2)
FREQ_R2_DECIMALS = 4


def note(ax, text):
    # white backing: points must not run through the label (POLAR-1 shifts reach +330 cm-1)
    ax.text(0.05, 0.94, text, transform=ax.transAxes, ha="left", va="top", fontsize=style.ANNOT_FS, color=style.INK,
            zorder=5, bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.8))  # fmt: skip


def num(v: float, nd: int) -> str:
    """Number for a panel label, with a typographic minus."""
    return f"{v:.{nd}f}".replace("-", "$-$")


def parity_row(axes, frames, xcol, ycol, lim, ticks, r2: int | None = None):
    """r2: decimals of the R2 (against y = x) under the MAE label; None = no R2."""
    for ax, (model, d) in zip(axes, frames):
        x, y = d[xcol], d[ycol]
        assert x.between(*lim).all() and y.between(*lim).all(), (model, xcol, x.min(), x.max(), y.min(), y.max())
        ax.plot(lim, lim, color=style.INK_MUTED, lw=0.8, ls="--", zorder=1)
        ax.plot(x, y, "o", color=style.model_color(model), zorder=2, **DOT)
        style.clean_axes(ax)
        mae = (y - x).abs().mean()
        note(ax, f"MAE {mae:.1f}" + (f"\n$R^2$ {num(ts.r2(x, y), r2)}" if r2 is not None else ""))
        ax.set_xlim(*lim), ax.set_ylim(*lim), ax.set_xticks(ticks), ax.set_yticks(ticks)
        ax.set_aspect("equal")
        print(f"{ycol:14s} {model:11s} n {len(d):4d} MAE {mae:5.1f} R2 {ts.r2(x, y):.4f}")


def titles(axes, models):
    for ax, m in zip(axes, models):
        ax.set_title(style.tex(ELAB[m]), color=style.model_color(m), fontsize=style.BASE_FS)


def main(out: Path) -> None:
    _, chi, fund, harm, _ = ts.load()
    E = ts.E_MODELS
    h = harm[harm.dipole_model == "mace_ml"]
    f = fund[(fund.dipole_model == "mace_ml") & ~fund.molecule.isin(ts.SYMTOP)].copy()
    f["dft_shift"], f["ml_shift"] = f.dft_anharmonic - f.dft_harmonic, f.ml_anharmonic - f.ml_harmonic
    x = chi[(chi.kind == "diagonal") & ~chi.molecule.isin(ts.SYMTOP)].dropna(subset=["dft_xmatrix", "ml_xmatrix"])
    style.apply_style()
    W = style.figwidth()

    # harmonic panels follow the harmonic MAE (Section 4.1), everything else the fixed chapter order
    for name, df, xc, yc, lab, order in (("regression_harmonic", h, "dft_harmonic", "ml_harmonic", "harmonic frequency", H_ORDER),
                                         ("regression_vpt2", f, "dft_anharmonic", "ml_anharmonic", "VPT2 fundamental", E)):
        fig, axes = plt.subplots(1, 4, figsize=(W, 2.15), sharex=True, sharey=True)
        parity_row(axes, [(m, df[df.energy_model == m]) for m in order], xc, yc, (0, 4200), [0, 2000, 4000],
                   r2=FREQ_R2_DECIMALS)
        titles(axes, order)
        axes[0].set_ylabel("MACE / cm$^{-1}$")
        fig.supxlabel(f"B3LYP {lab} / cm$^{{-1}}$", fontsize=style.BASE_FS)
        fig.tight_layout(w_pad=0.6)
        style.save(fig, out, name)

    both = pd.concat([x.dft_xmatrix, x.ml_xmatrix])
    clim = (np.floor(both.min() / 25) * 25 - 5, np.ceil(both.max() / 25) * 25 + 5)
    shift_row = ([(m, f[f.energy_model == m]) for m in E], "dft_shift", "ml_shift", (-480, 480), [-400, 0, 400])
    xii_row = ([(m, x[x.energy_model == m]) for m in E], "dft_xmatrix", "ml_xmatrix", clim,
               [t for t in (-200, -100, 0, 100) if clim[0] < t < clim[1]])  # fmt: skip
    # the two rows on their own (Section 4.2 cites them in separate paragraphs); same ranges
    for name, row, ylab, xlab in (("regression_shift", shift_row, "MACE shift / cm$^{-1}$", "B3LYP anharmonic shift / cm$^{-1}$"),
                                  ("regression_xii", xii_row, r"MACE $\chi_{ii}$ / cm$^{-1}$", r"B3LYP $\chi_{ii}$ / cm$^{-1}$")):  # fmt: skip
        fig, axes = plt.subplots(1, 4, figsize=(W, 2.15), sharex=True, sharey=True)
        parity_row(axes, *row, r2=2)
        titles(axes, E)
        axes[0].set_ylabel(ylab)
        fig.supxlabel(xlab, fontsize=style.BASE_FS)
        fig.tight_layout(w_pad=0.6)
        style.save(fig, out, name)

    fig, axes = plt.subplots(2, 4, figsize=(W, 3.9), sharex="row", sharey="row")
    parity_row(axes[0], *shift_row, r2=2)
    parity_row(axes[1], *xii_row, r2=2)
    titles(axes[0], E)
    axes[0, 0].set_ylabel("MACE shift / cm$^{-1}$"), axes[1, 0].set_ylabel(r"MACE $\chi_{ii}$ / cm$^{-1}$")
    for ax, lab in ((axes[0, 1], "B3LYP anharmonic shift / cm$^{-1}$"), (axes[1, 1], r"B3LYP $\chi_{ii}$ / cm$^{-1}$")):
        ax.set_xlabel(lab, x=1.1)
    fig.tight_layout(w_pad=0.6, h_pad=1.0)
    style.save(fig, out, "regression_anharmonicity")

    vis = harm[harm.dft_int >= ts.INT_MIN]
    lo, hi = 10 ** np.floor(np.log10(vis.ml_int.min())), 10 ** np.ceil(np.log10(max(vis.ml_int.max(), vis.dft_int.max())))
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.6), sharex=True, sharey=True)
    for ax, dm in zip(axes, ts.D_MODELS):
        d = vis[vis.dipole_model == dm]
        c = style.dipole_color(dm)
        xs = np.array([1, hi])
        ax.fill_between(xs, xs / 2, xs * 2, color=style.GRID, alpha=0.5, lw=0, zorder=0)
        ax.plot(xs, xs, color=style.INK_MUTED, lw=0.8, ls="--", zorder=1)
        ax.plot(d.dft_int, d.ml_int, "o", color=c, zorder=2, ms=2.2, mec="none", alpha=0.45)
        style.clean_axes(ax)
        within = ((d.ml_int / d.dft_int).between(0.5, 2)).mean() * 100
        note(ax, f"MAE {(d.ml_int - d.dft_int).abs().mean():.1f}\nwithin $2\\times$: {within:.0f} \\%")
        ax.set_xscale("log"), ax.set_yscale("log"), ax.set_xlim(1, hi), ax.set_ylim(lo, hi)
        ax.set_title(style.tex(DLAB[dm]), color=c, fontsize=style.BASE_FS)
        print(f"int harmonic   {dm:11s} n {len(d):4d} MAE {(d.ml_int - d.dft_int).abs().mean():5.1f} within2x {within:.1f}  ml min {d.ml_int.min():.3g}")
    axes[0].set_ylabel("MACE intensity / km mol$^{-1}$")
    axes[1].set_xlabel("B3LYP harmonic intensity / km mol$^{-1}$")
    fig.tight_layout(w_pad=0.6)
    style.save(fig, out, "regression_intensities")

    # checks for the prose of 4.1
    print("\nSO2 harmonic:")
    print(h[h.molecule == "sulfur_dioxide"][["energy_model", "dft_harmonic", "ml_harmonic"]].round(0).to_string(index=False))
    r = h[(h.region == "1600-2700")].assign(err=lambda d: d.ml_harmonic - d.dft_harmonic)
    print("\n1600-2700 harmonic, |err| > 40:")
    print(r[r.err.abs() > 40][["molecule", "energy_model", "dft_harmonic", "err"]].round(0).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out")
