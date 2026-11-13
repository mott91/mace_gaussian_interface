"""Residuals of the fundamentals against B3LYP.

  residuals_fundamentals   harmonic (top) and VPT2 (bottom), one column per energy model, fixed chapter order
  residuals_harmonic       harmonic row alone with its own y-range, panels in the harmonic MAE order (Section 4.1)

Uses the band tables of analysis/followup3/thesis_stats.py, so the points are the ones behind the
pooled tables of the results chapter (harmonic with NH3/PH3, VPT2 without). The MAE printed in
each panel is recomputed here and must equal the table value.

Run from any directory:  python analysis/results_figures/fig_residuals.py [output_dir]
"""

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]


def _load(name: str, path: Path):
    # by file path: importing the mace_gaussian package would write a .model.pt into the cwd
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ts = _load("thesis_stats", ROOT / "analysis" / "followup3" / "thesis_stats.py")
style = _load("thesis_style", ROOT / "mace_gaussian" / "analysis" / "thesis_style.py")

LABEL = {"mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23",
         "mace_polar": "MACE-POLAR-1"}  # fmt: skip
EDGES = [500, 1600, 2700, 3200]
YLIM = (-330, 500)
H_ORDER = ["mace_anicc", "mace_omol", "mace_polar", "mace_off"]  # harmonic MAE order of Table 4.1
H_YLIM = (-130, 200)  # data: -115.3 (PH3, MACE-OFF23) to +176.9 (P-H stretch of PH3, MACE-OFF23)


def main(out_dir: Path) -> None:
    _, _, fund, harm, _ = ts.load()
    f = fund[(fund.dipole_model == "mace_ml") & ~fund.molecule.isin(ts.SYMTOP)]
    h = harm[harm.dipole_model == "mace_ml"]
    rows = [("harmonic", h, "ml_harmonic", "dft_harmonic"), ("VPT2", f, "ml_anharmonic", "dft_anharmonic")]

    style.apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(style.figwidth(), 3.7), sharex=True, sharey=True)
    for r, (level, df, ml, ref) in enumerate(rows):
        for c, model in enumerate(ts.E_MODELS):
            ax, d = axes[r, c], df[df.energy_model == model]
            err = d[ml] - d[ref]
            assert err.min() > YLIM[0] and err.max() < YLIM[1], (level, model, err.min(), err.max())
            for e in EDGES:
                ax.axvline(e, color=style.GRID, lw=0.6, zorder=0)
            ax.axhline(0, color=style.INK_MUTED, lw=0.8, ls="--", zorder=1)
            ax.plot(d.dft_harmonic, err, "o", ms=2.6, mec="none", alpha=0.6, color=style.model_color(model), zorder=2)
            style.clean_axes(ax, grid=False)
            ax.text(0.04, 0.93, f"MAE {err.abs().mean():.1f}", transform=ax.transAxes, ha="left", va="top",
                    fontsize=style.ANNOT_FS, color=style.INK)  # fmt: skip
            if r == 0:
                ax.set_title(style.tex(LABEL[model]), color=style.model_color(model), fontsize=style.BASE_FS)
            print(f"{level:9s} {model:11s} n {len(d):4d}  MAE {err.abs().mean():5.1f}  min {err.min():7.1f}  max {err.max():6.1f}")
        axes[r, 0].set_ylabel(f"{level} error / cm$^{{-1}}$")
    axes[0, 0].set_ylim(*YLIM)
    axes[0, 0].set_xlim(0, 4000)
    axes[0, 0].set_xticks([0, 1500, 3000])
    axes[0, 0].set_yticks([-250, 0, 250, 500])
    fig.supxlabel("B3LYP harmonic frequency / cm$^{-1}$", fontsize=style.BASE_FS)
    fig.tight_layout(w_pad=0.6, h_pad=0.8)
    style.save(fig, out_dir, "residuals_fundamentals")



def harmonic(out_dir: Path) -> None:
    _, _, _, harm, _ = ts.load()
    h = harm[harm.dipole_model == "mace_ml"]
    style.apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(style.figwidth(), 2.3), sharex=True, sharey=True)
    for ax, model in zip(axes, H_ORDER):
        d = h[h.energy_model == model]
        err = d.ml_harmonic - d.dft_harmonic
        assert err.min() > H_YLIM[0] and err.max() < H_YLIM[1], (model, err.min(), err.max())
        for e in EDGES:
            ax.axvline(e, color=style.GRID, lw=0.6, zorder=0)
        ax.axhline(0, color=style.INK_MUTED, lw=0.8, ls="--", zorder=1)
        ax.plot(d.dft_harmonic, err, "o", ms=2.6, mec="none", alpha=0.6, color=style.model_color(model), zorder=2)
        style.clean_axes(ax, grid=False)
        ax.text(0.04, 0.95, f"MAE {err.abs().mean():.1f}", transform=ax.transAxes, ha="left", va="top",
                fontsize=style.ANNOT_FS, color=style.INK)  # fmt: skip
        ax.set_title(style.tex(LABEL[model]), color=style.model_color(model), fontsize=style.BASE_FS)
        print(f"harmonic  {model:11s} n {len(d):4d}  MAE {err.abs().mean():5.1f}  min {err.min():7.1f}  max {err.max():6.1f}")
    axes[0].set_ylabel("harmonic error / cm$^{-1}$")
    axes[0].set_ylim(*H_YLIM)
    axes[0].set_xlim(0, 4000)
    axes[0].set_xticks([0, 1500, 3000])
    axes[0].set_yticks([-100, 0, 100, 200])
    fig.supxlabel("B3LYP harmonic frequency / cm$^{-1}$", fontsize=style.BASE_FS)
    fig.tight_layout(w_pad=0.6)
    style.save(fig, out_dir, "residuals_harmonic")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out"
    main(out)
    harmonic(out)
