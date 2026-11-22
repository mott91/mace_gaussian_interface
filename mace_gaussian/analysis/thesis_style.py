"""Shared figure style for the thesis (Matplotlib).

Adopted from the user's ``figstyle.py`` (PR Numerical Methods protocol) so every
data figure in the thesis forms one visual system:

* text is rendered by LaTeX in Computer Modern, the document's own face;
* figures are authored at their true display width, so a PDF placed at
  ``width=\\textwidth`` shows its text at the real point size;
* colours come from one categorical palette, assigned in fixed order and never
  cycled. The five energy models each own one slot; DFT is ink, experiment is
  muted ink, dipole models reuse the first three slots inside a facet that is
  already labelled by energy model.

Geometry follows ``thesis/latex/main.tex``: A4, 2.5 cm margins, 6 mm binding
offset, 11 pt body. Use ``figwidth(0.75)`` for a figure included at
``width=0.75\\textwidth``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- geometry -------------------------------------------------------------
TEXTWIDTH = (21.0 - 2 * 2.5 - 0.6) / 2.54  # inches (A4 minus margins minus binding offset)


def figwidth(frac: float = 1.0) -> float:
    return TEXTWIDTH * frac


# --- palette --------------------------------------------------------------
# Categorical slots 1-5 (colourblind-validated on a light surface). Assign in
# order; never cycle, never reorder.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK = "#1a1a1a"  # reference curves and label text
INK_MUTED = "#6b6b6b"  # axes, ticks, guide lines, experiment
GRID = "#d8d8d8"

# Fixed slot per energy model (same assignment as analysis/palette.py, so the
# HTML report and the thesis agree).
MODEL_SLOT = {"mace_omol": 0, "mace_mp": 1, "mace_off": 2, "mace_polar": 3, "mace_anicc": 4}
MODEL_LABEL = {
    "mace_omol": "MACE-OMOL",
    "mace_mp": "MACE-MP",
    "mace_off": "MACE-OFF",
    "mace_polar": "MACE-POLAR",
    "mace_anicc": "MACE-ANI-cc",
}
# mace_mdp took espaloma's place in the 2026 campaign and shares its slot; the two never
# appear in one figure (espaloma only exists in the older results).
DIPOLE_SLOT = {"mace_ml": 0, "mace_polar1": 1, "mace_mdp": 2, "espaloma": 2}
DIPOLE_LABEL = {
    "mace_ml": "MACE4IR",
    "mace_polar1": "POLAR-1",
    "mace_mdp": "MDP",
    "espaloma": "espaloma",
}


# MACE-POLAR-1 in one red, as energy model and as dipole model (figures-ch4-final, item 1a):
# the yellow-orange of slot 4 was hard to see in print. Chosen from six candidate reds as the
# one with the largest smallest CAM02-UCS distance to the other colours in use (OFF23 green,
# ANICC pink, OMOL blue, B3LYP ink, muted ink) under simulated deuteranopia and protanopia
# (colorspacious and daltonlens, Machado 2009, severity 1): smallest distance to a model colour
# or the ink 22.6 (green, deuteranopia), to the muted grey of the experiment 19.7 (protanopia);
# lightness J' 44 against 65 for the green. analysis/results_figures/cvd_check.py
POLAR_RED = "#a93226"
COLOR_OVERRIDE = {"mace_polar": POLAR_RED, "mace_polar1": POLAR_RED}


def model_color(energy_model: str) -> str:
    return COLOR_OVERRIDE.get(energy_model) or SERIES[MODEL_SLOT.get(energy_model, 0)]


def dipole_color(dipole_model: str) -> str:
    return COLOR_OVERRIDE.get(dipole_model) or SERIES[DIPOLE_SLOT.get(dipole_model, 0)]


BASE_FS = 9.5  # against the document's 11 pt body
ANNOT_FS = 8.5


def apply_style() -> None:
    plt.rcParams.update(
        {
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "font.family": "serif",
            "font.size": BASE_FS,
            "axes.titlesize": BASE_FS,
            "axes.labelsize": BASE_FS,
            "axes.edgecolor": INK_MUTED,
            "axes.linewidth": 0.8,
            "axes.labelcolor": INK,
            "xtick.color": INK_MUTED,
            "ytick.color": INK_MUTED,
            "xtick.labelcolor": INK,
            "ytick.labelcolor": INK,
            "xtick.labelsize": BASE_FS,
            "ytick.labelsize": BASE_FS,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "legend.fontsize": BASE_FS,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "grid.alpha": 0.7,
            "lines.linewidth": 1.6,
            "savefig.facecolor": "white",
            "savefig.dpi": 320,
            "pdf.fonttype": 42,
        }
    )


def clean_axes(ax, grid: bool = True):
    """Recessive grid behind the marks, no top/right spines."""
    ax.set_axisbelow(True)
    ax.grid(grid)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return ax


def tex(s: str) -> str:
    """Escape a plain string for LaTeX text rendering."""
    return (
        str(s)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("#", r"\#")
    )


def tight(fig, **kw) -> None:
    """tight_layout as a layout engine, so that it runs again whenever ``save`` changes the width
    (``fig.tight_layout`` runs once and leaves no engine behind). Same keywords; the outer pad
    defaults to 0.4 font sizes instead of 1.08: ``save`` crops the outer margin anyway, and the
    pad also sets the gap to a shared x label (supxlabel)."""
    kw.setdefault("pad", 0.4)
    fig.set_layout_engine("tight", **kw)


PAD_PT = 1.0  # blank margin around the cropped figure, in PostScript points (bp)


def _cropped_width_pt(fig) -> float:
    fig.canvas.draw()  # runs the layout engine, if any
    return fig.get_tightbbox(fig.canvas.get_renderer()).width * 72 + 2 * PAD_PT


def fit_width(fig, frac: float = 1.0, tol_pt: float = 0.2) -> float:
    """Change the figure width until the tightly cropped PDF (pad PAD_PT) is ``frac`` text widths
    wide, so that LaTeX places it at ``width=frac\\textwidth`` without rescaling and the font
    sizes are the true ones. Secant steps: with a layout engine the cropped width grows one to
    one with the figure width, without one (fixed axes fractions) more slowly. Returns the
    cropped width in bp."""
    target = TEXTWIDTH * frac * 72
    w0, c0 = fig.get_figwidth(), _cropped_width_pt(fig)
    w1 = w0 + (target - c0) / 72
    for _ in range(12):
        if abs(c0 - target) < tol_pt:
            break
        fig.set_figwidth(w1)
        c1 = _cropped_width_pt(fig)
        slope = (c1 - c0) / ((w1 - w0) * 72) if w1 != w0 else 1.0
        w0, c0 = w1, c1
        step = (target - c0) / 72 / (slope if slope > 0.2 else 1.0)
        w1 = w0 + float(np.clip(step, -0.3 * w0, 0.3 * w0))
    assert abs(c0 - target) < tol_pt, (c0, target)
    return c0


def save(fig, out_dir: Path, name: str, png_dpi: int = 200, frac: float = 1.0) -> tuple[Path, Path]:
    """Write ``name.pdf`` (for LaTeX) and ``name.png`` (for the review gallery), both cropped to the
    drawn content plus PAD_PT, the PDF ``frac`` text widths wide (``fit_width``); ``frac=None``
    keeps the authored width."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = out_dir / f"{name}.pdf"
    png = out_dir / f"{name}.png"
    if frac is not None:
        fit_width(fig, frac)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=PAD_PT / 72)
    fig.savefig(png, dpi=png_dpi, bbox_inches="tight", pad_inches=PAD_PT / 72)
    plt.close(fig)
    return pdf, png
