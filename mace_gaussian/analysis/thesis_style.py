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
DIPOLE_SLOT = {"mace_ml": 0, "mace_polar1": 1, "espaloma": 2}
DIPOLE_LABEL = {"mace_ml": "MACE4IR", "mace_polar1": "POLAR-1", "espaloma": "espaloma"}


def model_color(energy_model: str) -> str:
    return SERIES[MODEL_SLOT.get(energy_model, 0)]


def dipole_color(dipole_model: str) -> str:
    return SERIES[DIPOLE_SLOT.get(dipole_model, 0)]


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


def save(fig, out_dir: Path, name: str, png_dpi: int = 200) -> tuple[Path, Path]:
    """Write ``name.pdf`` (for LaTeX) and ``name.png`` (for the review gallery)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = out_dir / f"{name}.pdf"
    png = out_dir / f"{name}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=png_dpi)
    plt.close(fig)
    return pdf, png
