"""Two figures that replace tables of the results chapter.

  molecules_vpt2       VPT2 MAE of every molecule and energy model, ordered by group
                       (analysis/followup3/out/freq_vpt2_per_molecule.csv)
  molecules_harmonic   the same at the harmonic level, NH3 and PH3 included, legend in the harmonic
                       MAE order (analysis/followup3/out/freq_harmonic_per_molecule.csv)
  time_split       share of the wall time per component, by molecule size and energy model
                   (analysis/ch5_timing/timing_runs.csv, time-weighted as in timing_breakdown.py)

python analysis/results_figures/fig_tables.py [output_dir]
"""

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("thesis_style", ROOT / "mace_gaussian" / "analysis" / "thesis_style.py")
style = importlib.util.module_from_spec(spec)
sys.modules["thesis_style"] = style
spec.loader.exec_module(style)

ELAB = {"mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23", "mace_polar": "MACE-POLAR-1"}
COL = {"ANI-cc": "mace_anicc", "OMOL": "mace_omol", "OFF": "mace_off", "POLAR": "mace_polar"}
GROUPS = [
    ("acids", [("formic_acid_trans", "t-formic"), ("formic_acid_cis", "c-formic"), ("acetic_acid", "acetic"),
               ("propionic_acid", "propionic"), ("butyric_acid", "butyric"), ("valeric_acid", "valeric"),
               ("caproic_acid", "caproic")]),
    ("amino acid", [("glycine", "glycine")]),
    ("alcohols", [("methanol", "methanol"), ("ethanol", "ethanol"), ("1-propanol", "1-propanol"),
                  ("isopropanol", "isopropanol"), ("1-butanol", "1-butanol"), ("tert-butanol", "t-butanol")]),
    ("small molecules", [("carbon_dioxide", "CO$_2$"), ("sulfur_dioxide", "SO$_2$"), ("hydrogen_peroxide", "H$_2$O$_2$"),
                         ("nitric_acid", "HNO$_3$"), ("carbonic_acid", "H$_2$CO$_3$")]),
    ("sugar", [("beta_glucose", r"$\beta$-glucose")]),
]
# harmonic level: NH3 and PH3 included, small molecules in the order of Table 3.3
GROUPS_H = [(g, ms) if g != "small molecules" else
            (g, [ms[0], ms[1], ("ammonia", "NH$_3$"), ("phosphine", "PH$_3$"), *ms[2:]]) for g, ms in GROUPS]  # fmt: skip
H_ORDER = ["mace_anicc", "mace_omol", "mace_polar", "mace_off"]  # harmonic MAE order of Table 4.1
# molecule names and group labels in ink, about 1 pt larger than before (figures-ch4-final)
NAME_FS, GROUP_FS = style.BASE_FS + 1.0, style.ANNOT_FS + 1.0


def molecules(out: Path, level: str = "vpt2") -> None:
    groups, models = (GROUPS, list(COL.values())) if level == "vpt2" else (GROUPS_H, H_ORDER)
    key = {m: k for k, m in COL.items()}
    df = pd.read_csv(ROOT / "analysis" / "followup3" / "out" / f"freq_{level}_per_molecule.csv", index_col=0)
    order = [m for _, ms in groups for m, _ in ms]
    assert set(order) == set(df.index), set(order) ^ set(df.index)
    fig, ax = plt.subplots(figsize=(style.figwidth(), 3.1))
    x = range(len(order))
    for model in models:
        y = pd.to_numeric(df.loc[order, key[model]], errors="coerce")
        assert y.max() < 115, (model, y.max())
        ax.plot(x, y, "o", ms=4.2, mec="none", alpha=0.85, color=style.model_color(model), label=ELAB[model])
        print(f"{level} {model:11s} {y.notna().sum()} molecules, max {y.max():.1f}")
    style.clean_axes(ax)
    pos = 0
    for name, ms in groups:
        if pos:
            ax.axvline(pos - 0.5, color=style.INK_MUTED, lw=0.6)
        ax.text(pos + (len(ms) - 1) / 2, 1.02, name if len(ms) > 1 else name.replace(" ", "\n"), transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=GROUP_FS, color=style.INK)
        pos += len(ms)
    ax.set_xticks(list(x), [lab for _, ms in groups for _, lab in ms], rotation=60, ha="right", rotation_mode="anchor",
                  fontsize=NAME_FS, color=style.INK)
    ax.set_xlim(-0.6, len(order) - 0.4), ax.set_ylim(0, 115)
    ax.set_ylabel(f"{'VPT2' if level == 'vpt2' else 'harmonic'} MAE / cm$^{{-1}}$")
    # legend directly above the group labels; borderpad 0: no blank band above the crop
    fig.legend(loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0), handletextpad=0.2, columnspacing=1.2,
               borderpad=0)  # fmt: skip
    style.tight(fig, pad=0.4)
    style.save(fig, out, f"molecules_{level}")


DIPOLE_GREY = "#9e9e9e"  # medium grey between INK_MUTED (Gaussian) and GRID (Python interface)


class _MultiPatch:
    """Legend key made of the four energy-model colours side by side."""

    def __init__(self, label: str):
        self._label = label

    def get_label(self) -> str:
        return self._label


class _MultiHandler:
    def legend_artist(self, legend, orig, fontsize, box):
        x0, y0, w, h = box.xdescent, box.ydescent, box.width, box.height
        n = len(ELAB)
        for i, m in enumerate(ELAB):
            box.add_artist(Rectangle((x0 + i * w / n, y0), w / n, h, color=style.model_color(m), lw=0,
                                     transform=box.get_transform()))  # fmt: skip
        return box


def time_split(out: Path) -> None:
    df = pd.read_csv(ROOT / "analysis" / "ch5_timing" / "timing_runs.csv")
    labels = ["3-4", "5-6", "8-12", "14-15", "17-24"]
    df["size"] = pd.cut(df.atoms, [0, 4, 6, 12, 15, 30], labels=labels)
    # energy model in the colour of the panel's own model (None below); dipole model, Gaussian and the Python
    # interface in three greys (figures-ch4-tweaks: the old blue and green repeated the OMOL-0 and OFF23 titles).
    # CIELAB L*: 45.4 (#6b6b6b), 65.0 (#9e9e9e), 86.6 (#d8d8d8)
    parts = [("energy_s", "energy model", None), ("dipole_s", "dipole model", DIPOLE_GREY),
             ("gaussian_own_s", "Gaussian", style.INK_MUTED), ("python_s", "Python interface", style.GRID)]
    t = df.groupby(["energy_model", "size"], observed=True)[[p for p, *_ in parts]].sum()
    t = t.div(t.sum(axis=1), axis=0) * 100
    print(t.loc["mace_anicc"].round(1).to_string())
    fig, axes = plt.subplots(1, 4, figsize=(style.figwidth(), 2.5), sharey=True)
    for ax, model in zip(axes, ELAB):
        bottom = pd.Series(0.0, index=labels)
        for col, _, color in parts:
            v = t.loc[model, col].reindex(labels).fillna(0)
            ax.bar(labels, v, bottom=bottom, color=color or style.model_color(model), width=0.75)
            bottom += v
        style.clean_axes(ax, grid=False)
        ax.set_title(style.tex(ELAB[model]), color=style.model_color(model), fontsize=style.BASE_FS)
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", labelrotation=60)
    axes[0].set_ylabel(r"share of wall time / \%")
    fig.supxlabel("number of atoms", fontsize=style.BASE_FS)
    # legend key of the energy model: the four panel colours side by side in one patch
    handles = [Patch(color=c, label=lab) if c else _MultiPatch(lab) for _, lab, c in parts]
    fig.legend(handles=handles, handler_map={_MultiPatch: _MultiHandler()}, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 1.0), handlelength=1.6, columnspacing=1.2, borderpad=0)  # fmt: skip
    style.tight(fig, pad=0.4, w_pad=0.6)
    style.save(fig, out, "time_split")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "analysis" / "results_figures" / "out"
    style.apply_style()
    molecules(out, "vpt2")
    molecules(out, "harmonic")
    time_split(out)
