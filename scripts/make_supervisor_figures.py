# ruff: noqa: E501  (embedded HTML/CSS template of the gallery)
#!/usr/bin/env python3
"""Supervisor figures: ML vs B3LYP per band type, plus chi and mode-alignment heatmaps.

Reads ``thesis/figures/supervisor/data/{bands,chi}.csv`` (written by
``mace_gaussian.analysis.band_table``; ``--rebuild`` regenerates them first) and writes
into ``thesis/figures/supervisor/``:

- ``grid_<mol>_<model>``: 3x3 per molecule and energy model. Columns fundamentals /
  overtones / combinations; rows frequency, intensity (one marker per dipole model) and
  anharmonic constants (x_ii under overtones, x_ij under combinations).
- ``grid_pooled_<model>``: the same, every molecule pooled (ethane left out: its B3LYP
  VPT2 baseline is broken, see thesis/TODO.md).
- ``*_xmatrix`` variants: chi row from Gaussian's deperturbed X matrix instead of band
  positions.
- ``overlap_<mol>``: eigenvector overlap ML vs B3LYP harmonic modes, every energy model.
- ``stats.csv`` and ``index.html`` (gallery with molecule / model selectors).

Usage::

    python scripts/make_supervisor_figures.py                     # everything
    python scripts/make_supervisor_figures.py --molecules water methanol --no-heatmaps
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, MaxNLocator, NullLocator  # noqa: E402

from mace_gaussian.analysis.thesis_style import (  # noqa: E402
    ANNOT_FS,
    DIPOLE_LABEL,
    INK,
    INK_MUTED,
    MODEL_LABEL,
    apply_style,
    clean_axes,
    dipole_color,
    figwidth,
    model_color,
    save,
    tex,
)

OUT = REPO / "thesis" / "figures" / "supervisor"
DATA = OUT / "data"
MODEL_ORDER = ["mace_omol", "mace_off", "mace_anicc", "mace_polar", "mace_mp"]
DIPOLE_ORDER = ["mace_ml", "mace_polar1", "espaloma"]
CHI_RUN_DIPOLE = "mace_ml"  # frequencies are dipole-independent; one run per energy model
POOL_EXCLUDE = {"ethane"}  # broken B3LYP VPT2 baseline (D3d symmetric-top treatment)
INT_FLOOR = 0.01  # km/mol; log axes, DFT intensities below this are left out
COLUMNS = [
    ("fundamental", "fundamentals"),
    ("overtone", "overtones"),
    ("combination", "combination bands"),
]
CM = r"cm$^{-1}$"


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def stats(dft: np.ndarray, ml: np.ndarray) -> dict:
    """R², slope, intercept (fit ML = a·DFT + b), MAE, mean signed error, max |error|."""
    dft, ml = np.asarray(dft, float), np.asarray(ml, float)
    n = len(dft)
    out = {"n": n}
    if n == 0:
        return out
    err = ml - dft
    out.update(
        mae=float(np.mean(np.abs(err))), mse=float(np.mean(err)), max=float(np.max(np.abs(err)))
    )
    if n >= 3 and np.ptp(dft) > 0:
        slope, intercept = np.polyfit(dft, ml, 1)
        r = np.corrcoef(dft, ml)[0, 1]
        out.update(slope=float(slope), intercept=float(intercept), r2=float(r * r))
    return out


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------


def _parity(ax, lo, hi, log=False):
    """Square parity panel in the thesis-figure style (dashed y = x, few ticks)."""
    clean_axes(ax)
    ax.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=0.9, ls=(0, (5, 3)), zorder=1)
    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")
        # at most ~3 labelled decades, so the exponents never collide
        decades = round(np.log10(hi / lo))
        step = max(1, int(np.ceil(decades / 3)))
        ticks = 10.0 ** np.arange(np.log10(lo), np.log10(hi) + 0.5, step)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator(ticks))
            axis.set_minor_locator(NullLocator())
    else:
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(MaxNLocator(nbins=3))
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")


def _ref_limits(ref, log=False, pad=0.08):
    """Axis range from the B3LYP (reference) values, so a few exploding ML values cannot
    squash the panel; ML points outside it are counted, not drawn."""
    ref = np.asarray(ref, float)
    if len(ref) == 0:
        return 0.0, 1.0
    lo, hi = ref.min(), ref.max()
    if log:
        return 10 ** np.floor(np.log10(lo)), 10 ** np.ceil(np.log10(hi))
    span = hi - lo or abs(hi) or 1.0
    return lo - pad * span, hi + pad * span


def _inside(x, y, lo, hi):
    return (x >= lo) & (x <= hi) & (y >= lo) & (y <= hi)


def _note(ax, mae=None, off=0, extra=""):
    lines = [f"MAE {mae:.0f}" if mae is not None and mae >= 10 else
             f"MAE {mae:.1f}" if mae is not None else ""]  # fmt: skip
    if off:
        lines.append(f"{off} off-scale")
    if extra:
        lines.append(extra)
    text = "\n".join(t for t in lines if t)
    if text:
        ax.text(0.96, 0.04, text, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=ANNOT_FS - 1, color=INK)  # fmt: skip


def _dots(ax, x, y, color, hollow=False):
    many = len(x) > 400  # pooled panels: smaller, lighter, so the density shows
    ms, alpha = (2.2, 0.5) if many else (3.2, 0.75)
    if hollow:
        ax.plot(x, y, "o", ms=ms, mfc="white", mec=color, mew=0.7, zorder=2)
    else:
        ax.plot(x, y, "o", ms=ms, color=color, mec="none", alpha=alpha, zorder=3)


def panel_frequency(ax, df, color):
    s = stats(df.dft_anharmonic, df.ml_anharmonic)
    if s["n"]:
        lo, hi = _ref_limits(df.dft_anharmonic)
        _parity(ax, lo, hi)
        ok = _inside(df.dft_anharmonic, df.ml_anharmonic, lo, hi)
        _dots(ax, df.dft_anharmonic[ok], df.ml_anharmonic[ok], color)
        _note(ax, s["mae"], int((~ok).sum()))
    else:
        _empty(ax)
    return s


def panel_intensity(ax, df):
    keep = df[(df.dft_intensity >= INT_FLOOR) & (df.ml_intensity > 0)]
    out, off = {}, 0
    if not len(keep):
        # e.g. overtones of centrosymmetric molecules: IR-forbidden, B3LYP gives exactly 0
        _empty(ax, f"all B3LYP intensities\nbelow {INT_FLOOR} km/mol" if len(df) else "no bands")
        return out
    lo, hi = _ref_limits(keep.dft_intensity, log=True)
    _parity(ax, lo, hi, log=True)
    for dip in DIPOLE_ORDER:
        d = keep[keep.dipole_model == dip]
        if not len(d):
            continue
        ok = _inside(d.dft_intensity, d.ml_intensity, lo, hi)
        off += int((~ok).sum())
        _dots(ax, d.dft_intensity[ok], d.ml_intensity[ok], dipole_color(dip))
        out[dip] = stats(np.log10(d.dft_intensity), np.log10(d.ml_intensity))
    _note(ax, None, off)
    return out


def panel_chi(ax, df, color, source):
    """source: 'bands' (resonance-shifted points hollow, stats on the rest) or 'xmatrix'."""
    if source == "bands":
        d_col, m_col = "dft_bands", "ml_bands"
        res = df.dft_resonant | df.ml_resonant
    else:
        df = df.dropna(subset=["dft_xmatrix", "ml_xmatrix"])
        d_col, m_col = "dft_xmatrix", "ml_xmatrix"
        res = pd.Series(False, index=df.index)
    clean = df[~res]
    s = stats(clean[d_col], clean[m_col])
    if not len(df):
        _empty(ax)
        return s
    lo, hi = _ref_limits(clean[d_col] if len(clean) else df[d_col])
    _parity(ax, lo, hi)
    ok = _inside(df[d_col], df[m_col], lo, hi)
    _dots(ax, df[d_col][res & ok], df[m_col][res & ok], color, hollow=True)
    _dots(ax, df[d_col][~res & ok], df[m_col][~res & ok], color)
    _note(ax, s.get("mae"), int((~ok).sum()))
    return s


def _empty(ax, text="no bands"):
    clean_axes(ax, grid=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center",
            color=INK_MUTED, fontsize=ANNOT_FS)  # fmt: skip


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def grid_figure(bands, chi, energy, heading, name, source="bands"):
    """3x3 in the style of the thesis figures: method colour for the single-series rows,
    dipole colours in the intensity row (as in the approved intensity figure)."""
    color = model_color(energy)
    w = figwidth(1.0)
    fig, axes = plt.subplots(3, 3, figsize=(w, w * 1.03))
    rows = []
    freq_df = bands[bands.dipole_model == CHI_RUN_DIPOLE]
    for c, (bt, col_label) in enumerate(COLUMNS):
        s = panel_frequency(axes[0, c], freq_df[freq_df.band_type == bt], color)
        rows.append(("frequency", bt, "", s))
        for dip, s in panel_intensity(axes[1, c], bands[bands.band_type == bt]).items():
            rows.append(("intensity_log10", bt, dip, s))
        axes[0, c].set_title(col_label, fontsize=ANNOT_FS, color=INK)
    for c, kind in ((1, "diagonal"), (2, "offdiagonal")):
        s = panel_chi(axes[2, c], chi[chi.kind == kind], color, source)
        rows.append((f"chi_{source}", kind, "", s))
    axes[2, 0].axis("off")
    how = "band positions" if source == "bands" else "X matrix"
    axes[2, 0].text(
        0.5, 0.5, f"$x_{{ii}}$ from overtones\n$x_{{ij}}$ from combinations\n({how})",
        transform=axes[2, 0].transAxes, ha="center", va="center", color=INK_MUTED,
        fontsize=ANNOT_FS,
    )  # fmt: skip
    axes[0, 0].set_ylabel(f"ML / {CM}")
    axes[0, 1].set_xlabel(f"B3LYP frequency / {CM}")
    axes[1, 0].set_ylabel(r"ML / km mol$^{-1}$")
    axes[1, 1].set_xlabel(r"B3LYP intensity / km mol$^{-1}$")
    axes[2, 1].set_ylabel(f"ML / {CM}")
    axes[2, 1].set_xlabel(f"B3LYP $x_{{ii}}$ / {CM}")
    axes[2, 2].set_xlabel(f"B3LYP $x_{{ij}}$ / {CM}")

    handles = [Line2D([], [], ls="none", label="dipole model:")] + [
        Line2D([], [], marker="o", ls="none", color=dipole_color(d), label=tex(DIPOLE_LABEL[d]))
        for d in DIPOLE_ORDER
    ]
    if source == "bands":
        handles.append(
            Line2D([], [], marker="o", ls="none", mfc="white", mec=INK_MUTED, mew=0.8,
                   label="resonance-shifted")
        )  # fmt: skip
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.965),
               handlelength=1.0, columnspacing=1.2, fontsize=ANNOT_FS)  # fmt: skip
    fig.suptitle(
        f"{heading} \\textperiodcentered{{}} {tex(MODEL_LABEL[energy])}",
        color=model_color(energy),
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.945), h_pad=0.6)
    save(fig, OUT, name, png_dpi=150)
    return rows


# ---------------------------------------------------------------------------
# Heatmaps (reuse the thesis-figure implementation, force_harmonic=True there)
# ---------------------------------------------------------------------------


def heatmap(mol: str):
    import make_thesis_figures as mtf

    mtf.OUT = OUT
    mtf.FIGURES.clear()
    mtf.fig_overlap_matrices(mol)
    src = OUT / f"overlap_matrices_{mol}"
    for ext in (".pdf", ".png"):
        p = src.with_suffix(ext)
        if p.exists():
            p.rename(OUT / f"overlap_{mol}{ext}")
    return (OUT / f"overlap_{mol}.png").exists()


# ---------------------------------------------------------------------------
# Gallery
# ---------------------------------------------------------------------------


def write_gallery(index: dict, stats_df: pd.DataFrame, pooled: list[str]):
    mols = sorted(index)
    models = [m for m in MODEL_ORDER if any(m in index[mol]["grids"] for mol in mols)]
    tbl = {}
    for (mol, energy), g in stats_df.groupby(["molecule", "energy_model"]):
        tbl[f"{mol}|{energy}"] = (
            g.drop(columns=["molecule", "energy_model"])
            .round(3)
            .to_html(index=False, na_rep="", classes="st", border=0)
        )
    payload = json.dumps({"index": index, "tables": tbl})
    opts = lambda xs, labels=None: "".join(  # noqa: E731
        f"<option value='{x}'>{html.escape((labels or {}).get(x, x))}</option>" for x in xs
    )
    page = f"""<!DOCTYPE html><html><head><meta charset='utf-8'><title>Supervisor figures</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#1a1a1a;background:#fff}}
h2{{border-bottom:1px solid #ddd;padding-bottom:4px;margin-top:2.5rem}}
select,label{{font-size:1rem;margin-right:1rem}} img{{max-width:100%;box-shadow:0 0 0 1px #eee;margin:1rem 0}}
table.st{{border-collapse:collapse;font-size:.8rem}} .st td,.st th{{padding:2px 8px;border-bottom:1px solid #eee;text-align:right}}
.muted{{color:#6b6b6b}}</style></head><body>
<h1>Supervisor figures</h1>
<p class='muted'>ML vs B3LYP-VPT2 by band type. Generated by <code>scripts/make_supervisor_figures.py</code>
from <code>data/bands.csv</code> and <code>data/chi.csv</code>. Bands are paired through the
eigenvector mode assignment. Intensity panels use log axes and leave out bands with DFT intensity
below {INT_FLOOR} km/mol; their statistics are on log10 intensities. Frequencies and anharmonic
constants come from the MACE4IR-dipole run of each energy model (they do not depend on the
dipole model).</p>
<h2>Per molecule</h2>
<label>molecule <select id='mol'>{opts(mols)}</select></label>
<label>energy model <select id='model'>{opts(models, MODEL_LABEL)}</select></label>
<label>χ from <select id='src'><option value='bands'>band positions</option><option value='xmatrix'>X matrix</option></select></label>
<div id='view'></div>
<h2>Pooled over all molecules</h2>
<p class='muted'>Ethane is left out (broken B3LYP VPT2 baseline).</p>
<label>energy model <select id='pmodel'>{opts(pooled, MODEL_LABEL)}</select></label>
<label>χ from <select id='psrc'><option value='bands'>band positions</option><option value='xmatrix'>X matrix</option></select></label>
<div id='pview'></div>
<h2>Mode alignment (eigenvector overlap, harmonic modes)</h2>
<label>molecule <select id='hmol'>{opts([m for m in mols if index[m]["heatmap"]])}</select></label>
<div id='hview'></div>
<script>
const D = {payload};
const $ = id => document.getElementById(id);
function show() {{
  const mol = $('mol').value, model = $('model').value, src = $('src').value;
  const ok = (D.index[mol].grids || []).includes(model);
  const name = 'grid_' + mol + '_' + model + (src === 'xmatrix' ? '_xmatrix' : '');
  $('view').innerHTML = ok ? `<img src="${{name}}.png"><p class='muted'>${{name}}.pdf</p>` +
    (D.tables[mol + '|' + model] || '') : '<p class=muted>no run for this model</p>';
}}
function pshow() {{
  const name = 'grid_pooled_' + $('pmodel').value + ($('psrc').value === 'xmatrix' ? '_xmatrix' : '');
  $('pview').innerHTML = `<img src="${{name}}.png">` + (D.tables['pooled|' + $('pmodel').value] || '');
}}
function hshow() {{ $('hview').innerHTML = `<img src="overlap_${{$('hmol').value}}.png">`; }}
['mol','model','src'].forEach(i => $(i).onchange = show);
['pmodel','psrc'].forEach(i => $(i).onchange = pshow);
$('hmol').onchange = hshow; show(); pshow(); hshow();
</script></body></html>"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--molecules", nargs="*", default=None)
    ap.add_argument("--rebuild", action="store_true", help="rebuild data/*.csv first")
    ap.add_argument("--no-heatmaps", action="store_true")
    args = ap.parse_args(argv)

    if args.rebuild or not (DATA / "bands.csv").exists():
        from mace_gaussian.analysis.band_table import write_tables

        avail = sorted(
            p.parent.name for p in (REPO / "analysis_results").glob("*/report_data.json")
        )
        write_tables(avail, DATA, results_base=REPO / "comparison_results",
                     analysis_base=REPO / "analysis_results")  # fmt: skip
    bands = pd.read_csv(DATA / "bands.csv")
    chi = pd.read_csv(DATA / "chi.csv")
    mols = args.molecules or sorted(bands.molecule.unique())
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)

    stat_rows, index = [], {}
    for mol in mols:
        index[mol] = {"grids": [], "heatmap": False}
        for energy in MODEL_ORDER:
            b = bands[(bands.molecule == mol) & (bands.energy_model == energy)]
            c = chi[(chi.molecule == mol) & (chi.energy_model == energy)]
            if b.empty:
                continue
            title = tex(mol)
            for src in ("bands", "xmatrix"):
                suffix = "" if src == "bands" else "_xmatrix"
                rows = grid_figure(b, c, energy, title, f"grid_{mol}_{energy}{suffix}", src)
                for qty, bt, dip, s in rows:
                    if src == "xmatrix" and not qty.startswith("chi"):
                        continue
                    stat_rows.append({"molecule": mol, "energy_model": energy, "quantity": qty,
                                      "band_type": bt, "dipole_model": dip, **s})  # fmt: skip
            index[mol]["grids"].append(energy)
            print(f"  {mol} {energy}", flush=True)
        if not args.no_heatmaps:
            index[mol]["heatmap"] = heatmap(mol)

    pooled = []

    # pooled: leave out ethane (broken baseline) and runs in a different conformer than B3LYP
    def _keep(df):
        wrong = df.get("same_conformer", pd.Series(dtype=object)).astype(str) == "False"
        return df[~df.molecule.isin(POOL_EXCLUDE) & ~wrong]

    pb, pc = _keep(bands), _keep(chi)
    for energy in MODEL_ORDER:
        b, c = pb[pb.energy_model == energy], pc[pc.energy_model == energy]
        if b.empty:
            continue
        title = "all molecules"
        for src in ("bands", "xmatrix"):
            suffix = "" if src == "bands" else "_xmatrix"
            rows = grid_figure(b, c, energy, title, f"grid_pooled_{energy}{suffix}", src)
            for qty, bt, dip, s in rows:
                if src == "xmatrix" and not qty.startswith("chi"):
                    continue
                stat_rows.append({"molecule": "pooled", "energy_model": energy, "quantity": qty,
                                  "band_type": bt, "dipole_model": dip, **s})  # fmt: skip
        pooled.append(energy)
        print(f"  pooled {energy}", flush=True)

    stats_df = pd.DataFrame(stat_rows)
    stats_df.to_csv(OUT / "stats.csv", index=False)
    write_gallery(index, stats_df, pooled)
    print(f"wrote {OUT / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
