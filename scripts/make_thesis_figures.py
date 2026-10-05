# ruff: noqa: RUF001, RUF002  (Greek nu in labels is intentional; long f-strings)
#!/usr/bin/env python3
"""Generate every thesis figure into thesis/figures/ plus a review gallery.

Reads only the analysis exports (``analysis_results/<mol>/report_data.json``,
``master_table.csv``) and ``comparison_results/*/results.json`` for timings;
never the HTML. Writes ``<name>.pdf`` (for ``\\includegraphics``), ``<name>.png``
(preview) and ``index.html`` (gallery, served by the local web server).

Usage::

    python scripts/make_thesis_figures.py                 # all figures
    python scripts/make_thesis_figures.py --only spectra   # one family
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import ScalarFormatter  # noqa: E402

from mace_gaussian.analysis.master_table import split_method  # noqa: E402
from mace_gaussian.analysis.thesis_style import (  # noqa: E402
    ANNOT_FS,
    DIPOLE_LABEL,
    INK,
    INK_MUTED,
    MODEL_LABEL,
    SERIES,
    apply_style,
    clean_axes,
    dipole_color,
    figwidth,
    model_color,
    save,
    tex,
)

ANALYSIS = REPO / "analysis_results"
COMPARISON = REPO / "comparison_results"
OUT = REPO / "thesis" / "figures"

MODEL_ORDER = ["mace_omol", "mace_off", "mace_anicc", "mace_polar", "mace_mp"]
DIPOLE_PREF = ["mace_ml", "mace_polar1", "mace_mdp", "espaloma"]
FWHM = 10.0  # cm^-1, Lorentzian, as in the report

FIGURES: list[dict] = []  # gallery entries


def register(name: str, width_frac: float, height_in: float, caption: str, section: str):
    FIGURES.append(
        {
            "name": name,
            "width": width_frac,
            "height": height_in,
            "caption": caption,
            "section": section,
        }
    )


# ---------------------------------------------------------------------------
# Data access
# ---------------------------------------------------------------------------


def load_report(mol: str) -> dict | None:
    p = ANALYSIS / mol / "report_data.json"
    return json.loads(p.read_text()) if p.exists() else None


def load_master_csv(mol: str) -> list[dict]:
    p = ANALYSIS / mol / "master_table.csv"
    if not p.exists():
        return []
    with p.open() as f:
        return list(csv.DictReader(f))


def representative_runs(report: dict) -> dict[str, dict]:
    """energy model -> the comparison with the preferred dipole model."""
    best: dict[str, tuple[int, dict]] = {}
    for c in report["comparisons"]:
        energy, dip = split_method(c["name"])
        rank = DIPOLE_PREF.index(dip) if dip in DIPOLE_PREF else 99
        if energy not in best or rank < best[energy][0]:
            best[energy] = (rank, c)
    return {k: v[1] for k, v in best.items()}


def broaden(freqs, ints, grid, fwhm=FWHM, scale=None):
    """Lorentzian-broadened spectrum. ``scale`` (the DFT peak) puts every trace on one
    shared scale so relative band strengths stay comparable; without it the trace is
    normalized to its own maximum. Returns (spectrum, peak)."""
    freqs = np.asarray(freqs, float)
    ints = np.asarray(ints, float)
    g = fwhm / 2.0
    out = np.zeros_like(grid)
    for f, i in zip(freqs, ints):
        out += i * g * g / ((grid - f) ** 2 + g * g)
    peak = float(out.max())
    denom = scale if scale else peak
    return (out / denom if denom > 0 else out), peak


def exp_on_grid(report: dict, grid, smooth_fwhm: float | None):
    e = report.get("experimental")
    if not e or not e.get("wavenumbers_cm"):
        return None
    order = np.argsort(e["wavenumbers_cm"])
    xp = np.asarray(e["wavenumbers_cm"], float)[order]
    fp = np.asarray(e["absorbance_normalized"], float)[order]
    y = np.interp(grid, xp, fp, left=0, right=0)
    if smooth_fwhm:
        step = float(grid[1] - grid[0])
        sigma = smooth_fwhm / 2.3548 / step
        half = int(4 * sigma) + 1
        x = np.arange(-half, half + 1)
        k = np.exp(-0.5 * (x / sigma) ** 2)
        y = np.convolve(y, k / k.sum(), mode="same")
    peak = y.max()
    return y / peak if peak > 0 else y


def band_labels(report: dict) -> list[tuple[str, float]]:
    seen: dict[str, float] = {}
    for r in report.get("master_table", {}).get("rows", []):
        e = r.get("experimental")
        if e and e.get("freq_cm") is not None:
            seen.setdefault(e["label"], float(e["freq_cm"]))
    out: list[tuple[str, float]] = []
    for label, f in sorted(seen.items(), key=lambda kv: kv[1]):
        if out and abs(f - out[-1][1]) <= 8:
            out[-1] = (out[-1][0] + "/" + label, (out[-1][1] + f) / 2)
        else:
            out.append((label, f))
    return out


def nu_tex(label: str) -> str:
    """'ν4/ν10' -> '$\\nu_{4}/\\nu_{10}$'."""
    parts = []
    for p in label.split("/"):
        p = p.replace("ν", "")
        parts.append(rf"\nu_{{{p}}}" if p.isdigit() else tex(p))
    return "$" + "/".join(parts) + "$"


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig_morse():
    """Morse vs harmonic well with the first vibrational levels (pedagogy, ch. 2)."""
    # HF-like parameters: D_e in cm^-1, a in 1/Angstrom, r_e in Angstrom
    De, a, re = 49000.0, 2.22, 0.917
    we = 4138.0
    wexe = we**2 / (4 * De)
    r = np.linspace(0.45, 2.2, 600)
    morse = De * (1 - np.exp(-a * (r - re))) ** 2
    harm = De * (a * (r - re)) ** 2

    fig, ax = plt.subplots(figsize=(figwidth(0.75), 3.4))
    clean_axes(ax)
    ax.plot(r, harm / 1000, color=INK_MUTED, lw=1.4, ls=(0, (5, 3)), label="harmonic")
    ax.plot(r, morse / 1000, color=INK, lw=2.0, label="Morse")
    for n in range(0, 6):
        e_h = we * (n + 0.5)
        e_m = we * (n + 0.5) - wexe * (n + 0.5) ** 2
        # turning points for the level lines
        rm = r[morse <= e_m]
        rh = r[harm <= e_h]
        if len(rm):
            ax.hlines(e_m / 1000, rm.min(), rm.max(), color=SERIES[0], lw=1.2)
            ax.annotate(
                rf"$v = {n}$",
                xy=(rm.max(), e_m / 1000),
                xytext=(4, 0),
                textcoords="offset points",
                fontsize=ANNOT_FS,
                color=SERIES[0],
                ha="left",
                va="center",
            )
        if len(rh):
            ax.hlines(e_h / 1000, rh.min(), rh.max(), color=INK_MUTED, lw=0.9, ls=(0, (2, 2)))
    ax.text(
        0.97,
        0.06,
        r"$E_v = \omega_e (v + \tfrac12) - \omega_e x_e (v + \tfrac12)^2$",
        transform=ax.transAxes,
        fontsize=ANNOT_FS,
        color=INK,
        ha="right",
        va="bottom",
    )
    ax.set_xlim(0.45, 2.2)
    ax.set_ylim(-1, 30)
    ax.set_xlabel(r"internuclear distance $r$ / \AA")
    ax.set_ylabel(r"energy / $10^{3}$ cm$^{-1}$")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, OUT, "morse_vs_harmonic")
    register(
        "morse_vs_harmonic",
        0.75,
        3.4,
        "Morse potential (ink) against its harmonic approximation (dashed) with the lowest "
        "vibrational levels; the anharmonic levels crowd together and overtones fall below "
        "integer multiples of the fundamental.",
        "2.1 Molecular vibrations",
    )


def fig_stacked_spectra(mol: str, smooth_small: bool = True):
    report = load_report(mol)
    if not report:
        return
    grid = np.arange(400.0, 4000.0, 1.0)
    reps = representative_runs(report)
    models = [m for m in MODEL_ORDER if m in reps]
    n_modes = len(report.get("master_table", {}).get("rows", []))
    exp = exp_on_grid(report, grid, 30.0 if (smooth_small and n_modes <= 9) else None)
    dft_spec = reps[models[0]]["spectrum_dft"]
    dft, dft_peak = broaden(dft_spec["frequencies_cm"], dft_spec["intensities"], grid)
    # Every calculated trace in units of the DFT peak; the lane height follows the
    # tallest one so a stronger spectrum gets more room instead of clipping.
    ml_traces = []
    for m in models:
        sp = reps[m]["spectrum_ml"]
        y, peak = broaden(sp["frequencies_cm"], sp["intensities"], grid, scale=dft_peak)
        ml_traces.append((m, y, peak / dft_peak if dft_peak else 1.0))
    offset = 1.1 * max([1.0, *(r for _m, _y, r in ml_traces)])
    n_traces = len(models) + 1 + (1 if exp is not None else 0)
    fig, ax = plt.subplots(figsize=(figwidth(1.0), 0.42 * n_traces + 1.1))
    clean_axes(ax, grid=False)
    level = 0.0
    baselines: list[float] = []

    def trace_label(text: str, color: str) -> None:
        # Outside the axes on the right (x = 400 is the right edge of the reversed axis).
        ax.annotate(
            text,
            xy=(400, level + 0.45),
            xytext=(4, 0),
            textcoords="offset points",
            fontsize=ANNOT_FS,
            color=color,
            va="center",
            ha="left",
            annotation_clip=False,
        )

    if exp is not None:
        ax.fill_between(grid, level, exp + level, color=INK_MUTED, alpha=0.25, lw=0)
        ax.plot(grid, exp + level, color=INK_MUTED, lw=0.9)
        trace_label("experiment (NIST)", INK_MUTED)
        baselines.append(level)
        level += offset
    ax.plot(grid, dft + level, color=INK, lw=1.4)
    trace_label("B3LYP (VPT2)", INK)
    baselines.append(level)
    for m, y, _ratio in ml_traces:
        level += offset
        ax.plot(grid, y + level, color=model_color(m), lw=1.4)
        trace_label(tex(MODEL_LABEL[m]), model_color(m))
        baselines.append(level)
    top = level + 1.05
    prev_f, tier = None, 0
    for label, f in band_labels(report):
        if not (grid[0] <= f <= grid[-1]):
            continue
        # neighbours closer than ~90 cm^-1 alternate between two label heights
        tier = (tier + 1) % 2 if prev_f is not None and abs(f - prev_f) < 90 else 0
        prev_f = f
        ax.vlines(f, -0.05, top + 0.55 * tier, color=INK_MUTED, lw=0.6, ls=(0, (1, 2)), zorder=0)
        ax.text(
            f,
            top + 0.04 + 0.55 * tier,
            nu_tex(label),
            fontsize=ANNOT_FS - 1,
            color=INK_MUTED,
            ha="center",
            va="bottom",
            rotation=90,
        )
    ax.set_xlim(4000, 400)
    ax.set_ylim(-0.05, top + 1.3)
    # y axis stays: one tick per trace baseline, the first one scaled 0 to 1
    ax.set_yticks(baselines)
    ax.set_yticklabels(["0"] + [""] * (len(baselines) - 1))
    ax.set_ylabel(r"absorbance (DFT peak $=1$)")
    ax.set_xlabel(r"wavenumber $\tilde{\nu}$ / cm$^{-1}$")
    fig.tight_layout()
    fig.subplots_adjust(right=0.80)
    save(fig, OUT, f"spectra_stacked_{mol}")
    register(
        f"spectra_stacked_{mol}",
        1.0,
        0.42 * n_traces + 1.1,
        f"{mol}: VPT2 spectra of the five energy models (MACE4IR dipoles) stacked above the "
        "B3LYP reference and the NIST gas-phase spectrum. All calculated traces share one scale "
        "(B3LYP peak = 1), so relative band strengths are comparable. The experimental trace is "
        "normalized to its own maximum (arbitrary units). "
        "Dotted lines: experimental band origins (Shimanouchi).",
        "5.1 / 5.6 spectra",
    )


def fig_per_mode_errors(mol: str):
    report = load_report(mol)
    if not report:
        return
    rows = report.get("master_table", {}).get("rows", [])
    if not rows:
        return
    reps = representative_runs(report)
    models = [m for m in MODEL_ORDER if m in reps]
    has_exp = any(r.get("experimental") for r in rows)
    panels = [("VPT2", "vpt2", "dft_vpt2"), ("harmonic", "harmonic", "dft_harmonic")]
    if has_exp:
        panels.append((r"VPT2 vs.\ experiment", "vpt2", "exp"))
    x = np.arange(len(rows))
    compact = len(rows) > 24
    labels = []
    for r in rows:
        e = r.get("experimental") or {}
        ref = r["dft_vpt2"] if r["dft_vpt2"] is not None else r["dft_harmonic"]
        if compact or not e:
            labels.append(f"{r['dft_mode']}" if compact else f"{r['dft_mode']} ({ref:.0f})")
        else:
            labels.append(nu_tex(e["label"]) + rf" {tex(e.get('description', ''))} ({ref:.0f})")

    fig, axes = plt.subplots(
        len(panels), 1, sharex=True, figsize=(figwidth(1.0), 1.55 * len(panels) + 1.0)
    )
    for ax, (title, key, ref_key) in zip(axes, panels):
        clean_axes(ax)
        ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
        if ref_key == "exp":
            y = [
                r["dft_vpt2"] - r["experimental"]["freq_cm"]
                if r.get("experimental") and r["dft_vpt2"] is not None
                else np.nan
                for r in rows
            ]
            ax.plot(x, y, marker="*", ms=7, color=INK, lw=0.8, label="B3LYP", zorder=4)
        for m in models:
            name = reps[m]["name"]
            ys, hollow = [], []
            for r in rows:
                c = r["methods"].get(name) or {}
                v = c.get(key)
                ref = (
                    r["experimental"]["freq_cm"]
                    if ref_key == "exp" and r.get("experimental")
                    else r.get(ref_key)
                )
                ys.append(v - ref if v is not None and ref is not None else np.nan)
                hollow.append(c.get("overlap") is not None and c["overlap"] < 0.7)
            ys = np.array(ys)
            hollow = np.array(hollow)
            ax.plot(x, ys, color=model_color(m), lw=1.6, alpha=0.9, zorder=3)
            ax.plot(
                x[~hollow],
                ys[~hollow],
                "o",
                ms=5.0 if compact else 6.5,
                color=model_color(m),
                mec="white",
                mew=0.8,
                zorder=4,
                label=tex(MODEL_LABEL[m]),
            )
            if hollow.any():
                ax.plot(
                    x[hollow],
                    ys[hollow],
                    "o",
                    ms=5.5,
                    mfc="white",
                    mec=model_color(m),
                    mew=1.2,
                    zorder=4,
                )
        ax.set_ylabel(title + "\n" + r"$\Delta\tilde{\nu}$ / cm$^{-1}$")
    axes[0].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.28),
        ncol=len(models) + 1,
        handlelength=1.2,
        columnspacing=1.2,
    )
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(
        labels,
        rotation=90 if compact else 40,
        ha="center" if compact else "right",
        fontsize=ANNOT_FS - (1 if compact else 0),
    )
    axes[-1].set_xlabel("DFT normal mode" + (" (index)" if compact else ""))
    fig.tight_layout()
    save(fig, OUT, f"per_mode_errors_{mol}")
    register(
        f"per_mode_errors_{mol}",
        1.0,
        1.55 * len(panels) + 1.0,
        f"{mol}: signed frequency error of every energy model per DFT normal mode (paired by "
        "eigenvector overlap; hollow markers below 0.7), for the VPT2 and harmonic frequencies "
        + ("and against the experimental band origins (B3LYP as stars)." if has_exp else "."),
        "5.2 per-mode errors",
    )


def _mae_table(mols: list[str]) -> tuple[list[str], np.ndarray, np.ndarray]:
    """rows = molecules, cols = MODEL_ORDER; MAE and bias of VPT2 fundamentals vs DFT."""
    mae = np.full((len(mols), len(MODEL_ORDER)), np.nan)
    bias = np.full_like(mae, np.nan)
    for i, mol in enumerate(mols):
        report = load_report(mol)
        if not report:
            continue
        reps = representative_runs(report)
        for j, m in enumerate(MODEL_ORDER):
            if m not in reps:
                continue
            name = reps[m]["name"]
            errs = []
            for r in report.get("master_table", {}).get("rows", []):
                c = r["methods"].get(name) or {}
                if c.get("vpt2") is not None and r["dft_vpt2"] is not None:
                    errs.append(c["vpt2"] - r["dft_vpt2"])
            if errs:
                mae[i, j] = np.mean(np.abs(errs))
                bias[i, j] = np.mean(errs)
    return mols, mae, bias


def fig_mae_heatmap(mols: list[str]):
    """Molecules as rows, energy models as columns; the whole figure is roughly square."""
    mols = [m for m in mols if load_report(m)]
    if not mols:
        return
    mols, mae, _bias = _mae_table(mols)
    keep = ~np.isnan(mae).all(axis=1)
    mols = [m for m, k in zip(mols, keep) if k]
    mae = mae[keep]
    side = figwidth(0.72)
    fig, ax = plt.subplots(figsize=(side + 0.9, side))
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("#e6e6e6")
    vmax = float(np.nanpercentile(mae, 95))
    im = ax.imshow(np.ma.masked_invalid(mae), cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i in range(mae.shape[0]):
        for j in range(mae.shape[1]):
            v = mae[i, j]
            if np.isnan(v):
                ax.text(j, i, "--", ha="center", va="center", fontsize=ANNOT_FS, color=INK_MUTED)
                continue
            ax.text(
                j,
                i,
                f"{v:.0f}",
                ha="center",
                va="center",
                fontsize=ANNOT_FS,
                color="white" if v > 0.6 * vmax else INK,
            )
    ax.set_xticks(range(len(MODEL_ORDER)))
    ax.set_xticklabels([tex(MODEL_LABEL[m].replace("MACE-", "")) for m in MODEL_ORDER])
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(mols)))
    ax.set_yticklabels([tex(m.replace("_", " ")) for m in mols])
    ax.tick_params(length=0)
    for side_ in ax.spines.values():
        side_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.05)
    cb.set_label(r"MAE / cm$^{-1}$")
    cb.outline.set_visible(False)
    fig.tight_layout()
    save(fig, OUT, "mae_heatmap")
    register(
        "mae_heatmap",
        0.72,
        side,
        "Mean absolute error of the VPT2 fundamentals against B3LYP per molecule and energy "
        "model (MACE4IR dipoles), in cm$^{-1}$.",
        "5.2 panel overview",
    )


def fig_intensity_scatter(mols: list[str]):
    """DFT vs ML intensities pooled over molecules, one facet per energy model."""
    data: dict[str, dict[str, list]] = {m: {} for m in MODEL_ORDER}
    for mol in mols:
        rows = load_master_csv(mol)
        if not rows:
            continue
        runs = sorted(
            {
                k[: -len("_intensity")]
                for k in rows[0]
                if k.endswith("_intensity") and k != "dft_intensity"
            }
        )
        for run in runs:
            energy, dip = split_method(run)
            if energy not in data or dip not in DIPOLE_PREF:
                continue
            pts = data[energy].setdefault(dip, [])
            for r in rows:
                try:
                    d, m = float(r["dft_intensity"]), float(r[f"{run}_intensity"])
                except ValueError:
                    continue
                if d >= 0.1 or m >= 0.1:
                    pts.append((max(d, 1e-3), max(m, 1e-3)))
    models = [m for m in MODEL_ORDER if data[m]]
    if not models:
        return
    fig, axes = plt.subplots(
        1,
        len(models),
        figsize=(figwidth(1.0), figwidth(1.0) / len(models) + 0.95),
        sharex=True,
        sharey=True,
    )
    axes = np.atleast_1d(axes)
    lo, hi = 1e-2, 1e3
    for ax, m in zip(axes, models):
        clean_axes(ax)
        ax.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=0.9, ls=(0, (5, 3)), zorder=1)
        for dip in DIPOLE_PREF:
            pts = np.array(data[m].get(dip, []))
            if len(pts):
                ax.plot(
                    pts[:, 0],
                    pts[:, 1],
                    "o",
                    ms=3.2,
                    color=dipole_color(dip),
                    mec="none",
                    alpha=0.75,
                    zorder=3,
                    label=tex(DIPOLE_LABEL[dip]),
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(tex(MODEL_LABEL[m]), fontsize=ANNOT_FS, color=model_color(m))
    axes[0].set_ylabel(r"ML / km mol$^{-1}$")
    axes[len(axes) // 2].set_xlabel(r"B3LYP / km mol$^{-1}$")
    handles = [Line2D([], [], ls="none", label="dipole model:")] + [
        Line2D([], [], marker="o", ls="none", color=dipole_color(d), label=tex(DIPOLE_LABEL[d]))
        for d in DIPOLE_PREF
        if any(data[m].get(d) for m in models)
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 1.0),
        handlelength=1.0,
        columnspacing=1.6,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, OUT, "intensity_scatter")
    # Single-model versions, half text width, same axes and colours.
    for m in models:
        fig1, ax = plt.subplots(figsize=(figwidth(0.5), figwidth(0.5) + 0.3))
        clean_axes(ax)
        ax.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=0.9, ls=(0, (5, 3)), zorder=1)
        for dip in DIPOLE_PREF:
            pts = np.array(data[m].get(dip, []))
            if len(pts):
                ax.plot(
                    pts[:, 0],
                    pts[:, 1],
                    "o",
                    ms=4.0,
                    color=dipole_color(dip),
                    mec="none",
                    alpha=0.75,
                    zorder=3,
                    label=tex(DIPOLE_LABEL[dip]),
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(r"B3LYP / km mol$^{-1}$")
        ax.set_ylabel(r"ML / km mol$^{-1}$")
        ax.set_title(tex(MODEL_LABEL[m]), color=model_color(m), pad=24)
        ax.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=3,
            handlelength=0.8,
            handletextpad=0.4,
            columnspacing=0.9,
            borderaxespad=0.2,
            fontsize=ANNOT_FS,
        )
        fig1.tight_layout()
        save(fig1, OUT, f"intensity_scatter_{m}")
        register(
            f"intensity_scatter_{m}",
            0.5,
            figwidth(0.5),
            f"IR intensities of the VPT2 fundamentals for {MODEL_LABEL[m]}, ML against B3LYP, "
            "pooled over the molecule panel; colour by dipole model.",
            "5.4 intensities (single)",
        )
    register(
        "intensity_scatter",
        1.0,
        figwidth(1.0) / len(models) + 0.7,
        "IR intensities of the VPT2 fundamentals, ML against B3LYP, pooled over the molecule "
        "panel; one facet per energy model, colour by dipole model. Bands below 0.1 km/mol in "
        "both are omitted.",
        "5.4 intensities",
    )


def fig_cost_scaling():
    """Wall time of the full VPT2 run vs number of atoms for the alkane ladder."""
    ladder = [
        "ethane",
        "propane",
        "butane",
        "pentane",
        "hexane",
        "heptane",
        "octane",
        "nonane",
        "decane",
    ]
    pts: dict[str, list[tuple[int, float]]] = {}
    for mol in ladder:
        xyz = REPO / "molecules" / f"{mol}.xyz"
        if not xyz.exists():
            continue
        natoms = int(xyz.read_text().splitlines()[0].strip())
        for p in (COMPARISON / mol).glob("*/results.json"):
            r = json.loads(p.read_text())
            if r.get("calculator_type") == "dft":
                t = (r.get("timing") or r.get("gaussian_timing") or {}).get(
                    "total_elapsed_s"
                ) or r.get("runtime_s")
                key = "dft"
            else:
                if not p.parent.name.endswith("_mace_ml"):
                    continue
                t = r.get("runtime_s")
                key = split_method(p.parent.name)[0]
            if t:
                pts.setdefault(key, []).append((natoms, float(t)))
    if not pts:
        return
    fig, ax = plt.subplots(figsize=(figwidth(0.75), 3.3))
    clean_axes(ax)
    labels_right: list[tuple[float, float, str]] = []
    for key in ["dft"] + [m for m in MODEL_ORDER if m in pts]:
        if key not in pts:
            continue
        arr = np.array(sorted(pts[key]))
        color = INK if key == "dft" else model_color(key)
        label = "B3LYP" if key == "dft" else MODEL_LABEL[key]
        ax.plot(
            arr[:, 0],
            arr[:, 1],
            "o-",
            color=color,
            lw=1.4,
            ms=4.5,
            mec="white",
            mew=0.8,
            label=tex(label),
        )
        if len(arr) >= 3:
            slope = np.polyfit(np.log(arr[:, 0]), np.log(arr[:, 1]), 1)[0]
            labels_right.append((float(np.log10(arr[-1, 1])), slope, color))
    # exponent labels in one column right of the data, pushed apart where lines end close
    x_lab = max(n for arr in pts.values() for n, _ in arr) * 1.05
    labels_right.sort()
    min_gap = 0.09  # decades
    for k in range(1, len(labels_right)):
        y_prev = labels_right[k - 1][0]
        if labels_right[k][0] - y_prev < min_gap:
            labels_right[k] = (y_prev + min_gap, *labels_right[k][1:])
    for y_log, slope, color in labels_right:
        ax.text(
            x_lab,
            10**y_log,
            rf"$\propto N^{{{slope:.1f}}}$",
            fontsize=ANNOT_FS - 1,
            color=color,
            va="center",
            ha="left",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    all_n = sorted({int(n) for arr in pts.values() for n, _ in arr})
    ax.set_xticks(all_n)
    ax.set_xticks([], minor=True)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xlim(right=x_lab * 1.28)
    ax.set_xlabel(r"number of atoms $N$")
    ax.set_ylabel("wall time of the VPT2 run / s")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=3,
        handlelength=1.4,
        columnspacing=1.0,
        borderaxespad=0.2,
        fontsize=ANNOT_FS,
    )
    fig.tight_layout()
    save(fig, OUT, "cost_scaling_alkanes")
    register(
        "cost_scaling_alkanes",
        0.75,
        3.3,
        "Wall time of a complete anharmonic (VPT2) run against molecule size for the n-alkane "
        "ladder, B3LYP on the cluster node versus the ML models on one GPU workstation; "
        "power-law exponents from a log-log fit.",
        "6.1 cost scaling",
    )


# ---------------------------------------------------------------------------
# Second batch (2026-09-18): pooled per-model views and the chapter 5.5 plot
# ---------------------------------------------------------------------------

B3LYP_SCALE = 0.961  # harmonic scaling factor for B3LYP/6-31G(d,p), CCCBDB-style


def _pooled_fundamentals(mols: list[str]) -> dict[str, list[dict]]:
    """model -> rows {mol, dft_harm, dft_vpt2, ml_harm, ml_vpt2, overlap, exp}."""
    out: dict[str, list[dict]] = {m: [] for m in MODEL_ORDER}
    for mol in mols:
        report = load_report(mol)
        if not report:
            continue
        reps = representative_runs(report)
        for m, comp in reps.items():
            if m not in out:
                continue
            for r in report.get("master_table", {}).get("rows", []):
                c = r["methods"].get(comp["name"]) or {}
                if c.get("vpt2") is None or r["dft_vpt2"] is None:
                    continue
                out[m].append(
                    {
                        "mol": mol,
                        "dft_harm": r["dft_harmonic"],
                        "dft_vpt2": r["dft_vpt2"],
                        "ml_harm": c.get("harmonic"),
                        "ml_vpt2": c["vpt2"],
                        "overlap": c.get("overlap"),
                        "exp": (r.get("experimental") or {}).get("freq_cm"),
                    }
                )
    return out


def _facets(n: int, height: float = 2.2):
    fig, axes = plt.subplots(1, n, figsize=(figwidth(1.0), height), sharey=True)
    return fig, np.atleast_1d(axes)


def fig_residuals(mols: list[str]):
    """Signed VPT2 error against DFT frequency, pooled, one facet per model."""
    data = _pooled_fundamentals(mols)
    models = [m for m in MODEL_ORDER if data[m]]
    if not models:
        return
    fig, axes = _facets(len(models), 3.0)
    allerr = np.concatenate([[r["ml_vpt2"] - r["dft_vpt2"] for r in data[m]] for m in models])
    lim = float(np.percentile(np.abs(allerr), 98)) * 1.1
    for ax, m in zip(axes, models):
        clean_axes(ax)
        rows = data[m]
        x = np.array([r["dft_vpt2"] for r in rows])
        y = np.array([r["ml_vpt2"] - r["dft_vpt2"] for r in rows])
        low = np.array([r["overlap"] is not None and r["overlap"] < 0.7 for r in rows])
        ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
        ax.plot(x[~low], y[~low], "o", ms=3.2, color=model_color(m), mec="none", alpha=0.75)
        if low.any():
            ax.plot(x[low], y[low], "o", ms=3.2, mfc="white", mec=model_color(m), mew=0.8)
        n_clip = int((np.abs(y) > lim).sum())
        ax.set_ylim(-lim, lim)
        ax.set_xlim(4000, 300)
        ax.set_xticks([3000, 1500])
        ax.set_title(tex(MODEL_LABEL[m]), fontsize=ANNOT_FS, color=model_color(m))
        ax.text(
            0.04,
            0.03,
            rf"MAE {np.mean(np.abs(y)):.0f}" + (f"\n{n_clip} off-scale" if n_clip else ""),
            transform=ax.transAxes,
            fontsize=ANNOT_FS - 1,
            color=INK,
            va="bottom",
        )
    axes[0].set_ylabel(r"$\tilde{\nu}_\mathrm{ML} - \tilde{\nu}_\mathrm{DFT}$ / cm$^{-1}$")
    axes[len(axes) // 2].set_xlabel(r"B3LYP VPT2 fundamental / cm$^{-1}$")
    fig.tight_layout()
    save(fig, OUT, "residuals_pooled")
    register(
        "residuals_pooled",
        1.0,
        3.0,
        "Signed error of the VPT2 fundamentals against the B3LYP frequency, pooled over the "
        "molecule panel; one facet per energy model, hollow markers for eigenvector overlap "
        "below 0.7. The bending region (below 1700 cm$^{-1}$) and the stretching region separate "
        "the models.",
        "5.2 residuals",
    )


def fig_anharmonicity_ratio(mols: list[str]):
    """ML anharmonic correction against the DFT one, pooled, one facet per model."""
    data = _pooled_fundamentals(mols)
    models = [m for m in MODEL_ORDER if data[m]]
    if not models:
        return
    fig, axes = _facets(len(models), figwidth(1.0) / len(models) + 0.9)
    lo, hi = -2.0, 12.0
    for ax, m in zip(axes, models):
        clean_axes(ax)
        rows = [r for r in data[m] if r["ml_harm"] and r["dft_harm"]]
        xd = np.array([(r["dft_harm"] - r["dft_vpt2"]) / r["dft_harm"] * 100 for r in rows])
        ym = np.array([(r["ml_harm"] - r["ml_vpt2"]) / r["ml_harm"] * 100 for r in rows])
        ax.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=0.9, ls=(0, (5, 3)), zorder=1)
        ax.plot(xd, ym, "o", ms=3.2, color=model_color(m), mec="none", alpha=0.75, zorder=3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(tex(MODEL_LABEL[m]), fontsize=ANNOT_FS, color=model_color(m))
        if len(xd) > 2:
            # R^2 about the diagonal y = x (not a fitted line): negative means the ML
            # corrections deviate from the DFT ones by more than the DFT corrections vary.
            r2 = 1 - np.sum((ym - xd) ** 2) / np.sum((xd - xd.mean()) ** 2)
            txt = rf"$R^2 = {r2:.2f}$" if r2 >= 0 else r"$R^2 < 0$"
            ax.text(0.04, 0.9, txt, transform=ax.transAxes, fontsize=ANNOT_FS - 1)
    axes[0].set_ylabel(r"ML correction / \%")
    axes[len(axes) // 2].set_xlabel(r"B3LYP anharmonic correction $(\omega - \nu)/\omega$ / \%")
    fig.tight_layout()
    save(fig, OUT, "anharmonicity_ratio_pooled")
    register(
        "anharmonicity_ratio_pooled",
        1.0,
        figwidth(1.0) / len(models) + 0.9,
        "Relative anharmonic correction of each fundamental, ML against B3LYP, pooled over the "
        "panel; a model on the diagonal carries the same anharmonicity as the DFT surface "
        "regardless of its harmonic error.",
        "5.5 anharmonicity",
    )


def _box(ax, pos, values, color, width=0.6):
    bp = ax.boxplot(
        [values],
        positions=[pos],
        widths=width,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color=INK, lw=1.2),
        whiskerprops=dict(color=color, lw=0.9),
        capprops=dict(color=color, lw=0.9),
        boxprops=dict(facecolor=color, alpha=0.25, edgecolor=color),
    )
    rng = np.random.default_rng(0)
    ax.plot(
        pos + rng.uniform(-0.18, 0.18, len(values)),
        values,
        "o",
        ms=2.0,
        color=color,
        mec="none",
        alpha=0.6,
        zorder=3,
    )
    return bp


def fig_error_distribution(mols: list[str]):
    data = _pooled_fundamentals(mols)
    models = [m for m in MODEL_ORDER if data[m]]
    if not models:
        return
    fig, ax = plt.subplots(figsize=(figwidth(0.75), 3.0))
    clean_axes(ax)
    ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
    for i, m in enumerate(models):
        err = np.array([r["ml_vpt2"] - r["dft_vpt2"] for r in data[m]])
        _box(ax, i, err, model_color(m))
        ax.text(i, ax.get_ylim()[0], "", fontsize=ANNOT_FS)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels([tex(MODEL_LABEL[m].replace("MACE-", "")) for m in models])
    lim = (
        float(
            np.percentile(
                np.abs(
                    np.concatenate(
                        [[r["ml_vpt2"] - r["dft_vpt2"] for r in data[m]] for m in models]
                    )
                ),
                98,
            )
        )
        * 1.1
    )
    ax.set_ylim(-lim, lim)
    for i, m in enumerate(models):
        err = np.array([r["ml_vpt2"] - r["dft_vpt2"] for r in data[m]])
        ax.text(
            i,
            -lim * 0.97,
            rf"MAE {np.mean(np.abs(err)):.0f}",
            ha="center",
            va="bottom",
            fontsize=ANNOT_FS - 1,
            color=INK,
        )
    ax.set_ylabel(r"$\tilde{\nu}_\mathrm{ML} - \tilde{\nu}_\mathrm{DFT}$ / cm$^{-1}$")
    fig.tight_layout()
    save(fig, OUT, "error_distribution")
    register(
        "error_distribution",
        0.75,
        3.0,
        "Distribution of the signed VPT2 fundamental errors per energy model over the whole "
        "panel (boxes: quartiles and median; points: individual modes; MAE annotated).",
        "5.2 error distribution",
    )


def _paired_bands(report: dict, comp: dict) -> dict[str, list[tuple[float, float]]]:
    """(dft, ml) frequency pairs by band type, pairing through the mode mapping."""
    mapping = {int(k) + 1: int(v) + 1 for k, v in (comp.get("mode_mapping") or {}).items()}
    dft = {
        i: f
        for i, f in zip(comp["spectrum_dft"]["mode_ids"], comp["spectrum_dft"]["frequencies_cm"])
    }
    out: dict[str, list] = {"fundamental": [], "overtone": [], "combination": []}
    for mid, f in zip(comp["spectrum_ml"]["mode_ids"], comp["spectrum_ml"]["frequencies_cm"]):
        kind = mid[0]
        if kind == "F":
            j = mapping.get(int(mid[1:]))
            key, typ = (f"F{j}" if j else None), "fundamental"
        elif kind == "O":
            j, rest = mid[1:].split("_", 1)  # rest: level, plus "_l±2" for degenerate modes
            jj = mapping.get(int(j))
            key, typ = (f"O{jj}_{rest}" if jj else None), "overtone"
        else:
            a, b = (mapping.get(int(t)) for t in mid[1:].split("_"))
            key, typ = (f"C{min(a, b)}_{max(a, b)}" if a and b else None), "combination"
        if key in dft:
            out[typ].append((dft[key], f))
    return out


def fig_band_type_errors(mols: list[str]):
    errs: dict[str, dict[str, list[float]]] = {
        m: {"fundamental": [], "overtone": [], "combination": []} for m in MODEL_ORDER
    }
    for mol in mols:
        report = load_report(mol)
        if not report:
            continue
        for m, comp in representative_runs(report).items():
            if m not in errs or not comp.get("mode_mapping"):
                continue
            for typ, pairs in _paired_bands(report, comp).items():
                errs[m][typ] += [ml - d for d, ml in pairs]
    models = [m for m in MODEL_ORDER if errs[m]["fundamental"]]
    if not models:
        return
    types = [
        ("fundamental", "fundamentals"),
        ("overtone", "overtones"),
        ("combination", "combination bands"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(figwidth(1.0), 2.8), sharey=True)
    for ax, (typ, label) in zip(axes, types):
        clean_axes(ax)
        ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
        for i, m in enumerate(models):
            v = np.array(errs[m][typ])
            if len(v):
                _box(ax, i, v, model_color(m))
                ax.text(
                    i,
                    0.02,
                    rf"{np.mean(np.abs(v)):.0f}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=ANNOT_FS - 1,
                    color=INK,
                )
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels(
            [tex(MODEL_LABEL[m].replace("MACE-", "")) for m in models], rotation=40, ha="right"
        )
        ax.set_title(label, fontsize=ANNOT_FS)
    allv = np.concatenate([np.array(errs[m][t]) for m in models for t, _ in types if errs[m][t]])
    lim = float(np.percentile(np.abs(allv), 97)) * 1.1
    axes[0].set_ylim(-lim, lim)
    axes[0].set_ylabel(r"$\tilde{\nu}_\mathrm{ML} - \tilde{\nu}_\mathrm{DFT}$ / cm$^{-1}$")
    fig.tight_layout()
    save(fig, OUT, "band_type_errors")
    register(
        "band_type_errors",
        1.0,
        2.8,
        "Signed errors against B3LYP by band type, pooled over the panel: fundamentals, first "
        "overtones and binary combination bands (paired through the eigenvector mapping of the "
        "underlying fundamentals). Numbers at the bottom: MAE in cm$^{-1}$.",
        "5.5 anharmonic content",
    )


def fig_pareto(mols: list[str]):
    pts: dict[str, list[tuple[float, float]]] = {m: [] for m in MODEL_ORDER}
    for mol in mols:
        report = load_report(mol)
        if not report:
            continue
        for m, comp in representative_runs(report).items():
            sp = comp.get("runtime", {}).get("speedup", 0.0)
            if m in pts and sp > 0:
                pts[m].append((sp, comp["metrics"]["mae_freq"]))
    models = [m for m in MODEL_ORDER if pts[m]]
    if not models:
        return
    fig, ax = plt.subplots(figsize=(figwidth(0.6), 2.9))
    clean_axes(ax)
    for m in models:
        arr = np.array(pts[m])
        ax.plot(arr[:, 0], arr[:, 1], "o", ms=5.5, color=model_color(m), mec="none", alpha=0.5)
        ax.plot(
            np.exp(np.mean(np.log(arr[:, 0]))),
            arr[:, 1].mean(),
            "D",
            ms=11,
            color=model_color(m),
            mec="white",
            mew=1.0,
            label=tex(MODEL_LABEL[m]),
            zorder=4,
        )
    ax.set_xscale("log")
    ax.margins(x=0.08)
    all_mae = np.concatenate([np.array(pts[m])[:, 1] for m in models])
    ax.set_ylim(0, float(np.percentile(all_mae, 92)) * 1.25)
    ax.set_xlabel("speed-up over the B3LYP run (Gaussian wall time)")
    ax.set_ylabel(r"MAE of VPT2 fundamentals / cm$^{-1}$")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=3,
        handlelength=1.0,
        columnspacing=0.9,
        borderaxespad=0.2,
        fontsize=ANNOT_FS,
    )
    fig.tight_layout()
    save(fig, OUT, "pareto_cost_accuracy")
    register(
        "pareto_cost_accuracy",
        0.6,
        2.9,
        "Accuracy against cost: MAE of the VPT2 fundamentals versus the wall-time speed-up over "
        "the B3LYP run, one small point per molecule and a diamond at each model's panel mean.",
        "6 cost vs accuracy",
    )


def fig_central_vs_experiment(mols: list[str]):
    """Section 5.5: scaled-harmonic B3LYP vs B3LYP VPT2 vs ML VPT2, all against experiment."""
    data = _pooled_fundamentals(mols)
    models = [m for m in MODEL_ORDER if any(r["exp"] for r in data[m])]
    if not models:
        return
    approaches: list[tuple[str, str, np.ndarray]] = []
    ref = data[models[0]]
    exp_rows = [r for r in ref if r["exp"] is not None and r["dft_harm"] is not None]
    approaches.append(
        (
            rf"B3LYP harmonic $\times {B3LYP_SCALE}$",
            INK_MUTED,
            np.array([r["dft_harm"] * B3LYP_SCALE - r["exp"] for r in exp_rows]),
        )
    )
    approaches.append(("B3LYP VPT2", INK, np.array([r["dft_vpt2"] - r["exp"] for r in exp_rows])))
    for m in models:
        rows = [r for r in data[m] if r["exp"] is not None]
        approaches.append(
            (
                tex(MODEL_LABEL[m].replace("MACE-", "")) + " VPT2",
                model_color(m),
                np.array([r["ml_vpt2"] - r["exp"] for r in rows]),
            )
        )
    fig, ax = plt.subplots(figsize=(figwidth(1.0), 3.3))
    clean_axes(ax)
    ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
    allv = np.concatenate([a[2] for a in approaches])
    lim = float(np.percentile(np.abs(allv), 97)) * 1.15
    for i, (_label, color, v) in enumerate(approaches):
        _box(ax, i, v, color)
        ax.text(
            i,
            -lim * 0.97,
            rf"{np.mean(np.abs(v)):.0f}",
            ha="center",
            va="bottom",
            fontsize=ANNOT_FS - 1,
            color=INK,
        )
    ax.set_ylim(-lim, lim)
    ax.set_xticks(range(len(approaches)))
    ax.set_xticklabels([a[0] for a in approaches], rotation=30, ha="right")
    ax.set_ylabel(r"$\tilde{\nu}_\mathrm{calc} - \tilde{\nu}_\mathrm{exp}$ / cm$^{-1}$")
    fig.tight_layout()
    save(fig, OUT, "central_vs_experiment")
    register(
        "central_vs_experiment",
        1.0,
        3.3,
        "The central comparison: error against the experimental band origin for scaled-harmonic "
        f"B3LYP (factor {B3LYP_SCALE}), B3LYP VPT2 and ML VPT2 with each energy model, over all "
        "fundamentals with a Shimanouchi assignment. Numbers: MAE in cm$^{-1}$.",
        "5.5 central plot",
    )


def fig_overlap_matrices(mol: str = "methane"):
    from mace_gaussian.analysis.mode_matching import (
        create_alignment_matrix,
        extract_mode_data_from_checkpoint,
    )

    base = COMPARISON / mol
    # local baselines write gaussian_dft.fchk, cluster ones <molecule>_freq_anharm.fchk
    dft_fchk = next(iter(sorted(base.glob("b3lyp*/*.fchk"))), None)
    if dft_fchk is None:
        return
    modes_dft, freqs_dft, _, _, _ = extract_mode_data_from_checkpoint(
        str(dft_fchk), force_harmonic=True
    )
    panels = []
    for m in MODEL_ORDER:
        f = base / f"{m}_mace_ml" / "gaussian_freq.fchk"
        if f.exists():
            modes_ml, freqs_ml, _, _, _ = extract_mode_data_from_checkpoint(
                str(f), force_harmonic=True
            )
            panels.append((m, np.abs(create_alignment_matrix(modes_ml, modes_dft)), freqs_ml))
    if not panels:
        return
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(figwidth(1.0), figwidth(1.0) / n + 0.8), sharey=True)
    axes = np.atleast_1d(axes)
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("ov", ["#ffffff", SERIES[0]])
    for ax, (m, mat, freqs_ml) in zip(axes, panels):
        # rows = B3LYP modes (shared by every panel), columns = this model's own modes
        im = ax.imshow(mat.T, cmap=cmap, vmin=0, vmax=1, aspect="equal", origin="lower")
        ax.set_title(tex(MODEL_LABEL[m]), fontsize=ANNOT_FS, color=model_color(m))
        ax.set_xticks(range(mat.shape[0]))
        ax.set_xticklabels([f"{f:.0f}" for f in freqs_ml], rotation=90, fontsize=ANNOT_FS - 2.5)
        ax.set_yticks(range(mat.shape[1]))
        ax.set_yticklabels([f"{f:.0f}" for f in freqs_dft], fontsize=ANNOT_FS - 2.5)
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(True)
            side.set_edgecolor(INK)
            side.set_linewidth(0.8)
    axes[0].set_ylabel(r"B3LYP harmonic mode / cm$^{-1}$")
    axes[len(axes) // 2].set_xlabel(r"ML harmonic mode / cm$^{-1}$")
    cb = fig.colorbar(im, ax=list(axes), pad=0.02, fraction=0.03)
    cb.set_label("overlap")
    cb.outline.set_visible(False)
    save(fig, OUT, f"overlap_matrices_{mol}")
    register(
        f"overlap_matrices_{mol}",
        1.0,
        figwidth(1.0) / n + 0.8,
        f"{mol}: mass-weighted eigenvector overlap between every ML and B3LYP harmonic normal "
        "mode. Degenerate sets appear as blocks; the Hungarian assignment picks one mode per "
        "row and column, the subspace overlap of a block is what is compared for degenerate modes.",
        "5.3 degeneracy",
    )


SUPERVISOR_EXAMPLES = ("water", "methanol", "formic_acid")


def fig_supervisor():
    """Register the figures written by make_supervisor_figures.py (not redrawn here)."""
    sup = OUT / "supervisor"
    w = figwidth(1.0)
    for m in MODEL_ORDER:
        if (sup / f"grid_pooled_{m}.png").exists():
            register(
                f"supervisor/grid_pooled_{m}",
                1.0,
                w * 1.03,
                f"{MODEL_LABEL[m]} against B3LYP-VPT2, all molecules pooled. "
                "Columns: fundamentals, overtones, combination bands. "
                "Rows: frequency, IR intensity (colour = dipole model), anharmonic constants "
                "x_ii and x_ij from band positions (hollow = resonance-shifted, not in the MAE).",
                "S1 supervisor: pooled by energy model",
            )
    for mol in SUPERVISOR_EXAMPLES:
        if (sup / f"grid_{mol}_mace_omol.png").exists():
            register(
                f"supervisor/grid_{mol}_mace_omol",
                1.0,
                w * 1.03,
                f"{mol}: MACE-OMOL against B3LYP-VPT2 by band type. Every molecule and model: "
                "supervisor/index.html.",
                "S2 supervisor: examples",
            )
    for mol in SUPERVISOR_EXAMPLES:
        if (sup / f"overlap_{mol}.png").exists():
            register(
                f"supervisor/overlap_{mol}",
                1.0,
                figwidth(1.0) / 5 + 0.8,
                f"{mol}: eigenvector overlap of every ML harmonic mode with every B3LYP mode "
                "(the scalar products behind the mode assignment).",
                "S3 supervisor: mode alignment",
            )


# ---------------------------------------------------------------------------
# Gallery
# ---------------------------------------------------------------------------


def write_gallery():
    by_section: dict[str, list[dict]] = {}
    for f in FIGURES:
        by_section.setdefault(f["section"], []).append(f)
    parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Thesis figures</title>",
        "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;"
        "padding:0 1rem;color:#1a1a1a}",
        "h2{border-bottom:1px solid #ddd;padding-bottom:4px;margin-top:2.5rem}",
        ".fig{margin:1.5rem 0;padding:1rem;border:1px solid #e5e7eb;border-radius:6px;"
        "background:#fff}",
        ".fig img{display:block;margin:0 auto;box-shadow:0 0 0 1px #eee}",
        ".meta{color:#6b6b6b;font-size:0.85em;margin:6px 0}.cap{font-size:0.95em}",
        "code{background:#f3f4f6;padding:1px 4px;border-radius:3px}</style></head><body>",
        "<h1>Thesis figures</h1><p>Generated by <code>scripts/make_thesis_figures.py</code>. "
        "Each PNG is shown at its print width (1 cm = 37.8 px at 96 dpi), so text size on screen "
        "matches the page roughly at 100% zoom. The PDF next to each is what LaTeX includes.</p>",
        "<p>Supervisor figures for every molecule and energy model, with full statistics "
        "tables: <a href='supervisor/index.html'>supervisor/index.html</a>.</p>"
        if (OUT / "supervisor" / "index.html").exists()
        else "",
    ]
    px_per_in = 96.0
    for section, figs in by_section.items():
        parts.append(f"<h2>{html.escape(section)}</h2>")
        for f in figs:
            w_in = figwidth(f["width"])
            parts.append(
                f"<div class='fig'><img src='{f['name']}.png' width='{w_in * px_per_in:.0f}'>"
                f"<div class='meta'><code>{f['name']}.pdf</code> · {f['width']:.2f} textwidth = "
                f"{w_in * 2.54:.1f} cm &times; {f['height'] * 2.54:.1f} cm</div>"
                f"<div class='cap'>{html.escape(f['caption'])}</div></div>"
            )
    parts.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(parts), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    global ANALYSIS, COMPARISON, OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--only", help="figure family: morse, spectra, permode, heatmap, intensity, cost"
    )
    ap.add_argument(
        "--molecules", nargs="*", default=None, help="molecules for per-molecule figures"
    )
    ap.add_argument(
        "--analysis-dir", type=Path, default=ANALYSIS, help="per-molecule analysis exports"
    )
    ap.add_argument(
        "--comparison-dir", type=Path, default=COMPARISON, help="raw results (timings, fchk)"
    )
    ap.add_argument("--out-dir", type=Path, default=OUT, help="where figures and gallery go")
    ap.add_argument(
        "--campaign",
        default=None,
        help="read campaigns/<NAME>/ and write campaigns/<NAME>/figures/ (overrides the dirs)",
    )
    ap.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="energy models to leave out of every figure, e.g. --exclude mace_mp",
    )
    args = ap.parse_args(argv)
    ANALYSIS, COMPARISON, OUT = args.analysis_dir, args.comparison_dir, args.out_dir
    if args.campaign is not None:
        from mace_gaussian.campaign import campaign_paths

        paths = campaign_paths(args.campaign)
        ANALYSIS, COMPARISON = REPO / paths.analysis, REPO / paths.comparison
        OUT = paths.figures
    if args.exclude:
        MODEL_ORDER[:] = [m for m in MODEL_ORDER if m not in set(args.exclude)]
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)

    available = sorted(p.parent.name for p in ANALYSIS.glob("*/report_data.json"))
    mols = args.molecules or available
    # The legacy thesis gallery shows five worked examples; any other gallery (a campaign)
    # gets the per-molecule figures (spectra, per-mode errors, overlap) for every molecule.
    legacy_gallery = OUT.resolve() == (REPO / "thesis" / "figures").resolve()
    examples = ("water", "methanol", "formaldehyde", "methane", "octane")
    per_mol = [m for m in examples if m in mols] if legacy_gallery else mols
    overlap_mols = [m for m in ("water",) if m in mols] if legacy_gallery else mols

    families = {
        "spectra": lambda: [fig_stacked_spectra(m) for m in per_mol],
        "permode": lambda: [fig_per_mode_errors(m) for m in per_mol],
        "heatmap": lambda: fig_mae_heatmap(mols),
        "intensity": lambda: fig_intensity_scatter(mols),
        "cost": fig_cost_scaling,
        "residuals": lambda: fig_residuals(mols),
        "anharm": lambda: fig_anharmonicity_ratio(mols),
        "distribution": lambda: fig_error_distribution(mols),
        "bandtype": lambda: fig_band_type_errors(mols),
        "pareto": lambda: fig_pareto(mols),
        "central": lambda: fig_central_vs_experiment(mols),
        "overlap": lambda: [fig_overlap_matrices(m) for m in overlap_mols],
        "supervisor": fig_supervisor,
    }
    for name, fn in families.items():
        if args.only and args.only != name:
            continue
        try:
            fn()
            print(f"  ok   {name}")
        except Exception as e:  # keep going; one broken figure must not block the gallery
            print(f"  FAIL {name}: {e!r}")
    if args.only:
        # the gallery lists only what this run registered; a partial run must not shrink it
        print("--only: figures redrawn, gallery page left as it was")
        return 0
    write_gallery()
    print(f"gallery: {OUT / 'index.html'}  ({len(FIGURES)} figures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
