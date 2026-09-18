"""
Batch Report Generator

Creates a multi-molecule HTML report aggregating accuracy metrics (R^2, RMSE)
across all molecules and calculator combinations. Reads existing comparison_results/
JSON files without re-computation.

Report sections:
  1. Accuracy Leaderboard (primary view)
  2. Summary Heatmap (RMSE or R^2 by combo x molecule)
  3. Box Plots (RMSE distribution per calculator combo)
  4. Size Scaling (RMSE vs atom count)
  5. Per-Molecule Spectrum Overlays
"""

import base64
import contextlib
import html
import json
import logging
import os
import tempfile
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ._shared_css import build_css
from .nist_fetcher import fetch_experimental_spectrum


def _esc(value: object) -> str:
    """HTML-escape a value for safe interpolation (T-23-01 mitigation)."""
    return html.escape(str(value), quote=True)

logger = logging.getLogger(__name__)

DPI = 300
FONT_FAMILY = "Arial, Helvetica, sans-serif"


def aggregate_results(
    results_dir: str = "comparison_results",
    analysis_dir: str = "analysis_results_harmonic",
) -> pd.DataFrame:
    """Walk comparison_results/ and collect R^2/RMSE per molecule x combo.

    Metrics are taken from the eigenvector-matched harmonic comparison CSVs written by
    ``run_analysis_harmonic.py`` (``<analysis_dir>/<molecule>/data/comparison_<combo>.csv``)
    when they exist, so the leaderboard agrees with the per-molecule reports (review
    finding M5). If no CSV exists for a combo the old sorted-frequency pairing is used
    and the row is marked ``pairing = "sorted"``.

    Parameters
    ----------
    results_dir : str
        Path to directory containing per-molecule subdirectories
        with results.json files.
    analysis_dir : str
        Base directory of the harmonic analysis output.

    Returns
    -------
    pd.DataFrame
        Columns: molecule, combo, r2, rmse, n_atoms, n_freqs, pairing, ...
    """
    results_path = Path(results_dir)
    rows: list[dict] = []
    empty = pd.DataFrame(
        columns=["molecule", "combo", "r2", "rmse", "n_atoms", "n_freqs", "pairing"]
    )

    if not results_path.is_dir():
        return empty

    for mol_dir in sorted(results_path.iterdir()):
        if not mol_dir.is_dir():
            continue
        molecule = mol_dir.name

        # Load DFT reference
        dft_json = mol_dir / "b3lyp_6-31Gdp" / "results.json"
        if not dft_json.exists():
            logger.debug("No DFT reference for %s, skipping", molecule)
            continue

        try:
            with dft_json.open() as f:
                dft_data = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Failed to load DFT results for %s: %s", molecule, e)
            continue

        dft_freqs = sorted(
            entry["freq_cm"] for entry in dft_data.get("frequencies", {}).get("harmonic", [])
        )
        if not dft_freqs:
            logger.debug("No DFT harmonic frequencies for %s", molecule)
            continue

        # Count atoms from optimized geometry
        n_atoms = _count_atoms(mol_dir)

        # Extract DFT timing and hardware
        dft_timing = dft_data.get("timing", dft_data.get("gaussian_timing", {}))
        dft_runtime = dft_data.get("runtime_s", 0)
        dft_gauss_s = dft_timing.get("total_elapsed_s", 0) if dft_timing else 0
        dft_hw = dft_data.get("hardware", {})
        dft_version = dft_data.get("version_info", {})
        dft_cpu = dft_hw.get("cpu_model", dft_version.get("cpu_model", ""))
        dft_node = dft_hw.get("node", "")

        # Process each ML combo
        for combo_dir in sorted(mol_dir.iterdir()):
            if not combo_dir.is_dir():
                continue
            combo_name = combo_dir.name
            if combo_name in ("b3lyp_6-31Gdp", "geometry_opt"):
                continue

            matched_csv = Path(analysis_dir) / molecule / "data" / f"comparison_{combo_name}.csv"
            row = _compute_combo_metrics(
                molecule, combo_name, combo_dir, dft_freqs, n_atoms, matched_csv=matched_csv
            )
            if row is not None:
                row["dft_runtime_s"] = dft_runtime
                row["dft_gaussian_s"] = dft_gauss_s
                row["dft_cpu"] = dft_cpu
                row["dft_node"] = dft_node
                row["speedup"] = (
                    dft_gauss_s / row["ml_runtime_s"]
                    if row.get("ml_runtime_s")
                    else 0
                )
                rows.append(row)

    return pd.DataFrame(rows) if rows else empty


def _metrics_from_matched_csv(csv_path: Path) -> tuple[float, float, int] | None:
    """(r2, rmse, n_pairs) from an eigenvector-matched comparison CSV.

    Uses the same conventions as the per-molecule report: Pearson r² on the matched
    pairs, pairs with an imaginary (negative) frequency on either side excluded.
    Returns None if the file is unreadable or has fewer than two usable pairs.
    """
    try:
        df = pd.read_csv(csv_path)
        dft = df["DFT_Frequency_cm"].to_numpy(dtype=float)
        ml = df["ML_Frequency_cm"].to_numpy(dtype=float)
    except (OSError, KeyError, ValueError) as e:
        logger.warning("Could not read matched CSV %s: %s", csv_path, e)
        return None
    real = (dft > 0) & (ml > 0)
    dft, ml = dft[real], ml[real]
    if len(dft) < 2:
        return None
    r = np.corrcoef(dft, ml)[0, 1]
    r2 = float(r * r) if np.isfinite(r) else 0.0
    rmse = float(np.sqrt(np.mean((dft - ml) ** 2)))
    return r2, rmse, len(dft)


def _compute_combo_metrics(
    molecule: str,
    combo_name: str,
    combo_dir: Path,
    dft_freqs: list[float],
    n_atoms: int,
    matched_csv: Path | None = None,
) -> dict | None:
    """Compute R^2 and RMSE for a single molecule/combo pair.

    Prefers the eigenvector-matched CSV (review finding M5); falls back to
    sorted-frequency pairing only when no CSV is available.
    """
    ml_json = combo_dir / "results.json"
    if not ml_json.exists():
        return None

    try:
        with ml_json.open() as f:
            ml_data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(
            "Failed to load ML results for %s/%s: %s",
            molecule,
            combo_name,
            e,
        )
        return None

    pairing = "sorted"
    n_pairs = 0
    if matched_csv is not None and matched_csv.exists():
        matched = _metrics_from_matched_csv(matched_csv)
        if matched is not None:
            r2, rmse, n_pairs = matched
            pairing = "eigenvector"

    if pairing == "sorted":
        logger.warning(
            "%s/%s: no matched comparison CSV, using sorted-frequency pairing "
            "(run run_analysis_harmonic.py first for mode-matched metrics)",
            molecule,
            combo_name,
        )
        ml_freqs = sorted(
            entry["freq_cm"] for entry in ml_data.get("frequencies", {}).get("harmonic", [])
        )

        if len(dft_freqs) != len(ml_freqs) or len(dft_freqs) == 0:
            logger.debug(
                "Frequency count mismatch for %s/%s: DFT=%d ML=%d",
                molecule,
                combo_name,
                len(dft_freqs),
                len(ml_freqs),
            )
            return None

        dft_arr = np.array(dft_freqs)
        ml_arr = np.array(ml_freqs)

        ss_res = np.sum((dft_arr - ml_arr) ** 2)
        ss_tot = np.sum((dft_arr - np.mean(dft_arr)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        rmse = float(np.sqrt(np.mean((dft_arr - ml_arr) ** 2)))
        n_pairs = len(dft_freqs)

    # Extract timing data
    ml_runtime = ml_data.get("runtime_s", 0)
    ml_gauss = ml_data.get("gaussian_timing", {})
    ml_gauss_s = ml_gauss.get("total_elapsed_s", 0) if ml_gauss else 0

    # Hardware
    version = ml_data.get("version_info", {})
    gpu = version.get("gpu_name", "")
    cpu = version.get("cpu_model", "")

    return {
        "molecule": molecule,
        "combo": combo_name,
        "r2": r2,
        "rmse": rmse,
        "n_atoms": n_atoms,
        "n_freqs": n_pairs,
        "pairing": pairing,
        "ml_runtime_s": ml_runtime,
        "ml_gaussian_s": ml_gauss_s,
        "ml_gpu": gpu,
        "ml_cpu": cpu,
    }


def _count_atoms(mol_dir: Path) -> int:
    """Count atoms from optimized.xyz or return 0 if unavailable."""
    opt_xyz = mol_dir / "geometry_opt" / "optimized.xyz"
    if opt_xyz.exists():
        try:
            from ase.io import read as ase_read

            atoms = ase_read(str(opt_xyz))
            return len(atoms)
        except Exception:
            # Fallback: parse first line of XYZ file
            try:
                with opt_xyz.open() as f:
                    return int(f.readline().strip())
            except (ValueError, OSError):
                pass
    return 0


def _save_figure(fig: plt.Figure) -> str:
    """Save matplotlib figure to a temporary PNG and return the path."""
    fd, tmp_path = tempfile.mkstemp(suffix=".png")
    os.close(fd)
    fig.savefig(tmp_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return tmp_path


def _plot_heatmap(df: pd.DataFrame, metric: str = "rmse") -> str:
    """Create summary heatmap of metric values (combo x molecule).

    Returns path to saved PNG file.
    """
    pivot = df.pivot_table(index="combo", columns="molecule", values=metric, aggfunc="first")

    if metric == "rmse":
        cmap = "YlOrRd"
        fmt = ".1f"
        title = "RMSE Heatmap (cm$^{-1}$) -- Lower is Better"
    else:
        cmap = "YlGn"
        fmt = ".3f"
        title = "R$^2$ Heatmap -- Higher is Better"

    n_combos = len(pivot.index)
    n_mols = len(pivot.columns)
    fig_width = max(8, 1.5 * n_mols + 3)
    fig_height = max(4, 0.6 * n_combos + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=fmt,
        cmap=cmap,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": metric.upper()},
    )
    ax.set_title(title, fontsize=14, fontfamily="sans-serif", pad=12)
    ax.set_ylabel("Calculator Combo")
    ax.set_xlabel("Molecule")
    plt.tight_layout()

    return _save_figure(fig)


def _plot_boxplots(df: pd.DataFrame) -> str:
    """Create box plots of RMSE distribution per calculator combo.

    Returns path to saved PNG file.
    """
    combo_order = df.groupby("combo")["rmse"].median().sort_values().index.tolist()

    height = max(4, 0.5 * len(combo_order) + 2)
    fig, ax = plt.subplots(figsize=(10, height))
    palette = sns.color_palette("colorblind", n_colors=len(combo_order))

    sns.boxplot(
        data=df,
        x="rmse",
        y="combo",
        order=combo_order,
        palette=palette,
        ax=ax,
        width=0.6,
    )
    sns.stripplot(
        data=df,
        x="rmse",
        y="combo",
        order=combo_order,
        color="0.3",
        size=5,
        alpha=0.6,
        ax=ax,
    )
    ax.set_title(
        "RMSE Distribution by Calculator Combo",
        fontsize=14,
        fontfamily="sans-serif",
    )
    ax.set_xlabel("RMSE (cm$^{-1}$)")
    ax.set_ylabel("")
    plt.tight_layout()

    return _save_figure(fig)


def _plot_size_scaling(df: pd.DataFrame) -> str:
    """Create RMSE vs atom count trend plot.

    Returns path to saved PNG file.
    """
    plot_df = df[df["n_atoms"] > 0].copy()

    fig, ax = plt.subplots(figsize=(10, 6))
    n_combos = df["combo"].nunique()
    palette = sns.color_palette("colorblind", n_colors=n_combos)

    if not plot_df.empty:
        sns.lineplot(
            data=plot_df,
            x="n_atoms",
            y="rmse",
            hue="combo",
            marker="o",
            palette=palette,
            ax=ax,
        )
        ax.legend(
            title="Calculator Combo",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No atom count data available",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )

    ax.set_title(
        "Size Scaling: RMSE vs Atom Count",
        fontsize=14,
        fontfamily="sans-serif",
    )
    ax.set_xlabel("Number of Atoms")
    ax.set_ylabel("RMSE (cm$^{-1}$)")
    plt.tight_layout()

    return _save_figure(fig)


def _plot_spectrum_overlay(molecule: str, results_dir: str) -> str | None:
    """Create per-molecule spectrum overlay (ML combos vs DFT).

    Returns path to saved PNG file, or None if DFT data is missing.
    """
    results_path = Path(results_dir)
    mol_dir = results_path / molecule

    # Load DFT reference
    dft_json = mol_dir / "b3lyp_6-31Gdp" / "results.json"
    if not dft_json.exists():
        return None

    try:
        with dft_json.open() as f:
            dft_data = json.load(f)
    except (json.JSONDecodeError, OSError):
        return None

    dft_harmonic = dft_data.get("frequencies", {}).get("harmonic", [])
    if not dft_harmonic:
        return None

    fig, ax = plt.subplots(figsize=(12, 5))
    palette = sns.color_palette("colorblind", n_colors=10)

    # Plot DFT as black stems
    dft_freqs = [e["freq_cm"] for e in dft_harmonic]
    dft_ints = [e.get("ir_intensity", 1.0) for e in dft_harmonic]
    ax.stem(
        dft_freqs,
        dft_ints,
        linefmt="k-",
        markerfmt="ko",
        basefmt="k-",
        label="DFT (B3LYP/6-31G(d,p))",
    )

    # Load experimental spectrum if cached (best-effort, per D-13)
    try:
        experimental = fetch_experimental_spectrum(molecule, cache_dir=mol_dir)
    except Exception:
        experimental = None

    if experimental is not None:
        # Filter to plot range
        mask = (experimental.wavenumbers >= 400) & (experimental.wavenumbers <= 4200)
        if np.sum(mask) > 0:
            from scipy.interpolate import interp1d

            exp_interp = interp1d(
                experimental.wavenumbers[mask],
                experimental.absorbance[mask],
                kind="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            # Scale experimental to DFT intensity range for visual comparison
            dft_int_max = max(dft_ints) if dft_ints else 1.0
            freq_grid = np.linspace(400, 4200, 1000)
            exp_on_grid = exp_interp(freq_grid)
            exp_max = np.max(exp_on_grid)
            exp_scaled = exp_on_grid * dft_int_max / exp_max if exp_max > 0 else exp_on_grid
            ax.plot(
                freq_grid,
                exp_scaled,
                color="#888888",
                linestyle="--",
                linewidth=0.8,
                label=f"Experimental ({experimental.source})",
                alpha=0.5,
                zorder=3,
            )

    # Plot each ML combo
    color_idx = 0
    for combo_dir in sorted(mol_dir.iterdir()):
        if not combo_dir.is_dir():
            continue
        if combo_dir.name in ("b3lyp_6-31Gdp", "geometry_opt"):
            continue

        ml_json = combo_dir / "results.json"
        if not ml_json.exists():
            continue

        try:
            with ml_json.open() as f:
                ml_data = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue

        ml_harmonic = ml_data.get("frequencies", {}).get("harmonic", [])
        if not ml_harmonic:
            continue

        ml_freqs = [e["freq_cm"] for e in ml_harmonic]
        ml_ints = [e.get("ir_intensity", 1.0) for e in ml_harmonic]

        color = palette[color_idx % len(palette)]
        markerline, stemlines, baseline = ax.stem(
            ml_freqs,
            ml_ints,
            linefmt="-",
            markerfmt="o",
            basefmt="-",
            label=combo_dir.name,
        )
        plt.setp(stemlines, color=color, alpha=0.7)
        plt.setp(markerline, color=color, markersize=4)
        plt.setp(baseline, visible=False)
        color_idx += 1

    ax.set_title(
        f"Spectrum Overlay: {molecule}",
        fontsize=14,
        fontfamily="sans-serif",
    )
    ax.set_xlabel("Frequency (cm$^{-1}$)")
    ax.set_ylabel("IR Intensity")
    ax.set_xlim(400, 4200)
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()

    return _save_figure(fig)


def _encode_plot(path: str) -> str:
    """Encode a plot file as base64 data URI."""
    with Path(path).open("rb") as f:
        encoded = base64.b64encode(f.read()).decode()
    return f"data:image/png;base64,{encoded}"


def _generate_html(df: pd.DataFrame, plot_paths: dict, results_dir: str) -> str:
    """Build self-contained HTML report with embedded base64 plots.

    Parameters
    ----------
    df : pd.DataFrame
        Aggregated results (molecule, combo, r2, rmse, n_atoms, n_freqs)
    plot_paths : dict
        Mapping of plot name to file path
    results_dir : str
        Path to comparison_results directory

    Returns
    -------
    str
        Complete HTML string
    """
    leaderboard_html = _build_leaderboard_html(df)

    # Embed plots as base64
    embedded = {}
    for name, path in plot_paths.items():
        if path is not None and Path(path).exists():
            embedded[name] = _encode_plot(path)

    # Per-molecule spectrum section
    spectrum_parts = []
    for name, data_uri in sorted(embedded.items()):
        if name.startswith("spectrum_"):
            mol_name = name.replace("spectrum_", "")
            mol_esc = _esc(mol_name)
            spectrum_parts.append(
                f'<div class="plot-card">'
                f"<h3>{mol_esc}</h3>"
                f'<img src="{data_uri}" '
                f'alt="Spectrum overlay for {mol_esc}">'
                f"</div>"
            )
    spectrum_section = (
        "\n".join(spectrum_parts)
        if spectrum_parts
        else "<p>No spectrum overlay data available.</p>"
    )

    # Embed summary plots
    heatmap_img = _embed_or_fallback(embedded, "heatmap_rmse", "RMSE Heatmap", "No heatmap data.")
    heatmap_r2_img = _embed_or_fallback(embedded, "heatmap_r2", "R2 Heatmap", "")
    boxplot_img = _embed_or_fallback(embedded, "boxplots", "Box Plots", "No box plot data.")
    size_scaling_img = _embed_or_fallback(
        embedded, "size_scaling", "Size Scaling", "No size scaling data."
    )

    generated_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    n_molecules = df["molecule"].nunique()
    n_combos = df["combo"].nunique()

    css = build_css()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MACE-Gaussian Batch Report</title>
{css}
</head>
<body>

<nav>
    <a href="#leaderboard">Leaderboard</a>
    <a href="#timing">Timing</a>
    <a href="#heatmap">Heatmap</a>
    <a href="#boxplots">Box Plots</a>
    <a href="#size-scaling">Size Scaling</a>
    <a href="#spectra">Per-Molecule Spectra</a>
</nav>

<div class="container">

<h1>MACE-Gaussian Batch Accuracy Report</h1>
<p class="subtitle">
  Aggregated comparison of ML calculator combinations
  across {n_molecules} molecules
</p>

<div class="stats-bar">
    <div class="stat-box">
        <div class="value">{n_molecules}</div>
        <div class="label">Molecules</div>
    </div>
    <div class="stat-box">
        <div class="value">{n_combos}</div>
        <div class="label">Calculator Combos</div>
    </div>
    <div class="stat-box">
        <div class="value">{len(df)}</div>
        <div class="label">Total Comparisons</div>
    </div>
</div>

<h2 id="leaderboard">Accuracy Leaderboard</h2>
<table>
<thead>
<tr>
    <th>Rank</th>
    <th>Calculator Combo</th>
    <th>Mean R&sup2;</th>
    <th>Mean RMSE (cm<sup>-1</sup>)</th>
    <th>Median RMSE (cm<sup>-1</sup>)</th>
    <th>N Molecules</th>
</tr>
</thead>
<tbody>
{leaderboard_html}
</tbody>
</table>

{_build_timing_html(df)}

<h2 id="heatmap">Summary Heatmap</h2>
<div class="plot-section">
    <div class="plot-card">
        <h3>RMSE by Calculator Combo and Molecule</h3>
        {heatmap_img}
    </div>
    <div class="plot-card">
        <h3>R&sup2; by Calculator Combo and Molecule</h3>
        {heatmap_r2_img}
    </div>
</div>

<h2 id="boxplots">Box Plots</h2>
<div class="plot-section">
    <div class="plot-card">
        <h3>RMSE Distribution per Calculator Combo</h3>
        {boxplot_img}
    </div>
</div>

<h2 id="size-scaling">Size Scaling</h2>
<div class="plot-section">
    <div class="plot-card">
        <h3>RMSE vs Number of Atoms</h3>
        {size_scaling_img}
    </div>
</div>

<h2 id="spectra">Per-Molecule Spectra</h2>
<div class="plot-section">
    {spectrum_section}
</div>

</div>

<footer>
    Generated by mace-gaussian report on {generated_date}
</footer>

</body>
</html>"""

    return html


def _build_timing_html(df: pd.DataFrame) -> str:
    """Build timing comparison table with hardware context."""
    if "ml_gaussian_s" not in df.columns:
        return ""

    has_timing = df["ml_gaussian_s"].sum() > 0 or df.get("dft_gaussian_s", pd.Series([0])).sum() > 0
    if not has_timing:
        return (
            '<h2 id="timing">Timing</h2>'
            "<p>No Gaussian timing data available yet.</p>"
        )

    rows = []
    for _, r in df.iterrows():
        ml_t = r.get("ml_gaussian_s", 0)
        dft_t = r.get("dft_gaussian_s", 0)
        speedup = r.get("speedup", 0)
        ml_hw = r.get("ml_gpu", "") or r.get("ml_cpu", "") or "—"
        dft_hw = r.get("dft_cpu", "") or "—"
        if r.get("dft_node"):
            dft_hw += f" ({r['dft_node']})"

        rows.append(
            f"<tr>"
            f"<td>{_esc(r['molecule'])}</td>"
            f"<td>{_esc(r['combo'])}</td>"
            f"<td>{ml_t:.1f}</td>"
            f"<td>{dft_t:.1f}</td>"
            f"<td>{speedup:.1f}x</td>"
            f"<td style='font-size:0.8em'>{_esc(ml_hw)}</td>"
            f"<td style='font-size:0.8em'>{_esc(dft_hw)}</td>"
            f"</tr>"
        )

    return f"""
<h2 id="timing">Timing &amp; Hardware</h2>
<p style="color: #666; font-size: 0.9em;">
    Gaussian wall-clock time (elapsed) as reported in log files.
    Speedup = DFT time / ML time. Hardware shown for context &mdash;
    speedup numbers are only meaningful when comparing similar hardware.
</p>
<table>
<thead>
<tr>
    <th>Molecule</th>
    <th>ML Combo</th>
    <th>ML (s)</th>
    <th>DFT (s)</th>
    <th>Speedup</th>
    <th>ML Hardware</th>
    <th>DFT Hardware</th>
</tr>
</thead>
<tbody>
{"".join(rows)}
</tbody>
</table>
"""


def _build_leaderboard_html(df: pd.DataFrame) -> str:
    """Build HTML table rows for the accuracy leaderboard."""
    leaderboard = (
        df.groupby("combo")
        .agg(
            mean_r2=("r2", "mean"),
            mean_rmse=("rmse", "mean"),
            median_rmse=("rmse", "median"),
            n_molecules=("molecule", "nunique"),
        )
        .sort_values("mean_rmse")
        .reset_index()
    )

    rows = []
    n_total = len(leaderboard)
    for rank, (_, row) in enumerate(leaderboard.iterrows(), 1):
        if rank == 1:
            cls = ' class="best"'
        elif rank == n_total:
            cls = ' class="worst"'
        else:
            cls = ""
        rows.append(
            f"<tr{cls}>"
            f"<td>{rank}</td>"
            f"<td>{_esc(row['combo'])}</td>"
            f"<td>{row['mean_r2']:.4f}</td>"
            f"<td>{row['mean_rmse']:.1f}</td>"
            f"<td>{row['median_rmse']:.1f}</td>"
            f"<td>{int(row['n_molecules'])}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def _embed_or_fallback(embedded: dict, key: str, alt: str, fallback_text: str) -> str:
    """Return an <img> tag or fallback <p> for an embedded plot."""
    if key in embedded:
        return f'<img src="{embedded[key]}" alt="{alt}">'
    if fallback_text:
        return f"<p>{fallback_text}</p>"
    return ""


def generate_batch_report(
    results_dir: str = "comparison_results",
    output_dir: str = "batch_report",
    analysis_dir: str = "analysis_results_harmonic",
) -> str:
    """Generate multi-molecule batch accuracy report.

    Parameters
    ----------
    results_dir : str
        Path to comparison_results directory.
    output_dir : str
        Output directory for the HTML report.
    analysis_dir : str
        Harmonic analysis output; its mode-matched CSVs supply the metrics (M5).

    Returns
    -------
    str
        Path to the generated batch_report.html file.

    Raises
    ------
    ValueError
        If no comparison results found in results_dir.
    """
    df = aggregate_results(results_dir, analysis_dir=analysis_dir)
    if df.empty:
        raise ValueError(f"No comparison results found in {results_dir}")

    logger.info(
        "Aggregated %d results across %d molecules and %d combos",
        len(df),
        df["molecule"].nunique(),
        df["combo"].nunique(),
    )

    # Generate all plots
    plot_paths: dict[str, str | None] = {}
    temp_files: list[str] = []

    try:
        # Heatmaps
        path = _plot_heatmap(df, metric="rmse")
        plot_paths["heatmap_rmse"] = path
        temp_files.append(path)

        path = _plot_heatmap(df, metric="r2")
        plot_paths["heatmap_r2"] = path
        temp_files.append(path)

        # Box plots
        path = _plot_boxplots(df)
        plot_paths["boxplots"] = path
        temp_files.append(path)

        # Size scaling
        path = _plot_size_scaling(df)
        plot_paths["size_scaling"] = path
        temp_files.append(path)

        # Per-molecule spectrum overlays
        for molecule in sorted(df["molecule"].unique()):
            path = _plot_spectrum_overlay(molecule, results_dir)
            if path is not None:
                plot_paths[f"spectrum_{molecule}"] = path
                temp_files.append(path)

        # Generate HTML
        html = _generate_html(df, plot_paths, results_dir)

        # Write output
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        report_file = output_path / "batch_report.html"
        report_file.write_text(html, encoding="utf-8")

        logger.info("Batch report written to %s", report_file)
        return str(report_file)

    finally:
        # Clean up temp plot files
        for tmp in temp_files:
            with contextlib.suppress(OSError):
                Path(tmp).unlink()
