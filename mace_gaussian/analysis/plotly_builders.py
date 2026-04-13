"""Pure-function Plotly figure builders for IR spectral analysis.

Phase 23 D-06: all figures built as pure functions (no I/O, no HTML).
Phase 23 D-08: these replace the old matplotlib plot functions for interactive reports.

Each builder takes already-computed numpy arrays and returns a ``plotly.graph_objects.Figure``.
Plotly is imported function-locally to keep CLI startup fast.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# Color conventions: DFT is the reference (black), ML is orange, experimental is grey
_DFT_COLOR = "#000000"
_ML_COLOR = "#DE8F05"
_EXP_COLOR = "rgba(128,128,128,0.15)"


def build_spectrum_figure(
    freq_grid: np.ndarray,
    dft_norm: np.ndarray,
    ml_norm: np.ndarray,
    ml_name: str,
    experimental_norm: np.ndarray | None = None,
    offset: float = 1.5,
) -> object:
    """Build an IR spectrum comparison figure (DFT vs ML, optional experimental).

    Parameters
    ----------
    freq_grid : np.ndarray
        Wavenumber grid (cm-1).
    dft_norm : np.ndarray
        Normalized DFT absorbance on ``freq_grid``.
    ml_norm : np.ndarray
        Normalized ML absorbance on ``freq_grid``.
    ml_name : str
        Name of the ML method for labels.
    experimental_norm : np.ndarray or None
        Normalized experimental absorbance on ``freq_grid``, or None.
    offset : float
        Vertical offset for ML trace (DFT above, ML below convention).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=freq_grid,
            y=dft_norm,
            name="DFT",
            line=dict(color=_DFT_COLOR, width=2),
            hovertemplate="%{x:.1f} cm\u207b\u00b9<br>%{y:.3f}<extra>DFT</extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=freq_grid,
            y=ml_norm + offset,
            name=f"ML ({ml_name})",
            line=dict(color=_ML_COLOR, width=2),
            hovertemplate="%{x:.1f} cm\u207b\u00b9<br>%{y:.3f}<extra>ML</extra>",
        )
    )
    if experimental_norm is not None:
        fig.add_trace(
            go.Scatter(
                x=freq_grid,
                y=experimental_norm,
                name="Experimental (NIST)",
                line=dict(color=_EXP_COLOR, width=1.5),
                hovertemplate=("%{x:.1f} cm\u207b\u00b9<br>%{y:.3f}<extra>Exp.</extra>"),
            )
        )
    fig.update_layout(
        title=f"IR Spectrum: {ml_name} vs DFT",
        xaxis_title="Wavenumber (cm\u207b\u00b9)",
        yaxis_title="Absorbance (normalized)",
        xaxis=dict(autorange="reversed"),
        hovermode="x unified",
        template="simple_white",
        height=500,
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig


def build_regression_figure(
    dft_freqs: np.ndarray,
    ml_freqs: np.ndarray,
    ml_name: str,
    mode_ids: list[str] | None = None,
) -> object:
    """Build a regression scatter plot (ML vs DFT frequencies).

    Includes a y=x perfect-agreement reference line.  When *mode_ids* is
    supplied, points are color-coded by type (fundamental / overtone /
    combination).

    Parameters
    ----------
    dft_freqs : np.ndarray
        DFT frequencies (cm-1).
    ml_freqs : np.ndarray
        ML frequencies (cm-1).
    ml_name : str
        Name of the ML method for axis labels.
    mode_ids : list[str] or None
        Matched mode IDs (e.g. ``["F1", "F2", "O1_2", "C1_2"]``).
        Used to color-code by type.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    all_vals = np.concatenate([dft_freqs, ml_freqs])
    lo = float(np.min(all_vals))
    hi = float(np.max(all_vals))
    margin = (hi - lo) * 0.05
    lo -= margin
    hi += margin

    # Overall R²
    n_shown = len(dft_freqs)
    ss_res = float(np.sum((ml_freqs - dft_freqs) ** 2))
    ss_tot = float(np.sum((dft_freqs - np.mean(dft_freqs)) ** 2))
    r2_all = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    _type_style = {
        "fundamental": {"color": _ML_COLOR, "symbol": "circle", "name": "Fundamental"},
        "overtone": {"color": "#0173B2", "symbol": "diamond", "name": "Overtone"},
        "combination": {"color": "#CC78BC", "symbol": "square", "name": "Combination"},
    }

    fig = go.Figure()

    if mode_ids is not None and len(mode_ids) == len(dft_freqs):
        # Group by type
        groups: dict[str, tuple[list, list]] = {}
        for i, mid in enumerate(mode_ids):
            if mid.startswith("O"):
                t = "overtone"
            elif mid.startswith("C"):
                t = "combination"
            else:
                t = "fundamental"
            groups.setdefault(t, ([], []))
            groups[t][0].append(dft_freqs[i])
            groups[t][1].append(ml_freqs[i])

        for t in ("fundamental", "overtone", "combination"):
            if t not in groups:
                continue
            dx, mx = groups[t]
            style = _type_style[t]
            fig.add_trace(
                go.Scatter(
                    x=dx,
                    y=mx,
                    mode="markers",
                    name=f"{style['name']} (n={len(dx)})",
                    marker=dict(
                        color=style["color"],
                        symbol=style["symbol"],
                        size=8,
                        line=dict(width=1, color="#333"),
                    ),
                    hovertemplate=(
                        "DFT: %{x:.1f} cm\u207b\u00b9<br>"
                        "ML: %{y:.1f} cm\u207b\u00b9"
                        f"<extra>{style['name']}</extra>"
                    ),
                )
            )
    else:
        fig.add_trace(
            go.Scatter(
                x=dft_freqs,
                y=ml_freqs,
                mode="markers",
                name=f"Frequencies (n={n_shown})",
                marker=dict(color=_ML_COLOR, size=8, line=dict(width=1, color="#333")),
                hovertemplate=(
                    "DFT: %{x:.1f} cm\u207b\u00b9<br>"
                    "ML: %{y:.1f} cm\u207b\u00b9<extra></extra>"
                ),
            )
        )
    # y=x perfect agreement line
    r2_str = f"{r2_all:.3f}" if not np.isnan(r2_all) else "N/A"
    fig.add_trace(
        go.Scatter(
            x=[lo, hi],
            y=[lo, hi],
            mode="lines",
            name=f"y=x | R\u00b2={r2_str}, n={n_shown}",
            line=dict(color="#888", dash="dash", width=1.5),
            hoverinfo="skip",
        )
    )
    fig.update_layout(
        title=f"Frequency Regression: {ml_name} vs DFT",
        xaxis_title="DFT frequencies (cm\u207b\u00b9)",
        yaxis_title="ML frequencies (cm\u207b\u00b9)",
        xaxis=dict(range=[lo, hi], constrain="domain"),
        yaxis=dict(range=[lo, hi], scaleanchor="x", scaleratio=1, constrain="domain"),
        template="simple_white",
        height=550,
        width=550,
        margin=dict(l=60, r=20, t=60, b=50),
        legend=dict(x=0.05, y=0.98, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0.8)"),
    )
    return fig


def build_intensity_regression_figure(
    dft_intensities: np.ndarray,
    ml_intensities: np.ndarray,
    ml_name: str,
    mode_ids: list[str] | None = None,
) -> object:
    """Build a regression scatter plot (ML vs DFT intensities).

    When *mode_ids* is supplied, points are color-coded by type
    (fundamental / overtone / combination).

    Parameters
    ----------
    dft_intensities : np.ndarray
        DFT IR intensities (km/mol).
    ml_intensities : np.ndarray
        ML IR intensities (km/mol).
    ml_name : str
        Name of the ML method for axis labels.
    mode_ids : list[str] or None
        Matched mode IDs for color-coding by type.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    all_vals = np.concatenate([dft_intensities, ml_intensities])
    lo = float(np.min(all_vals))
    hi = float(np.max(all_vals))
    margin = (hi - lo) * 0.05
    lo = max(0, lo - margin)
    hi += margin

    # Overall R²
    n_shown = len(dft_intensities)
    ss_res = float(np.sum((ml_intensities - dft_intensities) ** 2))
    ss_tot = float(np.sum((dft_intensities - np.mean(dft_intensities)) ** 2))
    r2_all = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    _type_style = {
        "fundamental": {"color": _ML_COLOR, "symbol": "circle", "name": "Fundamental"},
        "overtone": {"color": "#0173B2", "symbol": "diamond", "name": "Overtone"},
        "combination": {"color": "#CC78BC", "symbol": "square", "name": "Combination"},
    }

    fig = go.Figure()

    if mode_ids is not None and len(mode_ids) == len(dft_intensities):
        groups: dict[str, tuple[list, list]] = {}
        for i, mid in enumerate(mode_ids):
            if mid.startswith("O"):
                t = "overtone"
            elif mid.startswith("C"):
                t = "combination"
            else:
                t = "fundamental"
            groups.setdefault(t, ([], []))
            groups[t][0].append(dft_intensities[i])
            groups[t][1].append(ml_intensities[i])

        for t in ("fundamental", "overtone", "combination"):
            if t not in groups:
                continue
            dx, mx = groups[t]
            style = _type_style[t]
            fig.add_trace(
                go.Scatter(
                    x=dx,
                    y=mx,
                    mode="markers",
                    name=f"{style['name']} (n={len(dx)})",
                    marker=dict(
                        color=style["color"],
                        symbol=style["symbol"],
                        size=8,
                        line=dict(width=1, color="#333"),
                    ),
                    hovertemplate=(
                        "DFT: %{x:.1f} km/mol<br>"
                        "ML: %{y:.1f} km/mol"
                        f"<extra>{style['name']}</extra>"
                    ),
                )
            )
    else:
        fig.add_trace(
            go.Scatter(
                x=dft_intensities,
                y=ml_intensities,
                mode="markers",
                name=f"Intensities (n={n_shown})",
                marker=dict(color=_ML_COLOR, size=8, line=dict(width=1, color="#333")),
                hovertemplate=(
                    "DFT: %{x:.1f} km/mol<br>ML: %{y:.1f} km/mol<extra></extra>"
                ),
            )
        )
    r2_str = f"{r2_all:.3f}" if not np.isnan(r2_all) else "N/A"
    fig.add_trace(
        go.Scatter(
            x=[lo, hi],
            y=[lo, hi],
            mode="lines",
            name=f"y=x | R\u00b2={r2_str}, n={n_shown}",
            line=dict(color="#888", dash="dash", width=1.5),
            hoverinfo="skip",
        )
    )
    fig.update_layout(
        title=f"Intensity Regression: {ml_name} vs DFT",
        xaxis_title="DFT intensities (km/mol)",
        yaxis_title="ML intensities (km/mol)",
        xaxis=dict(range=[lo, hi], constrain="domain"),
        yaxis=dict(range=[lo, hi], scaleanchor="x", scaleratio=1, constrain="domain"),
        template="simple_white",
        height=550,
        width=550,
        margin=dict(l=60, r=20, t=60, b=50),
        legend=dict(x=0.05, y=0.98, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0.8)"),
    )
    return fig


def build_combined_spectrum_figure(
    freq_grid: np.ndarray,
    dft_norm: np.ndarray,
    ml_norms: dict[str, np.ndarray],
    experimental_norm: np.ndarray | None = None,
    offset: float = 1.5,
) -> object:
    """Build a combined spectrum figure with one trace per ML method plus DFT.

    ML traces are stacked with increasing vertical offset for visual separation.

    Parameters
    ----------
    freq_grid : np.ndarray
        Wavenumber grid (cm-1).
    dft_norm : np.ndarray
        Normalized DFT absorbance.
    ml_norms : dict[str, np.ndarray]
        Mapping of ML method name to normalized absorbance on ``freq_grid``.
    experimental_norm : np.ndarray or None
        Normalized experimental absorbance on ``freq_grid``, or None.
    offset : float
        Vertical spacing between stacked traces.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    # Cycle through distinguishable colors for multiple ML methods
    _ml_colors = [
        "#E6194B",  # red
        "#3CB44B",  # green
        "#4363D8",  # blue
        "#F58231",  # orange
        "#911EB4",  # purple
        "#42D4F4",  # cyan
        "#F032E6",  # magenta
        "#BFEF45",  # lime
        "#FABED4",  # pink
        "#9A6324",  # brown
    ]

    fig = go.Figure()

    # DFT trace at bottom (y=0)
    fig.add_trace(
        go.Scatter(
            x=freq_grid,
            y=dft_norm,
            name="DFT",
            line=dict(color=_DFT_COLOR, width=2),
            hovertemplate="%{x:.1f} cm\u207b\u00b9<br>%{y:.3f}<extra>DFT</extra>",
        )
    )

    # ML traces stacked above
    for i, (name, ml_norm) in enumerate(ml_norms.items()):
        color = _ml_colors[i % len(_ml_colors)]
        fig.add_trace(
            go.Scatter(
                x=freq_grid,
                y=ml_norm + offset * (i + 1),
                name=f"ML ({name})",
                line=dict(color=color, width=2),
                hovertemplate=(f"%{{x:.1f}} cm\u207b\u00b9<br>%{{y:.3f}}<extra>{name}</extra>"),
            )
        )

    # Experimental trace at y=0 (dashed)
    if experimental_norm is not None:
        fig.add_trace(
            go.Scatter(
                x=freq_grid,
                y=experimental_norm,
                name="Experimental (NIST)",
                line=dict(color=_EXP_COLOR, width=1.5),
                hovertemplate=("%{x:.1f} cm\u207b\u00b9<br>%{y:.3f}<extra>Exp.</extra>"),
            )
        )

    fig.update_layout(
        title="Combined IR Spectrum Comparison",
        xaxis_title="Wavenumber (cm\u207b\u00b9)",
        yaxis_title="Absorbance (normalized, stacked)",
        xaxis=dict(autorange="reversed"),
        hovermode="x unified",
        template="simple_white",
        height=600,
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig


def experimental_on_grid(
    experimental: Any | None,
    freq_grid: np.ndarray,
) -> np.ndarray | None:
    """Interpolate an ExperimentalSpectrum onto a frequency grid.

    Parameters
    ----------
    experimental : ExperimentalSpectrum or None
        Experimental spectrum with ``.wavenumbers`` and ``.absorbance`` arrays.
        If None, returns None.
    freq_grid : np.ndarray
        Target wavenumber grid for interpolation.

    Returns
    -------
    np.ndarray or None
        Normalized (0-1) absorbance on ``freq_grid``, or None.
    """
    if experimental is None:
        return None

    from scipy.interpolate import interp1d

    interp = interp1d(
        experimental.wavenumbers,
        experimental.absorbance,
        bounds_error=False,
        fill_value=0.0,
    )
    on_grid = interp(freq_grid)
    peak = float(np.max(on_grid))
    if peak > 0:
        on_grid = on_grid / peak
    return on_grid
