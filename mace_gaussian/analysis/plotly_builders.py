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
_EXP_COLOR = "rgba(128,128,128,0.45)"

# Eigenvector-overlap threshold below which a matched mode is rendered with a
# hollow marker — the matching algorithm placed it but the modes are physically
# different enough that the comparison should be read with caution.
LOW_OVERLAP_THRESHOLD = 0.7
# R\u00b2 is not shown below this many points (report review item 3, 2026-09-18)
MIN_N_FOR_R2 = 5

_TYPE_STYLE = {
    "fundamental": {"color": _ML_COLOR, "symbol": "circle", "name": "Fundamental"},
    "overtone":    {"color": "#0173B2", "symbol": "diamond", "name": "Overtone"},
    "combination": {"color": "#CC78BC", "symbol": "square", "name": "Combination"},
}


def _classify_mode_type(mode_id: str | None) -> str:
    if mode_id is None:
        return "fundamental"
    if mode_id.startswith("O"):
        return "overtone"
    if mode_id.startswith("C"):
        return "combination"
    return "fundamental"


def _add_regression_traces(
    fig: Any,
    x: np.ndarray,
    y: np.ndarray,
    mode_ids: list[str] | None,
    mode_overlaps: list[float | None] | None,
    *,
    unit_label: str,
    low_threshold: float = LOW_OVERLAP_THRESHOLD,
) -> int:
    """Add regression scatter traces grouped by (type, confidence).

    Low-overlap points (overlap < ``low_threshold``) render with hollow
    symbols and reduced opacity. Returns the count of low-overlap points so
    the caller can surface it in the report.
    """
    import plotly.graph_objects as go

    n = len(x)
    has_overlaps = mode_overlaps is not None and len(mode_overlaps) == n

    # Group key: (type, is_low_overlap) → (xs, ys, overlaps)
    groups: dict[tuple[str, bool], tuple[list, list, list]] = {}
    for i in range(n):
        mid = mode_ids[i] if mode_ids is not None and i < len(mode_ids) else None
        t = _classify_mode_type(mid)
        ovl = mode_overlaps[i] if has_overlaps else None
        is_low = has_overlaps and ovl is not None and ovl < low_threshold
        groups.setdefault((t, is_low), ([], [], []))
        groups[(t, is_low)][0].append(x[i])
        groups[(t, is_low)][1].append(y[i])
        groups[(t, is_low)][2].append(ovl)

    n_low = sum(len(v[0]) for k, v in groups.items() if k[1])

    for t in ("fundamental", "overtone", "combination"):
        for is_low in (False, True):
            key = (t, is_low)
            if key not in groups:
                continue
            xs, ys, ovls = groups[key]
            style = _TYPE_STYLE[t]
            symbol = style["symbol"] + ("-open" if is_low else "")
            label = f"{style['name']}{', low overlap' if is_low else ''} (n={len(xs)})"
            # Only fundamentals carry meaningful overlap values — show on hover.
            show_overlap = t == "fundamental" and any(o is not None for o in ovls)
            customdata = [[o if o is not None else float("nan")] for o in ovls] if show_overlap else None
            hover = f"DFT: %{{x:.1f}} {unit_label}<br>ML: %{{y:.1f}} {unit_label}"
            if show_overlap:
                hover += "<br>overlap: %{customdata[0]:.2f}"
            extra = f"{style['name']}{', low overlap' if is_low else ''}"

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="markers",
                    name=label,
                    legendgroup=t,
                    customdata=customdata,
                    marker=dict(
                        color=style["color"],
                        symbol=symbol,
                        size=8,
                        opacity=0.55 if is_low else 1.0,
                        line=dict(width=1.5 if is_low else 1, color="#333"),
                    ),
                    hovertemplate=hover + f"<extra>{extra}</extra>",
                )
            )
    return n_low


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
    mode_overlaps: list[float | None] | None = None,
) -> object:
    """Build a regression scatter plot (ML vs DFT frequencies).

    Includes a y=x perfect-agreement reference line.  When *mode_ids* is
    supplied, points are color-coded by type (fundamental / overtone /
    combination).  When *mode_overlaps* is supplied, fundamentals whose
    eigenvector overlap is below ``LOW_OVERLAP_THRESHOLD`` render with
    hollow markers — flagging that the matching algorithm may have paired
    physically different modes.

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
    mode_overlaps : list[float | None] or None
        Per-point eigenvector overlap, aligned with ``dft_freqs``.  ``None``
        for derived overtones / combinations.

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

    fig = go.Figure()

    if mode_ids is not None or mode_overlaps is not None:
        _add_regression_traces(
            fig, dft_freqs, ml_freqs, mode_ids, mode_overlaps,
            unit_label="cm\u207b\u00b9",
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
    r2_str = f"{r2_all:.3f}" if (not np.isnan(r2_all) and n_shown >= MIN_N_FOR_R2) else "n/a"
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
    mode_overlaps: list[float | None] | None = None,
) -> object:
    """Build a regression scatter plot (ML vs DFT intensities).

    When *mode_ids* is supplied, points are color-coded by type
    (fundamental / overtone / combination).  When *mode_overlaps* is
    supplied, fundamentals whose eigenvector overlap is below
    ``LOW_OVERLAP_THRESHOLD`` render with hollow markers.

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
    mode_overlaps : list[float | None] or None
        Per-point eigenvector overlap, aligned with ``dft_intensities``.

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

    fig = go.Figure()

    if mode_ids is not None or mode_overlaps is not None:
        _add_regression_traces(
            fig, dft_intensities, ml_intensities, mode_ids, mode_overlaps,
            unit_label="km/mol",
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
    r2_str = f"{r2_all:.3f}" if (not np.isnan(r2_all) and n_shown >= MIN_N_FOR_R2) else "n/a"
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


def build_residual_figure(
    dft_freqs: np.ndarray,
    ml_freqs: np.ndarray,
    ml_name: str,
    mode_ids: list[str] | None = None,
) -> object:
    """Build a residual plot: (ML - DFT) error vs DFT frequency.

    Reveals systematic bias (e.g. consistently blue-shifted regions).
    """
    import plotly.graph_objects as go

    errors = ml_freqs - dft_freqs

    _type_style = {
        "fundamental": {"color": _ML_COLOR, "symbol": "circle", "name": "Fundamental"},
        "overtone": {"color": "#0173B2", "symbol": "diamond", "name": "Overtone"},
        "combination": {"color": "#CC78BC", "symbol": "square", "name": "Combination"},
    }

    fig = go.Figure()

    if mode_ids is not None and len(mode_ids) == len(dft_freqs):
        groups: dict[str, tuple[list, list]] = {}
        for i, mid in enumerate(mode_ids):
            if mid.startswith("O"):
                t = "overtone"
            elif mid.startswith("C"):
                t = "combination"
            else:
                t = "fundamental"
            groups.setdefault(t, ([], []))
            groups[t][0].append(float(dft_freqs[i]))
            groups[t][1].append(float(errors[i]))

        for t in ("fundamental", "overtone", "combination"):
            if t not in groups:
                continue
            dx, ex = groups[t]
            style = _type_style[t]
            fig.add_trace(
                go.Scatter(
                    x=dx, y=ex, mode="markers",
                    name=f"{style['name']} (n={len(dx)})",
                    marker=dict(color=style["color"], symbol=style["symbol"], size=6,
                                line=dict(width=0.5, color="#333")),
                    hovertemplate=(
                        "DFT: %{x:.1f} cm\u207b\u00b9<br>"
                        "\u0394: %{y:.1f} cm\u207b\u00b9"
                        f"<extra>{style['name']}</extra>"
                    ),
                )
            )
    else:
        fig.add_trace(
            go.Scatter(
                x=dft_freqs, y=errors, mode="markers",
                name=f"Residuals (n={len(dft_freqs)})",
                marker=dict(color=_ML_COLOR, size=6, line=dict(width=0.5, color="#333")),
                hovertemplate="DFT: %{x:.1f}<br>\u0394: %{y:.1f}<extra></extra>",
            )
        )

    # Zero line
    fig.add_hline(y=0, line=dict(color="#888", dash="dash", width=1))

    mae = float(np.mean(np.abs(errors)))
    bias = float(np.mean(errors))
    fig.update_layout(
        title=f"Residuals: {ml_name} (bias={bias:+.1f}, MAE={mae:.1f} cm\u207b\u00b9)",
        xaxis_title="DFT frequency (cm\u207b\u00b9)",
        yaxis_title="ML \u2212 DFT (cm\u207b\u00b9)",
        template="simple_white",
        autosize=True,
        height=400,
        margin=dict(l=60, r=30, t=60, b=50),
        legend=dict(x=0.05, y=0.98, xanchor="left", yanchor="top",
                    bgcolor="rgba(255,255,255,0.8)"),
    )
    return fig


def build_error_histogram_figure(
    dft_freqs: np.ndarray,
    ml_freqs: np.ndarray,
    ml_name: str,
    mode_ids: list[str] | None = None,
) -> object:
    """Build a histogram of frequency errors (ML - DFT), stacked by mode type."""
    import plotly.graph_objects as go

    errors = ml_freqs - dft_freqs

    _type_colors = {
        "fundamental": _ML_COLOR,
        "overtone": "#0173B2",
        "combination": "#CC78BC",
    }

    fig = go.Figure()

    if mode_ids is not None and len(mode_ids) == len(errors):
        groups: dict[str, list[float]] = {}
        for i, mid in enumerate(mode_ids):
            if mid.startswith("O"):
                t = "overtone"
            elif mid.startswith("C"):
                t = "combination"
            else:
                t = "fundamental"
            groups.setdefault(t, [])
            groups[t].append(float(errors[i]))

        for t in ("fundamental", "overtone", "combination"):
            if t not in groups:
                continue
            fig.add_trace(
                go.Histogram(
                    x=groups[t], name=f"{t.capitalize()} (n={len(groups[t])})",
                    marker_color=_type_colors[t], opacity=0.7,
                )
            )
        fig.update_layout(barmode="overlay")
    else:
        fig.add_trace(
            go.Histogram(
                x=errors, name=f"All modes (n={len(errors)})",
                marker_color=_ML_COLOR, opacity=0.8,
            )
        )

    fig.add_vline(x=0, line=dict(color="#888", dash="dash", width=1))

    std = float(np.std(errors))
    # Symmetric x-axis so 0 is centred
    x_lim = float(np.max(np.abs(errors))) * 1.15
    fig.update_layout(
        title=f"Error Distribution: {ml_name} (\u03c3={std:.1f} cm\u207b\u00b9)",
        xaxis_title="ML \u2212 DFT (cm\u207b\u00b9)",
        xaxis=dict(range=[-x_lim, x_lim]),
        yaxis_title="Count",
        template="simple_white",
        autosize=True,
        height=400,
        margin=dict(l=60, r=30, t=60, b=50),
        legend=dict(x=0.95, y=0.98, xanchor="right", yanchor="top",
                    bgcolor="rgba(255,255,255,0.8)"),
    )
    return fig


def build_anharmonicity_ratio_figure(
    dft_harm_freqs: np.ndarray,
    dft_anharm_freqs: np.ndarray,
    ml_harm_freqs: np.ndarray,
    ml_anharm_freqs: np.ndarray,
    ml_name: str,
    labels: list[str] | None = None,
) -> object:
    """Anharmonic shift (VPT2 minus harmonic, in cm^-1) of every fundamental, ML vs DFT.

    Shows whether the ML surface has the same anharmonicity as the reference. The shift is
    given in cm^-1 and not in percent of the harmonic frequency: percent inflates the
    low-frequency modes. ``labels`` (one per mode) are shown when hovering over a point.
    """
    import plotly.graph_objects as go

    valid = np.isfinite(dft_harm_freqs) & np.isfinite(ml_harm_freqs)
    valid &= np.isfinite(dft_anharm_freqs) & np.isfinite(ml_anharm_freqs)

    if np.sum(valid) < 2:
        fig = go.Figure()
        fig.add_annotation(text="Insufficient data for the anharmonic shift",
                           xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    dft_shift = dft_anharm_freqs[valid] - dft_harm_freqs[valid]
    ml_shift = ml_anharm_freqs[valid] - ml_harm_freqs[valid]
    text = (
        [lab for lab, ok in zip(labels, valid) if ok]
        if labels is not None and len(labels) == len(valid)
        else [""] * int(np.sum(valid))
    )

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=dft_shift, y=ml_shift, mode="markers",
            name=f"Modes (n={int(np.sum(valid))})",
            text=text,
            marker=dict(color=_ML_COLOR, size=8, line=dict(width=1, color="#333")),
            hovertemplate=(
                "%{text}<br>DFT shift: %{x:+.1f} cm\u207b\u00b9"
                "<br>ML shift: %{y:+.1f} cm\u207b\u00b9<extra></extra>"
            ),
        )
    )

    all_r = np.concatenate([dft_shift, ml_shift])
    lo, hi = float(np.min(all_r)), float(np.max(all_r))
    margin = (hi - lo) * 0.05
    lo -= margin
    hi += margin
    mae = float(np.mean(np.abs(ml_shift - dft_shift)))

    fig.add_trace(
        go.Scatter(
            x=[lo, hi], y=[lo, hi], mode="lines",
            name=f"y=x | MAE {mae:.1f} cm\u207b\u00b9",
            line=dict(color="#888", dash="dash", width=1.5),
            hoverinfo="skip",
        )
    )

    fig.update_layout(
        title=f"Anharmonic shift: {ml_name} vs DFT",
        xaxis_title="DFT shift, VPT2 \u2212 harmonic (cm\u207b\u00b9)",
        yaxis_title="ML shift, VPT2 \u2212 harmonic (cm\u207b\u00b9)",
        xaxis=dict(range=[lo, hi], constrain="domain"),
        yaxis=dict(range=[lo, hi], scaleanchor="x", scaleratio=1, constrain="domain"),
        template="simple_white",
        width=400,
        height=400,
        margin=dict(l=60, r=30, t=60, b=50),
        legend=dict(x=0.95, y=0.05, xanchor="right", yanchor="bottom",
                    bgcolor="rgba(255,255,255,0.8)"),
    )
    return fig


def build_per_region_table(
    dft_freqs: np.ndarray,
    ml_freqs: np.ndarray,
    mode_ids: list[str] | None = None,
) -> dict[str, dict[str, float]]:
    """Accuracy broken down by wavenumber range and by band type.

    Returns an ordered dict of row_name -> {mae, rmse, bias, n}: neutral wavenumber
    ranges of the DFT frequency, then (when ``mode_ids`` is given) one row per band
    type. Not a Plotly figure. Report review item 2 (2026-09-18): the old rows were
    named after organic functional groups, which put water's bend overtone at
    3195 cm-1 into a "C-H stretch" row and all fundamentals into "Other".
    """
    ranges = [
        ("400\u20131000 cm\u207b\u00b9", 400, 1000),
        ("1000\u20132000 cm\u207b\u00b9", 1000, 2000),
        ("2000\u20133000 cm\u207b\u00b9", 2000, 3000),
        ("3000\u20134000 cm\u207b\u00b9", 3000, 4000),
        ("> 4000 cm\u207b\u00b9", 4000, float("inf")),
    ]
    dft_freqs = np.asarray(dft_freqs, dtype=float)
    ml_freqs = np.asarray(ml_freqs, dtype=float)
    errors = ml_freqs - dft_freqs

    def _stats(mask):
        if not np.any(mask):
            return None
        e = errors[mask]
        return {
            "mae": float(np.mean(np.abs(e))),
            "rmse": float(np.sqrt(np.mean(e**2))),
            "bias": float(np.mean(e)),
            "n": int(np.sum(mask)),
        }

    result: dict[str, dict[str, float]] = {}
    below = _stats(dft_freqs < 400)
    if below:
        result["< 400 cm\u207b\u00b9"] = below
    for name, lo, hi in ranges:
        row = _stats((dft_freqs >= lo) & (dft_freqs < hi))
        if row:
            result[name] = row
    if mode_ids is not None and len(mode_ids) == len(dft_freqs):
        prefixes = np.array([str(m)[:1] for m in mode_ids])
        for label, prefix in (
            ("Fundamentals", "F"),
            ("Overtones", "O"),
            ("Combination bands", "C"),
        ):
            row = _stats(prefixes == prefix)
            if row:
                result[label] = row
    return result


def build_pareto_figure(
    method_names: list[str],
    mae_values: list[float],
    speedup_values: list[float],
) -> object:
    """Build a cost-accuracy scatter: MAE vs speedup across methods.

    Duplicate (x, y) pairs get their labels alternated top/bottom so they
    don't pile on top of each other (e.g. methods sharing the same MACE
    frequency backbone have identical MAE_freq).
    """
    import plotly.graph_objects as go

    display_names = [n.replace("mace_", "") for n in method_names]

    # Pairs that share a MACE frequency backbone have identical MAE_freq;
    # alternate top/bottom label placement so names don't stack on each other.
    seen: dict[float, int] = {}
    positions: list[str] = []
    _alt = ("top right", "bottom right")
    for y in mae_values:
        key = round(float(y), 1)
        idx = seen.get(key, 0)
        positions.append(_alt[idx % len(_alt)])
        seen[key] = idx + 1

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=speedup_values,
            y=mae_values,
            mode="markers+text",
            text=display_names,
            textposition=positions,
            textfont=dict(size=10),
            cliponaxis=False,
            marker=dict(color=_ML_COLOR, size=10, line=dict(width=1, color="#333")),
            hovertemplate=(
                "%{text}<br>"
                "Speedup: %{x:.2f}\u00d7<br>"
                "MAE: %{y:.1f} cm\u207b\u00b9"
                "<extra></extra>"
            ),
        )
    )

    # Break-even reference: ML as fast as DFT
    if speedup_values and min(speedup_values) < 1.0 < max(speedup_values) * 1.05:
        fig.add_vline(
            x=1.0, line=dict(color="#888", dash="dash", width=1),
            annotation_text="ML = DFT cost", annotation_position="top",
            annotation_font=dict(size=10, color="#888"),
        )

    fig.update_layout(
        title="Cost vs Accuracy: Speedup vs Frequency MAE",
        xaxis_title="Speedup (\u00d7 vs DFT)",
        yaxis_title="Frequency MAE (cm\u207b\u00b9)",
        template="simple_white",
        autosize=True,
        height=420,
        margin=dict(l=70, r=40, t=60, b=55),
        showlegend=False,
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
