# ruff: noqa: RUF001  (typographic characters in user-facing strings are intentional)
"""The per-molecule HTML report: organized by energy model, with the per-mode error chart.

Introduced 2026-09-18 as "report v2" and made the default the same day; the previous
per-run layout (``html_report_generator.HTMLReportGenerator``) is still available
through ``MACE_GAUSSIAN_LEGACY_REPORT=1`` and its helpers are reused here.

What differs from the legacy layout:

1. One section per *energy model* (five), not per run (fifteen). The three dipole
   models of one energy model share every frequency, so the frequency story is told
   once and a "dipole models" panel compares only what they change: intensities.
2. A per-mode signed-error chart (ML minus DFT for VPT2 and harmonic, plus minus
   experiment where a band origin exists), all models on one axis: the visual of
   Hypothesis 1.
3. Combined spectrum: every trace normalized to its own maximum and stacked with a
   fixed offset, so nothing clips; range buttons for fingerprint and stretch regions.
4. Experimental band origins (Shimanouchi, via the master table) drawn as labelled
   ticks on every spectrum; the raw NIST trace is smoothed for small molecules.
5. Fixed color per energy model in every figure (``palette.py``); band type is
   carried by marker shape/color inside scatter plots, dipole model by symbol.
7. One overlap heatmap per energy model; merged summary table; grouped navigation.

Reuses the v1 generator's helpers (broadening, run-integrity table, anharmonicity
section, overtone tables) so the numbers are identical between the two reports.
"""

from __future__ import annotations

import base64
import html
import logging
from pathlib import Path
from typing import Any

import numpy as np

from .executive_summary import build_verdict, rank_methods
from .html_report_generator import HTMLReportGenerator
from .master_table import (
    _DIPOLE_PREFERENCE,
    build_master_table,
    display_name,
    render_master_table_html,
    split_method,
)
from .palette import (
    DFT_COLOR,
    EXP_COLOR,
    dipole_color,
    dipole_symbol,
    hex_to_rgba,
    model_color,
)
from .plotly_builders import (
    LOW_OVERLAP_THRESHOLD,
    build_error_histogram_figure,
    build_per_region_table,
    build_regression_figure,
    build_residual_figure,
    build_spectrum_figure,
)
from .report_data import export_report_data

logger = logging.getLogger(__name__)

# Smooth the raw NIST trace with this Gaussian FWHM (cm^-1) when the molecule has
# at most this many normal modes (rotational structure dominates small molecules).
EXP_SMOOTH_FWHM_CM = 30.0
EXP_SMOOTH_MAX_MODES = 9

_V2_CSS = """
<style>
.v2-banner{background:#fff7e6;border:1px solid #f0c36d;border-radius:6px;padding:8px 14px;
  margin:12px 0;font-size:0.9em}
.swatch{display:inline-block;width:14px;height:14px;border-radius:3px;margin-right:8px;
  vertical-align:-1px}
.model-section{border-left:6px solid var(--model-color,#888);padding-left:14px}
.model-section h2 small{font-weight:normal;color:#6b7280;font-size:0.6em;margin-left:8px}
.method-card{border-top:4px solid var(--model-color,#ddd)}
nav .nav-group{display:inline-block;margin:0 6px;padding:0 6px;border-left:1px solid #ccc}
.two-col{display:flex;gap:1rem;flex-wrap:wrap;overflow:hidden}
.two-col>div{flex:1 1 0;min-width:0;overflow:hidden}
.caption{color:#6b7280;font-size:0.85em;margin:4px 0 12px}
.method-cards{grid-template-columns:repeat(auto-fit,minmax(200px,1fr))}
.heatmap img{max-width:720px;margin:0 auto}
.side-by-side{display:flex;gap:1rem;flex-wrap:wrap;align-items:stretch}
.side-by-side>*{flex:1 1 320px;min-width:0;margin:1.5rem 0}
.plot-container{padding:0.75rem}
.two-col .plotly-graph-div{margin:0 auto}
.table-scroll{max-height:75vh;overflow:auto}
.table-scroll thead{position:sticky;top:0;background:#f6f8fa;z-index:1}
</style>
"""


def _esc(s: Any) -> str:
    return html.escape(str(s), quote=True)


class ReportV2Generator:
    """Build ``report.html`` (plus report_data.json and master table exports)."""

    def __init__(
        self,
        molecule_name: str,
        output_dir: Path,
        mode: str = "anharmonic",
        bandwidth_fwhm: float = 10.0,
        plotly_js: Any = "cdn",
    ) -> None:
        self.molecule_name = molecule_name
        self.output_dir = Path(output_dir)
        self.mode = mode
        self.plotly_js = plotly_js
        self._emitted = False
        # v1 helpers are reused for numbers and a few HTML blocks; its own plotly.js
        # emission is disabled so this report emits the script exactly once.
        self.v1 = HTMLReportGenerator(
            molecule_name, output_dir, mode=mode, plotly_js=plotly_js, bandwidth_fwhm=bandwidth_fwhm
        )
        self.v1._plotlyjs_emitted = True
        self._analyzer = self.v1._analyzer

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def generate(self, analysis_results: dict, filename: str = "report.html") -> Path:
        out_path = self.output_dir / filename
        out_path.parent.mkdir(parents=True, exist_ok=True)
        comparisons = analysis_results.get("comparisons") or []
        if not comparisons:
            out_path.write_text(
                "<!DOCTYPE html><html><body><h1>No comparisons available</h1></body></html>"
            )
            return out_path

        if "freq_grid" not in analysis_results:
            self.v1._enrich_with_broadened_spectra(analysis_results)
            self.v1._compute_and_attach_experimental_agreement(analysis_results)
        if not analysis_results.get("master_table"):
            analysis_results["master_table"] = build_master_table(analysis_results)
        table = analysis_results["master_table"]

        groups = self._group_by_energy_model(comparisons)
        for g in groups:
            for comp in g["runs"]:
                self._attach_matches(comp)

        ranked = rank_methods(comparisons)
        verdict = build_verdict(ranked, has_experimental=False)
        analysis_results["executive_summary"] = {"verdict": verdict, "ranked": ranked}
        exp_bands = self._experimental_bands(table)
        n_modes = len(table.get("rows") or [])
        exp_on = analysis_results.get("_exp_on_grid")
        exp_label = "Experimental (NIST)"
        if exp_on is not None and n_modes <= EXP_SMOOTH_MAX_MODES:
            exp_on = self._smooth(exp_on, float(self._analyzer.freq_step))
            exp_label = f"Experimental (NIST, smoothed {EXP_SMOOTH_FWHM_CM:.0f} cm⁻¹)"

        sections = [
            self.v1._create_head().replace("</head>", _V2_CSS + "</head>"),
            f"<header><h1>{_esc(self.molecule_name).upper()} — {self.mode} IR analysis "
            "</h1></header>",
            self._nav(groups),
            self._summary(groups, ranked, verdict),
            self._combined(analysis_results, groups, exp_on, exp_label, exp_bands),
            self._per_mode_errors(table, groups),
            render_master_table_html(table),
        ]
        for i, g in enumerate(groups, 1):
            sections.append(
                self._model_section(g, i, analysis_results, exp_on, exp_label, exp_bands)
            )
        if self.mode == "anharmonic":
            sections.append(
                self.v1._create_overtones_section({"comparisons": [g["rep"] for g in groups]})
            )
        sections.append(self._summary_table(groups))
        sections.append(self.v1._create_footer())

        html_out = "\n".join(sections) + "\n</body></html>\n"
        if html_out.count("cdn.plot.ly/plotly") > 1:
            raise RuntimeError("plotly.js referenced more than once")
        out_path.write_text(html_out, encoding="utf-8")
        export_report_data(analysis_results, self.output_dir / "report_data.json")
        return out_path

    # ------------------------------------------------------------------
    # Grouping and matching
    # ------------------------------------------------------------------

    @staticmethod
    def _group_by_energy_model(comparisons: list[dict]) -> list[dict]:
        groups: dict[str, dict] = {}
        for comp in comparisons:
            energy, dip = split_method(comp["name"])
            g = groups.setdefault(
                energy,
                {
                    "energy": energy,
                    "label": display_name(energy),
                    "color": model_color(energy),
                    "runs": [],
                },
            )
            comp["_dipole_model"] = dip
            g["runs"].append(comp)
        out = []
        for g in groups.values():
            g["runs"].sort(
                key=lambda c: _DIPOLE_PREFERENCE.index(c["_dipole_model"])
                if c["_dipole_model"] in _DIPOLE_PREFERENCE
                else 99
            )
            g["rep"] = g["runs"][0]
            out.append(g)
        out.sort(key=lambda g: g["rep"]["metrics"].mae_freq)
        return out

    def _attach_matches(self, comp: dict) -> None:
        """Matched DFT/ML arrays through the eigenvector mapping, cached on the comp."""
        if "_v2_match" in comp:
            return
        dft_f, ml_f, dft_i, ml_i, stats = self._analyzer.match_by_mode(
            comp["dft_spectrum"],
            comp["ml_spectrum"],
            mode_mapping=comp.get("mode_mapping"),
            mode_overlaps=comp.get("mode_overlaps"),
        )
        ids = stats.get("matched_mode_ids")
        overlaps = stats.get("matched_mode_overlaps")
        if len(dft_f) == 0:
            d, m = comp["dft_spectrum"], comp["ml_spectrum"]
            n = min(len(d.frequencies), len(m.frequencies))
            dft_f, ml_f = np.asarray(d.frequencies[:n], float), np.asarray(m.frequencies[:n], float)
            dft_i, ml_i = np.asarray(d.intensities[:n], float), np.asarray(m.intensities[:n], float)
            ids, overlaps = None, None
        cat_mae: dict[str, float | None] = {}
        if ids is not None and len(ids) == len(dft_f):
            for cat, prefix in (("fundamental", "F"), ("overtone", "O"), ("combination", "C")):
                mask = np.array([str(x).startswith(prefix) for x in ids])
                cat_mae[cat] = (
                    float(np.mean(np.abs(ml_f[mask] - dft_f[mask]))) if mask.any() else None
                )
        comp["_v2_match"] = {
            "dft_f": dft_f,
            "ml_f": ml_f,
            "dft_i": dft_i,
            "ml_i": ml_i,
            "ids": ids,
            "overlaps": overlaps,
            "cat_mae": cat_mae,
        }

    @staticmethod
    def _experimental_bands(table: dict, merge_within_cm: float = 8.0) -> list[tuple[str, float]]:
        """(label, band origin) per curated band; bands closer than ``merge_within_cm``
        share one tick (labels joined with a slash) so coincident bands do not overprint."""
        seen: dict[str, float] = {}
        for r in table.get("rows") or []:
            e = r.get("experimental")
            if e and e.get("freq_cm") is not None:
                seen.setdefault(e.get("label") or "", float(e["freq_cm"]))
        merged: list[tuple[str, float]] = []
        for label, freq in sorted(seen.items(), key=lambda kv: kv[1]):
            if merged and abs(freq - merged[-1][1]) <= merge_within_cm:
                prev_label, prev_freq = merged[-1]
                merged[-1] = (f"{prev_label}/{label}", (prev_freq + freq) / 2)
            else:
                merged.append((label, freq))
        return merged

    @staticmethod
    def _smooth(y: np.ndarray, step: float) -> np.ndarray:
        sigma = EXP_SMOOTH_FWHM_CM / 2.3548 / step
        half = int(4 * sigma) + 1
        x = np.arange(-half, half + 1)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        kernel /= kernel.sum()
        out = np.convolve(y, kernel, mode="same")
        peak = float(out.max()) if out.size else 0.0
        return out / peak if peak > 0 else out

    # ------------------------------------------------------------------
    # Plotly helpers
    # ------------------------------------------------------------------

    def _div(self, fig: Any, div_id: str, responsive: bool = True) -> str:
        import plotly.io as pio

        include: Any = self.plotly_js if not self._emitted else False
        self._emitted = True
        return pio.to_html(
            fig,
            include_plotlyjs=include,
            full_html=False,
            div_id=div_id,
            config={"displaylogo": False, "responsive": responsive},
        )

    @staticmethod
    def _add_band_ticks(fig: Any, exp_bands: list[tuple[str, float]], y_top: float) -> None:
        """Dotted band-origin lines with labels on two alternating heights, so
        neighbouring bands (methanol's CH3 modes, for instance) stay legible."""
        for i, (label, freq) in enumerate(exp_bands):
            fig.add_shape(
                type="line",
                x0=freq,
                x1=freq,
                y0=0,
                y1=y_top,
                line=dict(color=EXP_COLOR, dash="dot", width=1),
                layer="below",
            )
            fig.add_annotation(
                x=freq,
                y=y_top,
                text=label,
                showarrow=False,
                yanchor="bottom",
                yshift=0 if i % 2 == 0 else 13,
                font=dict(size=10, color=EXP_COLOR),
            )

    @staticmethod
    def _recolor(fig: Any, old: str, new: str) -> None:
        for tr in fig.data:
            line = getattr(tr, "line", None)
            if line is not None and line.color == old:
                tr.line.color = new
            marker = getattr(tr, "marker", None)
            if marker is not None and marker.color == old:
                tr.marker.color = new

    @staticmethod
    def _add_range_buttons(fig: Any, full_lo: float, full_hi: float) -> None:
        """Wavenumber range buttons; default view is the fundamental region."""
        ranges = [
            ("Fundamentals 400–4000", 400, 4000),
            ("Fingerprint 400–2000", 400, 2000),
            ("Stretch 2500–4000", 2500, 4000),
            ("Full incl. overtones", full_lo, full_hi),
        ]
        fig.update_layout(
            xaxis=dict(range=[4000, 400], autorange=False),
            title_x=0.0,
            margin=dict(t=90),
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    x=1.0,
                    xanchor="right",
                    y=1.0,
                    yanchor="bottom",
                    pad=dict(b=4),
                    buttons=[
                        dict(label=lbl, method="relayout", args=[{"xaxis.range": [hi, lo]}])
                        for lbl, lo, hi in ranges
                    ],
                )
            ],
        )

    # ------------------------------------------------------------------
    # Sections
    # ------------------------------------------------------------------

    def _nav(self, groups: list[dict]) -> str:
        links = [
            '<a href="#executive-summary">Summary</a>',
            '<a href="#combined">Combined</a>',
            '<a href="#per-mode-errors">Per-mode errors</a>',
            '<a href="#master-table">Master table</a>',
        ]
        models = "".join(
            f'<a href="#model-{i}"><span class="swatch" style="background:{g["color"]}"></span>'
            f"{_esc(g['label'])}</a>"
            for i, g in enumerate(groups, 1)
        )
        links.append(f'<span class="nav-group">{models}</span>')
        if self.mode == "anharmonic":
            links.append('<a href="#overtones">Overtones</a>')
        links.append('<a href="#summary-table">Table</a>')
        return f"<nav>{''.join(links)}</nav>"

    def _summary(self, groups: list[dict], ranked: list[dict], verdict: str) -> str:
        by_name = {r["name"]: r for r in ranked}
        cards = []
        for i, g in enumerate(groups):
            r = by_name.get(g["rep"]["name"], {})
            m = g["rep"]["metrics"]
            ov = g["rep"]["_v2_match"]["overlaps"] or []
            valid = [o for o in ov if o is not None]
            n_low = sum(1 for o in valid if o < LOW_OVERLAP_THRESHOLD)
            cls = "method-card best-method" if i == 0 else "method-card"
            runs = ", ".join(c["_dipole_model"] or "?" for c in g["runs"])
            cards.append(
                f'<div class="{cls}" style="--model-color:{g["color"]}">'
                f"<h3>{_esc(g['label'])}</h3>"
                f'<div class="metric"><span class="metric-label">MAE (freq)</span>'
                f'<span class="metric-value">{m.mae_freq:.1f} cm⁻¹</span></div>'
                f'<div class="metric"><span class="metric-label">RMSE (freq)</span>'
                f'<span class="metric-value">{m.rmse_freq:.1f} cm⁻¹</span></div>'
                f'<div class="metric"><span class="metric-label">Max error</span>'
                f'<span class="metric-value">{m.max_error_freq:.1f} cm⁻¹</span></div>'
                f'<div class="metric"><span class="metric-label">Low-overlap modes</span>'
                f'<span class="metric-value">{n_low}/{len(valid) or m.num_matched}</span></div>'
                f'<div class="metric"><span class="metric-label">Speedup</span>'
                f'<span class="metric-value">{r.get("speedup", 0.0):.1f}×</span></div>'
                f'<div class="metric"><span class="metric-label">Dipole runs</span>'
                f'<span class="metric-value">{_esc(runs)}</span></div>'
                "</div>"
            )
        return (
            '<section class="executive-summary" id="executive-summary">'
            f'<div class="verdict">{_esc(verdict)}</div>'
            '<p class="caption">Cards are energy models (frequencies); the dipole model only '
            "changes intensities and is compared inside each model section.</p>"
            f'<div class="method-cards">{"".join(cards)}</div></section>'
        )

    def _combined(
        self,
        analysis_results: dict,
        groups: list[dict],
        exp_on: np.ndarray | None,
        exp_label: str,
        exp_bands: list[tuple[str, float]],
    ) -> str:
        import plotly.graph_objects as go

        freq_grid = analysis_results["freq_grid"]
        # Shared scale: every ML trace in units of the DFT peak, so relative band
        # strengths stay comparable. The lane height follows the tallest trace, so a
        # model with stronger bands gets more room instead of clipping.
        ratios = [max(g["rep"].get("_ml_dft_peak_ratio", 1.0), 1.0) for g in groups]
        offset = 1.15 * max([1.0, *ratios])
        fig = go.Figure()
        level = 0.0
        if exp_on is not None:
            fig.add_trace(
                go.Scatter(
                    x=freq_grid,
                    y=exp_on,
                    name=exp_label,
                    line=dict(color=EXP_COLOR, width=1.5),
                    fill="tozeroy",
                    fillcolor=hex_to_rgba(EXP_COLOR, 0.15),
                    hovertemplate="%{x:.1f} cm⁻¹<br>%{y:.3f}<extra>Exp.</extra>",
                )
            )
            level += offset
        dft = groups[0]["rep"]["_dft_broadened"]
        fig.add_trace(
            go.Scatter(
                x=freq_grid,
                y=dft + level,
                name="DFT (B3LYP)",
                line=dict(color=DFT_COLOR, width=2),
                hovertemplate="%{x:.1f} cm⁻¹<br>%{y:.3f}<extra>DFT</extra>",
            )
        )
        for g in groups:
            level += offset
            ratio = g["rep"].get("_ml_dft_peak_ratio", 1.0)
            fig.add_trace(
                go.Scatter(
                    x=freq_grid,
                    y=g["rep"]["_ml_broadened"] * ratio + level,
                    name=(
                        f"{g['label']} (peak {ratio:.1f}\u00d7 DFT)"
                        if abs(ratio - 1.0) > 0.05
                        else g["label"]
                    ),
                    line=dict(color=g["color"], width=2),
                    hovertemplate=f"%{{x:.1f}} cm⁻¹<br>%{{y:.3f}}<extra>{g['label']}</extra>",
                )
            )
        y_top = level + 1.0
        self._add_band_ticks(fig, exp_bands, y_top)
        lo, hi = float(freq_grid[0]), float(freq_grid[-1])
        fig.update_layout(
            xaxis_title="Wavenumber (cm⁻¹)",
            yaxis=dict(
                title="Absorbance (DFT peak = 1, offset)",
                showticklabels=False,
                range=[0, y_top + 0.3],
            ),
            hovermode="x unified",
            template="simple_white",
            height=120 + 95 * (len(groups) + 2),
            margin=dict(l=60, r=20, t=50, b=50),
            legend=dict(orientation="h", y=-0.12, yanchor="top", x=0, xanchor="left"),
        )
        self._add_range_buttons(fig, lo, hi)
        div = self._div(fig, "combined-spectrum-v2")
        return (
            '<section class="comparison-section" id="combined"><h2>Combined spectra</h2>'
            f'<div class="plot-container">{div}</div>'
            '<p class="caption">All calculated traces share one scale (DFT peak = 1), so relative '
            "band strengths are comparable; the experimental trace is normalized to its own "
            "maximum (arbitrary units). Dotted ticks: experimental band origins "
            "(Shimanouchi 1972 via "
            "NIST WebBook). Lorentzian broadening "
            f"{self._analyzer.bandwidth_fwhm:.0f} cm⁻¹ FWHM.</p></section>"
        )

    def _per_mode_errors(self, table: dict, groups: list[dict]) -> str:
        from plotly.subplots import make_subplots

        rows = table.get("rows") or []
        anharm = table.get("mode") == "anharmonic"
        if not rows:
            return ""
        has_exp = any(r.get("experimental") for r in rows)
        panels = [("VPT2: ML − DFT", "vpt2", "dft_vpt2")] if anharm else []
        panels.append(("Harmonic: ML − DFT", "harmonic", "dft_harmonic"))
        if has_exp:
            panels.append(
                ("vs experiment: model − band origin", "vpt2" if anharm else "harmonic", "exp")
            )

        # With many modes the tick labels become a smear: show the mode index only
        # and keep the full label (assignment, DFT frequency) in the hover.
        compact = len(rows) > 24
        full_labels = []
        for r in rows:
            e = r.get("experimental") or {}
            ref = r["dft_vpt2"] if r["dft_vpt2"] is not None else r["dft_harmonic"]
            desc = e.get("description") or e.get("label") or ""
            full_labels.append(
                f"{r['dft_mode']}: {desc} ({ref:.0f})" if ref is not None else str(r["dft_mode"])
            )
        labels = [str(r["dft_mode"]) for r in rows] if compact else full_labels
        marker_size = 6 if compact else 9

        fig = make_subplots(
            rows=len(panels),
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            subplot_titles=[p[0] for p in panels],
        )
        for row_i, (_, key, ref_key) in enumerate(panels, 1):
            if ref_key == "exp":
                dft_key = "dft_vpt2" if anharm else "dft_harmonic"
                ys = [
                    (r[dft_key] - r["experimental"]["freq_cm"])
                    if r.get("experimental") and r.get(dft_key) is not None
                    else None
                    for r in rows
                ]
                fig.add_trace(
                    dict(
                        type="scatter",
                        x=labels,
                        y=ys,
                        name="DFT",
                        mode="markers+lines",
                        legendgroup="dft",
                        showlegend=True,
                        line=dict(color=DFT_COLOR, width=1),
                        text=full_labels,
                        marker=dict(color=DFT_COLOR, size=marker_size, symbol="star"),
                        hovertemplate="%{text}<br>DFT − exp: %{y:+.1f} cm⁻¹<extra></extra>",
                    ),
                    row=row_i,
                    col=1,
                )
            for g in groups:
                ys, symbols = [], []
                for r in rows:
                    c = r["methods"].get(g["rep"]["name"]) or {}
                    v = c.get(key)
                    ref = (
                        (r["experimental"]["freq_cm"] if r.get("experimental") else None)
                        if ref_key == "exp"
                        else r.get(ref_key)
                    )
                    ys.append(v - ref if v is not None and ref is not None else None)
                    ov = c.get("overlap")
                    symbols.append(
                        "circle-open" if ov is not None and ov < LOW_OVERLAP_THRESHOLD else "circle"
                    )
                fig.add_trace(
                    dict(
                        type="scatter",
                        x=labels,
                        y=ys,
                        name=g["label"],
                        mode="markers+lines",
                        legendgroup=g["energy"],
                        showlegend=(row_i == 1),
                        text=full_labels,
                        line=dict(color=g["color"], width=1),
                        marker=dict(
                            color=g["color"],
                            size=marker_size,
                            symbol=symbols,
                            line=dict(width=1.5, color=g["color"]),
                        ),
                        hovertemplate=(
                            f"%{{text}}<br>{g['label']}: %{{y:+.1f}} cm⁻¹<extra></extra>"
                        ),
                    ),
                    row=row_i,
                    col=1,
                )
            fig.add_hline(y=0, line=dict(color="#888", dash="dash", width=1), row=row_i, col=1)
            fig.update_yaxes(title_text="Δν (cm⁻¹)", row=row_i, col=1)
        fig.update_xaxes(
            tickangle=-90 if compact else -30,
            tickfont=dict(size=9 if compact else 11),
            title_text="DFT mode index" if compact else None,
            row=len(panels),
            col=1,
        )
        fig.update_layout(
            template="simple_white",
            height=240 * len(panels) + 150,
            margin=dict(l=70, r=20, t=70, b=120),
            legend=dict(orientation="h", y=1.06, yanchor="bottom", x=0, xanchor="left"),
            hovermode="x unified",
        )
        # The div id must differ from the section id, otherwise Plotly draws the
        # figure into the <section> element (first match of getElementById).
        div = self._div(fig, "per-mode-errors-fig")
        return (
            '<section class="comparison-section" id="per-mode-errors">'
            "<h2>Per-mode signed errors</h2>"
            '<p class="caption">One point per DFT normal mode and energy model, paired by '
            "eigenvector overlap (hollow marker: overlap below "
            f"{LOW_OVERLAP_THRESHOLD:.1f}). Mode labels use the Shimanouchi assignment where "
            "a band origin exists.</p>"
            f'<div class="plot-container">{div}</div></section>'
        )

    def _model_section(
        self,
        g: dict,
        index: int,
        analysis_results: dict,
        exp_on: np.ndarray | None,
        exp_label: str,
        exp_bands: list[tuple[str, float]],
    ) -> str:
        rep = g["rep"]
        mt = rep["_v2_match"]
        freq_grid = analysis_results["freq_grid"]
        color = g["color"]
        label = g["label"]

        ratio = rep.get("_ml_dft_peak_ratio", 1.0)
        spec_fig = build_spectrum_figure(
            freq_grid,
            rep["_dft_broadened"],
            rep["_ml_broadened"] * ratio,
            label,
            experimental_norm=exp_on,
            offset=1.15 * max(1.0, ratio),
        )
        self._recolor(spec_fig, "#DE8F05", color)
        for tr in spec_fig.data:
            if tr.name.startswith("Experimental"):
                tr.name = exp_label
        self._add_band_ticks(spec_fig, exp_bands, 2.6)
        spec_fig.update_layout(
            title=f"IR spectrum: {label} vs DFT",
            height=440,
            margin=dict(l=60, r=20, t=70, b=50),
            legend=dict(orientation="h", y=-0.15, yanchor="top", x=0, xanchor="left"),
        )
        self._add_range_buttons(spec_fig, float(freq_grid[0]), float(freq_grid[-1]))
        spec_div = self._div(spec_fig, f"spectrum-v2-{index}")

        reg_div = res_div = hist_div = anharm_div = ""
        row_h = 440
        if len(mt["dft_f"]) > 0:
            reg_fig = build_regression_figure(
                mt["dft_f"], mt["ml_f"], label, mode_ids=mt["ids"], mode_overlaps=mt["overlaps"]
            )
            reg_fig.update_layout(
                height=row_h,
                title_x=0.0,
                margin=dict(l=60, r=20, t=60, b=50),
                legend=dict(x=0.98, y=0.02, xanchor="right", yanchor="bottom"),
            )
            reg_div = self._div(reg_fig, f"regression-v2-{index}")
            res_fig = build_residual_figure(mt["dft_f"], mt["ml_f"], label, mode_ids=mt["ids"])
            res_fig.update_layout(height=row_h, title_x=0.0, margin=dict(l=60, r=20, t=60, b=50))
            res_div = self._div(res_fig, f"residual-v2-{index}")
            hist_fig = build_error_histogram_figure(
                mt["dft_f"], mt["ml_f"], label, mode_ids=mt["ids"]
            )
            hist_fig.update_layout(height=row_h, title_x=0.0, margin=dict(l=60, r=20, t=60, b=50))
            hist_div = self._div(hist_fig, f"hist-v2-{index}")
        if self.mode == "anharmonic" and rep.get("_ml_results") and rep.get("_dft_results"):
            anharm_fig = self._anharmonicity_figure(rep, label)
            if anharm_fig is not None:
                # v1 fixes this figure at 400 x 400 px and non-responsive; here it
                # fills its column like its neighbours.
                anharm_fig.update_layout(
                    width=row_h, height=row_h, title_x=0.0, margin=dict(l=60, r=20, t=60, b=50)
                )
                anharm_div = self._div(anharm_fig, f"anharm-v2-{index}", responsive=False)

        region_html = ""
        if len(mt["dft_f"]) > 0:
            region = build_per_region_table(mt["dft_f"], mt["ml_f"], mode_ids=mt["ids"])
            trs = "".join(
                f"<tr><td>{_esc(k)}</td><td>{v['n']}</td><td>{v['mae']:.1f}</td>"
                f"<td>{v['rmse']:.1f}</td><td>{v['bias']:+.1f}</td></tr>"
                for k, v in region.items()
            )
            region_html = (
                '<div class="stats-box"><h4>Accuracy by wavenumber range and band type</h4>'
                '<table class="summary-table"><thead><tr><th>Rows: DFT range, then band type</th>'
                "<th>n</th><th>MAE</th><th>RMSE</th><th>Bias</th></tr></thead>"
                f"<tbody>{trs}</tbody></table></div>"
            )

        heat_html = ""
        plots_dir = self.output_dir / "plots"
        files = (
            sorted(plots_dir.glob(f"mode_overlap_{rep['name']}_*.png"))
            if plots_dir.exists()
            else []
        )
        if files:
            b64 = base64.b64encode(files[0].read_bytes()).decode("ascii")
            heat_html = (
                '<div class="plot-container heatmap"><img src="data:image/png;base64,'
                f'{b64}" alt="Mode overlap heatmap ({_esc(label)})"></div>'
            )

        m = rep["metrics"]

        def stat(label_txt: str, value: str) -> str:
            return (
                f'<div class="stat-item"><div class="stat-label">{label_txt}</div>'
                f'<div class="stat-value">{value}</div></div>'
            )

        stats = (
            '<div class="stats-box"><h4>Frequencies (identical for all dipole runs)</h4>'
            '<div class="stats-grid">'
            + stat("MAE", f"{m.mae_freq:.2f} cm⁻¹")
            + stat("RMSE", f"{m.rmse_freq:.2f} cm⁻¹")
            + stat("Max error", f"{m.max_error_freq:.2f} cm⁻¹")
            + stat("Matched bands", f"{m.num_matched}/{m.num_matched + m.num_dft_only}")
            + self.v1._low_overlap_stat({"_matched_overlaps": mt["overlaps"]})
            + stat("Speedup", f"{rep.get('speedup', 0.0):.1f}×")
            + "</div></div>"
        )

        runs = ", ".join(_esc(c["name"]) for c in g["runs"])
        return (
            f'<section class="comparison-section model-section" id="model-{index}" '
            f'style="--model-color:{color}">'
            f'<h2><span class="swatch" style="background:{color}"></span>{_esc(label)}'
            f"<small>{_esc(g['energy'])} · runs: {runs}</small></h2>"
            f'<div class="side-by-side">{stats}{self.v1._create_run_integrity(rep)}</div>'
            f'<div class="plot-container">{spec_div}</div>'
            f'<div class="two-col"><div class="plot-container">{reg_div}</div>'
            f'<div class="plot-container">{res_div}</div></div>'
            f'<div class="two-col"><div class="plot-container">{anharm_div}</div>'
            f'<div class="plot-container">{hist_div}</div></div>'
            f"{region_html}"
            f"{self._dipole_panel(g, index)}"
            f"{heat_html}</section>"
        )

    @staticmethod
    def _anharmonicity_figure(rep: dict, label: str) -> Any | None:
        """Anharmonicity-ratio figure paired through the eigenvector mapping (as v1)."""
        from .analyze_spectra import gaussian_mode_to_checkpoint_index
        from .plotly_builders import build_anharmonicity_ratio_figure

        ml_anharm = rep["_ml_results"].get("frequencies", {}).get("anharmonic", [])
        dft_anharm = rep["_dft_results"].get("frequencies", {}).get("anharmonic", [])
        if not ml_anharm or not dft_anharm:
            return None
        ml_to_ckpt = gaussian_mode_to_checkpoint_index(ml_anharm)
        dft_to_ckpt = gaussian_mode_to_checkpoint_index(dft_anharm)
        ml_by = {ml_to_ckpt.get(m["mode"], m["mode"]) - 1: m for m in ml_anharm}
        dft_by = {dft_to_ckpt.get(m["mode"], m["mode"]) - 1: m for m in dft_anharm}
        mapping = rep.get("mode_mapping") or {i: i for i in ml_by}
        pairs = [
            (ml_by[i], dft_by[j]) for i, j in sorted(mapping.items()) if i in ml_by and j in dft_by
        ]
        if len(pairs) < 3:
            return None
        return build_anharmonicity_ratio_figure(
            np.array([d["freq_harmonic"] for _, d in pairs]),
            np.array([d["freq_cm"] for _, d in pairs]),
            np.array([m["freq_harmonic"] for m, _ in pairs]),
            np.array([m["freq_cm"] for m, _ in pairs]),
            label,
        )

    def _dipole_panel(self, g: dict, index: int) -> str:
        """Intensities are the only thing the dipole model changes: compare them here."""
        import plotly.graph_objects as go

        fig = go.Figure()
        rows = []
        all_vals = []
        for comp in g["runs"]:
            mt = comp["_v2_match"]
            dip = comp["_dipole_model"] or "?"
            mask = (mt["dft_i"] >= 0.1) | (mt["ml_i"] >= 0.1)
            x, y = mt["dft_i"][mask], mt["ml_i"][mask]
            if len(x):
                all_vals += list(x) + list(y)
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="markers",
                        name=dip,
                        marker=dict(
                            color=dipole_color(dip),
                            symbol=dipole_symbol(dip),
                            size=9,
                            line=dict(width=1, color="#333"),
                        ),
                        hovertemplate=f"DFT %{{x:.2f}}<br>{dip} %{{y:.2f}} km/mol<extra></extra>",
                    )
                )
            m = comp["metrics"]
            n_int = m.num_peaks - m.num_intensity_filtered
            rows.append(
                f"<tr><td>{_esc(dip)}</td><td>{m.mae_intensity:.2f}</td>"
                f"<td>{m.rmse_intensity:.2f}</td>"
                f"<td>{self.v1._r2_text(m.r2_intensity, n_int, 3)}</td>"
                f"<td>{comp.get('ml_runtime', 0.0):.1f} s</td></tr>"
            )
        if all_vals:
            lo, hi = max(min(all_vals), 1e-3), max(all_vals) * 1.1
            fig.add_trace(
                go.Scatter(
                    x=[lo, hi],
                    y=[lo, hi],
                    mode="lines",
                    name="y = x",
                    line=dict(color="#888", dash="dash", width=1),
                    hoverinfo="skip",
                )
            )
            log = hi / lo > 100
            axis = dict(type="log", dtick=1, exponentformat="power") if log else dict()
            fig.update_layout(
                title=f"IR intensities by dipole model: {g['label']}",
                xaxis=dict(title="DFT intensity (km/mol)", constrain="domain", **axis),
                yaxis=dict(
                    title="ML intensity (km/mol)",
                    scaleanchor="x",
                    scaleratio=1,
                    constrain="domain",
                    **axis,
                ),
                template="simple_white",
                height=440,
                margin=dict(l=60, r=20, t=60, b=50),
                legend=dict(
                    x=0.98,
                    y=0.02,
                    xanchor="right",
                    yanchor="bottom",
                    bgcolor="rgba(255,255,255,0.8)",
                ),
            )
            div = self._div(fig, f"dipole-v2-{index}")
        else:
            div = ""
        return (
            '<div class="side-by-side">'
            '<div class="stats-box"><h4>Dipole models (intensities only)</h4>'
            '<table class="summary-table"><thead><tr><th>Dipole model</th>'
            "<th>MAE (km/mol)</th>"
            "<th>RMSE (km/mol)</th><th>R²</th><th>Pipeline</th>"
            "</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table></div>"
            f'<div class="plot-container">{div}</div></div>'
        )

    def _summary_table(self, groups: list[dict]) -> str:
        dips = []
        for g in groups:
            for c in g["runs"]:
                if c["_dipole_model"] not in dips:
                    dips.append(c["_dipole_model"])
        head = (
            "<th>Model</th><th>MAE fund.</th><th>MAE overt.</th><th>MAE comb.</th><th>RMSE</th>"
            "<th>Max err.</th><th>Low overlap</th>"
            + "".join(f"<th>MAE int. ({_esc(d)})</th>" for d in dips)
            + "<th>Speedup</th>"
        )
        trs = []
        for g in groups:
            rep = g["rep"]
            m = rep["metrics"]
            cat = rep["_v2_match"]["cat_mae"]
            ov = [o for o in (rep["_v2_match"]["overlaps"] or []) if o is not None]
            n_low = sum(1 for o in ov if o < LOW_OVERLAP_THRESHOLD)
            by_dip = {c["_dipole_model"]: c for c in g["runs"]}
            cells = [
                f'<td><span class="swatch" style="background:{g["color"]}"></span>'
                f"{_esc(g['label'])}</td>",
                *(
                    f"<td>{cat.get(k):.1f}</td>" if cat.get(k) is not None else "<td>—</td>"
                    for k in ("fundamental", "overtone", "combination")
                ),
                f"<td>{m.rmse_freq:.1f}</td>",
                f"<td>{m.max_error_freq:.1f}</td>",
                f"<td>{n_low}/{len(ov)}</td>",
                *(
                    f"<td>{by_dip[d]['metrics'].mae_intensity:.1f}</td>"
                    if d in by_dip
                    else "<td>—</td>"
                    for d in dips
                ),
                f"<td>{rep.get('speedup', 0.0):.1f}×</td>",
            ]
            trs.append(f"<tr>{''.join(cells)}</tr>")
        return (
            '<section id="summary-table"><h2>Summary by energy model</h2>'
            '<p class="caption">Frequency columns in cm⁻¹ against DFT; intensity MAE in km/mol '
            "per dipole model. Replaces the separate summary table and category ranking of v1.</p>"
            f'<table class="data-table"><thead><tr>{head}</tr></thead>'
            f"<tbody>{''.join(trs)}</tbody></table></section>"
        )
