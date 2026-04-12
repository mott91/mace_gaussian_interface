"""HTML Report Generator -- Plotly-powered, mode-aware, executive-summary-first.

Phase 23 overhaul: replaces old matplotlib-base64 approach with interactive Plotly
figures, shared CSS from _shared_css, executive summary ranking, and structured
data export.  Mode flag ('harmonic' | 'anharmonic') controls overtones section.
"""

from __future__ import annotations

import base64
import html
import math
from pathlib import Path
from typing import Any, Literal

import numpy as np

from ._shared_css import build_css
from .analyze_spectra import SpectrumAnalyzer
from .executive_summary import (
    build_verdict,
    compute_experimental_agreement,
    rank_methods,
)
from .plotly_builders import (
    build_combined_spectrum_figure,
    build_regression_figure,
    build_spectrum_figure,
    experimental_on_grid,
)
from .report_data import export_report_data


class HTMLReportGenerator:
    """Generates HTML reports for spectral analysis with Plotly interactive figures."""

    def __init__(
        self,
        molecule_name: str,
        output_dir: Path,
        mode: Literal["harmonic", "anharmonic"] = "anharmonic",
        plotly_js: Any = "cdn",
        bandwidth_fwhm: float = 10.0,
        degenerate_groups: list | None = None,
    ):
        """Initialize report generator.

        Parameters
        ----------
        molecule_name : str
            Name of molecule.
        output_dir : Path
            Output directory containing plots and data.
        mode : str
            Report mode -- ``'harmonic'`` or ``'anharmonic'``.
        plotly_js : str or False
            How to include plotly.js: ``'cdn'``, ``'inline'``, or ``False``.
        bandwidth_fwhm : float
            Full width at half maximum for Lorentzian broadening (cm-1).
        degenerate_groups : list, optional
            Degenerate group info dicts.
        """
        self.molecule_name = molecule_name
        self.output_dir = Path(output_dir)
        self.mode = mode
        self.plotly_js: Any = plotly_js
        self.bandwidth_fwhm = bandwidth_fwhm
        self.degenerate_groups = degenerate_groups or []
        self._plotlyjs_emitted = False
        self._analyzer = SpectrumAnalyzer(bandwidth_fwhm=bandwidth_fwhm)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate_report(self, analysis_results: dict) -> None:
        """Generate complete HTML report.

        Parameters
        ----------
        analysis_results : dict
            Results from comparison workflow.
        """
        self._enrich_with_broadened_spectra(analysis_results)
        self._compute_and_attach_experimental_agreement(analysis_results)

        comparisons = analysis_results["comparisons"]
        ranked = rank_methods(comparisons)
        has_exp = analysis_results.get("experimental") is not None
        verdict = build_verdict(ranked, has_experimental=has_exp)
        analysis_results["executive_summary"] = {
            "verdict": verdict,
            "ranked": ranked,
        }

        sections = [
            self._create_head(),
            self._create_header(),
            self._create_navigation(comparisons),
            self._create_executive_summary(ranked, verdict),
            self._create_combined_plots(analysis_results),
        ]
        for i, comp in enumerate(comparisons, 1):
            sections.append(self._create_comparison_section(comp, i, analysis_results))
        sections.append(self._create_experimental_info_section(analysis_results))
        sections.append(self._create_summary_table(comparisons))
        if self.mode == "anharmonic":
            sections.append(self._create_overtones_section(analysis_results))
        sections.append(self._create_footer())

        html_out = "\n".join(sections) + "\n</body></html>\n"

        # Defensive emit-once check (T-23-05)
        cdn_count = html_out.count("cdn.plot.ly/plotly")
        if cdn_count > 1:
            raise RuntimeError(f"Plotly CDN referenced {cdn_count} times; expected at most 1.")

        out_path = self.output_dir / "report.html"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(html_out, encoding="utf-8")

        export_report_data(analysis_results, self.output_dir / "report_data.json")

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _esc(self, s: Any) -> str:
        """HTML-escape a value (T-23-01 mitigation)."""
        return html.escape(str(s), quote=True)

    @staticmethod
    def encode_image(image_path: Path) -> str:
        """Encode image to base64 for embedding."""
        with Path(image_path).open("rb") as f:
            encoded = base64.b64encode(f.read()).decode()
        return f"data:image/png;base64,{encoded}"

    def _fig_to_div(self, fig: Any, div_id: str) -> str:
        """Convert a Plotly figure to an HTML div string.

        Uses the emit-once pattern: first call includes plotly.js (per
        ``self.plotly_js``), subsequent calls use ``include_plotlyjs=False``.
        """
        import plotly.io as pio

        include: Any = self.plotly_js if not self._plotlyjs_emitted else False
        self._plotlyjs_emitted = True
        return pio.to_html(
            fig,
            include_plotlyjs=include,
            full_html=False,
            div_id=div_id,
            config={"displaylogo": False, "responsive": True},
        )

    # ------------------------------------------------------------------
    # Broadening / experimental enrichment
    # ------------------------------------------------------------------

    def _enrich_with_broadened_spectra(self, analysis_results: dict) -> None:
        """Broaden ML and DFT spectra onto the shared frequency grid."""
        freq_grid = self._analyzer.freq_grid
        analysis_results["freq_grid"] = freq_grid
        for comp in analysis_results["comparisons"]:
            ml_br = self._analyzer.broaden_spectrum(comp["ml_spectrum"])
            dft_br = self._analyzer.broaden_spectrum(comp["dft_spectrum"])
            comp["_ml_broadened"] = self._normalize(ml_br)
            comp["_dft_broadened"] = self._normalize(dft_br)

    @staticmethod
    def _normalize(arr: np.ndarray) -> np.ndarray:
        """Normalize an array to [0, 1]."""
        peak = float(np.max(arr)) if arr.size else 0.0
        return arr / peak if peak > 0 else arr

    def _compute_and_attach_experimental_agreement(self, analysis_results: dict) -> None:
        """Compute per-method experimental agreement and stash on comparisons."""
        exp = analysis_results.get("experimental")
        freq_grid = analysis_results["freq_grid"]
        exp_on = experimental_on_grid(exp, freq_grid) if exp is not None else None
        analysis_results["_exp_on_grid"] = exp_on
        for comp in analysis_results["comparisons"]:
            agree = compute_experimental_agreement(comp["_ml_broadened"], exp_on)
            comp["experimental_agreement"] = agree

    # ------------------------------------------------------------------
    # Section builders
    # ------------------------------------------------------------------

    def _create_head(self) -> str:
        """Create the ``<!DOCTYPE>`` preamble, ``<head>`` with shared CSS, and open ``<body>``.

        The ``<style>`` block is sourced from ``_shared_css.build_css()`` so that
        harmonic and anharmonic reports share byte-identical CSS.
        """
        mol = self._esc(self.molecule_name)
        return (
            '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
            '<meta charset="utf-8">\n'
            f"<title>{mol} -- {self.mode} report</title>\n"
            f"{build_css()}\n"
            "</head>\n<body>"
        )

    def _create_header(self) -> str:
        """Create the page header with molecule name and analysis mode.

        Molecule name is HTML-escaped to mitigate T-23-01.
        """
        mol = self._esc(self.molecule_name).upper()
        mode = self.mode.title()
        return f"<header>\n<h1>{mol} -- {mode} IR Analysis</h1>\n</header>"

    def _create_navigation(self, comparisons: list[dict]) -> str:
        """Create sticky navigation bar with anchor links.

        Links: executive summary, combined plots, per-method sections,
        summary table, and (if anharmonic) overtones section.
        """
        links = ['<a href="#executive-summary">Summary</a>']
        links.append('<a href="#combined">Combined</a>')
        for i, comp in enumerate(comparisons, 1):
            name = self._esc(comp["name"])
            links.append(f'<a href="#comparison-{i}">{name}</a>')
        links.append('<a href="#summary-table">Table</a>')
        if self.mode == "anharmonic":
            links.append('<a href="#overtones">Overtones</a>')
        return f"<nav>{''.join(links)}</nav>"

    def _create_executive_summary(self, ranked: list[dict], verdict: str) -> str:
        """Build the executive summary section (D-01, D-02).

        Renders a verdict line followed by one card per ranked method.
        The first card (best composite score) gets the ``best-method`` CSS class
        so it is visually highlighted.

        Parameters
        ----------
        ranked : list[dict]
            Methods sorted by composite score from :func:`rank_methods`.
        verdict : str
            Human-readable one-liner from :func:`build_verdict`.
        """
        cards: list[str] = []
        for idx, entry in enumerate(ranked):
            cls = "method-card best-method" if idx == 0 else "method-card"
            name = self._esc(entry["name"])
            exp_val = entry.get("experimental_agreement")
            exp_str = f"{exp_val:.2f}" if exp_val is not None else "\u2014"
            cards.append(
                f'<div class="{cls}">'
                f"<h3>{name}</h3>"
                f'<div class="metric">'
                f'<span class="metric-label">R\u00b2 (freq)</span>'
                f'<span class="metric-value">{entry["r2_freq"]:.3f}</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">R\u00b2 (intensity)</span>'
                f'<span class="metric-value">{entry["r2_intensity"]:.3f}</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">RMSE</span>'
                f'<span class="metric-value">{entry["rmse_freq"]:.1f} cm\u207b\u00b9</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">Speedup</span>'
                f'<span class="metric-value">{entry["speedup"]:.1f}\u00d7</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">Exp. agreement</span>'
                f'<span class="metric-value">{exp_str}</span></div>'
                "</div>"
            )
        return (
            '<section class="executive-summary" id="executive-summary">'
            f'<div class="verdict">{self._esc(verdict)}</div>'
            f'<div class="method-cards">{"".join(cards)}</div>'
            "</section>"
        )

    def _create_combined_plots(self, analysis_results: dict) -> str:
        """Build combined spectrum section with all ML methods overlaid.

        Uses ``build_combined_spectrum_figure`` from plotly_builders to produce
        a single stacked figure with DFT at bottom, ML methods above, and
        (optionally) the experimental NIST trace.
        """
        comparisons = analysis_results["comparisons"]
        freq_grid = analysis_results["freq_grid"]
        exp_on = analysis_results.get("_exp_on_grid")

        # Use first comparison's DFT broadened as the DFT reference
        dft_norm = comparisons[0]["_dft_broadened"]
        ml_norms: dict[str, np.ndarray] = {}
        for comp in comparisons:
            ml_norms[comp["name"]] = comp["_ml_broadened"]

        fig = build_combined_spectrum_figure(
            freq_grid, dft_norm, ml_norms, experimental_norm=exp_on
        )
        div = self._fig_to_div(fig, "combined-spectrum")
        return (
            '<section class="comparison-section" id="combined">'
            "<h2>Combined Spectrum Comparison</h2>"
            f'<div class="plot-container">{div}</div>'
            "</section>"
        )

    def _create_comparison_section(self, comp: dict, index: int, analysis_results: dict) -> str:
        """Build a per-method comparison section (D-11).

        Each section contains:
        - Interactive Plotly spectrum plot (ML vs DFT, optional experimental overlay)
        - Interactive Plotly regression scatter plot
        - Static PNG mode-overlap heatmap (if the file exists in output_dir/plots/)
        - Inline timing block (ML/DFT Gaussian elapsed, speedup)
        - Degenerate mode notes (if ``comp["deg_result"]`` has groups)
        - Per-method metrics mini-table

        Parameters
        ----------
        comp : dict
            Single comparison dict from ``analysis_results["comparisons"]``.
        index : int
            1-based index for HTML id attributes.
        analysis_results : dict
            Full analysis results (for freq_grid and experimental data).
        """
        freq_grid = analysis_results["freq_grid"]
        exp_on = analysis_results.get("_exp_on_grid")
        ml_name = comp["name"]
        ml_name_esc = self._esc(ml_name)

        # Plotly spectrum figure
        spec_fig = build_spectrum_figure(
            freq_grid,
            comp["_dft_broadened"],
            comp["_ml_broadened"],
            ml_name,
            experimental_norm=exp_on,
        )
        spec_div = self._fig_to_div(spec_fig, f"spectrum-{index}")

        # Plotly regression figure
        dft_freqs = np.asarray(comp["dft_spectrum"].frequencies, dtype=float)
        ml_freqs = np.asarray(comp["ml_spectrum"].frequencies, dtype=float)
        reg_fig = build_regression_figure(dft_freqs, ml_freqs, ml_name)
        reg_div = self._fig_to_div(reg_fig, f"regression-{index}")

        # Heatmap PNG (stays as static image)
        heatmap_html = ""
        plots_dir = self.output_dir / "plots"
        heatmap_files = (
            list(plots_dir.glob(f"mode_overlap_{ml_name}_*.png")) if plots_dir.exists() else []
        )
        if heatmap_files:
            heatmap_rel = f"plots/{heatmap_files[0].name}"
            heatmap_html = (
                f'<div class="plot-container">'
                f'<img src="{self._esc(heatmap_rel)}" alt="Mode overlap heatmap">'
                f"</div>"
            )

        # Timing block (D-09, D-11)
        ml_gauss_s = comp.get("ml_gaussian_timing", {}).get("total_elapsed_s", 0.0)
        dft_gauss_s = comp.get("dft_gaussian_timing", {}).get("total_elapsed_s", 0.0)
        speedup = comp.get("speedup", 0.0)
        timing_html = (
            '<div class="timing-block">'
            f"ML pipeline: <strong>{ml_gauss_s:.2f}s</strong> \u00b7 "
            f"DFT pipeline: <strong>{dft_gauss_s:.2f}s</strong> \u00b7 "
            f"Speedup: <strong>{speedup:.1f}\u00d7</strong>"
            "</div>"
        )

        # Degenerate notes
        deg_html = ""
        deg_result = comp.get("deg_result")
        if deg_result is not None:
            groups = getattr(deg_result, "groups", None) or []
            for g in groups:
                label = self._esc(g.get("label", ""))
                mult = g.get("multiplicity", 0)
                overlap = g.get("subspace_overlap", 0.0)
                deg_html += (
                    f'<div class="degenerate-note">'
                    f"Degenerate group {label}: {mult}-fold, "
                    f"subspace overlap {overlap:.3f}"
                    f"</div>"
                )

        # Per-method metrics mini table
        m = comp["metrics"]
        exp_agree = comp.get("experimental_agreement")
        exp_str = self._format_exp_agreement(exp_agree)
        metrics_table = (
            '<table class="data-table">'
            "<thead><tr>"
            "<th>R\u00b2 freq</th><th>R\u00b2 int</th>"
            "<th>RMSE</th><th>MAE</th>"
            "<th>Max err</th><th>Matched</th>"
            "<th>Speedup</th><th>Exp. agree</th>"
            "</tr></thead>"
            f"<tbody><tr>"
            f"<td>{m.r2_freq:.4f}</td>"
            f"<td>{m.r2_intensity:.4f}</td>"
            f"<td>{m.rmse_freq:.2f}</td>"
            f"<td>{m.mae_freq:.2f}</td>"
            f"<td>{m.max_error_freq:.2f}</td>"
            f"<td>{m.num_matched}/{m.num_matched + m.num_dft_only}</td>"
            f"<td>{speedup:.1f}\u00d7</td>"
            f"<td>{exp_str}</td>"
            f"</tr></tbody></table>"
        )

        return (
            f'<section class="comparison-section" id="comparison-{index}">'
            f"<h2>{ml_name_esc} vs DFT</h2>"
            f'<div class="plot-container">{spec_div}</div>'
            f'<div class="plot-container">{reg_div}</div>'
            f"{heatmap_html}"
            f"{timing_html}"
            f"{deg_html}"
            f"{metrics_table}"
            f"</section>"
        )

    def _create_summary_table(self, comparisons: list[dict]) -> str:
        """Build the overall summary comparison table.

        Columns: method name, R2 (freq), R2 (intensity), RMSE, speedup, and
        experimental agreement.  Method names are HTML-escaped.
        """
        if not comparisons:
            return (
                '<section id="summary-table">'
                "<h2>Summary Comparison</h2>"
                '<div class="warning-box">No comparisons available.</div>'
                "</section>"
            )

        rows: list[str] = []
        for comp in comparisons:
            m = comp["metrics"]
            name = self._esc(comp["name"])
            speedup = comp.get("speedup", 0.0)
            exp_agree = comp.get("experimental_agreement")
            exp_str = self._format_exp_agreement(exp_agree)
            rows.append(
                f"<tr>"
                f"<td>{name}</td>"
                f"<td>{m.r2_freq:.4f}</td>"
                f"<td>{m.r2_intensity:.4f}</td>"
                f"<td>{m.rmse_freq:.2f}</td>"
                f"<td>{speedup:.1f}\u00d7</td>"
                f"<td>{exp_str}</td>"
                f"</tr>"
            )
        return (
            '<section id="summary-table">'
            "<h2>Summary Comparison</h2>"
            '<table class="data-table">'
            "<thead><tr>"
            "<th>Method</th><th>R\u00b2 freq</th><th>R\u00b2 int</th>"
            "<th>RMSE</th><th>Speedup</th><th>Exp. agreement</th>"
            "</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>"
            "</section>"
        )

    @staticmethod
    def _format_exp_agreement(exp_agree: float | None) -> str:
        """Format experimental agreement for display.

        Returns the value formatted to 2 decimal places, or an em-dash
        when the value is None or NaN.
        """
        if exp_agree is None:
            return "\u2014"
        if isinstance(exp_agree, float) and math.isnan(exp_agree):
            return "\u2014"
        return f"{exp_agree:.2f}"

    def _create_experimental_info_section(self, analysis_results: dict) -> str:
        """Build the experimental data source info section.

        Shows NIST source, molecule name, CAS number, and data range
        when experimental data is available.  All user-facing strings
        are HTML-escaped (T-23-01).
        """
        experimental = analysis_results.get("experimental")
        if experimental is None:
            return ""

        source = self._esc(getattr(experimental, "source", ""))
        mol_name = self._esc(getattr(experimental, "molecule_name", ""))
        cas = self._esc(getattr(experimental, "cas_number", ""))

        # Wavenumber range
        wn = getattr(experimental, "wavenumbers", None)
        if wn is not None and len(wn) > 0:
            wn_range = f"{float(wn[0]):.0f} -- {float(wn[-1]):.0f} cm\u207b\u00b9"
        else:
            wn_range = "unknown"

        return (
            '<section class="comparison-section" id="experimental">'
            "<h2>Experimental Reference Data</h2>"
            '<div class="stats-box">'
            f"<p><strong>Source:</strong> {source}</p>"
            f"<p><strong>Molecule:</strong> {mol_name}</p>"
            f"<p><strong>CAS Number:</strong> {cas}</p>"
            f"<p><strong>Data range:</strong> {wn_range}</p>"
            "<p>Experimental IR spectrum overlaid as black dashed line "
            "on all spectrum plots above.</p>"
            "</div>"
            "</section>"
        )

    def _create_timing_hardware_section(self, comp: dict) -> str:
        """Build detailed timing and hardware context table for a comparison.

        This preserves the Phase 20 timing/hardware info that was in the
        original report generator.  Shows a two-row table (ML / DFT) with
        hardware device, pipeline time, and Gaussian wall time columns.

        Parameters
        ----------
        comp : dict
            Single comparison dict containing ``ml_hardware``, ``dft_hardware``,
            ``ml_gaussian_timing``, ``dft_gaussian_timing`` sub-dicts.

        Returns
        -------
        str
            HTML fragment, or empty string if no timing data is available.
        """
        ml_hw = comp.get("ml_hardware", {})
        dft_hw = comp.get("dft_hardware", {})
        ml_timing = comp.get("ml_gaussian_timing", {})
        dft_timing = comp.get("dft_gaussian_timing", {})

        # If no hardware or timing data at all, return empty
        has_data = any(
            [
                ml_hw.get("cpu"),
                ml_hw.get("gpu"),
                dft_hw.get("cpu"),
                dft_hw.get("node"),
                ml_timing.get("total_elapsed_s"),
                dft_timing.get("total_elapsed_s"),
            ]
        )
        if not has_data:
            return ""

        # ML row
        ml_gpu = ml_hw.get("gpu", "")
        ml_cpu = ml_hw.get("cpu", "")
        ml_device = self._esc(ml_gpu if ml_gpu else ml_cpu if ml_cpu else "\u2014")
        ml_gauss_s = ml_timing.get("total_elapsed_s", 0)
        ml_gauss = f"{ml_gauss_s:.1f} s" if ml_gauss_s else "\u2014"
        ml_pipeline = comp.get("ml_runtime", 0.0)

        # DFT row
        dft_cpu = dft_hw.get("cpu", "")
        dft_node = dft_hw.get("node", "")
        dft_cpus = dft_hw.get("cpus", "")
        dft_device = self._esc(dft_cpu if dft_cpu else "\u2014")
        if dft_node:
            dft_device += f" ({self._esc(dft_node)})"
        if dft_cpus:
            dft_device += f" x{self._esc(str(dft_cpus))}"
        dft_gauss_s = dft_timing.get("total_elapsed_s", 0)
        dft_gauss = f"{dft_gauss_s:.1f} s" if dft_gauss_s else "\u2014"
        dft_pipeline = comp.get("dft_runtime", 0.0)

        ml_name = self._esc(comp.get("name", "ML"))

        return (
            "<h3>Timing &amp; Hardware</h3>"
            '<table class="data-table">'
            "<thead><tr>"
            "<th>Method</th>"
            "<th>Hardware</th>"
            "<th>Pipeline Time</th>"
            "<th>Gaussian Wall Time</th>"
            "</tr></thead>"
            "<tbody>"
            f"<tr>"
            f"<td><strong>ML</strong> ({ml_name})</td>"
            f"<td>{ml_device}</td>"
            f"<td>{ml_pipeline:.1f} s</td>"
            f"<td>{ml_gauss}</td>"
            f"</tr>"
            f"<tr>"
            f"<td><strong>DFT</strong></td>"
            f"<td>{dft_device}</td>"
            f"<td>{dft_pipeline:.1f} s</td>"
            f"<td>{dft_gauss}</td>"
            f"</tr>"
            "</tbody></table>"
        )

    # ------------------------------------------------------------------
    # Overtones section (D-04 mode-aware)
    # ------------------------------------------------------------------

    # TODO(phase-follow-up): populate overtone/combination data from
    # analysis_workflow once upstream emits it as structured records in
    # comparison dicts.  As of Plan 23-04, overtone/combination band data
    # is embedded in SpectrumData labels but NOT surfaced as separate
    # comparison["overtones"] records -- see Step 0 grep survey.
    def _create_overtones_section(self, analysis_results: dict) -> str:
        """Render overtones section -- real data (Option A) or explicit stub (Option B).

        Option A (real data): If ``analysis_results["overtones"]`` or any
        ``comparison["overtones"]`` contains records, renders a
        ``<table class="overtone-table">`` with columns: mode_id,
        frequency, intensity, type.

        Option B (explicit stub): Renders an ``overtones-placeholder`` div
        explaining that the data is not yet surfaced by the pipeline.

        Both options produce grep-able markers for automated verification.
        The section heading always contains the word "overtone" so the
        harmonic-skip test can detect its presence/absence.
        """
        overtones_data = self._collect_overtones(analysis_results)
        if overtones_data:
            # Option A: real data
            rows = "".join(
                f"<tr><td>{self._esc(row.get('mode_id', ''))}</td>"
                f"<td>{row.get('freq_cm', 0):.1f}</td>"
                f"<td>{row.get('intensity', 0):.3f}</td>"
                f"<td>{self._esc(row.get('type', 'overtone'))}</td></tr>"
                for row in overtones_data
            )
            body = (
                '<table class="overtone-table data-table">'
                "<thead><tr><th>Mode</th><th>Freq (cm\u207b\u00b9)</th>"
                "<th>Intensity</th><th>Type</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>"
            )
        else:
            # Option B: explicit stub, grep-able
            body = (
                '<div class="overtones-placeholder">'
                "Placeholder -- anharmonic overtone/combination data "
                "not yet surfaced by the analysis pipeline. Populated "
                "in a follow-up phase."
                "</div>"
            )
        return (
            '<section class="comparison-section" id="overtones">'
            "<h2>Overtones and combination bands</h2>"
            f"{body}"
            "</section>"
        )

    def _collect_overtones(self, analysis_results: dict) -> list[dict]:
        """Collect overtone/combination band records from analysis_results.

        Checks ``analysis_results["overtones"]`` first (top-level), then
        iterates per-comparison ``comp["overtones"]``.  Returns an empty
        list if no data is found, which triggers the stub path in
        ``_create_overtones_section``.

        Returns
        -------
        list[dict]
            Each dict has keys: mode_id, freq_cm, intensity, type.
            Empty list if no overtone data is available.
        """
        # Check top-level first
        top = analysis_results.get("overtones")
        if top:
            return list(top)
        # Check per-comparison
        collected: list[dict] = []
        for comp in analysis_results.get("comparisons", []):
            ot = comp.get("overtones")
            if ot:
                collected.extend(ot)
        return collected

    def _create_footer(self) -> str:
        """Create the report footer with generation info."""
        return (
            f'<footer class="footer">'
            f"Generated by mace-gaussian {self._esc(self.mode)} analysis "
            f"\u00b7 Phase 23"
            f"</footer>"
        )
