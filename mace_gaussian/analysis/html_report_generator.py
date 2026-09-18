"""HTML Report Generator -- Plotly-powered, mode-aware, executive-summary-first.

Phase 23 overhaul: replaces old matplotlib-base64 approach with interactive Plotly
figures, shared CSS from _shared_css, executive summary ranking, and structured
data export.  Mode flag ('harmonic' | 'anharmonic') controls overtones section.
"""

from __future__ import annotations

import base64
import html
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
    build_anharmonicity_ratio_figure,
    build_combined_spectrum_figure,
    build_error_histogram_figure,
    build_intensity_regression_figure,
    build_pareto_figure,
    build_per_region_table,
    build_regression_figure,
    build_residual_figure,
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
        comparisons = analysis_results.get("comparisons") or []
        if not comparisons:
            # Review finding L11: run_full_analysis returns an empty comparison list
            # (with an "error" key) when no DFT baseline or no ML results are found;
            # write a short error page instead of failing on comparisons[0].
            self._write_error_report(analysis_results.get("error", "No comparisons available"))
            return

        self._enrich_with_broadened_spectra(analysis_results)
        self._compute_and_attach_experimental_agreement(analysis_results)
        # Sort comparisons: non-espaloma first (by MAE), then espaloma (by MAE)
        comparisons.sort(
            key=lambda c: ("espaloma" in c["name"].lower(), c["metrics"].mae_freq)
        )
        analysis_results["comparisons"] = comparisons

        ranked = rank_methods(comparisons)
        has_exp = analysis_results.get("experimental") is not None
        verdict = build_verdict(ranked, has_experimental=has_exp)
        analysis_results["executive_summary"] = {
            "verdict": verdict,
            "ranked": ranked,
        }

        # Mode count overview from first comparison's DFT spectrum
        mode_overview = ""
        if comparisons and "dft_spectrum" in comparisons[0]:
            from collections import Counter
            dft_labels = Counter(comparisons[0]["dft_spectrum"].labels)
            n_fund = dft_labels.get("fundamental", 0)
            n_ot = dft_labels.get("overtone", 0)
            n_cb = dft_labels.get("combination", 0)
            items = [
                f'<div class="stat-item"><div class="stat-label">Fundamentals</div>'
                f'<div class="stat-value">{n_fund}</div></div>'
            ]
            if self.mode == "anharmonic":
                items.append(
                    f'<div class="stat-item"><div class="stat-label">Overtones</div>'
                    f'<div class="stat-value">{n_ot}</div></div>'
                )
                items.append(
                    f'<div class="stat-item"><div class="stat-label">Combination Bands</div>'
                    f'<div class="stat-value">{n_cb}</div></div>'
                )
            items.append(
                f'<div class="stat-item"><div class="stat-label">Total Modes</div>'
                f'<div class="stat-value">{n_fund + n_ot + n_cb}</div></div>'
            )
            mode_overview = (
                '<div class="stats-grid" style="margin:16px 0">'
                + "".join(items)
                + "</div>"
            )

        sections = [
            self._create_head(),
            self._create_header(),
            self._create_navigation(comparisons),
            self._create_executive_summary(ranked, verdict),
            mode_overview,
            self._create_combined_plots(analysis_results),
        ]
        for i, comp in enumerate(comparisons, 1):
            sections.append(self._create_comparison_section(comp, i, analysis_results))
        sections.append(self._create_experimental_info_section(analysis_results))
        if self.mode == "anharmonic":
            sections.append(self._create_category_awards(comparisons))
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

    def _write_error_report(self, message: str) -> None:
        """Minimal report.html explaining why there is nothing to show (L11)."""
        html_out = "\n".join(
            [
                self._create_head(),
                self._create_header(),
                f'<section class="section"><h2>No comparisons</h2>'
                f"<p>{self._esc(message)}</p>"
                "<p>Check that <code>comparison_results/&lt;molecule&gt;/</code> contains a "
                "DFT directory with <code>calculator_type: dft</code> and at least one ML "
                "directory with <code>calculator_type: ml</code> in its results.json.</p>"
                "</section>",
                self._create_footer(),
                "</body></html>",
            ]
        )
        out_path = self.output_dir / "report.html"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(html_out, encoding="utf-8")

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _esc(self, s: Any) -> str:
        """HTML-escape a value (T-23-01 mitigation)."""
        return html.escape(str(s), quote=True)

    @staticmethod
    def _r2_text(r2: float, n: int | None, digits: int = 3) -> str:
        """R² for display, or 'n/a' when there are too few points for it to mean anything."""
        from .plotly_builders import MIN_N_FOR_R2

        if n is not None and n < MIN_N_FOR_R2:
            return f"n/a (n={n})"
        return f"{r2:.{digits}f}"

    def _create_run_integrity(self, comp: dict) -> str:
        """Was this run healthy? One row per check, ML and DFT side by side.

        Reads calculation_parameters written by workflow/dft_baseline since 2026-09-18
        (re-optimization, dipole fallbacks, VPT2 diagnostics). Older results.json
        files show "n/a". Report review item 4.
        """
        ml = (comp.get("_ml_results") or {}).get("calculation_parameters") or {}
        dft = (comp.get("_dft_results") or {}).get("calculation_parameters") or {}
        reopt = ml.get("reoptimization") or {}
        fb = ml.get("dipole_fallbacks") or {}
        vml = ml.get("vpt2_diagnostics") or {}
        vdft = dft.get("vpt2_diagnostics") or {}

        def cell(v, fmt=None, good=None):
            if v is None:
                return "<td>n/a</td>"
            txt = fmt(v) if fmt else str(v)
            cls = "" if good is None else (' class="metric-good"' if good(v) else ' class="metric-bad"')
            return f"<td{cls}>{self._esc(txt)}</td>"

        rows = [
            ("Geometry re-optimized on its own surface", cell(reopt.get("reoptimized_with")), "<td>B3LYP opt</td>"),
            ("Optimizer converged", cell(reopt.get("converged"), good=bool), "<td>n/a</td>"),
            ("Final max force (eV/Å)", cell(reopt.get("max_force_eV_A"), lambda v: f"{v:.1e}", lambda v: v < 1e-3), "<td>n/a</td>"),
            ("Max atom shift from start (Å)", cell(reopt.get("max_atom_shift_A"), lambda v: f"{v:.3f}", lambda v: v < 0.1), "<td>n/a</td>"),
            ("Dipole fallbacks (zeroed calls)", cell(fb.get("count"), good=lambda v: v == 0), "<td>0</td>"),
            ("Unreliable cubic force constants", cell(vml.get("unreliable_cubic"), good=lambda v: v == 0), cell(vdft.get("unreliable_cubic"), good=lambda v: v == 0)),
            ("Fermi resonance deperturbed", cell(vml.get("fermi_resonances")), cell(vdft.get("fermi_resonances"))),
            ("Darling-Dennison deperturbed", cell(vml.get("darling_dennison_resonances")), cell(vdft.get("darling_dennison_resonances"))),
        ]
        body = "".join(f"<tr><td>{self._esc(k)}</td>{a}{b}</tr>" for k, a, b in rows)
        return (
            '<div class="stats-box"><h4>Run integrity</h4>'
            '<table class="summary-table"><thead><tr><th>Check</th><th>ML</th><th>DFT</th></tr></thead>'
            f"<tbody>{body}</tbody></table></div>"
        )

    @staticmethod
    def encode_image(image_path: Path) -> str:
        """Encode image to base64 for embedding."""
        with Path(image_path).open("rb") as f:
            encoded = base64.b64encode(f.read()).decode()
        return f"data:image/png;base64,{encoded}"

    def _fig_to_div(self, fig: Any, div_id: str, responsive: bool = True) -> str:
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
            config={"displaylogo": False, "responsive": responsive},
        )

    # ------------------------------------------------------------------
    # Broadening / experimental enrichment
    # ------------------------------------------------------------------

    def _enrich_with_broadened_spectra(self, analysis_results: dict) -> None:
        """Broaden ML and DFT spectra onto the shared frequency grid.

        In anharmonic mode, extends the grid to cover the highest detected
        frequency (overtones/combinations can exceed 4000 cm-1) plus margin.
        """
        if self.mode == "anharmonic":
            max_freq = 4000.0
            for comp in analysis_results["comparisons"]:
                for spec in (comp["ml_spectrum"], comp["dft_spectrum"]):
                    if len(spec.frequencies) > 0:
                        max_freq = max(max_freq, float(np.max(spec.frequencies)))
            upper = max_freq + 200.0  # 200 cm-1 margin beyond highest band
            self._analyzer.freq_range = (400, upper)
            self._analyzer.freq_grid = np.arange(400, upper, self._analyzer.freq_step)

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
            cards.append(
                f'<div class="{cls}">'
                f"<h3>{name}</h3>"
                f'<div class="metric">'
                f'<span class="metric-label">MAE (freq)</span>'
                f'<span class="metric-value">{entry.get("mae_freq", float("nan")):.1f} cm\u207b\u00b9</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">MAE (intensity)</span>'
                f'<span class="metric-value">{entry.get("mae_intensity", float("nan")):.1f} km/mol</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">R\u00b2 (freq)</span>'
                f'<span class="metric-value">{self._r2_text(entry["r2_freq"], entry.get("n_freq"))}</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">RMSE</span>'
                f'<span class="metric-value">{entry["rmse_freq"]:.1f} cm\u207b\u00b9</span></div>'
                f'<div class="metric">'
                f'<span class="metric-label">Speedup</span>'
                f'<span class="metric-value">{entry["speedup"]:.1f}\u00d7</span></div>'
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

    def _create_pareto_section(self, comparisons: list[dict]) -> str:
        """Build a compact cost-vs-accuracy Pareto scatter (MAE vs speedup)."""
        names: list[str] = []
        maes: list[float] = []
        speedups: list[float] = []
        for comp in comparisons:
            speedup = comp.get("speedup", 0.0)
            mae = getattr(comp.get("metrics"), "mae_freq", None)
            if mae is None or speedup <= 0:
                continue
            names.append(comp["name"])
            maes.append(float(mae))
            speedups.append(float(speedup))
        if len(names) < 2:
            return ""

        fig = build_pareto_figure(names, maes, speedups)
        div = self._fig_to_div(fig, "pareto")
        return (
            '<section class="comparison-section" id="pareto">'
            "<h2>Cost vs Accuracy</h2>"
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

        # Get paired arrays for regression plots via mode matching;
        # fall back to raw spectrum arrays when mode IDs are absent.
        mode_mapping = comp.get("mode_mapping")
        mode_overlaps_dict = comp.get("mode_overlaps")
        matched_dft_freq, matched_ml_freq, matched_dft_int, matched_ml_int, match_stats = (
            self._analyzer.match_by_mode(
                comp["dft_spectrum"], comp["ml_spectrum"],
                mode_mapping=mode_mapping,
                mode_overlaps=mode_overlaps_dict,
            )
        )
        matched_ids = match_stats.get("matched_mode_ids")
        matched_overlaps = match_stats.get("matched_mode_overlaps")
        if len(matched_dft_freq) == 0:
            dft_s, ml_s = comp["dft_spectrum"], comp["ml_spectrum"]
            n = min(len(dft_s.frequencies), len(ml_s.frequencies))
            matched_dft_freq = np.asarray(dft_s.frequencies[:n], dtype=float)
            matched_ml_freq = np.asarray(ml_s.frequencies[:n], dtype=float)
            matched_dft_int = np.asarray(dft_s.intensities[:n], dtype=float)
            matched_ml_int = np.asarray(ml_s.intensities[:n], dtype=float)
            matched_ids = None
            matched_overlaps = None

        # Per-category MAE for category awards
        cat_mae: dict[str, float | None] = {"fundamental": None, "overtone": None, "combination": None}
        if matched_ids is not None and len(matched_ids) == len(matched_dft_freq):
            for cat, prefix in (("fundamental", "F"), ("overtone", "O"), ("combination", "C")):
                mask = [mid.startswith(prefix) for mid in matched_ids]
                if any(mask):
                    m = np.array(mask)
                    cat_mae[cat] = float(np.mean(np.abs(matched_ml_freq[m] - matched_dft_freq[m])))
        comp["_cat_mae"] = cat_mae
        comp["_matched_overlaps"] = matched_overlaps

        # Frequency regression (left) — skip if no matched modes
        reg_div = ""
        if len(matched_dft_freq) > 0:
            reg_fig = build_regression_figure(
                matched_dft_freq, matched_ml_freq, ml_name,
                mode_ids=matched_ids,
                mode_overlaps=matched_overlaps,
            )
            reg_div = self._fig_to_div(reg_fig, f"regression-{index}")

        # Intensity regression (right) — filter symmetry-forbidden modes
        # Keep a mode if either DFT or ML intensity >= threshold
        int_reg_div = ""
        if len(matched_dft_int) > 0:
            int_mask = (matched_dft_int >= 0.1) | (matched_ml_int >= 0.1)
            if np.sum(int_mask) > 1:
                int_ids = [matched_ids[i] for i, m in enumerate(int_mask) if m] if matched_ids else None
                int_overlaps = (
                    [matched_overlaps[i] for i, m in enumerate(int_mask) if m]
                    if matched_overlaps else None
                )
                int_fig = build_intensity_regression_figure(
                    matched_dft_int[int_mask], matched_ml_int[int_mask], ml_name,
                    mode_ids=int_ids,
                    mode_overlaps=int_overlaps,
                )
                int_reg_div = self._fig_to_div(int_fig, f"int-regression-{index}")

        # Residual plot + error histogram
        residual_row = ""
        res_div = ""
        hist_div = ""
        if len(matched_dft_freq) > 0:
            res_fig = build_residual_figure(
                matched_dft_freq, matched_ml_freq, ml_name, mode_ids=matched_ids
            )
            res_div = self._fig_to_div(res_fig, f"residual-{index}")
            hist_fig = build_error_histogram_figure(
                matched_dft_freq, matched_ml_freq, ml_name, mode_ids=matched_ids
            )
            hist_div = self._fig_to_div(hist_fig, f"error-hist-{index}")
            residual_row = (
                '<div style="display:flex;gap:1rem;flex-wrap:wrap;overflow:hidden">'
                f'<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">{res_div}</div>'
                f'<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">{hist_div}</div>'
                "</div>"
            )

        # Per-region accuracy breakdown
        region_html = ""
        if len(matched_dft_freq) > 0:
            region_data = build_per_region_table(
                matched_dft_freq, matched_ml_freq, mode_ids=matched_ids
            )
            if region_data:
                rows = ""
                for rname, rd in region_data.items():
                    bias_sign = "+" if rd["bias"] >= 0 else ""
                    rows += (
                        f"<tr><td>{self._esc(rname)}</td>"
                        f"<td>{rd['n']}</td>"
                        f"<td>{rd['mae']:.1f}</td>"
                        f"<td>{rd['rmse']:.1f}</td>"
                        f"<td>{bias_sign}{rd['bias']:.1f}</td></tr>"
                    )
                region_html = (
                    '<div class="stats-box">'
                    "<h4>Accuracy by wavenumber range and band type</h4>"
                    '<table class="summary-table"><thead>'
                    "<tr><th>Rows: DFT wavenumber range, then band type</th><th>n</th><th>MAE</th>"
                    "<th>RMSE</th><th>Bias</th></tr></thead>"
                    f"<tbody>{rows}</tbody></table></div>"
                )

        # Anharmonicity ratio (only in anharmonic mode)
        anharm_div = ""
        if self.mode == "anharmonic":
            ml_results = comp.get("_ml_results")
            dft_results = comp.get("_dft_results")
            if ml_results and dft_results:
                anharm_div = self._build_anharmonicity_section(
                    ml_results, dft_results, ml_name, index, mode_mapping=mode_mapping
                )

        # Heatmap PNG (stays as static image)
        heatmap_html = ""
        plots_dir = self.output_dir / "plots"
        heatmap_files = (
            list(plots_dir.glob(f"mode_overlap_{ml_name}_*.png")) if plots_dir.exists() else []
        )
        if heatmap_files:
            import base64

            img_data = heatmap_files[0].read_bytes()
            b64 = base64.b64encode(img_data).decode("ascii")
            heatmap_html = (
                f'<div class="plot-container">'
                f'<img src="data:image/png;base64,{b64}" alt="Mode overlap heatmap">'
                f"</div>"
            )

        # Timing — prefer gaussian_timing, fall back to runtime_s
        ml_gauss_s = comp.get("ml_gaussian_timing", {}).get("total_elapsed_s", 0.0)
        if ml_gauss_s == 0.0:
            ml_gauss_s = comp.get("ml_runtime", 0.0)
        dft_gauss_s = comp.get("dft_gaussian_timing", {}).get("total_elapsed_s", 0.0)
        if dft_gauss_s == 0.0:
            dft_gauss_s = comp.get("dft_runtime", 0.0)
        speedup = comp.get("speedup", 0.0)

        # Degenerate notes — compact single-line summary
        deg_html = ""
        deg_result = comp.get("deg_result")
        if deg_result is not None:
            groups = getattr(deg_result, "groups", None) or []
            if groups:
                items = []
                for g in groups:
                    label = self._esc(getattr(g, "symmetry_label", ""))
                    mult = getattr(g, "multiplicity", 0)
                    overlap = getattr(g, "subspace_overlap", 0.0)
                    tag = f"{label} " if label else ""
                    items.append(f"{tag}{mult}-fold ({overlap:.2f})")
                deg_html = (
                    f'<div class="degenerate-note">'
                    f"{len(groups)} degenerate group{'s' if len(groups) != 1 else ''}: "
                    + ", ".join(items)
                    + "</div>"
                )

        # Per-method metrics — punchy stat boxes with timing integrated
        m = comp["metrics"]

        r2_class = (
            "metric-good" if m.r2_freq > 0.95
            else "metric-warning" if m.r2_freq > 0.90
            else "metric-bad"
        )
        mae_class = (
            "metric-good" if m.mae_freq < 10
            else "metric-warning" if m.mae_freq < 20
            else "metric-bad"
        )

        metrics_table = (
            '<div class="stats-box">'
            "<h4>Statistical Metrics</h4>"
            '<div class="stats-grid">'
            f'<div class="stat-item">'
            f'<div class="stat-label">R\u00b2 (Frequency)</div>'
            f'<div class="stat-value {r2_class}">{self._r2_text(m.r2_freq, m.num_peaks, 4)}</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">R\u00b2 (Intensity)</div>'
            f'<div class="stat-value">{self._r2_text(m.r2_intensity, m.num_peaks - m.num_intensity_filtered, 4)}</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">MAE</div>'
            f'<div class="stat-value {mae_class}">{m.mae_freq:.2f} cm\u207b\u00b9</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">RMSE</div>'
            f'<div class="stat-value">{m.rmse_freq:.2f} cm\u207b\u00b9</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">MAE (Intensity)</div>'
            f'<div class="stat-value">{m.mae_intensity:.2f} km/mol</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">RMSE (Intensity)</div>'
            f'<div class="stat-value">{m.rmse_intensity:.2f} km/mol</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">Max Error</div>'
            f'<div class="stat-value">{m.max_error_freq:.2f} cm\u207b\u00b9</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">Matched Modes</div>'
            f'<div class="stat-value">{m.num_matched}/{m.num_matched + m.num_dft_only}</div></div>'
            f"{self._low_overlap_stat(comp)}"
            f'<div class="stat-item">'
            f'<div class="stat-label">ML Pipeline</div>'
            f'<div class="stat-value">{ml_gauss_s:.1f}s</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">DFT Pipeline</div>'
            f'<div class="stat-value">{dft_gauss_s:.1f}s</div></div>'
            f'<div class="stat-item">'
            f'<div class="stat-label">Speedup</div>'
            f'<div class="stat-value">{speedup:.1f}\u00d7</div></div>'
            "</div></div>"
        )

        # Side-by-side regression container (skip entirely if no matched modes)
        reg_row = ""
        if reg_div or int_reg_div:
            reg_left = (
                f'<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">'
                f"{reg_div}</div>"
            ) if reg_div else ""
            reg_right = (
                '<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">'
                f"{int_reg_div}</div>"
            ) if int_reg_div else ""
            reg_row = (
                '<div style="display:flex;gap:1rem;flex-wrap:wrap;overflow:hidden">'
                f"{reg_left}{reg_right}"
                "</div>"
            )

        # Residual gets full width; anharmonicity (square) pairs with error histogram
        if anharm_div and hist_div:
            anharm_hist_row = (
                '<div style="display:flex;gap:1rem;flex-wrap:wrap;overflow:hidden">'
                f'<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">{anharm_div}</div>'
                f'<div class="plot-container" style="flex:1 1 0;min-width:0;overflow:hidden">{hist_div}</div>'
                "</div>"
            )
            residual_full = (
                f'<div class="plot-container">{res_div}</div>' if res_div else ""
            )
        else:
            anharm_hist_row = f'<div class="plot-container">{anharm_div}</div>' if anharm_div else ""
            residual_full = residual_row

        return (
            f'<section class="comparison-section" id="comparison-{index}">'
            f"<h2>{ml_name_esc} vs DFT</h2>"
            f"{metrics_table}"
            f"{self._create_run_integrity(comp)}"
            f"{deg_html}"
            f'<div class="plot-container">{spec_div}</div>'
            f"{reg_row}"
            f"{residual_full}"
            f"{anharm_hist_row}"
            f"{region_html}"
            f"{heatmap_html}"
            f"</section>"
        )

    def _build_anharmonicity_section(
        self,
        ml_results: dict,
        dft_results: dict,
        ml_name: str,
        index: int,
        mode_mapping: dict[int, int] | None = None,
    ) -> str:
        """Build anharmonicity ratio plot from raw results dicts.

        Pairs modes through the eigenvector mapping (ML checkpoint index -> DFT
        checkpoint index), translating Gaussian's symmetry-block mode numbers to
        checkpoint indices first. Review finding H2b: this used to pair ML Mode(n)
        with DFT Mode(n) by Gaussian's number and ignore the mapping.
        """
        from .analyze_spectra import gaussian_mode_to_checkpoint_index

        try:
            ml_anharm = ml_results.get("frequencies", {}).get("anharmonic", [])
            dft_anharm = dft_results.get("frequencies", {}).get("anharmonic", [])
            if not ml_anharm or not dft_anharm:
                return ""

            ml_to_ckpt = gaussian_mode_to_checkpoint_index(ml_anharm)
            dft_to_ckpt = gaussian_mode_to_checkpoint_index(dft_anharm)
            # checkpoint index (0-based) -> row
            ml_by_ckpt = {ml_to_ckpt.get(m["mode"], m["mode"]) - 1: m for m in ml_anharm}
            dft_by_ckpt = {dft_to_ckpt.get(m["mode"], m["mode"]) - 1: m for m in dft_anharm}
            mapping = mode_mapping if mode_mapping else {i: i for i in ml_by_ckpt}
            pairs = [
                (ml_by_ckpt[i], dft_by_ckpt[j])
                for i, j in sorted(mapping.items())
                if i in ml_by_ckpt and j in dft_by_ckpt
            ]
            if len(pairs) < 3:
                return ""

            dft_harm = np.array([d["freq_harmonic"] for _, d in pairs])
            dft_anh = np.array([d["freq_cm"] for _, d in pairs])
            ml_harm = np.array([m["freq_harmonic"] for m, _ in pairs])
            ml_anh = np.array([m["freq_cm"] for m, _ in pairs])

            fig = build_anharmonicity_ratio_figure(dft_harm, dft_anh, ml_harm, ml_anh, ml_name)
            return self._fig_to_div(fig, f"anharm-ratio-{index}", responsive=False)
        except Exception:
            return ""

    def _create_summary_table(self, comparisons: list[dict]) -> str:
        """Build the overall summary comparison table.

        Columns: method name, R2 (freq), R2 (intensity), RMSE, and speedup.
        Method names are HTML-escaped.
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
            rows.append(
                f"<tr>"
                f"<td>{name}</td>"
                f"<td>{m.mae_freq:.2f}</td>"
                f"<td>{self._r2_text(m.r2_freq, m.num_peaks, 4)}</td>"
                f"<td>{self._r2_text(m.r2_intensity, m.num_peaks - m.num_intensity_filtered, 4)}</td>"
                f"<td>{m.rmse_freq:.2f}</td>"
                f"<td>{speedup:.1f}\u00d7</td>"
                f"</tr>"
            )
        return (
            '<section id="summary-table">'
            "<h2>Summary Comparison</h2>"
            '<table class="data-table">'
            "<thead><tr>"
            "<th>Method</th><th>MAE</th><th>R\u00b2 freq</th><th>R\u00b2 int</th>"
            "<th>RMSE</th><th>Speedup</th>"
            "</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>"
            "</section>"
        )

    @staticmethod
    def _low_overlap_stat(comp: dict) -> str:
        """Render a stat-item showing the count of low-overlap fundamentals.

        Reads the eigenvector-overlap list cached on ``comp`` by the
        regression-figure section.  Empty string when no overlap data is
        available so the box layout collapses gracefully.
        """
        from .plotly_builders import LOW_OVERLAP_THRESHOLD

        overlaps = comp.get("_matched_overlaps")
        if not overlaps:
            return ""
        valid = [o for o in overlaps if o is not None]
        if not valid:
            return ""
        n_low = sum(1 for o in valid if o < LOW_OVERLAP_THRESHOLD)
        cls = "metric-warning" if n_low > 0 else "metric-good"
        return (
            '<div class="stat-item">'
            f'<div class="stat-label">Low-overlap (&lt; {LOW_OVERLAP_THRESHOLD:.1f})</div>'
            f'<div class="stat-value {cls}">{n_low}/{len(valid)}</div></div>'
        )

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
            # Option A: real data — DFT vs ML comparison table
            def _fmt(v, fmt=".1f"):
                return f"{v:{fmt}}" if v is not None else "\u2014"

            # Group by method
            by_method: dict[str, list[dict]] = {}
            for row in overtones_data:
                by_method.setdefault(row["method"], []).append(row)

            _MAX_TABLE_ROWS = 15  # Show full table up to this many rows

            tables = []
            for method, rows in by_method.items():

                # Compute errors for rows that have both DFT and ML
                paired = [r for r in rows if r["dft_freq"] is not None and r["ml_freq"] is not None]
                n_ot = sum(1 for r in rows if r["type"] == "overtone")
                n_cb = sum(1 for r in rows if r["type"] == "combination")
                errors = [abs(r["ml_freq"] - r["dft_freq"]) for r in paired]
                mae = sum(errors) / len(errors) if errors else 0.0

                summary = (
                    f"<h3>{self._esc(method)} vs DFT</h3>"
                    f'<div class="stats-grid" style="margin-bottom:12px">'
                    f'<div class="stat-item"><div class="stat-label">Overtones</div>'
                    f'<div class="stat-value">{n_ot}</div></div>'
                    f'<div class="stat-item"><div class="stat-label">Combinations</div>'
                    f'<div class="stat-value">{n_cb}</div></div>'
                    f'<div class="stat-item"><div class="stat-label">MAE</div>'
                    f'<div class="stat-value">{mae:.1f} cm\u207b\u00b9</div></div>'
                    f"</div>"
                )

                # Show full table for small molecules, top-N worst for large
                if len(rows) <= _MAX_TABLE_ROWS:
                    display_rows = rows
                    table_note = ""
                else:
                    # Sort by absolute error, show worst N
                    paired.sort(key=lambda r: abs(r["ml_freq"] - r["dft_freq"]), reverse=True)
                    display_rows = paired[:_MAX_TABLE_ROWS]
                    table_note = (
                        f'<p style="color:#6b7280;font-size:0.85em;margin-top:4px">'
                        f"Showing {_MAX_TABLE_ROWS} largest errors out of {len(rows)} total entries.</p>"
                    )

                def _delta(r):
                    if r["dft_freq"] is not None and r["ml_freq"] is not None:
                        return f"{r['ml_freq'] - r['dft_freq']:+.1f}"
                    return "\u2014"

                trs = "".join(
                    f"<tr><td>{self._esc(r['mode_id'])}</td>"
                    f"<td>{self._esc(r['type'])}</td>"
                    f"<td>{_fmt(r['dft_freq'])}</td>"
                    f"<td>{_fmt(r['ml_freq'])}</td>"
                    f"<td>{_delta(r)}</td>"
                    f"<td>{_fmt(r['dft_int'], '.3f')}</td>"
                    f"<td>{_fmt(r['ml_int'], '.3f')}</td></tr>"
                    for r in display_rows
                )
                tables.append(
                    summary
                    + f'<table class="overtone-table data-table">'
                    f"<thead><tr><th>Mode</th><th>Type</th>"
                    f"<th>DFT freq (cm\u207b\u00b9)</th>"
                    f"<th>ML freq (cm\u207b\u00b9)</th>"
                    f"<th>\u0394freq (cm\u207b\u00b9)</th>"
                    f"<th>DFT int (km/mol)</th>"
                    f"<th>ML int (km/mol)</th></tr></thead>"
                    f"<tbody>{trs}</tbody></table>"
                    + table_note
                )
            body = "\n".join(tables)
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
        """Collect overtone/combination band records from SpectrumData objects.

        Extracts entries with ``label == "overtone"`` or
        ``label == "combination"`` from both DFT and ML ``SpectrumData``
        in each comparison, then matches them by ``mode_id`` to produce
        side-by-side DFT vs ML rows.

        Returns
        -------
        list[dict]
            Each dict has keys: method, mode_id, dft_freq, dft_int,
            ml_freq, ml_int, type.  Empty list if no overtone data.
        """
        collected: list[dict] = []
        for comp in analysis_results.get("comparisons", []):
            dft_spec = comp.get("dft_spectrum")
            ml_spec = comp.get("ml_spectrum")
            if dft_spec is None or ml_spec is None:
                continue

            method = comp.get("name", "ML")

            # Build {mode_id: (freq, intensity)} for non-fundamental entries
            def _extract_ot(spec):
                out = {}
                for i, label in enumerate(spec.labels):
                    if label in ("overtone", "combination"):
                        out[spec.mode_ids[i]] = (
                            float(spec.frequencies[i]),
                            float(spec.intensities[i]),
                            label,
                        )
                return out

            dft_ot = _extract_ot(dft_spec)
            ml_ot = _extract_ot(ml_spec)

            # All mode_ids from both sides, sorted
            all_ids = sorted(set(dft_ot) | set(ml_ot))
            for mid in all_ids:
                d = dft_ot.get(mid)
                m = ml_ot.get(mid)
                collected.append(
                    {
                        "method": method,
                        "mode_id": mid,
                        "dft_freq": d[0] if d else None,
                        "dft_int": d[1] if d else None,
                        "ml_freq": m[0] if m else None,
                        "ml_int": m[1] if m else None,
                        "type": (d or m)[2],
                    }
                )
        return collected

    def _create_category_awards(self, comparisons: list[dict]) -> str:
        """Render a category awards section crowning the best ML method per mode type.

        Shows best MAE for fundamentals, overtones, and combination bands
        separately, so the user can see which model excels at each.
        """
        categories = [
            ("fundamental", "Fundamentals"),
            ("overtone", "Overtones"),
            ("combination", "Combination Bands"),
        ]

        cards = []
        for cat_key, cat_label in categories:
            # Collect (method_name, mae) for this category
            entries = []
            for comp in comparisons:
                cat_mae = comp.get("_cat_mae", {})
                mae = cat_mae.get(cat_key)
                if mae is not None:
                    entries.append((comp["name"], mae))

            if not entries:
                cards.append(
                    f'<div class="award-card">'
                    f"<h3>{self._esc(cat_label)}</h3>"
                    f'<div class="award-no-data">No data</div>'
                    f"</div>"
                )
                continue

            entries.sort(key=lambda x: x[1])
            best_name, best_mae = entries[0]

            rows = "".join(
                f'<tr class="{"award-winner" if i == 0 else ""}">'
                f"<td>{i + 1}</td>"
                f"<td>{self._esc(name)}</td>"
                f"<td>{mae:.1f}</td></tr>"
                for i, (name, mae) in enumerate(entries)
            )

            cards.append(
                f'<div class="award-card">'
                f"<h3>{self._esc(cat_label)}</h3>"
                f'<div class="award-winner-name">{self._esc(best_name)}</div>'
                f'<div class="award-winner-mae">MAE: {best_mae:.1f} cm\u207b\u00b9</div>'
                f'<table class="data-table award-table">'
                f"<thead><tr><th>#</th><th>Method</th>"
                f"<th>MAE (cm\u207b\u00b9)</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>"
                f"</div>"
            )

        return (
            '<section class="comparison-section" id="per-category-accuracy">'
            "<h2>Per-Category Accuracy Ranking</h2>"
            '<div class="awards-grid">'
            + "\n".join(cards)
            + "</div></section>"
        )

    def _create_footer(self) -> str:
        """Create the report footer with generation info."""
        return (
            f'<footer class="footer">'
            f"Generated by mace-gaussian {self._esc(self.mode)} analysis "
            f"\u00b7 Phase 23"
            f"</footer>"
        )
