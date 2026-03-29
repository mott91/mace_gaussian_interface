"""
HTML Report Generator

Creates comprehensive, beautiful HTML reports with embedded plots,
tables, and statistics for IR spectral analysis.
"""

import base64
from datetime import datetime
from pathlib import Path

import pandas as pd


class HTMLReportGenerator:
    """Generates HTML reports for spectral analysis"""

    def __init__(self, molecule_name: str, output_dir: Path, bandwidth_fwhm: float = 10.0):
        """
        Initialize report generator

        Parameters
        ----------
        molecule_name : str
            Name of molecule
        output_dir : Path
            Output directory containing plots and data
        bandwidth_fwhm : float
            Full width at half maximum used for Lorentzian broadening (cm-1)
        """
        self.molecule_name = molecule_name
        self.output_dir = Path(output_dir)
        self.bandwidth_fwhm = bandwidth_fwhm
        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"

    def encode_image(self, image_path: Path) -> str:
        """Encode image to base64 for embedding"""
        with Path(image_path).open("rb") as f:
            encoded = base64.b64encode(f.read()).decode()
        return f"data:image/png;base64,{encoded}"

    def create_css(self) -> str:
        """Create CSS styling"""
        return """
        <style>
            * {
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }

            body {
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                line-height: 1.6;
                color: #1a1a1a;
                background: #f5f5f5;
                padding: 0;
            }

            .container {
                max-width: 100%;
                margin: 0;
                background: white;
                min-height: 100vh;
            }

            header {
                background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
                color: white;
                padding: 3rem 4rem;
                border-bottom: 3px solid #3498db;
            }

            header h1 {
                font-size: 2.2em;
                margin-bottom: 0.5rem;
                font-weight: 600;
                letter-spacing: -0.5px;
            }

            header .subtitle {
                font-size: 1.1em;
                opacity: 0.85;
                font-weight: 300;
            }

            nav {
                background: white;
                padding: 1rem 4rem;
                border-bottom: 1px solid #e5e7eb;
                position: sticky;
                top: 0;
                z-index: 100;
                box-shadow: 0 1px 3px rgba(0,0,0,0.05);
            }

            nav a {
                color: #374151;
                text-decoration: none;
                margin-right: 2rem;
                font-weight: 500;
                font-size: 0.95rem;
                transition: all 0.2s;
                padding-bottom: 0.25rem;
                border-bottom: 2px solid transparent;
            }

            nav a:hover {
                color: #3498db;
                border-bottom-color: #3498db;
            }

            .content {
                padding: 3rem 4rem;
                max-width: 1600px;
                margin: 0 auto;
            }

            section {
                margin-bottom: 4rem;
            }

            h2 {
                color: #1a1a1a;
                font-size: 1.75em;
                margin-bottom: 1.5rem;
                padding-bottom: 0.75rem;
                border-bottom: 2px solid #e5e7eb;
                font-weight: 600;
                letter-spacing: -0.3px;
            }

            h3 {
                color: #374151;
                font-size: 1.3em;
                margin: 2rem 0 1rem 0;
                font-weight: 600;
            }

            .summary-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 1.5rem;
                margin: 2rem 0;
            }

            .summary-card {
                background: white;
                border: 1px solid #e5e7eb;
                padding: 1.5rem;
                border-radius: 8px;
                transition: all 0.2s;
            }

            .summary-card:hover {
                box-shadow: 0 4px 12px rgba(0,0,0,0.08);
                transform: translateY(-2px);
            }

            .summary-card .value {
                font-size: 2.2em;
                font-weight: 600;
                margin: 0.5rem 0;
                color: #1a1a1a;
            }

            .summary-card .label {
                font-size: 0.85rem;
                color: #6b7280;
                text-transform: uppercase;
                letter-spacing: 0.5px;
                font-weight: 500;
            }

            .comparison-section {
                background: #fafafa;
                padding: 2rem;
                border-radius: 8px;
                margin-bottom: 2.5rem;
                border: 1px solid #e5e7eb;
            }

            .plot-container {
                background: white;
                padding: 1.5rem;
                border-radius: 6px;
                border: 1px solid #e5e7eb;
                margin: 1.5rem 0;
            }

            .plot-container img {
                width: 100%;
                height: auto;
                display: block;
                border-radius: 4px;
            }

            table {
                width: 100%;
                border-collapse: collapse;
                margin: 1.5rem 0;
                background: white;
                border: 1px solid #e5e7eb;
                border-radius: 6px;
                overflow: hidden;
                font-size: 0.9rem;
            }

            thead {
                background: #f9fafb;
                border-bottom: 2px solid #e5e7eb;
            }

            th, td {
                padding: 0.875rem 1rem;
                text-align: left;
                border-bottom: 1px solid #f3f4f6;
            }

            th {
                font-weight: 600;
                color: #374151;
                font-size: 0.85rem;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }

            td {
                color: #4b5563;
            }

            tbody tr:hover {
                background: #f9fafb;
            }

            .metric-good {
                color: #059669;
                font-weight: 600;
            }

            .metric-warning {
                color: #d97706;
                font-weight: 600;
            }

            .metric-bad {
                color: #dc2626;
                font-weight: 600;
            }

            .stats-box {
                background: white;
                padding: 1.5rem;
                border-radius: 6px;
                border: 1px solid #e5e7eb;
                margin: 1.5rem 0;
            }

            .stats-box h4 {
                color: #374151;
                margin-bottom: 1rem;
                font-weight: 600;
                font-size: 1.1rem;
            }

            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 1rem;
            }

            .stat-item {
                padding: 1rem;
                background: #f9fafb;
                border-radius: 6px;
                border: 1px solid #f3f4f6;
            }

            .stat-label {
                font-size: 0.8rem;
                color: #6b7280;
                text-transform: uppercase;
                letter-spacing: 0.5px;
                margin-bottom: 0.25rem;
                font-weight: 500;
            }

            .stat-value {
                font-size: 1.5em;
                font-weight: 600;
                color: #1a1a1a;
            }

            .warning-box {
                background: #fef3c7;
                border-left: 3px solid #f59e0b;
                padding: 1rem;
                margin: 1.5rem 0;
                border-radius: 4px;
                color: #92400e;
            }

            footer {
                background: #f9fafb;
                padding: 2rem;
                text-align: center;
                color: #6b7280;
                border-top: 1px solid #e5e7eb;
                font-size: 0.9rem;
            }

            .timestamp {
                font-style: italic;
                color: #9ca3af;
                margin-top: 0.5rem;
            }
        </style>
        """

    def create_header(self) -> str:
        """Create HTML header"""
        return f"""
        <header>
            <h1>IR Spectral Analysis Report</h1>
            <div class="subtitle">{self.molecule_name.upper()}</div>
            <div class="subtitle">ML vs DFT Comparison</div>
        </header>
        """

    def create_navigation(self, comparisons: list[dict]) -> str:
        """Create navigation menu"""
        nav_items = ['<a href="#overview">Overview</a>', '<a href="#combined">Combined</a>']
        for i, comp in enumerate(comparisons, 1):
            nav_items.append(f'<a href="#comp{i}">{comp["name"]}</a>')
        nav_items.append('<a href="#summary">Summary</a>')

        return f"<nav>{' '.join(nav_items)}</nav>"

    def create_mode_overlap_section(self) -> str:
        """Create section with mode overlap heatmaps"""
        # Find all mode overlap heatmap files
        mode_overlap_files = list(self.plots_dir.glob("mode_overlap_*.png"))

        if not mode_overlap_files:
            return ""

        # Build heatmap sections
        heatmap_sections = []
        for heatmap_file in sorted(mode_overlap_files):
            filename = heatmap_file.name
            # Extract method names from filename (e.g., mode_overlap_mace_mp_espaloma_vs_wb97xd.png)
            parts = filename.replace("mode_overlap_", "").replace(".png", "").split("_vs_")
            ml_method = parts[0] if len(parts) > 0 else "ML"
            dft_method = parts[1] if len(parts) > 1 else "DFT"

            heatmap_sections.append(f"""
            <h3>Mode Overlap: {ml_method.upper().replace("_", " ")} vs {dft_method.upper()}</h3>
            <div class="plot-container">
                <img src="plots/{filename}" alt="Mode overlap" style="width: 100%;">
            </div>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 30px;">
                Heatmap shows the overlap (dot product) between vibrational modes.
                Red = high overlap (modes match),
                Blue = weak overlap, White = no overlap (modes are orthogonal).
            </p>
            """)

        heatmaps_html = "\n".join(heatmap_sections)

        return f"""
        <section id="mode-overlap" class="comparison-section">
            <h2>Mode Matching Analysis</h2>
            <p style="color: #666; margin-bottom: 20px;">
                Vibrational modes can change order between DFT and ML calculations.
                Mode matching via normal mode overlap ensures we compare the correct physical modes.
            </p>

            <div style="padding: 15px; background: #fff3cd; border-left: 4px solid #ffc107;">
                <strong>Important:</strong> Off-diagonal red spots indicate mode reordering.
                Perfect overlap = 1.00 (dark red), orthogonal modes = 0.00 (white).
            </div>

            {heatmaps_html}
        </section>
        """

    def create_combined_plots_section(self) -> str:
        """Create section with combined plots showing all ML methods"""
        combined_spectrum = self.plots_dir / "spectrum_combined.png"
        combined_spectrum_extended = self.plots_dir / "spectrum_combined_extended.png"
        combined_regression = self.plots_dir / "regression_combined.png"

        # Check if combined plots exist
        if not combined_spectrum.exists() or not combined_regression.exists():
            return ""

        # Encode images as base64
        spectrum_img = self.encode_image(combined_spectrum)
        regression_img = self.encode_image(combined_regression)

        # Build extended spectrum section if it exists
        extended_section = ""
        if combined_spectrum_extended.exists():
            extended_img = self.encode_image(combined_spectrum_extended)
            extended_section = f"""
            <h3>Extended IR Spectrum (400-8000 cm⁻¹, includes overtones)</h3>
            <div class="plot-container">
                <img src="{extended_img}" alt="Extended spectrum" style="width: 100%;">
            </div>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 30px;">
                Extended frequency range shows fundamental modes, overtones, and combination bands.
            </p>
            """

        return f"""
        <section id="combined" class="comparison-section">
            <h2>Combined Analysis: All Methods</h2>
            <p style="color: #666; margin-bottom: 20px;">
                Direct comparison of all ML methods against the DFT baseline in a single view.
            </p>

            <h3>Combined IR Spectrum (Fundamentals)</h3>
            <div class="plot-container">
                <img src="{spectrum_img}" alt="Combined spectrum" style="width: 100%;">
            </div>

            {extended_section}

            <h3>Combined Regression Analysis</h3>
            <div class="plot-container">
                <img src="{regression_img}" alt="Combined regression" style="width: 100%;">
            </div>

            <p style="margin-top: 20px; padding: 15px; background: #f8f9fa;">
                <strong>Note:</strong> Combined plots compare ML methods vs DFT reference.
                R² values of "N/A" indicate insufficient data points (N&lt;3).
            </p>
        </section>
        """

    def create_overview(self, comparisons: list[dict]) -> str:
        """Create overview section"""
        # FIXED: Check list length instead of truthy value
        if len(comparisons) == 0:
            return """
            <section id="overview">
                <h2>Overview</h2>
                <div class="warning-box">
                    <strong>No successful comparisons found.</strong>
                    <p>All ML calculators either failed or produced no matching peaks.</p>
                </div>
            </section>
            """

        num_comparisons = len(comparisons)
        avg_mae = sum(c["metrics"].mae_freq for c in comparisons) / num_comparisons
        best_r2 = max(c["metrics"].r2_freq for c in comparisons)
        avg_speedup = sum(c["speedup"] for c in comparisons) / num_comparisons

        return f"""
        <section id="overview">
            <h2>Overview</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="label">ML Calculators Tested</div>
                    <div class="value">{num_comparisons}</div>
                </div>
                <div class="summary-card">
                    <div class="label">Average MAE</div>
                    <div class="value">{avg_mae:.1f}</div>
                    <div class="label">cm^-1</div>
                </div>
                <div class="summary-card">
                    <div class="label">Best R^2</div>
                    <div class="value">{best_r2:.4f}</div>
                </div>
                <div class="summary-card">
                    <div class="label">Avg Speedup</div>
                    <div class="value">{avg_speedup:.1f}x</div>
                </div>
            </div>

            <p style="margin-top: 30px; font-size: 1.1em; line-height: 1.8;">
                This report compares machine learning (ML) force fields
                against density functional theory (DFT) anharmonic calculations for IR spectroscopy.
                All spectra include anharmonic fundamentals, overtones, and combination bands.
                Broadening was applied using Lorentzian line shapes
                with {self.bandwidth_fwhm} cm<sup>-1</sup> FWHM.
                Modes with DFT IR intensity below 0.1 km/mol were excluded
                from intensity regression metrics.
            </p>
        </section>
        """

    def create_comparison_section(
        self, comparison: dict, index: int, dft_method: str = "DFT"
    ) -> str:
        """Create detailed comparison section for one ML calculator"""
        m = comparison["metrics"]
        ml_name = comparison["name"]

        # Determine quality coloring
        r2_class = (
            "metric-good"
            if m.r2_freq > 0.95
            else "metric-warning"
            if m.r2_freq > 0.90
            else "metric-bad"
        )
        mae_class = (
            "metric-good"
            if m.mae_freq < 10
            else "metric-warning"
            if m.mae_freq < 20
            else "metric-bad"
        )

        # Find corresponding mode overlap heatmap
        mode_overlap_files = list(self.plots_dir.glob(f"mode_overlap_{ml_name}_*.png"))
        heatmap_section = ""
        if mode_overlap_files:
            heatmap_file = mode_overlap_files[0]
            heatmap_img = self.encode_image(heatmap_file)
            heatmap_section = f"""
            <h3>Mode Overlap Matrix</h3>
            <div class="plot-container">
                <img src="{heatmap_img}" alt="Mode overlap heatmap" style="width: 100%;">
            </div>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 20px;">
                Heatmap shows the eigenvector overlap between vibrational modes.
                Red = high overlap, White = no overlap. Perfect overlap = 1.00, orthogonal = 0.00.
                Off-diagonal red spots indicate mode reordering between ML and DFT.
            </p>
            """

        # Read comparison table - FIXED: Handle empty DataFrames properly
        table_path = self.data_dir / comparison["table_file"]

        # FIXED: Use proper pandas empty check
        try:
            df = pd.read_csv(table_path)
            if not df.empty:
                # Custom formatter for Mode_Overlap column
                formatters = {}
                if "Mode_Overlap" in df.columns:
                    formatters["Mode_Overlap"] = lambda x: f"{x:.3f}" if pd.notna(x) else "N/A"

                table_html = df.head(20).to_html(
                    index=False,
                    float_format=lambda x: f"{x:.2f}",
                    formatters=formatters,
                    classes="data-table",
                    na_rep="N/A",
                )
            else:
                table_html = (
                    '<p class="warning-box">No matched peaks found for this calculator.</p>'
                )
        except Exception as e:
            table_html = f'<p class="warning-box">Error loading comparison table: {e}</p>'

        spectrum_img = self.encode_image(self.plots_dir / comparison["spectrum_plot"])
        regression_img = self.encode_image(self.plots_dir / comparison["regression_plot"])
        return f"""
        <section id="comp{index}" class="comparison-section">
            <h2>{comparison["name"]}</h2>

            {heatmap_section}

            <div class="stats-box">
                <h4>Statistical Metrics</h4>
                <div class="stats-grid">
                    <div class="stat-item">
                        <div class="stat-label">R^2 (Frequency)</div>
                        <div class="stat-value {r2_class}">{m.r2_freq:.4f}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">MAE</div>
                        <div class="stat-value {mae_class}">{m.mae_freq:.2f} cm^-1</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">RMSE</div>
                        <div class="stat-value">{m.rmse_freq:.2f} cm^-1</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Max Error</div>
                        <div class="stat-value">{m.max_error_freq:.2f} cm^-1</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Matched Modes</div>
                        <div class="stat-value">{m.num_matched}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Match Rate</div>
                        <div class="stat-value">{m.match_rate:.1%}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Missing Modes (DFT only)</div>
                        <div class="stat-value">{m.num_dft_only}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Spurious Modes (ML only)</div>
                        <div class="stat-value">{m.num_ml_only}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Speedup</div>
                        <div class="stat-value">{comparison["speedup"]:.1f}x</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">ML Runtime</div>
                        <div class="stat-value">{comparison["ml_runtime"]:.1f} s</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">DFT Runtime</div>
                        <div class="stat-value">{comparison["dft_runtime"]:.1f} s</div>
                    </div>
                </div>
            </div>

            <h3>Spectral Comparison</h3>
            <div class="plot-container">
                <img src="{spectrum_img}" alt="Spectrum">
            </div>

            <h3>Regression Analysis</h3>
            <div class="plot-container">
                <img src="{regression_img}" alt="Regression">
            </div>

            <h3>Detailed Frequency Comparison</h3>
            <p>First 20 matched peaks (full table: data/{comparison["table_file"]})</p>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
                <strong>Mode_Overlap</strong>: eigenvector dot product (1.0=perfect, 0.0=orth.).
            </p>
            {table_html}
        </section>
        """

    def create_summary_table(self, comparisons: list[dict]) -> str:
        """Create summary comparison table"""
        # FIXED: Check list length
        if len(comparisons) == 0:
            return """
            <section id="summary">
                <h2>Summary Comparison</h2>
                <div class="warning-box">
                    <p>No comparisons available to summarize.</p>
                </div>
            </section>
            """

        rows = []
        for comp in comparisons:
            m = comp["metrics"]
            rows.append(f"""
            <tr>
                <td>{comp["name"]}</td>
                <td>{m.r2_freq:.4f}</td>
                <td>{m.mae_freq:.2f}</td>
                <td>{m.rmse_freq:.2f}</td>
                <td>{m.max_error_freq:.2f}</td>
                <td>{m.num_matched}/{m.num_matched + m.num_dft_only}</td>
                <td>{m.match_rate:.1%}</td>
                <td>{comp["speedup"]:.1f}x</td>
            </tr>
            """)

        return f"""
        <section id="summary">
            <h2>Summary Comparison</h2>
            <table>
                <thead>
                    <tr>
                        <th>Calculator</th>
                        <th>R^2</th>
                        <th>MAE (cm^-1)</th>
                        <th>RMSE (cm^-1)</th>
                        <th>Max Error (cm^-1)</th>
                        <th>Matched Modes</th>
                        <th>Match Rate</th>
                        <th>Speedup</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(rows)}
                </tbody>
            </table>
        </section>
        """

    def create_footer(self) -> str:
        """Create footer"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
        <footer>
            <p>Generated by IR Spectral Analysis Framework</p>
            <p class="timestamp">Report generated: {timestamp}</p>
        </footer>
        """

    def generate_report(self, analysis_results: dict):
        """
        Generate complete HTML report

        Parameters
        ----------
        analysis_results : dict
            Results from comparison workflow
        """
        comparisons = analysis_results["comparisons"]

        html_parts = [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="UTF-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
            f"<title>IR Analysis: {self.molecule_name}</title>",
            self.create_css(),
            "</head>",
            "<body>",
            '<div class="container">',
            self.create_header(),
            self.create_navigation(comparisons),
            '<div class="content">',
            self.create_overview(comparisons),
            self.create_combined_plots_section(),  # Add combined plots after overview
        ]

        # Add comparison sections (now includes mode overlap heatmaps within each section)
        for i, comp in enumerate(comparisons, 1):
            html_parts.append(self.create_comparison_section(comp, i))

        # Add summary and footer
        html_parts.extend(
            [
                self.create_summary_table(comparisons),
                "</div>",  # content
                self.create_footer(),
                "</div>",  # container
                "</body>",
                "</html>",
            ]
        )

        # Write to file
        output_path = self.output_dir / "report.html"
        with output_path.open("w") as f:
            f.write("\n".join(html_parts))

        print(f"\n{'=' * 60}")
        print("HTML REPORT GENERATED")
        print(f"{'=' * 60}")
        print(f"Location: {output_path}")
        print(f"Open in browser: file://{output_path.absolute()}")
        print(f"{'=' * 60}\n")
