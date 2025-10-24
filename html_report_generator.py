"""
HTML Report Generator

Creates comprehensive, beautiful HTML reports with embedded plots,
tables, and statistics for IR spectral analysis.
"""

import base64
from pathlib import Path
from typing import Dict, List
import pandas as pd
from datetime import datetime


class HTMLReportGenerator:
    """Generates HTML reports for spectral analysis"""
    
    def __init__(self, molecule_name: str, output_dir: Path):
        """
        Initialize report generator
        
        Parameters
        ----------
        molecule_name : str
            Name of molecule
        output_dir : Path
            Output directory containing plots and data
        """
        self.molecule_name = molecule_name
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"
    
    def encode_image(self, image_path: Path) -> str:
        """Encode image to base64 for embedding"""
        with open(image_path, 'rb') as f:
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
                font-family: 'Helvetica Neue', Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                padding: 20px;
            }
            
            .container {
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                border-radius: 10px;
                box-shadow: 0 10px 40px rgba(0,0,0,0.3);
                overflow: hidden;
            }
            
            header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 40px;
                text-align: center;
            }
            
            header h1 {
                font-size: 2.5em;
                margin-bottom: 10px;
                text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
            }
            
            header .subtitle {
                font-size: 1.2em;
                opacity: 0.9;
            }
            
            nav {
                background: #f8f9fa;
                padding: 15px 40px;
                border-bottom: 2px solid #e0e0e0;
                position: sticky;
                top: 0;
                z-index: 100;
            }
            
            nav a {
                color: #667eea;
                text-decoration: none;
                margin-right: 25px;
                font-weight: 600;
                transition: color 0.3s;
            }
            
            nav a:hover {
                color: #764ba2;
            }
            
            .content {
                padding: 40px;
            }
            
            section {
                margin-bottom: 60px;
            }
            
            h2 {
                color: #667eea;
                font-size: 2em;
                margin-bottom: 20px;
                padding-bottom: 10px;
                border-bottom: 3px solid #667eea;
            }
            
            h3 {
                color: #764ba2;
                font-size: 1.5em;
                margin: 30px 0 15px 0;
            }
            
            .summary-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px;
                margin: 30px 0;
            }
            
            .summary-card {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 25px;
                border-radius: 10px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                transition: transform 0.3s;
            }
            
            .summary-card:hover {
                transform: translateY(-5px);
            }
            
            .summary-card .value {
                font-size: 2.5em;
                font-weight: bold;
                margin: 10px 0;
            }
            
            .summary-card .label {
                font-size: 0.9em;
                opacity: 0.9;
            }
            
            .comparison-section {
                background: #f8f9fa;
                padding: 30px;
                border-radius: 10px;
                margin: 20px 0;
            }
            
            .plot-container {
                text-align: center;
                margin: 30px 0;
                background: white;
                padding: 20px;
                border-radius: 10px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            
            .plot-container img {
                max-width: 100%;
                height: auto;
                border-radius: 5px;
            }
            
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            
            thead {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
            }
            
            th, td {
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #e0e0e0;
            }
            
            tbody tr:hover {
                background: #f8f9fa;
            }
            
            .metric-good {
                color: #28a745;
                font-weight: bold;
            }
            
            .metric-warning {
                color: #ffc107;
                font-weight: bold;
            }
            
            .metric-bad {
                color: #dc3545;
                font-weight: bold;
            }
            
            .stats-box {
                background: white;
                padding: 20px;
                border-radius: 10px;
                border-left: 4px solid #667eea;
                margin: 20px 0;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            
            .stats-box h4 {
                color: #667eea;
                margin-bottom: 15px;
            }
            
            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
            }
            
            .stat-item {
                padding: 10px;
                background: #f8f9fa;
                border-radius: 5px;
            }
            
            .stat-label {
                font-size: 0.9em;
                color: #666;
            }
            
            .stat-value {
                font-size: 1.4em;
                font-weight: bold;
                color: #667eea;
            }
            
            footer {
                background: #f8f9fa;
                padding: 30px;
                text-align: center;
                color: #666;
                border-top: 2px solid #e0e0e0;
            }
            
            .timestamp {
                font-style: italic;
                color: #999;
                margin-top: 10px;
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
    
    def create_navigation(self, comparisons: List[Dict]) -> str:
        """Create navigation menu"""
        nav_items = ['<a href="#overview">Overview</a>']
        for i, comp in enumerate(comparisons, 1):
            nav_items.append(f'<a href="#comp{i}">{comp["name"]}</a>')
        nav_items.append('<a href="#summary">Summary</a>')
        
        return f'<nav>{" ".join(nav_items)}</nav>'
    
    def create_overview(self, comparisons: List[Dict]) -> str:
        """Create overview section"""
        num_comparisons = len(comparisons)
        avg_mae = sum(c['metrics'].mae_freq for c in comparisons) / num_comparisons
        best_r2 = max(c['metrics'].r2_freq for c in comparisons)
        avg_speedup = sum(c['speedup'] for c in comparisons) / num_comparisons
        
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
                This report presents a comprehensive comparison of machine learning (ML) force fields 
                against density functional theory (DFT) anharmonic calculations for IR spectroscopy. 
                All spectra include anharmonic fundamentals, overtones, and combination bands. 
                Broadening was applied using Gaussian convolution with 8 cm^-1 FWHM.
            </p>
        </section>
        """
    
    def create_comparison_section(self, comparison: Dict, index: int) -> str:
        """Create detailed comparison section for one ML calculator"""
        m = comparison['metrics']
        
        # Determine quality coloring
        r2_class = 'metric-good' if m.r2_freq > 0.95 else 'metric-warning' if m.r2_freq > 0.90 else 'metric-bad'
        mae_class = 'metric-good' if m.mae_freq < 10 else 'metric-warning' if m.mae_freq < 20 else 'metric-bad'
        
        # Read comparison table
        table_path = self.data_dir / comparison['table_file']
        df = pd.read_csv(table_path)
        table_html = df.head(20).to_html(
            index=False, 
            float_format=lambda x: f'{x:.2f}',
            classes='data-table'
        )
        
        return f"""
        <section id="comp{index}" class="comparison-section">
            <h2>{comparison['name']}</h2>
            
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
                        <div class="stat-label">Matched Peaks</div>
                        <div class="stat-value">{m.num_peaks}</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">Speedup</div>
                        <div class="stat-value">{comparison['speedup']:.1f}x</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">ML Runtime</div>
                        <div class="stat-value">{comparison['ml_runtime']:.1f} s</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-label">DFT Runtime</div>
                        <div class="stat-value">{comparison['dft_runtime']:.1f} s</div>
                    </div>
                </div>
            </div>
            
            <h3>Spectral Comparison</h3>
            <div class="plot-container">
                <img src="plots/{comparison['spectrum_plot']}" alt="Spectrum comparison">
            </div>
            
            <h3>Regression Analysis</h3>
            <div class="plot-container">
                <img src="plots/{comparison['regression_plot']}" alt="Regression plot">
            </div>
            
            <h3>Detailed Frequency Comparison</h3>
            <p>Showing first 20 matched peaks (full table available in data/{comparison['table_file']})</p>
            {table_html}
        </section>
        """
    
    def create_summary_table(self, comparisons: List[Dict]) -> str:
        """Create summary comparison table"""
        rows = []
        for comp in comparisons:
            m = comp['metrics']
            rows.append(f"""
            <tr>
                <td>{comp['name']}</td>
                <td>{m.r2_freq:.4f}</td>
                <td>{m.mae_freq:.2f}</td>
                <td>{m.rmse_freq:.2f}</td>
                <td>{m.max_error_freq:.2f}</td>
                <td>{m.num_peaks}</td>
                <td>{comp['speedup']:.1f}x</td>
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
                        <th>Peaks</th>
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
    
    def generate_report(self, analysis_results: Dict):
        """
        Generate complete HTML report
        
        Parameters
        ----------
        analysis_results : dict
            Results from comparison workflow
        """
        comparisons = analysis_results['comparisons']
        
        html_parts = [
            '<!DOCTYPE html>',
            '<html lang="en">',
            '<head>',
            '<meta charset="UTF-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
            f'<title>IR Analysis: {self.molecule_name}</title>',
            self.create_css(),
            '</head>',
            '<body>',
            '<div class="container">',
            self.create_header(),
            self.create_navigation(comparisons),
            '<div class="content">',
            self.create_overview(comparisons),
        ]
        
        # Add comparison sections
        for i, comp in enumerate(comparisons, 1):
            html_parts.append(self.create_comparison_section(comp, i))
        
        # Add summary and footer
        html_parts.extend([
            self.create_summary_table(comparisons),
            '</div>',  # content
            self.create_footer(),
            '</div>',  # container
            '</body>',
            '</html>'
        ])
        
        # Write to file
        output_path = self.output_dir / "report.html"
        with open(output_path, 'w') as f:
            f.write('\n'.join(html_parts))
        
        print(f"\n{'='*60}")
        print(f"HTML REPORT GENERATED")
        print(f"{'='*60}")
        print(f"Location: {output_path}")
        print(f"Open in browser: file://{output_path.absolute()}")
        print(f"{'='*60}\n")
