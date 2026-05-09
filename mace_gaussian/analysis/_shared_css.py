"""Shared CSS for html_report_generator.py and batch_report.py.

Phase 23 D-04: harmonic and anharmonic reports must share identical styling.
Both report generators must import ``build_css()`` from this module -- no inline CSS.
"""


def build_css() -> str:
    """Return the full ``<style>...</style>`` block for all report HTML."""
    return """<style>
    /* ---- Typography & layout base ---- */
    body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        max-width: 95%;
        margin: 0 auto;
        padding: 24px 40px;
        color: #1a1a1a;
        background: #fafafa;
        line-height: 1.5;
    }
    h1, h2, h3 { color: #0173B2; margin-top: 2em; }
    h1 { font-size: 2em; border-bottom: 3px solid #0173B2; padding-bottom: 8px; }
    h2 { font-size: 1.5em; border-bottom: 1px solid #ddd; padding-bottom: 4px; }

    /* ---- Header (from html_report_generator) ---- */
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
        color: white;
        border-bottom: none;
    }
    header .subtitle {
        font-size: 1.1em;
        opacity: 0.85;
        font-weight: 300;
    }

    /* ---- Navigation (sticky) ---- */
    nav {
        background: white;
        padding: 1rem 4rem;
        border-bottom: 1px solid #e5e7eb;
        position: sticky;
        top: 0;
        z-index: 100;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        display: flex;
        flex-wrap: wrap;
        gap: 0.5rem 0;
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

    /* ---- Content wrapper ---- */
    .content {
        padding: 3rem 4rem;
        max-width: 100%;
        margin: 0 auto;
    }
    section { margin-bottom: 4rem; }

    /* ---- Executive summary cards (Phase 23 D-01, D-02) ---- */
    .executive-summary {
        background: #fff;
        padding: 24px;
        border-radius: 8px;
        margin: 24px 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    .verdict {
        font-size: 1.15em;
        color: #0173B2;
        font-weight: 600;
        margin-bottom: 16px;
    }
    .method-cards {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
        gap: 16px;
    }
    .method-card {
        background: #f6f8fa;
        border: 1px solid #e1e4e8;
        border-radius: 6px;
        padding: 16px;
    }
    .method-card.best-method {
        border: 2px solid #DE8F05;
        background: #fffbf2;
    }
    .method-card h3 {
        margin-top: 0;
        color: #1a1a1a;
        font-size: 1.1em;
    }
    .method-card .metric {
        display: flex;
        justify-content: space-between;
        padding: 4px 0;
        font-size: 0.95em;
    }
    .method-card .metric-label { color: #666; }
    .method-card .metric-value {
        font-weight: 600;
        font-family: "SF Mono", Menlo, monospace;
    }

    /* ---- Summary grid (from html_report_generator) ---- */
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

    /* ---- Comparison sections (D-11) ---- */
    .comparison-section {
        background: #fff;
        padding: 24px;
        border-radius: 8px;
        margin: 24px 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    .comparison-section .plot-container { margin: 16px 0; }
    .plotly-graph-div { display: block; }

    /* ---- Plot containers ---- */
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

    /* ---- Tables ---- */
    .data-table {
        border-collapse: collapse;
        width: 100%;
        margin: 12px 0;
        font-size: 0.9em;
    }
    .data-table th, .data-table td {
        text-align: right;
        padding: 8px 12px;
        border-bottom: 1px solid #eee;
    }
    .data-table th {
        background: #f6f8fa;
        font-weight: 600;
    }
    .data-table tr:hover { background: #f6f8fa; }

    /* ---- Generic table styling (from html_report_generator) ---- */
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
    td { color: #4b5563; }
    tbody tr:hover { background: #f9fafb; }

    /* ---- Metric coloring ---- */
    .metric-good { color: #059669; font-weight: 600; }
    .metric-warning { color: #d97706; font-weight: 600; }
    .metric-bad { color: #dc2626; font-weight: 600; }

    /* ---- Stats boxes ---- */
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

    /* ---- Timing + hardware inline (D-09) ---- */
    .timing-block {
        background: #f6f8fa;
        padding: 12px 16px;
        border-radius: 6px;
        margin: 12px 0;
        font-size: 0.9em;
    }
    .timing-block strong { color: #0173B2; }

    /* ---- Degenerate mode notes ---- */
    .degenerate-note {
        background: #fff4e6;
        border-left: 3px solid #DE8F05;
        padding: 8px 12px;
        margin: 8px 0;
        font-size: 0.9em;
    }

    /* ---- Warning boxes ---- */
    .warning-box {
        background: #fef3c7;
        border-left: 3px solid #f59e0b;
        padding: 1rem;
        margin: 1.5rem 0;
        border-radius: 4px;
        color: #92400e;
    }

    /* ---- Batch report specific (from batch_report.py) ---- */
    .container {
        max-width: 1400px;
        margin: 0 auto;
        padding: 24px;
    }
    .subtitle { color: #7f8c8d; margin-bottom: 24px; }
    .plot-section { margin: 24px 0; }
    .plot-section img {
        max-width: 100%;
        height: auto;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    .plot-card {
        background: white;
        padding: 16px;
        margin: 16px 0;
        border-radius: 8px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    .plot-card h3 { margin-bottom: 12px; color: #2c3e50; }
    .plot-card img { max-width: 100%; height: auto; }
    .stats-bar {
        display: flex;
        gap: 24px;
        margin: 16px 0;
        flex-wrap: wrap;
    }
    .stat-box {
        background: white;
        padding: 16px 24px;
        border-radius: 8px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        text-align: center;
    }
    .stat-box .value {
        font-size: 28px;
        font-weight: 700;
        color: #3498db;
    }
    .stat-box .label { font-size: 13px; color: #7f8c8d; }
    tr.best td { background: #d4edda; }
    tr.worst td { background: #f8d7da; }

    /* ---- Category awards ---- */
    .awards-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
        gap: 16px;
        margin: 16px 0;
    }
    .award-card {
        background: #f6f8fa;
        border: 1px solid #e1e4e8;
        border-radius: 6px;
        padding: 16px;
        text-align: center;
    }
    .award-card h3 {
        margin-top: 0;
        color: #0173B2;
        font-size: 1.1em;
    }
    .award-winner-name {
        font-size: 1.2em;
        font-weight: 700;
        color: #1a1a1a;
        margin: 8px 0 4px;
    }
    .award-winner-mae {
        font-size: 0.95em;
        color: #059669;
        font-weight: 600;
        margin-bottom: 12px;
    }
    .award-no-data {
        color: #9ca3af;
        font-style: italic;
        padding: 16px 0;
    }
    .award-table { font-size: 0.85em; }
    .award-table td, .award-table th { text-align: center; }
    tr.award-winner td {
        background: #fffbf2;
        font-weight: 600;
    }

    /* ---- Footer ---- */
    .footer, footer {
        margin-top: 48px;
        padding-top: 16px;
        border-top: 1px solid #ddd;
        color: #666;
        font-size: 0.85em;
        text-align: center;
    }
    .timestamp {
        font-style: italic;
        color: #9ca3af;
        margin-top: 0.5rem;
    }
</style>"""
