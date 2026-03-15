"""Frequency range coverage analysis for ML vs DFT spectral comparison.

Bins matched frequency pairs into chemically meaningful regions and computes
per-region error metrics (RMSE, MAE, mean % error, mode count). Aggregates
results across molecules for heatmap and bar chart visualization.
"""

import base64
import json
import logging
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")  # Non-interactive backend for headless rendering

logger = logging.getLogger(__name__)

# 7 chemically meaningful frequency regions covering 400-4000 cm^-1.
# Each tuple: (low, high, label, description)
FREQUENCY_REGIONS: list[tuple[int, int, str, str]] = [
    (400, 700, "400-700", "Skeletal bending, torsions, heavy-atom stretches"),
    (700, 1000, "700-1000", "C-X stretches, aromatic C-H out-of-plane bends"),
    (1000, 1300, "1000-1300", "C-O/C-C stretches, CH bending"),
    (1300, 1500, "1300-1500", "CH2/CH3 deformations, ring stretches"),
    (1500, 1800, "1500-1800", "C=C, C=O, C=N double bond stretches"),
    (1800, 2800, "1800-2800", "Silent zone / triple bonds"),
    (2800, 4000, "2800-4000", "C-H, N-H, O-H stretches"),
]


class CoverageAnalyzer:
    """Analyze ML vs DFT frequency accuracy across spectral regions.

    Parameters
    ----------
    molecules : list[str]
        Molecule names to analyze.
    results_base : Path or str
        Base directory containing per-molecule harmonic analysis results.
        Default: ``analysis_results_harmonic``.
    """

    def __init__(
        self,
        molecules: list[str],
        results_base: Path | str = Path("analysis_results_harmonic"),
    ) -> None:
        self.molecules = molecules
        self.results_base = Path(results_base)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _bin_frequencies(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add a 'region' column by binning DFT_Frequency_cm into spectral regions.

        Frequencies outside 400-4000 cm^-1 receive NaN (silently excluded).
        Uses right=False so that bin [low, high) includes the left boundary.
        """
        bins = [r[0] for r in FREQUENCY_REGIONS] + [FREQUENCY_REGIONS[-1][1]]
        labels = [r[2] for r in FREQUENCY_REGIONS]
        df = df.copy()
        df["region"] = pd.cut(
            df["DFT_Frequency_cm"],
            bins=bins,
            labels=labels,
            right=False,
        )
        return df

    def _compute_region_metrics(self, df: pd.DataFrame) -> dict:
        """Compute error metrics for a single-region DataFrame.

        Returns
        -------
        dict with keys: rmse, mae, mean_pct_error, mode_count.
        Empty DataFrame returns NaN metrics with mode_count=0.
        """
        n = len(df)
        if n == 0:
            return {
                "rmse": float("nan"),
                "mae": float("nan"),
                "mean_pct_error": float("nan"),
                "mode_count": 0,
            }
        diff = df["ML_Frequency_cm"] - df["DFT_Frequency_cm"]
        rmse = float(np.sqrt(np.mean(diff**2)))
        mae = float(np.mean(np.abs(diff)))
        mean_pct_error = float(np.mean(np.abs(df["Percent_Error"])))
        return {
            "rmse": rmse,
            "mae": mae,
            "mean_pct_error": mean_pct_error,
            "mode_count": int(n),
        }

    def _discover_calc_combos(self, molecule: str) -> list[str]:
        """Discover calculator combo names from comparison CSV files.

        Returns list of combo names (e.g. ``mace_mp_espaloma``).
        If the data directory does not exist, logs a warning and returns [].
        """
        data_dir = self.results_base / molecule / "data"
        if not data_dir.is_dir():
            logger.warning(
                "No harmonic results directory for '%s': %s", molecule, data_dir
            )
            return []
        csvs = sorted(data_dir.glob("comparison_*.csv"))
        combos = []
        for csv_path in csvs:
            name = csv_path.stem  # comparison_mace_mp_espaloma
            combo = name.removeprefix("comparison_")
            combos.append(combo)
        return combos

    def _load_molecule_data(
        self, molecule: str, calc_combo: str
    ) -> pd.DataFrame:
        """Load a comparison CSV and add region binning."""
        csv_path = (
            self.results_base / molecule / "data" / f"comparison_{calc_combo}.csv"
        )
        df = pd.read_csv(csv_path)
        return self._bin_frequencies(df)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def analyze(self) -> dict:
        """Run coverage analysis across all molecules and calculator combos.

        Returns
        -------
        dict
            Nested structure: ``{calc_combo: {molecule: {region_label: metrics_dict}}}``
        """
        results: dict = {}
        region_labels = [r[2] for r in FREQUENCY_REGIONS]

        for molecule in self.molecules:
            combos = self._discover_calc_combos(molecule)
            for combo in combos:
                if combo not in results:
                    results[combo] = {}

                df = self._load_molecule_data(molecule, combo)
                mol_metrics: dict = {}
                for label in region_labels:
                    region_df = df[df["region"] == label]
                    mol_metrics[label] = self._compute_region_metrics(region_df)

                results[combo][molecule] = mol_metrics

        return results

    def build_heatmap_data(
        self, metrics: dict, calc_combo: str
    ) -> pd.DataFrame:
        """Build a molecules x regions DataFrame of RMSE values for heatmap.

        Parameters
        ----------
        metrics : dict
            Output of :meth:`analyze`.
        calc_combo : str
            Calculator combo key to extract.

        Returns
        -------
        pd.DataFrame
            Rows = molecules, columns = region labels, values = RMSE.
            Missing molecule-combo pairs produce NaN rows.
        """
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        combo_data = metrics.get(calc_combo, {})

        rows = {}
        for molecule in self.molecules:
            if molecule in combo_data:
                rows[molecule] = {
                    label: combo_data[molecule][label]["rmse"]
                    for label in region_labels
                }
            else:
                rows[molecule] = {label: float("nan") for label in region_labels}

        return pd.DataFrame.from_dict(rows, orient="index", columns=region_labels)

    def build_barchart_data(
        self, metrics: dict, calc_combo: str
    ) -> tuple[dict, dict, dict]:
        """Compute per-region mean/std/count of RMSE across molecules.

        Parameters
        ----------
        metrics : dict
            Output of :meth:`analyze`.
        calc_combo : str
            Calculator combo key to extract.

        Returns
        -------
        tuple[dict, dict, dict]
            (region_means, region_stds, region_counts) — each keyed by region label.
            NaN RMSE values are excluded from mean/std computation.
        """
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        combo_data = metrics.get(calc_combo, {})

        region_means: dict[str, float] = {}
        region_stds: dict[str, float] = {}
        region_counts: dict[str, int] = {}

        for label in region_labels:
            values = []
            for molecule in combo_data:
                rmse = combo_data[molecule][label]["rmse"]
                if not np.isnan(rmse):
                    values.append(rmse)
            if values:
                region_means[label] = float(np.mean(values))
                region_stds[label] = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
                region_counts[label] = len(values)
            else:
                region_means[label] = float("nan")
                region_stds[label] = float("nan")
                region_counts[label] = 0

        return region_means, region_stds, region_counts

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------

    def plot_heatmap(
        self,
        heatmap_df: pd.DataFrame,
        calc_combo: str,
        save_path: Path,
    ) -> Path:
        """Render a molecules x regions RMSE heatmap to PNG.

        Empty cells (NaN) are masked and annotated with "--" in gray.
        """
        save_path = Path(save_path)
        fig, ax = plt.subplots(
            figsize=(12, max(4, len(heatmap_df) * 0.6 + 2))
        )
        mask = heatmap_df.isna()
        sns.heatmap(
            heatmap_df,
            annot=True,
            fmt=".1f",
            cmap="YlOrRd",
            mask=mask,
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": r"RMSE (cm$^{-1}$)"},
            ax=ax,
        )
        # Annotate masked (empty) cells with "--"
        for i in range(heatmap_df.shape[0]):
            for j in range(heatmap_df.shape[1]):
                if mask.iloc[i, j]:
                    ax.text(
                        j + 0.5,
                        i + 0.5,
                        "--",
                        ha="center",
                        va="center",
                        color="gray",
                        fontsize=10,
                    )
        ax.set_title(f"Frequency Region Coverage: {calc_combo}", fontsize=13)
        ax.set_xlabel("Frequency Region (cm$^{-1}$)")
        ax.set_ylabel("Molecule")
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return save_path

    def plot_barchart(
        self,
        means: dict,
        stds: dict,
        counts: dict,
        calc_combo: str,
        save_path: Path,
    ) -> Path:
        """Render a mean-RMSE-by-region bar chart to PNG."""
        save_path = Path(save_path)
        region_labels = list(means.keys())
        vals = [means[k] for k in region_labels]
        errs = [stds[k] if not math.isnan(stds.get(k, float("nan"))) else 0 for k in region_labels]
        cnts = [counts.get(k, 0) for k in region_labels]

        colors = sns.color_palette("colorblind", len(region_labels))
        fig, ax = plt.subplots(figsize=(10, 5))
        x = range(len(region_labels))
        bars = ax.bar(x, vals, yerr=errs, capsize=4, color=colors, edgecolor="white")

        # Annotate mode counts above each bar
        for idx, bar in enumerate(bars):
            if not math.isnan(vals[idx]):
                y_top = vals[idx] + (errs[idx] if errs[idx] else 0)
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    y_top + 0.5,
                    f"n={cnts[idx]}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        ax.set_xticks(list(x))
        ax.set_xticklabels(region_labels, rotation=45, ha="right")
        ax.set_ylabel(r"RMSE (cm$^{-1}$)")
        ax.set_title(f"Mean RMSE by Frequency Region: {calc_combo}", fontsize=13)
        fig.tight_layout()
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return save_path

    # ------------------------------------------------------------------
    # JSON output
    # ------------------------------------------------------------------

    def save_metrics_json(self, metrics: dict, save_path: Path) -> Path:
        """Serialize metrics to JSON, converting NaN to null."""
        save_path = Path(save_path)

        def _nan_to_none(obj):
            if isinstance(obj, float) and math.isnan(obj):
                return None
            if isinstance(obj, dict):
                return {k: _nan_to_none(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [_nan_to_none(v) for v in obj]
            return obj

        clean = _nan_to_none(metrics)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        save_path.write_text(json.dumps(clean, indent=2))
        return save_path

    # ------------------------------------------------------------------
    # HTML report
    # ------------------------------------------------------------------

    def generate_html_report(
        self,
        output_dir: Path,
        metrics: dict,
        plot_paths: dict,
    ) -> Path:
        """Generate a self-contained HTML report with embedded images.

        Parameters
        ----------
        output_dir : Path
            Directory to write ``coverage_report.html``.
        metrics : dict
            Output of :meth:`analyze`.
        plot_paths : dict
            ``{calc_combo: {"heatmap": Path, "barchart": Path}}``.
        """
        output_dir = Path(output_dir)
        html_path = output_dir / "coverage_report.html"

        sections = []
        for combo in sorted(metrics.keys()):
            paths = plot_paths.get(combo, {})
            heatmap_b64 = self._encode_image(paths["heatmap"]) if "heatmap" in paths else ""
            barchart_b64 = self._encode_image(paths["barchart"]) if "barchart" in paths else ""

            # Summary table rows
            region_labels = [r[2] for r in FREQUENCY_REGIONS]
            combo_data = metrics[combo]
            table_rows = ""
            for label in region_labels:
                rmse_vals, mae_vals, pct_vals, total_modes = [], [], [], 0
                for mol in combo_data:
                    m = combo_data[mol][label]
                    if m["mode_count"] > 0:
                        rmse_vals.append(m["rmse"])
                        mae_vals.append(m["mae"])
                        pct_vals.append(m["mean_pct_error"])
                    total_modes += m["mode_count"]
                mean_rmse = f"{np.mean(rmse_vals):.1f}" if rmse_vals else "--"
                mean_mae = f"{np.mean(mae_vals):.1f}" if mae_vals else "--"
                mean_pct = f"{np.mean(pct_vals):.1f}%" if pct_vals else "--"
                table_rows += (
                    f"<tr><td>{label}</td><td>{mean_rmse}</td>"
                    f"<td>{mean_mae}</td><td>{mean_pct}</td>"
                    f"<td>{total_modes}</td></tr>\n"
                )

            heatmap_img = (
                f'<div class="plot-container"><img src="{heatmap_b64}" alt="Heatmap"></div>'
                if heatmap_b64
                else ""
            )
            barchart_img = (
                f'<div class="plot-container"><img src="{barchart_b64}" alt="Bar chart"></div>'
                if barchart_b64
                else ""
            )

            sections.append(f"""
            <section>
                <h2>{combo}</h2>
                {heatmap_img}
                {barchart_img}
                <h3>Summary Table</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Region</th>
                            <th>Mean RMSE</th>
                            <th>Mean MAE</th>
                            <th>Mean % Error</th>
                            <th>Total Modes</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </section>
            """)

        now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Frequency Range Coverage Analysis</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        line-height: 1.6; color: #1a1a1a; background: #f5f5f5;
    }}
    .container {{ max-width: 100%; background: white; min-height: 100vh; }}
    header {{
        background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        color: white; padding: 3rem 4rem; border-bottom: 3px solid #3498db;
    }}
    header h1 {{ font-size: 2.2em; margin-bottom: 0.5rem; font-weight: 600; }}
    header .subtitle {{ font-size: 1.1em; opacity: 0.85; font-weight: 300; }}
    .content {{ padding: 3rem 4rem; max-width: 1600px; margin: 0 auto; }}
    section {{ margin-bottom: 4rem; }}
    h2 {{
        color: #1a1a1a; font-size: 1.75em; margin-bottom: 1.5rem;
        padding-bottom: 0.75rem; border-bottom: 2px solid #e5e7eb; font-weight: 600;
    }}
    h3 {{ color: #374151; font-size: 1.3em; margin: 2rem 0 1rem 0; font-weight: 600; }}
    .plot-container {{
        background: white; padding: 1.5rem; border-radius: 6px;
        border: 1px solid #e5e7eb; margin: 1.5rem 0;
    }}
    .plot-container img {{ max-width: 100%; height: auto; display: block; margin: 0 auto; }}
    table {{
        width: 100%; border-collapse: collapse; margin: 1rem 0;
        font-size: 0.95rem;
    }}
    th, td {{ padding: 0.75rem 1rem; text-align: left; border-bottom: 1px solid #e5e7eb; }}
    th {{ background: #f9fafb; font-weight: 600; color: #374151; }}
    tr:hover {{ background: #f9fafb; }}
    footer {{
        padding: 2rem 4rem; color: #6b7280; font-size: 0.85rem;
        border-top: 1px solid #e5e7eb; text-align: center;
    }}
</style>
</head>
<body>
<div class="container">
    <header>
        <h1>Frequency Range Coverage Analysis</h1>
        <div class="subtitle">ML vs DFT accuracy across spectral regions
        &mdash; Molecules: {', '.join(self.molecules)}</div>
    </header>
    <div class="content">
        {"".join(sections)}
    </div>
    <footer>
        Generated {now} &mdash; MACE-Gaussian Coverage Analysis
    </footer>
</div>
</body>
</html>"""

        html_path.write_text(html)
        return html_path

    @staticmethod
    def _encode_image(image_path: Path) -> str:
        """Encode a PNG to a base64 data URI."""
        with Path(image_path).open("rb") as f:
            encoded = base64.b64encode(f.read()).decode()
        return f"data:image/png;base64,{encoded}"

    # ------------------------------------------------------------------
    # Orchestrator
    # ------------------------------------------------------------------

    def run(self, output_dir: Path = Path("coverage_analysis")) -> Path:
        """Run full coverage analysis: analyze, plot, report.

        Creates per-combo subdirectories with heatmap/barchart PNGs and
        JSON metrics, plus a combined HTML report.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        metrics = self.analyze()
        plot_paths: dict[str, dict[str, Path]] = {}

        for combo in sorted(metrics.keys()):
            combo_dir = output_dir / combo
            combo_dir.mkdir(parents=True, exist_ok=True)

            # Heatmap
            heatmap_df = self.build_heatmap_data(metrics, combo)
            heatmap_path = self.plot_heatmap(
                heatmap_df, combo, combo_dir / f"heatmap_{combo}.png"
            )

            # Bar chart
            means, stds, counts = self.build_barchart_data(metrics, combo)
            barchart_path = self.plot_barchart(
                means, stds, counts, combo, combo_dir / f"barchart_{combo}.png"
            )

            plot_paths[combo] = {"heatmap": heatmap_path, "barchart": barchart_path}

            # Per-combo JSON
            self.save_metrics_json(
                {combo: metrics[combo]}, combo_dir / f"metrics_{combo}.json"
            )

        # Combined JSON
        self.save_metrics_json(metrics, output_dir / "all_metrics.json")

        # HTML report
        self.generate_html_report(output_dir, metrics, plot_paths)

        # Summary to stdout
        print("Coverage analysis complete:")
        print(f"  Molecules: {', '.join(self.molecules)}")
        print(f"  Calculator combos: {', '.join(sorted(metrics.keys()))}")
        print(f"  Output: {output_dir}/")

        return output_dir
