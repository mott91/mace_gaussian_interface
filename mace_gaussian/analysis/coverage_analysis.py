"""Frequency range coverage analysis for ML vs DFT spectral comparison.

Bins matched frequency pairs into chemically meaningful regions and computes
per-region error metrics (RMSE, MAE, mean % error, mode count). Aggregates
results across molecules for heatmap and bar chart visualization.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

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
