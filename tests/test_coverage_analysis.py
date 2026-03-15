"""Tests for frequency range coverage analysis."""

import json
import math
from pathlib import Path

import pandas as pd
import pytest

from mace_gaussian.analysis.coverage_analysis import (
    FREQUENCY_REGIONS,
    CoverageAnalyzer,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def analyzer():
    """CoverageAnalyzer with no molecules (for unit-level helpers)."""
    return CoverageAnalyzer(molecules=[], results_base="/nonexistent")


@pytest.fixture
def sample_df():
    """DataFrame mimicking a comparison CSV with known values."""
    return pd.DataFrame({
        "DFT_Frequency_cm": [500.0, 1600.0, 3500.0],
        "ML_Frequency_cm": [510.0, 1620.0, 3490.0],
        "Freq_Difference_cm": [10.0, 20.0, -10.0],
        "Percent_Error": [2.0, 1.25, -0.286],
        "Mode_Overlap": [0.99, 0.98, 0.97],
        "DFT_Intensity_km_mol": [10.0, 20.0, 30.0],
        "ML_Intensity_km_mol": [11.0, 22.0, 28.0],
        "Intensity_Difference": [1.0, 2.0, -2.0],
    })


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

class TestFrequencyRegions:
    def test_all_seven_regions_defined(self):
        assert len(FREQUENCY_REGIONS) == 7

    def test_regions_cover_400_to_4000(self):
        lows = [r[0] for r in FREQUENCY_REGIONS]
        highs = [r[1] for r in FREQUENCY_REGIONS]
        assert min(lows) == 400
        assert max(highs) == 4000

    def test_regions_contiguous(self):
        """Each region's high == next region's low."""
        for i in range(len(FREQUENCY_REGIONS) - 1):
            assert FREQUENCY_REGIONS[i][1] == FREQUENCY_REGIONS[i + 1][0]


class TestRegionBinning:
    def test_region_binning(self, analyzer, sample_df):
        result = analyzer._bin_frequencies(sample_df)
        assert "region" in result.columns
        regions = result["region"].tolist()
        assert regions[0] == "400-700"
        assert regions[1] == "1500-1800"
        assert regions[2] == "2800-4000"

    def test_out_of_range(self, analyzer):
        df = pd.DataFrame({
            "DFT_Frequency_cm": [31.0, 500.0, 4500.0],
            "ML_Frequency_cm": [30.0, 510.0, 4490.0],
            "Freq_Difference_cm": [-1.0, 10.0, -10.0],
            "Percent_Error": [-3.2, 2.0, -0.22],
            "Mode_Overlap": [0.9, 0.99, 0.95],
            "DFT_Intensity_km_mol": [1.0, 10.0, 50.0],
            "ML_Intensity_km_mol": [1.1, 11.0, 48.0],
            "Intensity_Difference": [0.1, 1.0, -2.0],
        })
        result = analyzer._bin_frequencies(df)
        # Only 500.0 should have a valid region
        valid = result.dropna(subset=["region"])
        assert len(valid) == 1
        assert valid.iloc[0]["DFT_Frequency_cm"] == 500.0


class TestRegionMetrics:
    def test_region_metrics(self, analyzer):
        df = pd.DataFrame({
            "DFT_Frequency_cm": [100.0, 200.0, 300.0],
            "ML_Frequency_cm": [110.0, 220.0, 290.0],
            "Freq_Difference_cm": [10.0, 20.0, -10.0],
            "Percent_Error": [10.0, 10.0, -3.33],
        })
        metrics = analyzer._compute_region_metrics(df)
        # RMSE = sqrt((100+400+100)/3) = sqrt(200) ~ 14.14
        assert metrics["mode_count"] == 3
        assert abs(metrics["rmse"] - math.sqrt(200)) < 0.01
        # MAE = (10+20+10)/3 ~ 13.33
        assert abs(metrics["mae"] - 40 / 3) < 0.01
        # mean_pct_error = mean(|10|, |10|, |-3.33|) = (10+10+3.33)/3 ~ 7.78
        assert abs(metrics["mean_pct_error"] - (10 + 10 + 3.33) / 3) < 0.1

    def test_empty_region(self, analyzer):
        df = pd.DataFrame({
            "DFT_Frequency_cm": pd.Series([], dtype=float),
            "ML_Frequency_cm": pd.Series([], dtype=float),
            "Freq_Difference_cm": pd.Series([], dtype=float),
            "Percent_Error": pd.Series([], dtype=float),
        })
        metrics = analyzer._compute_region_metrics(df)
        assert metrics["mode_count"] == 0
        assert math.isnan(metrics["rmse"])
        assert math.isnan(metrics["mae"])
        assert math.isnan(metrics["mean_pct_error"])

    def test_percent_error_uses_absolute(self, analyzer):
        df = pd.DataFrame({
            "DFT_Frequency_cm": [100.0, 200.0],
            "ML_Frequency_cm": [110.0, 190.0],
            "Freq_Difference_cm": [10.0, -10.0],
            "Percent_Error": [10.0, -5.0],
        })
        metrics = analyzer._compute_region_metrics(df)
        # mean_pct_error = mean(|10|, |-5|) = 7.5
        assert abs(metrics["mean_pct_error"] - 7.5) < 0.01


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

def _write_csv(path: Path, rows: list[dict]) -> None:
    """Helper: write a comparison CSV from row dicts."""
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


class TestIntegration:
    """Integration tests using synthetic CSV files on disk."""

    @pytest.fixture
    def synth_dir(self, tmp_path):
        """Create a temp directory with synthetic molecule comparison CSVs.

        molecule_a has two calc combos (calc_x, calc_y).
        molecule_b has one calc combo (calc_x).
        """
        # molecule_a / calc_x: modes in 400-700 and 2800-4000
        _write_csv(
            tmp_path / "molecule_a" / "data" / "comparison_calc_x.csv",
            [
                {
                    "DFT_Frequency_cm": 500.0,
                    "ML_Frequency_cm": 520.0,
                    "Freq_Difference_cm": 20.0,
                    "Percent_Error": 4.0,
                    "Mode_Overlap": 0.99,
                    "DFT_Intensity_km_mol": 10.0,
                    "ML_Intensity_km_mol": 11.0,
                    "Intensity_Difference": 1.0,
                },
                {
                    "DFT_Frequency_cm": 3000.0,
                    "ML_Frequency_cm": 3050.0,
                    "Freq_Difference_cm": 50.0,
                    "Percent_Error": 1.67,
                    "Mode_Overlap": 0.95,
                    "DFT_Intensity_km_mol": 100.0,
                    "ML_Intensity_km_mol": 95.0,
                    "Intensity_Difference": -5.0,
                },
            ],
        )
        # molecule_a / calc_y: mode only in 1500-1800
        _write_csv(
            tmp_path / "molecule_a" / "data" / "comparison_calc_y.csv",
            [
                {
                    "DFT_Frequency_cm": 1700.0,
                    "ML_Frequency_cm": 1680.0,
                    "Freq_Difference_cm": -20.0,
                    "Percent_Error": -1.18,
                    "Mode_Overlap": 0.98,
                    "DFT_Intensity_km_mol": 50.0,
                    "ML_Intensity_km_mol": 48.0,
                    "Intensity_Difference": -2.0,
                },
            ],
        )
        # molecule_b / calc_x: mode in 400-700
        _write_csv(
            tmp_path / "molecule_b" / "data" / "comparison_calc_x.csv",
            [
                {
                    "DFT_Frequency_cm": 600.0,
                    "ML_Frequency_cm": 610.0,
                    "Freq_Difference_cm": 10.0,
                    "Percent_Error": 1.67,
                    "Mode_Overlap": 0.97,
                    "DFT_Intensity_km_mol": 20.0,
                    "ML_Intensity_km_mol": 21.0,
                    "Intensity_Difference": 1.0,
                },
            ],
        )
        return tmp_path

    def test_analyze_nested_structure(self, synth_dir):
        analyzer = CoverageAnalyzer(
            molecules=["molecule_a", "molecule_b"],
            results_base=synth_dir,
        )
        metrics = analyzer.analyze()

        # calc_x has both molecules
        assert "calc_x" in metrics
        assert "molecule_a" in metrics["calc_x"]
        assert "molecule_b" in metrics["calc_x"]

        # calc_y only has molecule_a
        assert "calc_y" in metrics
        assert "molecule_a" in metrics["calc_y"]
        assert "molecule_b" not in metrics["calc_y"]

        # Check specific metric: molecule_a calc_x region 400-700 has 1 mode, RMSE=20
        m = metrics["calc_x"]["molecule_a"]["400-700"]
        assert m["mode_count"] == 1
        assert abs(m["rmse"] - 20.0) < 0.01

        # molecule_b calc_x region 400-700 has 1 mode, RMSE=10
        m2 = metrics["calc_x"]["molecule_b"]["400-700"]
        assert m2["mode_count"] == 1
        assert abs(m2["rmse"] - 10.0) < 0.01

    def test_build_heatmap_data(self, synth_dir):
        analyzer = CoverageAnalyzer(
            molecules=["molecule_a", "molecule_b"],
            results_base=synth_dir,
        )
        metrics = analyzer.analyze()
        heatmap = analyzer.build_heatmap_data(metrics, "calc_x")

        assert list(heatmap.index) == ["molecule_a", "molecule_b"]
        assert len(heatmap.columns) == 7

        # molecule_a 400-700 RMSE=20.0
        assert abs(heatmap.loc["molecule_a", "400-700"] - 20.0) < 0.01
        # molecule_b 2800-4000 is NaN (no modes there)
        assert math.isnan(heatmap.loc["molecule_b", "2800-4000"])

    def test_build_barchart_data(self, synth_dir):
        analyzer = CoverageAnalyzer(
            molecules=["molecule_a", "molecule_b"],
            results_base=synth_dir,
        )
        metrics = analyzer.analyze()
        means, stds, counts = analyzer.build_barchart_data(metrics, "calc_x")

        # 400-700: molecule_a RMSE=20, molecule_b RMSE=10 -> mean=15, std~7.07
        assert abs(means["400-700"] - 15.0) < 0.01
        assert abs(stds["400-700"] - math.sqrt(50)) < 0.01
        assert counts["400-700"] == 2

        # 2800-4000: only molecule_a has RMSE=50 -> mean=50, std=0
        assert abs(means["2800-4000"] - 50.0) < 0.01
        assert stds["2800-4000"] == 0.0
        assert counts["2800-4000"] == 1

        # 1000-1300: no molecules have modes -> NaN mean
        assert math.isnan(means["1000-1300"])
        assert counts["1000-1300"] == 0

    def test_missing_molecule_warns_and_skips(self, synth_dir, caplog):
        analyzer = CoverageAnalyzer(
            molecules=["molecule_a", "nonexistent_mol"],
            results_base=synth_dir,
        )
        with caplog.at_level("WARNING"):
            metrics = analyzer.analyze()

        # Should have warned about nonexistent_mol
        assert any("nonexistent_mol" in rec.message for rec in caplog.records)

        # calc_x should still have molecule_a data
        assert "calc_x" in metrics
        assert "molecule_a" in metrics["calc_x"]
        # nonexistent_mol should not appear in any combo
        for combo in metrics:
            assert "nonexistent_mol" not in metrics[combo]


# ---------------------------------------------------------------------------
# Smoke tests for visualization + output
# ---------------------------------------------------------------------------

class TestHeatmapGeneration:
    def test_heatmap_creates_png(self, tmp_path):
        """Synthetic heatmap_df with NaN renders to a non-empty PNG."""
        analyzer = CoverageAnalyzer(molecules=["mol_a", "mol_b"], results_base="/nonexistent")
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        data = {
            "mol_a": {label: (i + 1) * 5.0 for i, label in enumerate(region_labels)},
            "mol_b": {label: float("nan") if i % 2 == 0 else (i + 2) * 3.0
                      for i, label in enumerate(region_labels)},
        }
        heatmap_df = pd.DataFrame.from_dict(data, orient="index", columns=region_labels)
        out = tmp_path / "heatmap.png"
        result = analyzer.plot_heatmap(heatmap_df, "test_combo", out)
        assert result == out
        assert out.exists()
        assert out.stat().st_size > 0

    def test_single_molecule_heatmap(self, tmp_path):
        """Heatmap works with a single molecule row."""
        analyzer = CoverageAnalyzer(molecules=["only"], results_base="/nonexistent")
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        data = {"only": {label: 10.0 for label in region_labels}}
        heatmap_df = pd.DataFrame.from_dict(data, orient="index", columns=region_labels)
        out = tmp_path / "heatmap_single.png"
        analyzer.plot_heatmap(heatmap_df, "combo", out)
        assert out.exists()


class TestBarchartGeneration:
    def test_barchart_creates_png(self, tmp_path):
        """Bar chart with means/stds/counts renders to a non-empty PNG."""
        analyzer = CoverageAnalyzer(molecules=[], results_base="/nonexistent")
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        means = {label: (i + 1) * 8.0 for i, label in enumerate(region_labels)}
        stds = {label: (i + 1) * 1.5 for i, label in enumerate(region_labels)}
        counts = {label: (i + 1) * 2 for i, label in enumerate(region_labels)}
        out = tmp_path / "barchart.png"
        result = analyzer.plot_barchart(means, stds, counts, "test_combo", out)
        assert result == out
        assert out.exists()
        assert out.stat().st_size > 0

    def test_barchart_handles_nan(self, tmp_path):
        """Bar chart handles NaN mean values without crashing."""
        analyzer = CoverageAnalyzer(molecules=[], results_base="/nonexistent")
        region_labels = [r[2] for r in FREQUENCY_REGIONS]
        means = {label: float("nan") if i == 0 else 10.0
                 for i, label in enumerate(region_labels)}
        stds = {label: float("nan") if i == 0 else 2.0
                for i, label in enumerate(region_labels)}
        counts = {label: 0 if i == 0 else 3 for i, label in enumerate(region_labels)}
        out = tmp_path / "barchart_nan.png"
        analyzer.plot_barchart(means, stds, counts, "combo", out)
        assert out.exists()


class TestJsonSerialization:
    def test_nan_converted_to_null(self, tmp_path):
        """NaN values become null in JSON output."""
        analyzer = CoverageAnalyzer(molecules=[], results_base="/nonexistent")
        metrics = {
            "combo": {
                "mol": {
                    "400-700": {"rmse": 10.0, "mae": 8.0, "mean_pct_error": 5.0, "mode_count": 3},
                    "1800-2800": {"rmse": float("nan"), "mae": float("nan"),
                                  "mean_pct_error": float("nan"), "mode_count": 0},
                }
            }
        }
        out = tmp_path / "metrics.json"
        result = analyzer.save_metrics_json(metrics, out)
        assert result == out
        assert out.exists()

        loaded = json.loads(out.read_text())
        assert loaded["combo"]["mol"]["400-700"]["rmse"] == 10.0
        assert loaded["combo"]["mol"]["1800-2800"]["rmse"] is None
        assert loaded["combo"]["mol"]["1800-2800"]["mode_count"] == 0

    def test_json_is_valid(self, tmp_path):
        """Output is valid JSON with indent=2."""
        analyzer = CoverageAnalyzer(molecules=[], results_base="/nonexistent")
        metrics = {"a": {"b": {"c": {"rmse": 1.0, "mode_count": 1}}}}
        out = tmp_path / "test.json"
        analyzer.save_metrics_json(metrics, out)
        text = out.read_text()
        loaded = json.loads(text)
        assert loaded == metrics
        # Check indentation
        assert "  " in text
