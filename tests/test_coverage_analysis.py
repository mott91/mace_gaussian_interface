"""Tests for frequency range coverage analysis."""

import math

import numpy as np
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
