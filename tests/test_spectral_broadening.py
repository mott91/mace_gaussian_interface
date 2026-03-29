"""Tests for spectral broadening and intensity filtering (Phase 18, Plans 01-02)."""

import numpy as np
import pytest

from mace_gaussian.analysis.analyze_spectra import SpectrumAnalyzer, SpectrumData


class TestIntensityFiltering:
    """Tests for intensity filtering in calculate_metrics (SPEC-02)."""

    def test_intensity_filter_threshold(self):
        """Intensity regression excludes modes with DFT intensity < 0.1 km/mol (SPEC-02, D-03)."""
        analyzer = SpectrumAnalyzer(freq_range=(400, 4000), bandwidth_fwhm=8.0)
        dft = SpectrumData(
            frequencies=np.array([500.0, 1000.0, 1500.0, 2000.0, 3000.0]),
            intensities=np.array([0.05, 50.0, 0.01, 100.0, 75.0]),  # 2 below 0.1
            labels=["fundamental"] * 5,
            mode_ids=["F1", "F2", "F3", "F4", "F5"],
        )
        ml = SpectrumData(
            frequencies=np.array([505.0, 1005.0, 1505.0, 2005.0, 3005.0]),
            intensities=np.array([0.1, 55.0, 0.02, 105.0, 80.0]),
            labels=["fundamental"] * 5,
            mode_ids=["F1", "F2", "F3", "F4", "F5"],
        )
        metrics = analyzer.calculate_metrics(ml, dft)
        assert metrics.num_intensity_filtered == 2
        # R2 should be computed from 3 modes only (F2, F4, F5)
        assert metrics.r2_intensity > 0.0

    def test_frequency_metrics_unfiltered(self):
        """Frequency metrics include ALL modes regardless of intensity (SPEC-02, D-04)."""
        analyzer = SpectrumAnalyzer(freq_range=(400, 4000), bandwidth_fwhm=8.0)
        dft = SpectrumData(
            frequencies=np.array([500.0, 1000.0, 1500.0, 2000.0, 3000.0]),
            intensities=np.array([0.05, 50.0, 0.01, 100.0, 75.0]),
            labels=["fundamental"] * 5,
            mode_ids=["F1", "F2", "F3", "F4", "F5"],
        )
        ml = SpectrumData(
            frequencies=np.array([505.0, 1005.0, 1505.0, 2005.0, 3005.0]),
            intensities=np.array([0.1, 55.0, 0.02, 105.0, 80.0]),
            labels=["fundamental"] * 5,
            mode_ids=["F1", "F2", "F3", "F4", "F5"],
        )
        metrics = analyzer.calculate_metrics(ml, dft)
        # All 5 modes included in frequency regression (D-04: no intensity filtering on freq)
        assert metrics.num_peaks == 5
        assert metrics.num_matched == 5
        # Each frequency pair differs by exactly 5.0 cm-1, so MAE should be 5.0
        assert metrics.mae_freq == pytest.approx(5.0, rel=0.01)

    def test_intensity_filter_all_below(self):
        """When all DFT intensities < 0.1, r2_intensity is 0.0 and no error raised (SPEC-02)."""
        analyzer = SpectrumAnalyzer(freq_range=(400, 4000), bandwidth_fwhm=8.0)
        dft = SpectrumData(
            frequencies=np.array([500.0, 1000.0, 1500.0]),
            intensities=np.array([0.01, 0.05, 0.09]),
            labels=["fundamental"] * 3,
            mode_ids=["F1", "F2", "F3"],
        )
        ml = SpectrumData(
            frequencies=np.array([505.0, 1005.0, 1505.0]),
            intensities=np.array([0.1, 0.2, 0.3]),
            labels=["fundamental"] * 3,
            mode_ids=["F1", "F2", "F3"],
        )
        metrics = analyzer.calculate_metrics(ml, dft)
        assert metrics.r2_intensity == 0.0
        assert metrics.num_intensity_filtered == 3
