import numpy as np
import pytest

from mace_gaussian.analysis.analyze_spectra import SpectrumAnalyzer, SpectrumData


def test_lorentzian_peak_height():
    analyzer = SpectrumAnalyzer(freq_range=(1000, 2000), bandwidth_fwhm=10.0, freq_step=0.5)
    spectrum = SpectrumData(
        frequencies=np.array([1500.0]),
        intensities=np.array([100.0]),
        labels=["fundamental"],
        mode_ids=["F1"],
    )
    broadened = analyzer.broaden_spectrum(spectrum)
    peak_idx = np.argmin(np.abs(analyzer.freq_grid - 1500.0))
    assert broadened[peak_idx] == pytest.approx(100.0, rel=0.01)


def test_lorentzian_fwhm():
    analyzer = SpectrumAnalyzer(freq_range=(1000, 2000), bandwidth_fwhm=10.0, freq_step=0.5)
    spectrum = SpectrumData(
        frequencies=np.array([1500.0]),
        intensities=np.array([100.0]),
        labels=["fundamental"],
        mode_ids=["F1"],
    )
    broadened = analyzer.broaden_spectrum(spectrum)
    # At half-FWHM distance (gamma = 5.0), value should be half of peak
    hm_idx = np.argmin(np.abs(analyzer.freq_grid - 1505.0))
    assert broadened[hm_idx] == pytest.approx(50.0, rel=0.01)


def test_default_fwhm():
    analyzer = SpectrumAnalyzer()
    assert analyzer.bandwidth_fwhm == 10.0
