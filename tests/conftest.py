"""Shared test fixtures and configuration for mace-gaussian test suite."""

from pathlib import Path

import numpy as np
import pytest

from mace_gaussian.analysis.analyze_spectra import ComparisonMetrics, SpectrumData
from mace_gaussian.analysis.nist_fetcher import ExperimentalSpectrum

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir():
    """Return path to test fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def water_dft_log():
    """Path to water DFT B3LYP/6-31G(d,p) Gaussian log file."""
    return str(FIXTURES_DIR / "water" / "dft_b3lyp.log")


@pytest.fixture
def water_ml_log():
    """Path to water MACE-MP/Espaloma ML Gaussian log file."""
    return str(FIXTURES_DIR / "water" / "ml_mace_mp_esp.log")


@pytest.fixture
def water_dft_fchk():
    """Path to water DFT B3LYP/6-31G(d,p) formatted checkpoint file."""
    return FIXTURES_DIR / "water" / "dft_b3lyp.fchk"


@pytest.fixture
def water_ml_fchk():
    """Path to water MACE-MP/Espaloma ML formatted checkpoint file."""
    return str(FIXTURES_DIR / "water" / "ml_mace_mp_esp.fchk")


@pytest.fixture
def ch4_dft_log():
    """Path to CH4 DFT B3LYP/6-31G(d,p) Gaussian log file."""
    return str(FIXTURES_DIR / "CH4_ase" / "dft_b3lyp.log")


@pytest.fixture
def ch4_ml_log():
    """Path to CH4 MACE-MP/Espaloma ML Gaussian log file."""
    return str(FIXTURES_DIR / "CH4_ase" / "ml_mace_mp_esp.log")


@pytest.fixture
def ch4_dft_fchk():
    """Path to CH4 DFT formatted checkpoint file."""
    return str(FIXTURES_DIR / "CH4_ase" / "dft_b3lyp.fchk")


@pytest.fixture
def ch4_ml_fchk():
    """Path to CH4 MACE-MP/Espaloma ML formatted checkpoint file."""
    return str(FIXTURES_DIR / "CH4_ase" / "ml_mace_mp_esp.fchk")


@pytest.fixture
def acoh_ml_log():
    """Path to acetic acid MACE-MP/Espaloma ML log (demonstrates parsing bug)."""
    return str(FIXTURES_DIR / "acoh" / "ml_mace_mp_esp.log")


@pytest.fixture
def water_results_json():
    """Path to water ML reference results.json."""
    return FIXTURES_DIR / "water" / "results.json"


@pytest.fixture
def ch4_results_json():
    """Path to CH4 ML reference results.json."""
    return FIXTURES_DIR / "CH4_ase" / "results.json"


# Phase 23 fixtures — shared by report generator tests


@pytest.fixture
def fake_metrics() -> ComparisonMetrics:
    """Fake ComparisonMetrics for report tests."""
    return ComparisonMetrics(
        mae_freq=12.3,
        rmse_freq=18.5,
        r2_freq=0.987,
        slope_freq=0.99,
        intercept_freq=5.0,
        mae_intensity=15.2,
        r2_intensity=0.91,
        max_error_freq=42.1,
        num_peaks=9,
        num_matched=9,
        num_dft_only=0,
        num_ml_only=0,
        match_rate=1.0,
        num_intensity_filtered=1,
    )


@pytest.fixture
def fake_spectrum() -> SpectrumData:
    """Fake SpectrumData for report tests."""
    freqs = np.array([1600.0, 3700.0, 3800.0])
    return SpectrumData(
        frequencies=freqs,
        intensities=np.array([700.0, 370.0, 700.0]),
        labels=["fundamental", "fundamental", "fundamental"],
        mode_ids=["F1", "F2", "F3"],
    )


@pytest.fixture
def fake_experimental() -> ExperimentalSpectrum:
    """Fake ExperimentalSpectrum for report tests."""
    wn = np.linspace(400, 4000, 100)
    return ExperimentalSpectrum(
        wavenumbers=wn,
        absorbance=np.exp(-((wn - 3700) ** 2) / 500),
        source="NIST WebBook (fixture)",
        molecule_name="water",
        cas_number="7732-18-5",
    )


@pytest.fixture
def fake_analysis_results(fake_metrics, fake_spectrum, fake_experimental, tmp_path) -> dict:
    """Fake analysis results dict matching the canonical comparison shape."""
    return {
        "molecule": "water",
        "mode": "anharmonic",
        "bandwidth_fwhm": 10.0,
        "output_dir": str(tmp_path),
        "comparisons": [
            {
                "name": "mace_off_espaloma",
                "metrics": fake_metrics,
                "ml_spectrum": fake_spectrum,
                "dft_spectrum": fake_spectrum,
                "ml_runtime": 13.4,
                "dft_runtime": 60.0,
                "speedup": 4.5,
                "ml_gaussian_timing": {"total_elapsed_s": 11.5},
                "dft_gaussian_timing": {"total_elapsed_s": 58.0},
                "ml_hardware": {
                    "cpu": "i7-6800K",
                    "gpu": "RTX 2070 SUPER",
                    "ram_gb": 62.7,
                },
                "dft_hardware": {
                    "cpu": "Xeon Gold 6254",
                    "node": "rune03",
                    "cpus": "48",
                },
                "spectrum_plot": "spectrum_mace_off_espaloma.png",
                "regression_plot": "regression_mace_off_espaloma.png",
                "table_file": "comparison_mace_off_espaloma.csv",
                "mode_mapping": None,
                "deg_result": None,
            }
        ],
        "experimental": fake_experimental,
    }
