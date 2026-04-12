"""Shared test fixtures and configuration for mace-gaussian test suite."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from mace_gaussian.analysis.analyze_spectra import ComparisonMetrics, SpectrumData

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@dataclass
class FakeExperimentalSpectrum:
    """Stand-in for ExperimentalSpectrum until nist_fetcher lands."""

    source: str
    molecule_name: str
    cas_number: str
    wavenumbers: np.ndarray
    absorbance: np.ndarray


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
def fake_metrics_b():
    """A second ComparisonMetrics with slightly worse values."""
    return ComparisonMetrics(
        mae_freq=25.0,
        rmse_freq=30.0,
        r2_freq=0.980,
        slope_freq=0.98,
        intercept_freq=10.0,
        mae_intensity=0.2,
        r2_intensity=0.70,
        max_error_freq=55.0,
        num_peaks=3,
        num_matched=3,
        num_dft_only=0,
        num_ml_only=0,
        match_rate=1.0,
    )


@pytest.fixture
def fake_experimental():
    """A FakeExperimentalSpectrum for water."""
    return FakeExperimentalSpectrum(
        source="NIST",
        molecule_name="water",
        cas_number="7732-18-5",
        wavenumbers=np.linspace(400, 4000, 100),
        absorbance=np.random.default_rng(42).random(100),
    )


@pytest.fixture
def fake_analysis_results(fake_metrics, fake_metrics_b, fake_spectrum, fake_experimental):
    """Full analysis_results dict matching the shape consumed by report_data.export_report_data."""
    return {
        "molecule": "water",
        "mode": "anharmonic",
        "bandwidth_fwhm": 10.0,
        "output_dir": "/tmp/test_output",
        "experimental": fake_experimental,
        "comparisons": [
            {
                "name": "mace_off_espaloma",
                "metrics": fake_metrics,
                "ml_peaks": 3,
                "dft_peaks": 3,
                "ml_runtime": 12.5,
                "dft_runtime": 120.0,
                "speedup": 9.6,
                "spectrum_plot": "spectrum.png",
                "regression_plot": "regression.png",
                "table_file": "table.csv",
                "comparison_df": None,
                "ml_spectrum": fake_spectrum,
                "dft_spectrum": fake_spectrum,
                "mode_mapping": None,
                "deg_result": None,
                "experimental": fake_experimental,
                "experimental_agreement": 0.92,
                "ml_hardware": {"cpu": "Intel i7", "gpu": "RTX 3090", "ram_gb": 32.0},
                "dft_hardware": {"cpu": "Xeon E5", "node": "rune03", "cpus": "16"},
                "ml_gaussian_timing": {"total_elapsed_s": 10.0},
                "dft_gaussian_timing": {"total_elapsed_s": 100.0},
            },
            {
                "name": "mace_anicc_mace_ml",
                "metrics": fake_metrics_b,
                "ml_peaks": 3,
                "dft_peaks": 3,
                "ml_runtime": 15.0,
                "dft_runtime": 120.0,
                "speedup": 8.0,
                "spectrum_plot": "spectrum2.png",
                "regression_plot": "regression2.png",
                "table_file": "table2.csv",
                "comparison_df": None,
                "ml_spectrum": fake_spectrum,
                "dft_spectrum": fake_spectrum,
                "mode_mapping": None,
                "deg_result": None,
                "experimental": fake_experimental,
                "experimental_agreement": 0.85,
                "ml_hardware": {"cpu": "Intel i7", "gpu": "RTX 3090", "ram_gb": 32.0},
                "dft_hardware": {"cpu": "Xeon E5", "node": "rune03", "cpus": "16"},
                "ml_gaussian_timing": {"total_elapsed_s": 13.0},
                "dft_gaussian_timing": {"total_elapsed_s": 100.0},
            },
        ],
        "executive_summary": None,
    }
