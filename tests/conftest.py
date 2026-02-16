"""Shared test fixtures and configuration for mace-gaussian test suite."""

import pytest
from pathlib import Path

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
    return str(FIXTURES_DIR / "water" / "dft_b3lyp.fchk")


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
