"""Tests for physical constants and unit conversion utilities.

Verifies CODATA 2018 values and roundtrip conversion accuracy.
"""

import pytest

from mace_gaussian.utils.units import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HARTREE,
    HARTREE_TO_EV,
    angstrom_to_bohr,
    bohr_to_angstrom,
    ev_to_hartree,
    hartree_to_ev,
)


class TestConstants:
    """Verify physical constants match CODATA 2018 recommended values."""

    def test_bohr_to_angstrom_codata_2018(self):
        assert pytest.approx(0.529177210903, rel=1e-10) == BOHR_TO_ANGSTROM

    def test_angstrom_to_bohr_inverse(self):
        assert pytest.approx(1.0 / 0.529177210903, rel=1e-10) == ANGSTROM_TO_BOHR

    def test_hartree_to_ev_codata_2018(self):
        assert pytest.approx(27.211386245988, rel=1e-10) == HARTREE_TO_EV

    def test_ev_to_hartree_inverse(self):
        assert pytest.approx(1.0 / 27.211386245988, rel=1e-10) == EV_TO_HARTREE


class TestConversionFunctions:
    """Verify conversion functions produce correct results."""

    def test_hartree_to_ev(self):
        assert hartree_to_ev(1.0) == pytest.approx(27.211386245988, rel=1e-10)

    def test_ev_to_hartree(self):
        assert ev_to_hartree(27.211386245988) == pytest.approx(1.0, rel=1e-10)

    def test_bohr_to_angstrom(self):
        assert bohr_to_angstrom(1.0) == pytest.approx(0.529177210903, rel=1e-10)

    def test_angstrom_to_bohr(self):
        assert angstrom_to_bohr(1.0) == pytest.approx(1.8897261246257702, rel=1e-10)

    def test_roundtrip_energy(self):
        """Converting Hartree -> eV -> Hartree should recover the original value."""
        original = 13.6
        assert hartree_to_ev(ev_to_hartree(original)) == pytest.approx(original, rel=1e-12)

    def test_roundtrip_distance(self):
        """Converting Bohr -> Angstrom -> Bohr should recover the original value."""
        assert angstrom_to_bohr(bohr_to_angstrom(1.0)) == pytest.approx(1.0, rel=1e-12)
