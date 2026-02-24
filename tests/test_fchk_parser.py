"""Unit tests for fchk_parser.py.

Tests FCHK section parsing, mode extraction, coordinate conversion, and error handling.
Uses pre-committed .fchk fixture files (water and CH4) -- no Gaussian dependency.

Requirement: TEST-03
"""

import numpy as np
import pytest

from mace_gaussian.gaussian.fchk import extract_modes_from_fchk, parse_fchk_section


class TestExtractModesWater:
    """Test mode extraction from water DFT .fchk file."""

    def test_extract_modes_water(self, water_dft_fchk):
        """Water fchk yields 3 atoms, 3 modes, correct shapes and mass values."""
        modes, frequencies, coords, masses, n_atoms = extract_modes_from_fchk(
            water_dft_fchk, force_harmonic=True
        )

        # Returns tuple of 5 elements
        assert n_atoms == 3  # water: O, H, H

        # Modes: 3 vibrational modes for 3 atoms, xyz displacements
        assert modes.shape == (3, 3, 3)

        # Vib-E2 section has 14 values per mode (freq + thermodynamic properties)
        # Only first n_modes entries are actual frequencies
        assert frequencies.shape[0] > 0
        assert np.all(frequencies[:3] > 0), "First 3 frequencies should be positive"

        # Coordinates: 3 atoms, xyz
        assert coords.shape == (3, 3)

        # Masses: 3 atoms
        assert masses.shape == (3,)

        # Oxygen mass ~15.999, hydrogen ~1.008
        assert masses[0] == pytest.approx(15.9949, abs=0.01), "First atom should be oxygen"
        assert masses[1] == pytest.approx(1.008, abs=0.01), "Second atom should be hydrogen"
        assert masses[2] == pytest.approx(1.008, abs=0.01), "Third atom should be hydrogen"


class TestExtractModesCH4:
    """Test mode extraction from CH4 DFT .fchk file."""

    def test_extract_modes_ch4(self, ch4_dft_fchk):
        """CH4 fchk yields 5 atoms, 9 modes (3*5-6), correct shapes."""
        modes, frequencies, coords, masses, n_atoms = extract_modes_from_fchk(
            ch4_dft_fchk, force_harmonic=True
        )

        assert n_atoms == 5  # CH4: C + 4H
        assert modes.shape[0] == 9  # 3*5-6 = 9 vibrational modes
        assert modes.shape == (9, 5, 3)
        assert np.all(frequencies[:9] > 0), "First 9 frequencies should be positive"


class TestModeNormalization:
    """Test that modes extracted from .fchk files are properly normalized."""

    def test_modes_are_normalized(self, water_dft_fchk):
        """Gaussian .fchk normal modes should be normalized to unit length."""
        modes, _, _, _, _ = extract_modes_from_fchk(water_dft_fchk, force_harmonic=True)

        for i in range(modes.shape[0]):
            norm = np.linalg.norm(modes[i].flatten())
            assert norm == pytest.approx(1.0, abs=0.01), f"Mode {i} norm is {norm}, expected ~1.0"


class TestParseFchkSection:
    """Test low-level section parsing from .fchk content."""

    def test_parse_fchk_section_real(self, water_dft_fchk):
        """Parsing 'Real atomic weights' yields correct masses for water."""
        with open(water_dft_fchk) as f:
            content = f.read()

        masses = parse_fchk_section(content, "Real atomic weights", "R")
        assert len(masses) == 3
        assert masses[0] == pytest.approx(15.9949, abs=0.01)  # Oxygen
        assert masses[1] == pytest.approx(1.008, abs=0.01)  # Hydrogen
        assert masses[2] == pytest.approx(1.008, abs=0.01)  # Hydrogen

    def test_parse_fchk_section_integer(self, water_dft_fchk):
        """Parsing 'Atomic numbers' yields [8, 1, 1] for water (O, H, H)."""
        with open(water_dft_fchk) as f:
            content = f.read()

        atomic_numbers = parse_fchk_section(content, "Atomic numbers", "I")
        np.testing.assert_array_equal(atomic_numbers, [8, 1, 1])

    def test_parse_fchk_section_missing(self, water_dft_fchk):
        """Requesting a nonexistent section raises ValueError."""
        with open(water_dft_fchk) as f:
            content = f.read()

        with pytest.raises(ValueError, match="not found"):
            parse_fchk_section(content, "Nonexistent Section", "R")


class TestErrorHandling:
    """Test error paths in fchk_parser."""

    def test_extract_modes_file_not_found(self):
        """Requesting a nonexistent .fchk file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            extract_modes_from_fchk("nonexistent.fchk")


class TestForceHarmonicFlag:
    """Test the force_harmonic parameter behavior."""

    def test_force_harmonic_flag(self, water_dft_fchk):
        """force_harmonic=True and False both succeed; if anharmonic data exists, results differ."""
        modes_h, freqs_h, _, _, n_h = extract_modes_from_fchk(water_dft_fchk, force_harmonic=True)
        modes_d, freqs_d, _, _, n_d = extract_modes_from_fchk(water_dft_fchk, force_harmonic=False)

        # Both calls must succeed
        assert n_h == 3
        assert n_d == 3

        # Water DFT fchk has anharmonic sections, so default (False) should use them.
        # Anharmonic has 9 modes (fund + overtones + combos) vs 3 harmonic modes.
        assert modes_h.shape[0] == 3, "Harmonic should have 3 fundamental modes"
        assert modes_d.shape[0] == 9, "Anharmonic should have 9 modes (fund+overtone+combo)"

        # Frequency arrays should differ in shape since mode counts differ
        assert freqs_h.shape != freqs_d.shape


class TestCoordinateConversion:
    """Test that coordinates are properly converted from Bohr to Angstrom."""

    def test_coordinates_bohr_to_angstrom(self, water_dft_fchk):
        """O-H bond distance should be ~0.96 Angstrom, not ~1.8 Bohr."""
        _, _, coords, _, _ = extract_modes_from_fchk(water_dft_fchk, force_harmonic=True)

        # Compute O-H distance (first atom is O, second is H)
        oh_distance = np.linalg.norm(coords[0] - coords[1])

        # O-H bond in water is ~0.96 A (DFT-optimized geometry may vary slightly)
        assert oh_distance == pytest.approx(0.96, abs=0.1), (
            f"O-H distance {oh_distance:.4f} A should be ~0.96 A (not Bohr)"
        )

        # Sanity check: Bohr distances would be ~1.8, so must be < 1.2
        assert oh_distance < 1.2, "Distance looks like Bohr, not Angstrom"
