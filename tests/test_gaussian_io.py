"""Tests for Gaussian I/O: write_gaussian_output() polarizability handling.

Covers:
- Default (None) polarizability writes zeros -- backward compatible
- Non-zero polarizability array writes correct values in D-notation
- Polarizability block appears in correct position (after gradient, before dipole derivatives)
- Values appear with D-exponent notation
"""

import tempfile
from pathlib import Path

import numpy as np

from mace_gaussian.gaussian.io import write_gaussian_output

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_and_read(
    natoms: int = 3,
    deriv: int = 1,
    polarizability=None,
) -> list[str]:
    """Write a Gaussian output file and return its lines."""
    energy = -75.0  # eV (will be converted internally)
    gradient = np.zeros((natoms, 3))
    dipole = np.array([0.1, 0.2, 0.3])
    dipole_derivatives = np.zeros((3 * natoms, 3))
    hessian = None

    with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
        tmp_path = f.name

    try:
        write_gaussian_output(
            tmp_path,
            natoms,
            energy,
            gradient,
            dipole,
            dipole_derivatives,
            hessian,
            deriv,
            polarizability=polarizability,
        )
        return Path(tmp_path).read_text().splitlines()
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _parse_d_notation_line(line: str) -> list[float]:
    """Parse a line of D-notation floats (3 values, each 20 chars wide)."""
    values = []
    for i in range(3):
        token = line[i * 20 : (i + 1) * 20].strip()
        values.append(float(token.replace("D", "E")))
    return values


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestWriteGaussianOutputPolarizability:
    """Tests for polarizability kwarg in write_gaussian_output."""

    def test_default_none_writes_zeros(self):
        """Calling without polarizability arg writes np.zeros(6)."""
        lines = _write_and_read(polarizability=None)
        natoms = 3
        # Polarizability is 2 lines after gradient block:
        # Line 0: energy+dipole, Lines 1..natoms: gradient, then 2 polarizability lines
        pol_line1 = _parse_d_notation_line(lines[1 + natoms])
        pol_line2 = _parse_d_notation_line(lines[1 + natoms + 1])
        all_pol = pol_line1 + pol_line2
        np.testing.assert_allclose(all_pol, [0.0] * 6, atol=1e-15)

    def test_nonzero_polarizability_written(self):
        """Passing non-zero array writes those exact values."""
        pol_input = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 3.0])
        lines = _write_and_read(polarizability=pol_input)
        natoms = 3
        pol_line1 = _parse_d_notation_line(lines[1 + natoms])
        pol_line2 = _parse_d_notation_line(lines[1 + natoms + 1])
        all_pol = pol_line1 + pol_line2
        np.testing.assert_allclose(all_pol, [1.0, 0.0, 2.0, 0.0, 0.0, 3.0], rtol=1e-10)

    def test_d_notation_used(self):
        """Polarizability values use D-exponent notation (not E)."""
        pol_input = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 3.0])
        lines = _write_and_read(polarizability=pol_input)
        natoms = 3
        pol_line = lines[1 + natoms]
        assert "D" in pol_line
        assert "E" not in pol_line

    def test_polarizability_position_after_gradient_before_derivatives(self):
        """Polarizability block is after gradient and before dipole derivatives."""
        pol_input = np.array([1.5, 0.0, 2.5, 0.0, 0.0, 3.5])
        lines = _write_and_read(polarizability=pol_input)
        natoms = 3

        # Structure: 1 energy line + natoms gradient lines + 2 polarizability lines
        # + 3*natoms dipole derivative lines
        expected_total = 1 + natoms + 2 + 3 * natoms
        assert len(lines) == expected_total

        # Verify polarizability values in correct position
        pol_line1 = _parse_d_notation_line(lines[1 + natoms])
        assert abs(pol_line1[0] - 1.5) < 1e-10  # first polarizability component

    def test_backward_compatible_no_kwarg(self):
        """Calling write_gaussian_output without polarizability kwarg still works."""
        energy = -75.0
        gradient = np.zeros((3, 3))
        dipole = np.array([0.1, 0.2, 0.3])
        dipole_derivatives = np.zeros((9, 3))

        with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as f:
            tmp_path = f.name

        try:
            # Call WITHOUT the polarizability kwarg - must not raise
            write_gaussian_output(
                tmp_path, 3, energy, gradient, dipole, dipole_derivatives, None, 1
            )
            lines = Path(tmp_path).read_text().splitlines()
            assert len(lines) > 0  # file was written
        finally:
            Path(tmp_path).unlink(missing_ok=True)
