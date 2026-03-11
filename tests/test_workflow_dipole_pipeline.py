"""Tests for the 4-tuple calculate_dipole_properties() return and dalpha_dr threading.

Tests cover:
- calculate_dipole_properties() returns a 4-tuple
  (dipole, derivatives, charges, polarizability_voigt6)
- When dipole_calc has calculate_polarizability, result is converted to Voigt-6 in Bohr^3
- When dipole_calc lacks calculate_polarizability, polarizability_voigt6 is np.zeros(6)
- dalpha_dr is retrieved from dipole_calc._last_dalpha_dr after calculate_dipole_derivatives()
- run_next_calculation() passes polarizability to write_gaussian_output as keyword arg
- run_next_calculation() stores dalpha_dr in mol.info["dalpha_dr"]
- Exception fallback returns 4-tuple with np.zeros(6) for polarizability

All tests are GPU-free (mocked calculators).
"""

from unittest.mock import MagicMock, patch

import numpy as np
from ase import Atoms

from mace_gaussian.utils.units import BOHR_TO_ANGSTROM

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ANGSTROM3_TO_BOHR3 = (1.0 / BOHR_TO_ANGSTROM) ** 3


def _water() -> Atoms:
    """Build a minimal water molecule."""
    return Atoms(
        symbols=["O", "H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    )


def _mock_dipole_calc_with_polarizability(natoms: int = 3):
    """Create a mock dipole calculator WITH calculate_polarizability."""
    calc = MagicMock()
    calc.name = "mock_mace_dipole"
    calc.calculate_dipole.return_value = (np.array([0.1, 0.2, 0.3]), np.ones(natoms))
    calc.calculate_dipole_derivatives.return_value = np.ones((3 * natoms, 3))

    # Polarizability: identity-like (3,3) in Angstrom^3
    polar_3x3 = np.array([[1.0, 0.1, 0.0], [0.1, 2.0, 0.0], [0.0, 0.0, 3.0]])
    calc.calculate_polarizability.return_value = polar_3x3

    # Side-effect: _last_dalpha_dr set after derivatives call
    calc._last_dalpha_dr = np.ones((3 * natoms, 3, 3))
    return calc


def _mock_dipole_calc_without_polarizability(natoms: int = 3):
    """Create a mock dipole calculator WITHOUT calculate_polarizability (e.g. Espaloma)."""
    calc = MagicMock()
    calc.name = "mock_espaloma"
    calc.calculate_dipole.return_value = (np.array([0.1, 0.2, 0.3]), np.ones(natoms))
    calc.calculate_dipole_derivatives.return_value = np.ones((3 * natoms, 3))

    # Remove calculate_polarizability so hasattr returns False
    del calc.calculate_polarizability
    # Remove _last_dalpha_dr
    del calc._last_dalpha_dr
    return calc


# ---------------------------------------------------------------------------
# Tests for calculate_dipole_properties
# ---------------------------------------------------------------------------


class TestCalculateDipoleProperties4Tuple:
    """Tests for the 4-tuple return from calculate_dipole_properties."""

    def test_returns_4_tuple(self):
        """calculate_dipole_properties returns exactly 4 elements."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = _mock_dipole_calc_with_polarizability()

        result = calculate_dipole_properties(atoms, calc, deriv=1, calculate_derivatives=True)
        assert len(result) == 4

    def test_fourth_element_is_polarizability_voigt6(self):
        """Fourth return element is shape (6,) polarizability in Bohr^3."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = _mock_dipole_calc_with_polarizability()

        _dipole, _derivs, _charges, pol = calculate_dipole_properties(
            atoms, calc, deriv=1, calculate_derivatives=True
        )
        assert pol.shape == (6,)

    def test_polarizability_converted_to_bohr3(self):
        """Polarizability values are converted from Angstrom^3 to Bohr^3."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = _mock_dipole_calc_with_polarizability()

        # The mock returns identity-like: [[1,0.1,0],[0.1,2,0],[0,0,3]] in Angstrom^3
        _dipole, _derivs, _charges, pol = calculate_dipole_properties(
            atoms, calc, deriv=1, calculate_derivatives=True
        )

        # Voigt order: [p00, p01, p11, p02, p12, p22]
        expected = np.array([1.0, 0.1, 2.0, 0.0, 0.0, 3.0]) * ANGSTROM3_TO_BOHR3
        np.testing.assert_allclose(pol, expected, rtol=1e-10)

    def test_no_polarizability_attribute_returns_zeros(self):
        """When dipole_calc lacks calculate_polarizability, returns np.zeros(6)."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = _mock_dipole_calc_without_polarizability()

        _dipole, _derivs, _charges, pol = calculate_dipole_properties(
            atoms, calc, deriv=1, calculate_derivatives=True
        )
        np.testing.assert_array_equal(pol, np.zeros(6))

    def test_exception_fallback_returns_4_tuple(self):
        """Exception handler returns 4-tuple with zeros for polarizability."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = MagicMock()
        calc.name = "broken_calc"
        calc.calculate_dipole.side_effect = RuntimeError("boom")

        result = calculate_dipole_properties(atoms, calc, deriv=1, calculate_derivatives=True)
        assert len(result) == 4
        np.testing.assert_array_equal(result[3], np.zeros(6))

    def test_dalpha_dr_retrieved_from_calc(self):
        """After derivatives, _last_dalpha_dr is accessible from the calculator."""
        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = _water()
        calc = _mock_dipole_calc_with_polarizability()

        calculate_dipole_properties(atoms, calc, deriv=1, calculate_derivatives=True)

        # After call, _last_dalpha_dr should still be on the calc
        assert calc._last_dalpha_dr is not None
        assert calc._last_dalpha_dr.shape == (9, 3, 3)


# ---------------------------------------------------------------------------
# Tests for run_next_calculation pipeline threading
# ---------------------------------------------------------------------------


class TestRunNextCalculationPipeline:
    """Tests for polarizability kwarg and dalpha_dr storage in run_next_calculation."""

    def _setup_mocks(self):
        """Set up common mocks for run_next_calculation tests."""
        mol = _water()
        mol.info["charge"] = 0.0
        mol.info["spin"] = 1.0
        mock_calc = MagicMock()

        # Mock energy/forces
        mock_calc.get_potential_energy = MagicMock(return_value=-75.0)
        mock_calc.get_forces = MagicMock(return_value=np.zeros((3, 3)))
        mol.calc = mock_calc

        return mol, mock_calc

    @patch("mace_gaussian.workflow.write_gaussian_output")
    @patch("mace_gaussian.workflow.parse_gaussian_input")
    @patch("mace_gaussian.workflow.dipole_factory")
    def test_polarizability_passed_to_write_gaussian_output(
        self, mock_factory, mock_parse, mock_write
    ):
        """run_next_calculation passes polarizability= kwarg to write_gaussian_output."""
        mol, mock_calc = self._setup_mocks()

        # Setup parse return
        mock_parse.return_value = (3, 1, 0, 1, np.zeros((3, 3)), ["O", "H", "H"])

        # Setup dipole factory
        dipole_calc = _mock_dipole_calc_with_polarizability()
        mock_factory.get_calculator.return_value = dipole_calc

        from mace_gaussian.workflow import run_next_calculation

        run_next_calculation(mol, "infile|outfile", mock_calc, dipole_method="mace_ml")

        # Check write_gaussian_output was called with polarizability keyword
        mock_write.assert_called_once()
        call_kwargs = mock_write.call_args
        # polarizability should be passed as keyword arg
        assert "polarizability" in call_kwargs.kwargs

    @patch("mace_gaussian.workflow.write_gaussian_output")
    @patch("mace_gaussian.workflow.parse_gaussian_input")
    @patch("mace_gaussian.workflow.dipole_factory")
    def test_dalpha_dr_stored_in_mol_info(self, mock_factory, mock_parse, mock_write):
        """run_next_calculation stores dalpha_dr in mol.info['dalpha_dr']."""
        mol, mock_calc = self._setup_mocks()

        mock_parse.return_value = (3, 1, 0, 1, np.zeros((3, 3)), ["O", "H", "H"])

        dipole_calc = _mock_dipole_calc_with_polarizability()
        mock_factory.get_calculator.return_value = dipole_calc

        from mace_gaussian.workflow import run_next_calculation

        run_next_calculation(mol, "infile|outfile", mock_calc, dipole_method="mace_ml")

        assert "dalpha_dr" in mol.info
        assert mol.info["dalpha_dr"].shape == (9, 3, 3)

    @patch("mace_gaussian.workflow.write_gaussian_output")
    @patch("mace_gaussian.workflow.parse_gaussian_input")
    @patch("mace_gaussian.workflow.dipole_factory")
    def test_no_dalpha_dr_when_not_available(self, mock_factory, mock_parse, mock_write):
        """When _last_dalpha_dr is not set, mol.info does not get dalpha_dr."""
        mol, mock_calc = self._setup_mocks()

        mock_parse.return_value = (3, 1, 0, 1, np.zeros((3, 3)), ["O", "H", "H"])

        dipole_calc = _mock_dipole_calc_without_polarizability()
        mock_factory.get_calculator.return_value = dipole_calc

        from mace_gaussian.workflow import run_next_calculation

        run_next_calculation(mol, "infile|outfile", mock_calc, dipole_method="espaloma")

        assert "dalpha_dr" not in mol.info
