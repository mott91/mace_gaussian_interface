"""Abstract base class for dipole calculators."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod

import numpy as np

from ..utils.units import BOHR_TO_ANGSTROM

logger = logging.getLogger(__name__)


class DipoleCalculatorBase(ABC):
    """Abstract base class for dipole calculators.

    Unit contract (Gaussian external interface expects atomic units):
    - ``calculate_dipole`` returns the dipole in e*Bohr.
    - ``calculate_dipole_derivatives`` returns d(mu)/dr in e (= e*Bohr/Bohr).
    """

    def __init__(self, name: str):
        self.name = name
        self.available = None
        self._check_availability()

    @abstractmethod
    def _check_availability(self) -> bool:
        """Check if this calculator is available"""
        pass

    @abstractmethod
    def calculate_dipole(self, atoms, **kwargs) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Calculate dipole moment and partial charges
        Returns: (dipole_vector, partial_charges)
        """
        pass

    def calculate_dipole_derivatives(self, atoms, displacement=0.01, **kwargs) -> np.ndarray:
        """Calculate dipole derivatives numerically.

        Central differences of ``calculate_dipole`` (e*Bohr) with respect to
        Cartesian displacements in Angstrom, converted to atomic units (e)
        before returning. ``displacement`` is in Angstrom.
        """
        natoms = len(atoms)
        dipole_derivatives = np.zeros((3 * natoms, 3))
        base_positions = atoms.get_positions().copy()

        try:
            for i in range(natoms):
                for j in range(3):  # x, y, z directions
                    # Positive displacement
                    pos_disp = base_positions.copy()
                    pos_disp[i, j] += displacement
                    atoms_temp = atoms.copy()
                    atoms_temp.set_positions(pos_disp)
                    dipole_pos, _ = self.calculate_dipole(atoms_temp, **kwargs)

                    # Negative displacement
                    pos_disp[i, j] -= 2 * displacement
                    atoms_temp.set_positions(pos_disp)
                    dipole_neg, _ = self.calculate_dipole(atoms_temp, **kwargs)

                    # Central difference derivative
                    dipole_deriv = (dipole_pos - dipole_neg) / (2 * displacement)
                    dipole_derivatives[3 * i + j, :] = dipole_deriv

        except Exception as e:
            logger.warning(f"Dipole derivative calculation failed: {e}")

        finally:
            # Restore original positions
            atoms.set_positions(base_positions)

        # (e*Bohr)/Angstrom -> e*Bohr/Bohr = e (atomic units, as Gaussian expects)
        return dipole_derivatives * BOHR_TO_ANGSTROM
