"""Abstract base class for dipole calculators."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


class DipoleCalculatorBase(ABC):
    """Abstract base class for dipole calculators"""

    def __init__(self, name: str):
        self.name = name
        self.available = None
        self._check_availability()

    @abstractmethod
    def _check_availability(self) -> bool:
        """Check if this calculator is available"""
        pass

    @abstractmethod
    def calculate_dipole(self, atoms, **kwargs) -> tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Calculate dipole moment and partial charges
        Returns: (dipole_vector, partial_charges)
        """
        pass

    def calculate_dipole_derivatives(self, atoms, displacement=0.01, **kwargs) -> np.ndarray:
        """Calculate dipole derivatives numerically"""
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

        return dipole_derivatives
