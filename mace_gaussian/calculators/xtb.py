"""xTB-based dipole calculator."""

import logging

from ..utils.units import BOHR_TO_ANGSTROM
from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)


class XTBDipoleCalculator(DipoleCalculatorBase):
    """xTB-based dipole calculator"""

    def __init__(self):
        super().__init__("xtb")

    def _check_availability(self):
        try:
            from xtb.ase.calculator import XTB  # noqa: F401

            self.available = True
            logger.info("\u2713 xTB dipole calculator available")
        except Exception as e:  # L4
            self.available = False
            logger.warning(f"\u2717 xTB dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole using xTB"""
        from xtb.ase.calculator import XTB

        atoms_copy = atoms.copy()
        atoms_copy.calc = XTB(method="GFN2-xTB")

        # ASE convention: get_dipole_moment() returns e*Angstrom
        dipole_moment = atoms_copy.get_dipole_moment()

        # Get partial charges
        partial_charges = atoms_copy.calc.get_charges(atoms_copy)

        # Convert e*Angstrom -> e*Bohr (base-class unit contract)
        dipole_moment = dipole_moment / BOHR_TO_ANGSTROM

        logger.debug(f"xTB dipole: {dipole_moment} e*Bohr")
        return dipole_moment, partial_charges
