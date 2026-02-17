"""xTB-based dipole calculator."""

import logging

from calculators.base import DipoleCalculatorBase

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
        except ImportError as e:
            self.available = False
            logger.warning(f"\u2717 xTB dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole using xTB"""
        from xtb.ase.calculator import XTB

        atoms_copy = atoms.copy()
        atoms_copy.calc = XTB(method="GFN2-xTB")

        # Get dipole moment (in e*Bohr from xTB)
        dipole_moment = atoms_copy.get_dipole_moment()

        # Get partial charges
        partial_charges = atoms_copy.calc.get_charges(atoms_copy)

        logger.debug(f"xTB dipole: {dipole_moment}")
        return dipole_moment, partial_charges
