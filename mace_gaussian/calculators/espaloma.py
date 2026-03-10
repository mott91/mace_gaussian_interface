"""Espaloma-charge based dipole calculator."""

import logging

import numpy as np

from ..utils.units import BOHR_TO_ANGSTROM
from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)


class EspalomaDipoleCalculator(DipoleCalculatorBase):
    """Espaloma-charge based dipole calculator"""

    def __init__(self):
        super().__init__("espaloma")

    def _check_availability(self):
        try:
            import espaloma_charge
            from rdkit import Chem

            # Test basic functionality
            test_mol = Chem.MolFromSmiles("N")
            espaloma_charge.charge(test_mol)

            self.available = True
            logger.info("\u2713 Espaloma-charge dipole calculator available and tested")
        except ImportError as e:
            self.available = False
            logger.warning(f"\u2717 Espaloma-charge dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole using espaloma partial charges"""
        import espaloma_charge
        import torch
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds

        # Create RDKit molecule from ASE atoms
        mol = Chem.RWMol()
        for symbol in atoms.get_chemical_symbols():
            atom = Chem.Atom(symbol)
            mol.AddAtom(atom)

        # Set coordinates
        conf = Chem.Conformer(len(atoms))
        positions = atoms.get_positions()
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, (float(pos[0]), float(pos[1]), float(pos[2])))
        mol.AddConformer(conf)

        # Determine bonds
        rdDetermineBonds.DetermineBonds(mol, charge=0)

        # Patch espaloma to use float32
        # Temporarily set torch default dtype to float32
        original_dtype = torch.get_default_dtype()
        torch.set_default_dtype(torch.float32)

        try:
            # Get partial charges from espaloma
            charges = espaloma_charge.charge(mol)
        finally:
            # Restore original dtype
            torch.set_default_dtype(original_dtype)

        # Calculate dipole moment: mu = Sum q_i * r_i
        dipole = np.dot(charges, positions)

        # Convert from e*Angstrom to e*Bohr (Gaussian units)
        dipole = dipole / BOHR_TO_ANGSTROM

        logger.debug(f"Espaloma dipole: {dipole}, charges sum: {np.sum(charges):.6f}")
        return dipole, charges
