"""MACE-MDP dipole and polarizability calculator.

MACE-MDP (SPICE-alpha, wB97M-D3(BJ)/def2-TZVPPD) predicts the molecular dipole and the
full polarizability tensor for H, C, N, O, P, S, F, Cl, Br, I. It gives no energies or
forces, so it is a dipole-side model only.

mace 0.3.16 has no ``mace_mdp()`` loader yet (it was merged after that release), but the
model file loads through the generic ``MACECalculator`` with
``model_type="DipolePolarizabilityMACE"``.
"""

from __future__ import annotations

import logging
import os
import urllib.request
from pathlib import Path

import numpy as np
import torch

from ..utils.units import BOHR_TO_ANGSTROM
from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)

MODEL_URL = "https://raw.githubusercontent.com/Nilsgoe/MACE-MDP/main/models/MACE-MDP.model"
MODEL_PATH = Path(os.getenv("MACE_MDP_MODEL_PATH", Path.home() / ".cache/mace/MACE-MDP.model"))
SUPPORTED_ELEMENTS = frozenset({"H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"})

# The model returns the polarizability in e*Angstrom^2/V. Dividing by 4*pi*eps0
# (= 1/14.3996 e/(V*Angstrom)) gives the polarizability volume in Angstrom^3.
E_ANGSTROM2_PER_V_TO_ANGSTROM3 = 14.399645


class MACEMDPDipoleCalculator(DipoleCalculatorBase):
    """Molecular dipole, dipole derivatives and polarizability from MACE-MDP."""

    def __init__(self, device: str | None = None):
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.calc = None
        super().__init__("mace_mdp")

    def _check_availability(self):
        try:
            from mace.calculators import MACECalculator  # noqa: F401

            self.available = True
            logger.info("✓ MACE-MDP dipole calculator available")
        except Exception as e:  # L4
            self.available = False
            logger.warning(f"✗ MACE-MDP dipole calculator failed: {e}")

    def _ensure_calculator(self):
        if self.calc is None:
            from mace.calculators import MACECalculator

            if not MODEL_PATH.exists():
                MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
                logger.info("Downloading MACE-MDP model to %s", MODEL_PATH)
                urllib.request.urlretrieve(MODEL_URL, MODEL_PATH)
            self.calc = MACECalculator(
                model_paths=str(MODEL_PATH),
                model_type="DipolePolarizabilityMACE",
                device=self.device,
                default_dtype="float64",
            )
            logger.info("Loaded MACE-MDP on %s", self.device)

    def _prepared(self, atoms):
        if not self.available:
            raise RuntimeError("MACE-MDP dipole calculator not available")
        unsupported = sorted(set(atoms.get_chemical_symbols()) - SUPPORTED_ELEMENTS)
        if unsupported:
            raise ValueError(f"MACE-MDP does not support these elements: {unsupported}")
        self._ensure_calculator()
        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc
        return atoms_copy

    def calculate_dipole(self, atoms, **kwargs):
        dipole = np.asarray(self._prepared(atoms).get_dipole_moment()).reshape(3)
        # MACE-MDP outputs e*Angstrom -> convert to e*Bohr (base-class contract)
        return dipole / BOHR_TO_ANGSTROM, None

    def calculate_dipole_derivatives(self, atoms, **kwargs) -> np.ndarray:
        """Dipole derivatives dmu/dr by autograd, shape (3*N_atoms, 3), atomic units e."""
        atoms_copy = self._prepared(atoms)
        dmu_dr, _ = self.calc.get_dielectric_derivatives(atoms_copy)
        # mace returns (3, N, 3) = (dipole component, atom, Cartesian), in e*Angstrom/Angstrom = e.
        # Base-class contract: shape (3*N, 3), atom-major Cartesian-minor; dipole is the column.
        return np.asarray(dmu_dr).reshape(3, len(atoms), 3).transpose(1, 2, 0).reshape(-1, 3)

    def polarizability_angstrom3(self, atoms) -> np.ndarray:
        """Static polarizability in Angstrom^3, shape (3, 3).

        Not named ``calculate_polarizability`` on purpose: the workflow hands that to
        Gaussian, and the campaign runs IR only. Rename it to switch Raman input on.
        """
        atoms_copy = self._prepared(atoms)
        alpha = np.asarray(self.calc.get_property("polarizability", atoms_copy)).reshape(3, 3)
        return alpha * E_ANGSTROM2_PER_V_TO_ANGSTROM3
