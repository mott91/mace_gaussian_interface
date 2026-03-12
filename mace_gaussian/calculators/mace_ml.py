"""MACE ML-based dipole calculator."""

import logging
import os
from pathlib import Path

from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)

# Default path for the MACE dipole model (MACE4IR foundation model, 78 elements)
DEFAULT_MACE_DIPOLE_MODEL = os.getenv(
    "MACE_DIPOLE_MODEL_PATH",
    str(
        Path.home()
        / "mace_gaussian"
        / "mace4ir_models"
        / "pretrained_models"
        / "model_1_dipole.model"
    ),
)


class MACEMLDipoleCalculator(DipoleCalculatorBase):
    """MACE ML-based dipole calculator"""

    def __init__(self, model_path=None):
        # Use provided path, or fall back to configured default
        self.model_path = model_path if model_path is not None else DEFAULT_MACE_DIPOLE_MODEL
        self.mace_calc = None

        super().__init__("mace_ml")

    def _check_availability(self):
        try:
            from .mace_loader import MACEDipoleCalculator

            # Check if model file exists
            if not Path(self.model_path).exists():
                raise FileNotFoundError(f"MACE dipole model not found at: {self.model_path}")

            self.mace_calc = MACEDipoleCalculator(self.model_path)
            self.available = True
            logger.info(f"\u2713 MACE ML dipole calculator available (model: {self.model_path})")
        except (ImportError, FileNotFoundError) as e:
            self.available = False
            logger.warning(f"\u2717 MACE ML dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        if not self.available:
            raise RuntimeError("MACE ML dipole calculator not available")
        return self.mace_calc.calculate_dipole(atoms, **kwargs)

    def calculate_dipole_derivatives(self, atoms, **kwargs):
        if not self.available:
            raise RuntimeError("MACE ML dipole calculator not available")
        return self.mace_calc.calculate_dipole_derivatives(atoms, **kwargs)

    def calculate_polarizability(self, atoms):
        if not self.available:
            raise RuntimeError("MACE ML dipole calculator not available")
        return self.mace_calc.calculate_polarizability(atoms)
