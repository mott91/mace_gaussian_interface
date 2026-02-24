"""Factory for managing different dipole calculators."""

from __future__ import annotations

import logging

from calculators.base import DipoleCalculatorBase
from calculators.espaloma import EspalomaDipoleCalculator
from calculators.mace_ml import MACEMLDipoleCalculator
from calculators.xtb import XTBDipoleCalculator

logger = logging.getLogger(__name__)


class DipoleCalculatorFactory:
    """Factory for managing different dipole calculators"""

    def __init__(self):
        self.calculators = {}
        self.preferred_order = ["mace_ml", "espaloma", "xtb"]
        self._register_calculators()

    def _register_calculators(self):
        """Register all available calculators"""
        calculators = [
            EspalomaDipoleCalculator(),
            XTBDipoleCalculator(),
            MACEMLDipoleCalculator(),
        ]

        for calc in calculators:
            self.calculators[calc.name] = calc

    def get_calculator(self, method: str = "auto") -> DipoleCalculatorBase:
        """Get dipole calculator by name or auto-select best available"""

        if method == "auto":
            # Return first available calculator in preferred order
            for method_name in self.preferred_order:
                calc = self.calculators.get(method_name)
                if calc and calc.available:
                    logger.info(f"Auto-selected dipole calculator: {method_name}")
                    return calc

            raise RuntimeError("No dipole calculators available")

        calc = self.calculators.get(method)
        if not calc:
            raise ValueError(f"Unknown dipole calculator: {method}")

        if not calc.available:
            raise RuntimeError(f"Dipole calculator '{method}' not available")

        return calc

    def list_available(self) -> dict[str, bool]:
        """List all calculators and their availability"""
        return {name: calc.available for name, calc in self.calculators.items()}


# Global dipole calculator factory
dipole_factory = DipoleCalculatorFactory()
