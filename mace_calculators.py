import types
import importlib
import sys
import numpy as np
from ase.calculators.calculator import Calculator
import logging

logger = logging.getLogger(__name__)

def fake_module_from_real(module_name):
    real_mod = importlib.import_module(module_name)
    fake_mod = types.ModuleType(module_name)
    for attr in dir(real_mod):
        if not attr.startswith("__"):
            setattr(fake_mod, attr, getattr(real_mod, attr))
    return fake_mod

# Backup original MACE module (only once)
original_mace_module = sys.modules.get("mace.modules.models", None)

def load_standard_mace_calculator(*args, **kwargs):
    if original_mace_module:
        sys.modules["mace.modules.models"] = original_mace_module
    else:
        sys.modules["mace.modules.models"] = fake_module_from_real("mace.modules.models")
    importlib.reload(importlib.import_module("mace.calculators.mace"))
    from mace.calculators.mace import MACECalculator
    return MACECalculator(*args, **kwargs)

def load_dipole_mace_calculator(*args, **kwargs):
    sys.modules["mace.modules.models"] = fake_module_from_real("mace_dipole_core.modules.models")
    importlib.reload(importlib.import_module("mace_dipole_core.calculators.mace"))
    from mace_dipole_core.calculators.mace import MACECalculator_dipole
    return MACECalculator_dipole(*args, **kwargs)

def cleanup_mace_modules():
    """Clean up monkey-patched modules"""
    sys.modules.pop("mace.modules.models", None)
    importlib.invalidate_caches()

class MACEDipoleCalculator:
    """Simple wrapper for MACE dipole calculator"""
    
    def __init__(self, model_path, device="cuda"):
        self.model_path = model_path
        self.device = device
        self.calc = None
        
    def _ensure_calculator(self):
        if self.calc is None:
            self.calc = load_dipole_mace_calculator(
                model_paths=self.model_path,
                device=self.device,
                model_type="DipolePolarizabilityMACE",
                default_dtype="float64"
            )
    
    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole moment using MACE dipole model"""
        self._ensure_calculator()
        
        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc
        
        # Get dipole moment (already in correct units from MACE)
        dipole_moment = atoms_copy.get_dipole_moment()
        
        logger.debug(f"MACE dipole: {dipole_moment}")
        return dipole_moment, None  # No partial charges for now
