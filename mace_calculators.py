import types
import importlib
import sys
import numpy as np
import torch
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
        
        try:
            # Get dipole moment - this may return various formats
            dipole_moment = atoms_copy.get_dipole_moment()
            
            logger.debug(f"DEBUG: dipole_moment type: {type(dipole_moment)}")
            logger.debug(f"DEBUG: dipole_moment value: {dipole_moment}")
            
            # Handle different return formats
            if isinstance(dipole_moment, tuple):
                # If it's a tuple, take the first element
                logger.debug("Dipole returned as tuple, extracting first element")
                dipole_moment = dipole_moment[0]
            
            # Convert torch tensor to numpy if needed
            if isinstance(dipole_moment, torch.Tensor):
                dipole_moment = dipole_moment.detach().cpu().numpy()
            
            # Ensure it's a 1D numpy array
            dipole_moment = np.atleast_1d(np.array(dipole_moment))
            
            # If it's the wrong shape, try to fix it
            if dipole_moment.shape == (1, 3):
                dipole_moment = dipole_moment.flatten()
            elif dipole_moment.shape != (3,):
                logger.warning(f"Unexpected dipole shape: {dipole_moment.shape}, attempting to reshape")
                dipole_moment = dipole_moment.reshape(-1)[:3]  # Take first 3 elements
            
            logger.debug(f"MACE dipole (final): {dipole_moment}, shape: {dipole_moment.shape}")
            
            # Dipole should already be in correct units (e*Bohr) from MACE
            return dipole_moment, None  # No partial charges for now
            
        except Exception as e:
            logger.error(f"Error in MACE dipole calculation: {e}")
            logger.error("Full traceback:", exc_info=True)
            raise
