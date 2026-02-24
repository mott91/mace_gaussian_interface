"""Utility package for MACE-Gaussian interface."""

from .exceptions import (
    CUDANotAvailableWarning,
    GaussianParseError,
    InputValidationError,
    MaceGaussianError,
    PrerequisiteError,
)
from .results import ResultsManager
from .units import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HARTREE,
    HARTREE_TO_EV,
    angstrom_to_bohr,
    bohr_to_angstrom,
    ev_to_hartree,
    hartree_to_ev,
)

__all__ = [
    "ResultsManager",
    "ANGSTROM_TO_BOHR",
    "BOHR_TO_ANGSTROM",
    "CUDANotAvailableWarning",
    "EV_TO_HARTREE",
    "GaussianParseError",
    "HARTREE_TO_EV",
    "InputValidationError",
    "MaceGaussianError",
    "PrerequisiteError",
    "angstrom_to_bohr",
    "bohr_to_angstrom",
    "ev_to_hartree",
    "hartree_to_ev",
]
