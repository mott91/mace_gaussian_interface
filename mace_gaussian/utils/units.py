"""Physical constants and unit conversion utilities.

All constants use CODATA 2018 recommended values.
Internal units: eV (energy), Angstrom (distance), cm-1 (frequency).
"""

# --- Physical Constants (CODATA 2018) ---
BOHR_TO_ANGSTROM: float = 0.529177210903
ANGSTROM_TO_BOHR: float = 1.0 / BOHR_TO_ANGSTROM

HARTREE_TO_EV: float = 27.211386245988
EV_TO_HARTREE: float = 1.0 / HARTREE_TO_EV


# --- Convenience Functions ---
def hartree_to_ev(energy_hartree: float) -> float:
    """Convert energy from Hartree to eV."""
    return energy_hartree * HARTREE_TO_EV


def ev_to_hartree(energy_ev: float) -> float:
    """Convert energy from eV to Hartree."""
    return energy_ev * EV_TO_HARTREE


def bohr_to_angstrom(distance_bohr: float) -> float:
    """Convert distance from Bohr to Angstrom."""
    return distance_bohr * BOHR_TO_ANGSTROM


def angstrom_to_bohr(distance_angstrom: float) -> float:
    """Convert distance from Angstrom to Bohr."""
    return distance_angstrom * ANGSTROM_TO_BOHR
