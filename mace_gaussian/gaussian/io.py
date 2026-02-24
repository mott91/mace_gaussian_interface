"""Gaussian I/O functions: parse input files, write output files, generate .gjf inputs."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import numpy as np
from ase.data import chemical_symbols

from ..utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV

# Default paths - can be overridden by environment variables
DEFAULT_HELPER_SCRIPT = os.getenv(
    "MACE_HELPER_SCRIPT_PATH",
    str(Path(__file__).parent.parent / "gm_helper.py"),  # Use relative path by default
)


def parse_gaussian_input(infile: str) -> tuple[int, int, int, int, np.ndarray, list]:
    """
    Parse Gaussian external calculation input file.

    Args:
        infile: Path to Gaussian input file

    Returns:
        Tuple of (natoms, deriv, charge, spin, coordinates, atomnames)
        - natoms: Number of atoms
        - deriv: Derivative level (0=energy, 1=gradient, 2=hessian)
        - charge: Molecular charge
        - spin: Spin multiplicity
        - coordinates: Numpy array of shape (natoms, 3) in Angstroms
        - atomnames: List of element symbols
    """
    with open(infile) as f:
        lines = f.readlines()

    # Extract system info from header line
    header = lines[0].split()
    natoms = int(header[0])
    deriv = int(header[1])
    charge = int(header[2])
    spin = int(header[3])

    # Initialize arrays
    coordinates = np.zeros((natoms, 3))
    atomnames = []

    # Parse atomic data (atomic number + coordinates in Bohr)
    for i, line in enumerate(lines[1 : natoms + 1]):
        elements = line.split()

        # Convert atomic number to element symbol
        atomic_num = int(elements[0])
        atomnames.append(chemical_symbols[atomic_num])

        # Convert coordinates from Bohr to Angstrom
        coordinates[i] = BOHR_TO_ANGSTROM * np.array(
            [float(elements[1]), float(elements[2]), float(elements[3])]
        )

    return natoms, deriv, charge, spin, coordinates, atomnames


def write_gaussian_output(
    outfile: str,
    natoms: int,
    energy: float,
    gradient: np.ndarray,
    dipole: np.ndarray,
    dipole_derivatives: np.ndarray,
    hessian: Optional[np.ndarray],
    deriv: int,
):
    """
    Write results to Gaussian external calculation output file.

    Args:
        outfile: Path to output file
        natoms: Number of atoms
        energy: Energy in eV
        gradient: Gradient in eV/Angstrom, shape (natoms, 3)
        dipole: Dipole moment in e*Bohr, shape (3,)
        dipole_derivatives: Dipole derivatives in e, shape (3*natoms, 3)
        hessian: Hessian in Hartree/Bohr^2, shape (3*natoms, 3*natoms), or None
        deriv: Derivative level (0, 1, or 2)
    """
    # Convert energy from eV to Hartree
    energy_hartree = energy / HARTREE_TO_EV

    # Convert gradient from eV/Angstrom to Hartree/Bohr
    gradient_hartree_bohr = gradient * BOHR_TO_ANGSTROM / HARTREE_TO_EV

    # Polarizability (not implemented, set to zero)
    polarizability = np.zeros(6)

    with open(outfile, "w") as f:
        # Write energy and dipole (Fortran format with 'D' exponent)
        line = f"{energy_hartree:20.12E}{dipole[0]:20.12E}{dipole[1]:20.12E}{dipole[2]:20.12E}"
        f.write(line.replace("E", "D") + "\n")

        # Write gradient (forces)
        for i in range(natoms):
            line = f"{gradient_hartree_bohr[i, 0]:20.12E}{gradient_hartree_bohr[i, 1]:20.12E}{gradient_hartree_bohr[i, 2]:20.12E}"
            f.write(line.replace("E", "D") + "\n")

        # Write polarizability (2 lines, 3 components each)
        line = f"{polarizability[0]:20.12E}{polarizability[1]:20.12E}{polarizability[2]:20.12E}"
        f.write(line.replace("E", "D") + "\n")
        line = f"{polarizability[3]:20.12E}{polarizability[4]:20.12E}{polarizability[5]:20.12E}"
        f.write(line.replace("E", "D") + "\n")

        # Write dipole derivatives (3*natoms lines, 3 components each)
        for i in range(3 * natoms):
            line = f"{dipole_derivatives[i, 0]:20.12E}{dipole_derivatives[i, 1]:20.12E}{dipole_derivatives[i, 2]:20.12E}"
            f.write(line.replace("E", "D") + "\n")

        # Write Hessian if second derivatives requested
        if deriv >= 2 and hessian is not None:
            count = 0
            for i in range(3 * natoms):
                for j in range(i + 1):  # Lower triangle including diagonal
                    line = f"{hessian[i, j]:20.12E}"
                    f.write(line.replace("E", "D"))
                    count += 1

                    if count % 3 == 0:  # Three entries per line
                        f.write("\n")

            # Ensure file ends with newline
            if count % 3 != 0:
                f.write("\n")


def ase_to_gjf(
    atoms,
    filename="molecule.gjf",
    route=None,
    title="Gaussian input generated from ASE",
    charge=0,
    multiplicity=1,
):
    # Use default route with configured helper script path if none provided
    if route is None:
        route = f'# freq (anharm)\n# external="{DEFAULT_HELPER_SCRIPT}"'

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()  # Angstrom by ASE convention
    link0 = f"%chk={filename[:-3]}chk\n%mem=2GB\n%NProcShared=2"
    with open(filename, "w") as f:
        f.write(f"{link0}\n")
        f.write(f"{route}\n\n")
        f.write(f"{title}\n\n")
        f.write(f"{charge} {multiplicity}\n")
        for s, pos in zip(symbols, positions):
            f.write(f"{s:2s} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}\n")
        f.write("\n")  # blank line after coordinates
