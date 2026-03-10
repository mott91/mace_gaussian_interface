#!/usr/bin/env python3
"""
Formatted Checkpoint File (.fchk) Parser

This module lives at gaussian/fchk.py (relocated from top-level fchk_parser.py).
The top-level fchk_parser.py is kept intact for now; callers are updated in plan 05.

Parses Gaussian formatted checkpoint files to extract:
- Vibrational modes (normal mode displacement vectors)
- Atomic coordinates
- Atomic masses
- Frequencies

The .fchk format is much cleaner than parsing .log files.
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Optional

import numpy as np

from ..utils.units import BOHR_TO_ANGSTROM

logger = logging.getLogger(__name__)


def convert_chk_to_fchk(chk_file: str, fchk_file: Optional[str] = None) -> str:
    """
    Convert Gaussian checkpoint file (.chk) to formatted checkpoint (.fchk).

    Uses the formchk utility from Gaussian.

    Parameters
    ----------
    chk_file : str
        Path to .chk file
    fchk_file : str, optional
        Path to output .fchk file. If None, uses same name with .fchk extension

    Returns
    -------
    fchk_file : str
        Path to created .fchk file

    Raises
    ------
    RuntimeError
        If formchk command fails
    """
    chk_path = Path(chk_file)
    if not chk_path.exists():
        raise FileNotFoundError(f"Checkpoint file not found: {chk_file}")

    if fchk_file is None:
        fchk_file = str(chk_path.with_suffix(".fchk"))

    logger.info(f"Converting {chk_file} to {fchk_file}")

    try:
        result = subprocess.run(
            ["formchk", chk_file, fchk_file], capture_output=True, text=True, check=True, timeout=60
        )
        logger.debug(f"formchk output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"formchk failed: {e.stderr}") from e
    except FileNotFoundError as e:
        raise RuntimeError("formchk command not found. Is Gaussian loaded?") from e

    if not Path(fchk_file).exists():
        raise RuntimeError(f"formchk did not create output file: {fchk_file}")

    return fchk_file


def parse_fchk_section(content: str, section_name: str, data_type: str = "R") -> np.ndarray:
    """
    Parse a section from .fchk file.

    .fchk format:
    Section Name                              Type      N=           12
       value1  value2  value3  value4  value5

    Parameters
    ----------
    content : str
        Full .fchk file content
    section_name : str
        Name of section to extract (e.g., 'Current cartesian coordinates')
    data_type : str
        'R' for real numbers, 'I' for integers

    Returns
    -------
    data : np.ndarray
        Extracted values
    """
    # Find the section header
    pattern = rf"^{re.escape(section_name)}\s+{data_type}\s+N=\s*(\d+)"
    match = re.search(pattern, content, re.MULTILINE)

    if not match:
        raise ValueError(f"Section '{section_name}' not found in .fchk file")

    n_values = int(match.group(1))
    start_pos = match.end()

    # Extract the data lines following the header
    # Data comes in lines of up to 5 values (for reals) or more (for integers)
    values = []
    lines = content[start_pos:].split("\n")

    for line in lines:
        if not line.strip():
            continue  # Skip empty lines, don't break
        # Check if this is a new section starting (sections typically start at column 1)
        if re.match(r"^[A-Z]", line):
            break  # New section starting
        try:
            # Parse values from line
            line_values = line.split()
            values.extend([float(v) if data_type == "R" else int(v) for v in line_values])
            # Stop if we've collected enough values
            if len(values) >= n_values:
                break
        except ValueError:
            break  # Can't parse, must be end of section

    if len(values) < n_values:
        raise ValueError(f"Expected {n_values} values for '{section_name}', got {len(values)}")

    return np.array(values[:n_values])


def extract_modes_from_fchk(
    fchk_file: str, force_harmonic: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Extract vibrational mode data from .fchk file.

    Parameters
    ----------
    fchk_file : str
        Path to .fchk file
    force_harmonic : bool, optional
        If True, force extraction of harmonic frequencies/modes only,
        skipping anharmonic sections even if available. Default: False.
        Use this for mode overlap calculations where harmonic eigenvectors
        are needed.

    Returns
    -------
    modes : np.ndarray
        Normal mode displacement vectors, shape (n_modes, n_atoms, 3)
    frequencies : np.ndarray
        Vibrational frequencies in cm^-1, shape (n_modes,)
    coords : np.ndarray
        Atomic coordinates in Angstroms, shape (n_atoms, 3)
    masses : np.ndarray
        Atomic masses in amu, shape (n_atoms,)
    n_atoms : int
        Number of atoms
    """
    fchk_path = Path(fchk_file)
    if not fchk_path.exists():
        raise FileNotFoundError(f"Formatted checkpoint file not found: {fchk_file}")

    with fchk_path.open() as f:
        content = f.read()

    # Extract number of atoms
    natoms_match = re.search(r"^Number of atoms\s+I\s+(\d+)", content, re.MULTILINE)
    if not natoms_match:
        raise ValueError("Could not find 'Number of atoms' in .fchk file")
    n_atoms = int(natoms_match.group(1))

    logger.info(f"Found {n_atoms} atoms in .fchk file")

    # Extract atomic masses
    masses = parse_fchk_section(content, "Real atomic weights", "R")
    if len(masses) != n_atoms:
        raise ValueError(f"Expected {n_atoms} masses, got {len(masses)}")

    # Extract current coordinates (in Bohr)
    coords_bohr = parse_fchk_section(content, "Current cartesian coordinates", "R")
    if len(coords_bohr) != n_atoms * 3:
        raise ValueError(f"Expected {n_atoms * 3} coordinates, got {len(coords_bohr)}")

    # Convert Bohr to Angstroms
    coords = coords_bohr * BOHR_TO_ANGSTROM
    coords = coords.reshape(n_atoms, 3)

    # Extract vibrational frequencies and modes
    # Prefer anharmonic sections if available (these contain actual vibrational modes only)
    # UNLESS force_harmonic is True (for mode overlap calculations)
    frequencies = None
    vib_modes = None

    # Try anharmonic sections first (UNLESS force_harmonic=True)
    if not force_harmonic:
        try:
            frequencies = parse_fchk_section(content, "Anharmonic Vib-E2", "R")
            # Anharmonic E2 section has multiple values per mode, extract first N
            # where N = number of modes (check Anharmonic Number of Normal Modes)
            nmodes_match = re.search(
                r"^Anharmonic Number of Normal Modes\s+I\s+(\d+)", content, re.MULTILINE
            )
            if nmodes_match:
                n_modes = int(nmodes_match.group(1))
                frequencies = frequencies[:n_modes]
            else:
                n_modes = None
            logger.info(f"Found {len(frequencies)} anharmonic vibrational frequencies")
        except ValueError:
            logger.debug("No anharmonic frequencies found, using harmonic")
            n_modes = None

        try:
            vib_modes = parse_fchk_section(content, "Anharmonic Vib-Modes", "R")
            logger.info("Using anharmonic vibrational modes")
        except ValueError:
            logger.debug("No anharmonic modes found, using harmonic")
    else:
        logger.info("force_harmonic=True, skipping anharmonic sections")
        n_modes = None

    # Fall back to harmonic if anharmonic not available
    if frequencies is None:
        try:
            frequencies = parse_fchk_section(content, "Vib-E2", "R")
            logger.info(f"Found {len(frequencies)} harmonic vibrational frequencies")
        except ValueError:
            logger.warning("Could not find frequencies")

    if vib_modes is None:
        try:
            vib_modes = parse_fchk_section(content, "Vib-Modes", "R")
            logger.info("Using harmonic vibrational modes")
        except ValueError:
            raise ValueError(
                "Could not find 'Vib-Modes' or 'Anharmonic Vib-Modes' section in .fchk file"
            ) from None

    # Determine number of modes if not already found
    if n_modes is None:
        n_modes = len(vib_modes) // (n_atoms * 3)
        logger.info(f"Inferred {n_modes} modes from Vib-Modes section")

    expected_values = n_modes * n_atoms * 3
    if len(vib_modes) != expected_values:
        raise ValueError(f"Expected {expected_values} mode values, got {len(vib_modes)}")

    # Reshape modes: stored as all values for mode 1, then mode 2, etc.
    # Each mode has n_atoms * 3 values (x,y,z for each atom)
    modes = vib_modes.reshape(n_modes, n_atoms, 3)

    # If frequencies weren't found, set to None
    if frequencies is None:
        frequencies = np.zeros(n_modes)
        logger.warning("Frequencies not available in .fchk file")

    return modes, frequencies, coords, masses, n_atoms


def get_fchk_from_chk(chk_file: str) -> str:
    """
    Get .fchk file from .chk file, converting if necessary.

    Parameters
    ----------
    chk_file : str
        Path to .chk file

    Returns
    -------
    fchk_file : str
        Path to .fchk file (existing or newly created)
    """
    chk_path = Path(chk_file)
    fchk_path = chk_path.with_suffix(".fchk")

    # If .fchk already exists and is newer than .chk, use it
    if fchk_path.exists() and fchk_path.stat().st_mtime > chk_path.stat().st_mtime:
        logger.info(f"Using existing .fchk file: {fchk_path}")
        return str(fchk_path)

    # Convert .chk to .fchk
    logger.info("Converting .chk to .fchk")
    return convert_chk_to_fchk(str(chk_path))


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) < 2:
        print("Usage: python fchk.py <file.chk or file.fchk>")
        sys.exit(1)

    input_file = sys.argv[1]
    input_path = Path(input_file)

    if input_path.suffix == ".chk":
        fchk_file = get_fchk_from_chk(input_file)
    elif input_path.suffix == ".fchk":
        fchk_file = input_file
    else:
        print(f"Unknown file type: {input_path.suffix}")
        sys.exit(1)

    print(f"\nParsing {fchk_file}...")
    modes, frequencies, coords, masses, n_atoms = extract_modes_from_fchk(fchk_file)

    print(f"\n{'=' * 60}")
    print("FCHK FILE SUMMARY")
    print(f"{'=' * 60}")
    print(f"Number of atoms: {n_atoms}")
    print(f"Number of modes: {modes.shape[0]}")
    print("\nFirst 10 frequencies (cm^-1):")
    for i, freq in enumerate(frequencies[:10]):
        print(f"  Mode {i + 1:3d}: {freq:10.2f} cm^-1")
    print(f"{'=' * 60}")
