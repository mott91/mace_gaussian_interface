"""
NIST Experimental Spectrum Fetcher

Fetches, caches, and parses experimental IR spectra from NIST WebBook.
Uses nistchempy for NIST search/download and jcamp for JCAMP-DX parsing.

Design decisions:
- Gas-phase spectra only (D-03)
- Best-effort: never raises, returns None on any failure (D-04)
- Cache raw JDX in comparison_results/{molecule}/experimental/ (D-05, D-06)
- All spectra normalized to [0, 1] absorbance (D-08)
"""

import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ExperimentalSpectrum:
    """Parsed experimental IR spectrum from NIST."""

    wavenumbers: np.ndarray  # cm-1, sorted ascending
    absorbance: np.ndarray  # normalized to [0, 1]
    source: str  # e.g., "NIST WebBook"
    molecule_name: str
    cas_number: str


def fetch_experimental_spectrum(
    molecule_name: str,
    cache_dir: Path,
) -> ExperimentalSpectrum | None:
    """Fetch and cache experimental IR spectrum from NIST WebBook.

    Parameters
    ----------
    molecule_name : str
        Molecule name to search (e.g., 'water', 'methanol')
    cache_dir : Path
        Directory for caching (typically comparison_results/{molecule}).
        JDX file is stored at cache_dir/experimental/nist_ir.jdx

    Returns
    -------
    ExperimentalSpectrum or None
        Parsed spectrum, or None if not found / fetch failed.
        Never raises exceptions -- all errors are logged and return None.
    """
    try:
        cache_path = cache_dir / "experimental" / "nist_ir.jdx"

        # Check cache first (D-06)
        if cache_path.exists():
            logger.info(f"Loading cached experimental spectrum: {cache_path}")
            return _parse_jdx_file(cache_path, molecule_name=molecule_name)

        # Fetch from NIST
        import nistchempy as nist

        logger.info(f"Searching NIST WebBook for '{molecule_name}'...")
        search = nist.run_search(molecule_name, "name")

        if not search.success or search.num_compounds == 0:
            logger.warning(f"No NIST compounds found for '{molecule_name}'")
            return None

        if search.num_compounds > 1:
            logger.info(f"Multiple NIST compounds found for '{molecule_name}', using first match")

        compound = search.compounds[0]
        compound.get_ir_spectra()

        if not compound.ir_specs:
            logger.warning(f"No IR spectra available on NIST for '{molecule_name}'")
            return None

        # Find first gas-phase spectrum (D-03)
        for spec in compound.ir_specs:
            jdx_text = spec.jdx_text

            # Write JDX to temp file for parsing
            with tempfile.NamedTemporaryFile(mode="w", suffix=".jdx", delete=False) as f:
                f.write(jdx_text)
                tmp_path = Path(f.name)

            try:
                import jcamp

                data = jcamp.jcamp_readfile(str(tmp_path))
            finally:
                tmp_path.unlink(missing_ok=True)

            # Check for gas phase
            state = data.get("state", "").lower()
            if "gas" not in state:
                logger.debug(f"Skipping non-gas-phase spectrum (state='{state}')")
                continue

            x = np.asarray(data["x"], dtype=np.float64)
            y = np.asarray(data["y"], dtype=np.float64)
            yunits = data.get("yunits", "").upper()

            # Convert transmittance to absorbance
            if "TRANSMITTANCE" in yunits:
                if np.max(y) > 2:  # percent transmittance
                    y = np.clip(y, 0.01, 100)
                    y = 2.0 - np.log10(y)
                else:  # fractional transmittance [0, 1]
                    y = np.clip(y, 1e-6, 1.0)
                    y = -np.log10(y)

            # Sort by wavenumber ascending
            sort_idx = np.argsort(x)
            x = x[sort_idx]
            y = y[sort_idx]

            # Normalize to [0, 1] (D-08)
            y_max = np.max(y)
            if y_max > 0:
                y = y / y_max

            # Cache raw JDX (D-05)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(jdx_text)
            logger.info(f"Cached experimental spectrum to {cache_path}")

            return ExperimentalSpectrum(
                wavenumbers=x,
                absorbance=y,
                source="NIST WebBook",
                molecule_name=compound.name,
                cas_number=compound.ID,
            )

        logger.warning(f"No gas-phase IR spectrum found on NIST for '{molecule_name}'")
        return None

    except Exception as e:
        logger.error(f"Failed to fetch experimental spectrum for '{molecule_name}': {e}")
        return None


def _parse_jdx_file(
    jdx_path: Path,
    molecule_name: str = "",
) -> ExperimentalSpectrum | None:
    """Parse a cached JCAMP-DX file into an ExperimentalSpectrum.

    Parameters
    ----------
    jdx_path : Path
        Path to the .jdx file
    molecule_name : str
        Fallback molecule name if not in JDX metadata

    Returns
    -------
    ExperimentalSpectrum or None
        Parsed spectrum, or None on parse failure.
    """
    try:
        import jcamp

        data = jcamp.jcamp_readfile(str(jdx_path))

        x = np.asarray(data["x"], dtype=np.float64)
        y = np.asarray(data["y"], dtype=np.float64)
        yunits = data.get("yunits", "").upper()

        # Convert transmittance to absorbance
        if "TRANSMITTANCE" in yunits:
            if np.max(y) > 2:  # percent transmittance
                y = np.clip(y, 0.01, 100)
                y = 2.0 - np.log10(y)
            else:  # fractional transmittance [0, 1]
                y = np.clip(y, 1e-6, 1.0)
                y = -np.log10(y)

        # Sort by wavenumber ascending
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]

        # Normalize to [0, 1] (D-08)
        y_max = np.max(y)
        if y_max > 0:
            y = y / y_max

        # Extract metadata from JDX
        title = data.get("title", molecule_name) or molecule_name
        cas = data.get("cas registry no", "unknown")

        return ExperimentalSpectrum(
            wavenumbers=x,
            absorbance=y,
            source="NIST WebBook",
            molecule_name=title,
            cas_number=cas,
        )

    except Exception as e:
        logger.error(f"Failed to parse JDX file {jdx_path}: {e}")
        return None
