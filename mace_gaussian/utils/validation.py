"""Prerequisite checks, input validation, device detection, and version metadata.

Pure functions for validating the runtime environment before calculations begin.
All functions use stdlib or existing project dependencies only.
"""

from __future__ import annotations

import logging
import platform
import shutil
import sys
from pathlib import Path
from typing import Any

from utils.exceptions import InputValidationError, PrerequisiteError

logger = logging.getLogger(__name__)


def check_gaussian_available() -> str:
    """Check that Gaussian 16 (g16) is available on PATH.

    Returns
    -------
    str
        Path to g16 executable.

    Raises
    ------
    PrerequisiteError
        If g16 is not found on PATH.
    """
    path = shutil.which("g16")
    if path is None:
        raise PrerequisiteError(
            "Gaussian 16 (g16) not found on PATH. "
            "Try: module load gaussian"
        )
    return path


def check_formchk_available() -> str:
    """Check that formchk is available on PATH.

    Returns
    -------
    str
        Path to formchk executable.

    Raises
    ------
    PrerequisiteError
        If formchk is not found on PATH.
    """
    path = shutil.which("formchk")
    if path is None:
        raise PrerequisiteError(
            "formchk not found on PATH. "
            "It is typically installed alongside Gaussian 16."
        )
    return path


def check_dipole_model(model_path: str) -> Path:
    """Check that the dipole model file exists.

    Parameters
    ----------
    model_path : str
        Path to the dipole model file.

    Returns
    -------
    Path
        Validated path to the model file.

    Raises
    ------
    PrerequisiteError
        If the model file does not exist.
    """
    p = Path(model_path)
    if not p.exists():
        raise PrerequisiteError(
            f"Dipole model not found: {model_path}. "
            "Set MACE_DIPOLE_MODEL_PATH environment variable or check the path."
        )
    return p


def check_helper_script(script_path: str) -> Path:
    """Check that the Gaussian helper script exists.

    Parameters
    ----------
    script_path : str
        Path to the helper script.

    Returns
    -------
    Path
        Validated path to the script.

    Raises
    ------
    PrerequisiteError
        If the script does not exist.
    """
    p = Path(script_path)
    if not p.exists():
        raise PrerequisiteError(
            f"Helper script not found: {script_path}. "
            "Set MACE_HELPER_SCRIPT_PATH environment variable or check the path."
        )
    return p


def validate_xyz_file(file_path: str) -> dict:
    """Validate an XYZ molecular structure file.

    Parameters
    ----------
    file_path : str
        Path to the XYZ file.

    Returns
    -------
    dict
        Dictionary with 'path', 'n_atoms', and 'symbols' keys.

    Raises
    ------
    InputValidationError
        If the file is missing, has wrong extension, or cannot be parsed.
    """
    p = Path(file_path)

    if not p.exists():
        raise InputValidationError(f"XYZ file not found: {file_path}")

    if p.suffix.lower() != ".xyz":
        raise InputValidationError(
            f"Expected .xyz extension, got '{p.suffix}': {file_path}"
        )

    try:
        from ase.io import read

        atoms = read(str(p))
    except Exception as exc:
        raise InputValidationError(
            f"Failed to parse XYZ file: {file_path}: {exc}"
        ) from exc

    n_atoms = len(atoms)
    if n_atoms > 200:
        logger.warning(
            "Large molecule (%d atoms) in %s -- calculation may be slow",
            n_atoms,
            file_path,
        )

    return {
        "path": str(p.resolve()),
        "n_atoms": n_atoms,
        "symbols": list(atoms.get_chemical_symbols()),
    }


def validate_all_prerequisites(
    check_gaussian: bool = True,
    check_formchk_tool: bool = True,
    dipole_model_path: str | None = None,
    helper_script_path: str | None = None,
) -> dict:
    """Run all prerequisite checks and return a summary.

    Parameters
    ----------
    check_gaussian : bool
        Whether to check for g16.
    check_formchk_tool : bool
        Whether to check for formchk.
    dipole_model_path : str, optional
        Path to dipole model file.
    helper_script_path : str, optional
        Path to helper script.

    Returns
    -------
    dict
        Summary of check results.

    Raises
    ------
    PrerequisiteError
        On the first failing check.
    """
    results: dict[str, Any] = {}

    if check_gaussian:
        results["g16_path"] = check_gaussian_available()
    if check_formchk_tool:
        results["formchk_path"] = check_formchk_available()
    if dipole_model_path is not None:
        results["dipole_model"] = str(check_dipole_model(dipole_model_path))
    if helper_script_path is not None:
        results["helper_script"] = str(check_helper_script(helper_script_path))

    return results


def detect_device() -> str:
    """Detect whether CUDA is available and return the device string.

    Returns
    -------
    str
        'cuda' if GPU is available, 'cpu' otherwise.
    """
    try:
        import torch

        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            logger.info("CUDA available: %s", device_name)
            return "cuda"
        else:
            logger.warning("CUDA not available, falling back to CPU")
            return "cpu"
    except ImportError:
        logger.warning("PyTorch not installed, falling back to CPU")
        return "cpu"


def collect_version_metadata() -> dict:
    """Collect version information for reproducibility metadata.

    Returns
    -------
    dict
        JSON-serializable dictionary with version info including
        tool_version, python_version, platform, and package versions.
    """
    metadata: dict[str, Any] = {
        "tool_version": "0.2.0",
        "python_version": sys.version,
        "platform": platform.platform(),
    }

    # Try to get version from installed package
    try:
        from importlib.metadata import version

        metadata["tool_version"] = version("mace-gaussian-interface")
    except Exception:
        pass  # Keep fallback "0.2.0"

    # Collect package versions
    try:
        import torch

        metadata["torch_version"] = torch.__version__
        metadata["cuda_available"] = torch.cuda.is_available()
        if torch.cuda.is_available():
            metadata["cuda_version"] = torch.version.cuda or "unknown"
            metadata["gpu_name"] = torch.cuda.get_device_name(0)
    except ImportError:
        metadata["torch_version"] = "not installed"
        metadata["cuda_available"] = False

    try:
        import mace

        metadata["mace_version"] = getattr(mace, "__version__", "unknown")
    except ImportError:
        metadata["mace_version"] = "not installed"

    try:
        import ase

        metadata["ase_version"] = ase.__version__
    except ImportError:
        metadata["ase_version"] = "not installed"

    return metadata
