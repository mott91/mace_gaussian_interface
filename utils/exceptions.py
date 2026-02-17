"""
Custom exception hierarchy for MACE-Gaussian interface.

Provides specific error types for different failure modes, enabling
structured error handling and informative error messages throughout
the pipeline.
"""


class MaceGaussianError(Exception):
    """Base exception for all MACE-Gaussian interface errors.

    Use this as a catch-all when you want to handle any error
    originating from the MACE-Gaussian pipeline.
    """


class PrerequisiteError(MaceGaussianError):
    """A required external tool or file is missing.

    Raised when a prerequisite check fails, e.g., Gaussian 16 (g16)
    is not on PATH, the dipole model file does not exist, or formchk
    is unavailable.
    """


class GaussianParseError(MaceGaussianError):
    """Failed to parse expected data from Gaussian output.

    Raised when a Gaussian .log or .fchk file does not contain the
    expected sections, has malformed data, or yields unexpected
    structure (e.g., missing frequency table, truncated output).
    """


class InputValidationError(MaceGaussianError):
    """Input file or parameter failed validation.

    Raised when user-supplied input is invalid, e.g., an XYZ file
    that does not exist, has the wrong extension, cannot be parsed
    by ASE, or contains out-of-range parameters.
    """


class CUDANotAvailableWarning(UserWarning):
    """CUDA was requested but is not available.

    Issued as a warning (not an exception) when the user expects
    GPU acceleration but torch.cuda.is_available() returns False.
    The system will fall back to CPU.
    """
