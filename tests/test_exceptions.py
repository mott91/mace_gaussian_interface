"""Tests for the exception hierarchy in exceptions.py.

Verifies class relationships, inheritance chains, and message preservation.
"""

import pytest

from utils.exceptions import (
    CUDANotAvailableWarning,
    GaussianParseError,
    InputValidationError,
    MaceGaussianError,
    PrerequisiteError,
)


class TestExceptionHierarchy:
    """Verify exception class relationships."""

    def test_base_is_exception(self):
        """MaceGaussianError inherits from Exception."""
        assert issubclass(MaceGaussianError, Exception)

    @pytest.mark.parametrize(
        "exc_class",
        [PrerequisiteError, GaussianParseError, InputValidationError],
    )
    def test_subclass_of_base(self, exc_class):
        """All domain exceptions are subclasses of MaceGaussianError."""
        assert issubclass(exc_class, MaceGaussianError)

    def test_cuda_warning_is_user_warning(self):
        """CUDANotAvailableWarning inherits from UserWarning, not MaceGaussianError."""
        assert issubclass(CUDANotAvailableWarning, UserWarning)
        assert not issubclass(CUDANotAvailableWarning, MaceGaussianError)


class TestExceptionUsage:
    """Verify exceptions can be raised, caught, and carry messages."""

    def test_raise_and_catch_by_base(self):
        """Raising PrerequisiteError should be catchable as MaceGaussianError."""
        with pytest.raises(MaceGaussianError):
            raise PrerequisiteError("g16 not found")

    def test_raise_and_catch_by_exception(self):
        """All domain exceptions are catchable as plain Exception."""
        with pytest.raises(Exception):
            raise GaussianParseError("bad output")

    def test_message_preserved_prerequisite(self):
        """PrerequisiteError preserves its message string."""
        msg = "Gaussian 16 not found on PATH"
        with pytest.raises(PrerequisiteError, match=msg):
            raise PrerequisiteError(msg)

    def test_message_preserved_parse(self):
        """GaussianParseError preserves its message string."""
        msg = "Failed to parse frequency table"
        with pytest.raises(GaussianParseError, match=msg):
            raise GaussianParseError(msg)

    def test_message_preserved_input(self):
        """InputValidationError preserves its message string."""
        msg = "XYZ file not found"
        with pytest.raises(InputValidationError, match=msg):
            raise InputValidationError(msg)

    def test_cuda_warning_can_be_issued(self):
        """CUDANotAvailableWarning can be issued via warnings.warn."""
        import warnings

        with pytest.warns(CUDANotAvailableWarning, match="no GPU"):
            warnings.warn("no GPU", CUDANotAvailableWarning)
