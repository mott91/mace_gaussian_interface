"""Tests for validation.py -- prerequisite checks, XYZ validation, device detection, metadata.

Uses unittest.mock to isolate from real system state (g16 availability, CUDA, etc.).
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from mace_gaussian.utils.exceptions import InputValidationError, PrerequisiteError
from mace_gaussian.utils.validation import (
    check_dipole_model,
    check_formchk_available,
    check_gaussian_available,
    check_helper_script,
    collect_version_metadata,
    detect_device,
    validate_all_prerequisites,
    validate_xyz_file,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

WATER_XYZ = """\
3
water
O  0.0000  0.0000  0.1173
H  0.0000  0.7572 -0.4692
H  0.0000 -0.7572 -0.4692
"""


@pytest.fixture()
def water_xyz(tmp_path: Path) -> Path:
    """Write a minimal water XYZ file and return its path."""
    p = tmp_path / "water.xyz"
    p.write_text(WATER_XYZ)
    return p


# ---------------------------------------------------------------------------
# check_gaussian_available
# ---------------------------------------------------------------------------


class TestCheckGaussianAvailable:
    """Tests for check_gaussian_available."""

    @patch("mace_gaussian.utils.validation.shutil.which", return_value="/usr/bin/g16")
    def test_found(self, mock_which):
        """Returns path when g16 is on PATH."""
        result = check_gaussian_available()
        assert result == "/usr/bin/g16"
        mock_which.assert_called_once_with("g16")

    @patch("mace_gaussian.utils.validation.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        """Raises PrerequisiteError with 'g16' in message when missing."""
        with pytest.raises(PrerequisiteError, match="g16"):
            check_gaussian_available()


# ---------------------------------------------------------------------------
# check_formchk_available
# ---------------------------------------------------------------------------


class TestCheckFormchkAvailable:
    """Tests for check_formchk_available."""

    @patch("mace_gaussian.utils.validation.shutil.which", return_value="/usr/bin/formchk")
    def test_found(self, mock_which):
        """Returns path when formchk is on PATH."""
        assert check_formchk_available() == "/usr/bin/formchk"
        mock_which.assert_called_once_with("formchk")

    @patch("mace_gaussian.utils.validation.shutil.which", return_value=None)
    def test_not_found(self, mock_which):
        """Raises PrerequisiteError when formchk is missing."""
        with pytest.raises(PrerequisiteError, match="formchk"):
            check_formchk_available()


# ---------------------------------------------------------------------------
# check_dipole_model
# ---------------------------------------------------------------------------


class TestCheckDipoleModel:
    """Tests for check_dipole_model."""

    def test_existing_file(self, tmp_path: Path):
        """Returns Path for an existing model file."""
        model = tmp_path / "model_1.model"
        model.write_text("dummy")
        result = check_dipole_model(str(model))
        assert isinstance(result, Path)
        assert result.exists()

    def test_missing_file(self):
        """Raises PrerequisiteError with path in message for missing file."""
        bad_path = "/nonexistent/model.model"
        with pytest.raises(PrerequisiteError, match=bad_path):
            check_dipole_model(bad_path)


# ---------------------------------------------------------------------------
# check_helper_script
# ---------------------------------------------------------------------------


class TestCheckHelperScript:
    """Tests for check_helper_script."""

    def test_existing_script(self, tmp_path: Path):
        """Returns Path for an existing script."""
        script = tmp_path / "gm_helper.py"
        script.write_text("# helper")
        result = check_helper_script(str(script))
        assert isinstance(result, Path)
        assert result.exists()

    def test_missing_script(self):
        """Raises PrerequisiteError for missing script."""
        bad_path = "/nonexistent/helper.py"
        with pytest.raises(PrerequisiteError, match=bad_path):
            check_helper_script(bad_path)


# ---------------------------------------------------------------------------
# validate_xyz_file
# ---------------------------------------------------------------------------


class TestValidateXyzFile:
    """Tests for validate_xyz_file."""

    def test_valid_water(self, water_xyz: Path):
        """Valid water XYZ returns dict with n_atoms=3 and correct symbols."""
        result = validate_xyz_file(str(water_xyz))
        assert result["n_atoms"] == 3
        assert result["symbols"] == ["O", "H", "H"]
        assert "path" in result

    def test_nonexistent_file(self):
        """Raises InputValidationError for missing file."""
        with pytest.raises(InputValidationError, match="not found"):
            validate_xyz_file("/nonexistent/molecule.xyz")

    def test_wrong_extension(self, tmp_path: Path):
        """Raises InputValidationError for non-.xyz extension."""
        txt_file = tmp_path / "molecule.txt"
        txt_file.write_text(WATER_XYZ)
        with pytest.raises(InputValidationError, match="extension"):
            validate_xyz_file(str(txt_file))

    def test_unparseable_content(self, tmp_path: Path):
        """Raises InputValidationError for garbage XYZ content."""
        bad_xyz = tmp_path / "bad.xyz"
        bad_xyz.write_text("this is not valid xyz data\ngarbage\n")
        with pytest.raises(InputValidationError, match="parse"):
            validate_xyz_file(str(bad_xyz))


# ---------------------------------------------------------------------------
# validate_all_prerequisites
# ---------------------------------------------------------------------------


class TestValidateAllPrerequisites:
    """Tests for validate_all_prerequisites."""

    @patch("mace_gaussian.utils.validation.shutil.which", return_value="/usr/bin/g16")
    def test_all_pass(self, mock_which, tmp_path: Path):
        """Returns summary dict when all checks pass."""
        # Create temp files for model and script
        model = tmp_path / "model.model"
        model.write_text("dummy")
        script = tmp_path / "helper.py"
        script.write_text("# helper")

        # Mock which to return paths for both g16 and formchk
        mock_which.side_effect = ["/usr/bin/g16", "/usr/bin/formchk"]

        result = validate_all_prerequisites(
            check_gaussian=True,
            check_formchk_tool=True,
            dipole_model_path=str(model),
            helper_script_path=str(script),
        )

        assert "g16_path" in result
        assert "formchk_path" in result

    @patch("mace_gaussian.utils.validation.shutil.which", return_value=None)
    def test_raises_on_first_failure(self, mock_which):
        """Raises PrerequisiteError on first failing check (g16 missing)."""
        with pytest.raises(PrerequisiteError, match="g16"):
            validate_all_prerequisites(check_gaussian=True, check_formchk_tool=False)


# ---------------------------------------------------------------------------
# detect_device
# ---------------------------------------------------------------------------


class TestDetectDevice:
    """Tests for detect_device."""

    @patch("torch.cuda.get_device_name", return_value="NVIDIA RTX 2070")
    @patch("torch.cuda.is_available", return_value=True)
    def test_cuda_available(self, mock_avail, mock_name):
        """Returns 'cuda' when GPU is available."""
        assert detect_device() == "cuda"

    @patch("torch.cuda.is_available", return_value=False)
    def test_cpu_fallback(self, mock_avail):
        """Returns 'cpu' when CUDA is not available."""
        assert detect_device() == "cpu"


# ---------------------------------------------------------------------------
# collect_version_metadata
# ---------------------------------------------------------------------------


class TestCollectVersionMetadata:
    """Tests for collect_version_metadata."""

    def test_required_keys(self):
        """Returned dict has tool_version, python_version, platform."""
        meta = collect_version_metadata()
        assert "tool_version" in meta
        assert "python_version" in meta
        assert "platform" in meta

    def test_torch_version_present(self):
        """Returned dict has torch_version key."""
        meta = collect_version_metadata()
        assert "torch_version" in meta

    def test_json_serializable(self):
        """All values must be JSON-serializable (no numpy, torch types)."""
        meta = collect_version_metadata()
        # This will raise TypeError if any value is not serializable
        serialized = json.dumps(meta)
        assert isinstance(serialized, str)
