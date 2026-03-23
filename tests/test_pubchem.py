"""Tests for PubChem 3D structure fetching."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Minimal valid SDF for a water molecule (3 atoms)
WATER_SDF = """\

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1173 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
"""


def _mock_response(status_code, text):
    """Create a mock requests.Response."""
    resp = MagicMock()
    resp.status_code = status_code
    resp.text = text
    resp.raise_for_status = MagicMock()
    if status_code >= 400:
        resp.raise_for_status.side_effect = Exception(f"HTTP {status_code}")
    return resp


class TestFetch3DStructure:
    """Tests for fetch_3d_structure function."""

    @patch("mace_gaussian.pubchem.requests.get")
    def test_successful_fetch_writes_xyz(self, mock_get, tmp_path):
        """fetch_3d_structure downloads SDF, converts to XYZ, returns path."""
        from mace_gaussian.pubchem import fetch_3d_structure

        mock_get.return_value = _mock_response(200, WATER_SDF)

        result = fetch_3d_structure("water", tmp_path, force=False)

        assert result == tmp_path / "water.xyz"
        assert result.exists()
        # Verify the URL contains the molecule name
        call_url = mock_get.call_args[0][0]
        assert "water" in call_url
        assert "record/SDF" in call_url
        # Verify XYZ file has valid content (at least header + atoms)
        lines = result.read_text().strip().split("\n")
        atom_count = int(lines[0].strip())
        assert atom_count == 3

    @patch("mace_gaussian.pubchem.requests.get")
    def test_unknown_molecule_raises_valueerror(self, mock_get, tmp_path):
        """Unknown molecule name returns 404 with 'No CID found'."""
        from mace_gaussian.pubchem import fetch_3d_structure

        mock_get.return_value = _mock_response(
            404,
            "Status: 404\nCode: PUGREST.NotFound\nMessage: No CID found\nDetail: none",
        )

        with pytest.raises(ValueError, match="not found on PubChem"):
            fetch_3d_structure("nonexistent_xyz_molecule_12345", tmp_path)

    @patch("mace_gaussian.pubchem.requests.get")
    def test_no_3d_conformer_raises_valueerror(self, mock_get, tmp_path):
        """Molecule exists but has no 3D conformer."""
        from mace_gaussian.pubchem import fetch_3d_structure

        mock_get.return_value = _mock_response(
            404,
            "Status: 404\nCode: PUGREST.NotFound\n"
            "Message: No records found for the given CID(s)\nDetail: none",
        )

        with pytest.raises(ValueError, match="No 3D structure"):
            fetch_3d_structure("buckminsterfullerene", tmp_path)

    def test_file_exists_raises_error(self, tmp_path):
        """Existing file raises FileExistsError when force=False."""
        from mace_gaussian.pubchem import fetch_3d_structure

        existing = tmp_path / "aspirin.xyz"
        existing.write_text("existing content")

        with pytest.raises(FileExistsError, match="already exists"):
            fetch_3d_structure("aspirin", tmp_path, force=False)

    @patch("mace_gaussian.pubchem.requests.get")
    def test_force_overwrites_existing(self, mock_get, tmp_path):
        """force=True overwrites existing file."""
        from mace_gaussian.pubchem import fetch_3d_structure

        existing = tmp_path / "water.xyz"
        existing.write_text("old content")

        mock_get.return_value = _mock_response(200, WATER_SDF)

        result = fetch_3d_structure("water", tmp_path, force=True)

        assert result.exists()
        content = result.read_text()
        assert content != "old content"

    @patch("mace_gaussian.pubchem.requests.get")
    def test_network_timeout_raises(self, mock_get, tmp_path):
        """Network timeout propagates as requests exception."""
        import requests as req

        from mace_gaussian.pubchem import fetch_3d_structure

        mock_get.side_effect = req.ConnectionError("Connection timed out")

        with pytest.raises(req.ConnectionError):
            fetch_3d_structure("aspirin", tmp_path)
