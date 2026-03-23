"""Tests for mace_gaussian.batch module -- manifest-driven batch runner."""

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mace_gaussian.batch import (
    STATUS_COMPLETE,
    STATUS_FAILED,
    STATUS_PENDING,
    _combination_key,
    load_manifest,
    parse_batch_file,
    run_batch,
    save_manifest,
)


def _write_xyz(path: Path, name: str = "water") -> Path:
    """Write a minimal 3-atom XYZ file for testing."""
    xyz = path / f"{name}.xyz"
    xyz.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
    return xyz


def _write_batch_file(path: Path, xyz_paths: list[Path]) -> Path:
    """Write a molecules.txt file listing XYZ paths."""
    batch = path / "molecules.txt"
    batch.write_text("\n".join(str(p) for p in xyz_paths) + "\n")
    return batch


class TestLoadManifest:
    def test_returns_empty_skeleton_when_file_missing(self, tmp_path):
        result = load_manifest(tmp_path / "nonexistent.json")
        assert result == {"molecules": {}}

    def test_loads_existing_manifest(self, tmp_path):
        manifest_path = tmp_path / "batch_manifest.json"
        data = {"molecules": {"water": {"geometry_opt": "complete"}}}
        manifest_path.write_text(json.dumps(data))
        result = load_manifest(manifest_path)
        assert result == data


class TestSaveManifest:
    def test_roundtrip(self, tmp_path):
        manifest_path = tmp_path / "batch_manifest.json"
        data = {"molecules": {"water": {"geometry_opt": "complete", "combinations": {}}}}
        save_manifest(data, manifest_path)
        loaded = load_manifest(manifest_path)
        assert loaded == data

    def test_atomic_write_uses_os_replace(self, tmp_path):
        """save_manifest must use atomic write (tempfile + os.replace)."""
        manifest_path = tmp_path / "batch_manifest.json"
        data = {"molecules": {}}
        with patch("mace_gaussian.batch.os.replace", wraps=os.replace) as mock_replace:
            save_manifest(data, manifest_path)
            mock_replace.assert_called_once()
            # The second arg should be the final path
            assert mock_replace.call_args[0][1] == str(manifest_path)


class TestParseBatchFile:
    def test_reads_molecule_paths(self, tmp_path):
        xyz1 = _write_xyz(tmp_path, "water")
        xyz2 = _write_xyz(tmp_path, "methane")
        batch = _write_batch_file(tmp_path, [xyz1, xyz2])
        result = parse_batch_file(batch)
        assert len(result) == 2
        assert result[0] == xyz1
        assert result[1] == xyz2

    def test_skips_blank_lines_and_comments(self, tmp_path):
        xyz1 = _write_xyz(tmp_path, "water")
        batch = tmp_path / "molecules.txt"
        batch.write_text(f"# This is a comment\n\n{xyz1}\n\n# Another comment\n")
        result = parse_batch_file(batch)
        assert len(result) == 1
        assert result[0] == xyz1

    def test_resolves_relative_paths(self, tmp_path, monkeypatch):
        xyz = _write_xyz(tmp_path, "water")
        monkeypatch.chdir(tmp_path)
        batch = tmp_path / "molecules.txt"
        batch.write_text("water.xyz\n")
        result = parse_batch_file(batch)
        assert len(result) == 1
        assert result[0].is_absolute()
        assert result[0].exists()

    def test_raises_on_missing_file(self, tmp_path):
        batch = tmp_path / "molecules.txt"
        batch.write_text("nonexistent.xyz\n")
        with pytest.raises(FileNotFoundError):
            parse_batch_file(batch)


class TestCombinationKey:
    def test_format(self):
        assert _combination_key("mace_mp", "espaloma") == "mace_mp_espaloma"


def _make_mock_atoms():
    """Create a mock ASE Atoms object."""
    mock = MagicMock()
    mock.info = {"charge": 0.0, "spin": 1.0}
    return mock


class TestRunBatch:
    """Tests for run_batch with mocked workflow functions."""

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_basic_two_molecules(self, mock_geom_opt, mock_freq_calc, mock_read, tmp_path):
        """run_batch with 2 molecules and 1x1 combo calls geom opt and freq calc per molecule."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True

        xyz1 = _write_xyz(tmp_path, "water")
        xyz2 = _write_xyz(tmp_path, "methane")
        batch_file = _write_batch_file(tmp_path, [xyz1, xyz2])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        summary = run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        assert mock_geom_opt.call_count == 2
        assert mock_freq_calc.call_count == 2
        assert summary["complete"] == 2
        assert summary["failed"] == 0

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_restart_skips_complete(self, mock_geom_opt, mock_freq_calc, mock_read, tmp_path):
        """run_batch skips combinations marked 'complete' in existing manifest."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True

        xyz1 = _write_xyz(tmp_path, "water")
        batch_file = _write_batch_file(tmp_path, [xyz1])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        # Pre-populate manifest with one complete combination
        manifest_path = Path(output_dir) / "batch_manifest.json"
        manifest = {
            "molecules": {
                "water": {
                    "xyz_path": str(xyz1),
                    "geometry_opt": STATUS_COMPLETE,
                    "dft_baseline": "skipped",
                    "combinations": {
                        "mace_mp_espaloma": {"status": STATUS_COMPLETE, "runtime_s": 10.0},
                    },
                }
            }
        }
        manifest_path.write_text(json.dumps(manifest))

        # Also create the optimized geometry so it can be loaded
        opt_dir = Path(output_dir) / "water" / "geometry_opt"
        opt_dir.mkdir(parents=True, exist_ok=True)
        (opt_dir / "optimized.xyz").write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
        )

        summary = run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        # Should skip the already-complete combination
        mock_freq_calc.assert_not_called()
        assert summary["skipped"] == 1
        assert summary["complete"] == 0

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_failure_isolation(self, mock_geom_opt, mock_freq_calc, mock_read, tmp_path):
        """Failed combination is recorded but batch continues to next molecule."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms

        # First call fails, second succeeds
        mock_freq_calc.side_effect = [Exception("GPU OOM"), True]

        xyz1 = _write_xyz(tmp_path, "water")
        xyz2 = _write_xyz(tmp_path, "methane")
        batch_file = _write_batch_file(tmp_path, [xyz1, xyz2])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        summary = run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        assert summary["failed"] == 1
        assert summary["complete"] == 1

        # Check manifest records the error
        manifest = load_manifest(Path(output_dir) / "batch_manifest.json")
        water_combos = manifest["molecules"]["water"]["combinations"]
        assert water_combos["mace_mp_espaloma"]["status"] == STATUS_FAILED
        assert "GPU OOM" in water_combos["mace_mp_espaloma"]["error"]

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_dft_baselines")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_skip_dft_baseline_true(
        self, mock_geom_opt, mock_freq_calc, mock_dft, mock_read, tmp_path
    ):
        """skip_dft_baseline=True should not call run_dft_baselines."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True

        xyz1 = _write_xyz(tmp_path, "water")
        batch_file = _write_batch_file(tmp_path, [xyz1])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        mock_dft.assert_not_called()

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_dft_baselines")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_skip_dft_baseline_false(
        self, mock_geom_opt, mock_freq_calc, mock_dft, mock_read, tmp_path
    ):
        """skip_dft_baseline=False should call run_dft_baselines once per molecule."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True
        mock_dft.return_value = {"b3lyp": True}

        xyz1 = _write_xyz(tmp_path, "water")
        batch_file = _write_batch_file(tmp_path, [xyz1])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=False,
            output_dir=output_dir,
        )

        mock_dft.assert_called_once()

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_returns_summary_dict(self, mock_geom_opt, mock_freq_calc, mock_read, tmp_path):
        """run_batch returns summary dict with complete, failed, skipped counts."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True

        xyz1 = _write_xyz(tmp_path, "water")
        batch_file = _write_batch_file(tmp_path, [xyz1])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        summary = run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        assert "complete" in summary
        assert "failed" in summary
        assert "skipped" in summary
        assert isinstance(summary["complete"], int)
        assert isinstance(summary["failed"], int)
        assert isinstance(summary["skipped"], int)

    @patch("mace_gaussian.batch.read")
    @patch("mace_gaussian.workflow.run_frequency_calculation")
    @patch("mace_gaussian.workflow.run_geometry_optimization")
    def test_manifest_written_after_run(self, mock_geom_opt, mock_freq_calc, mock_read, tmp_path):
        """Manifest file exists after run_batch completes."""
        mock_atoms = _make_mock_atoms()
        mock_geom_opt.return_value = mock_atoms
        mock_read.return_value = mock_atoms
        mock_freq_calc.return_value = True

        xyz1 = _write_xyz(tmp_path, "water")
        batch_file = _write_batch_file(tmp_path, [xyz1])
        output_dir = str(tmp_path / "results")
        os.makedirs(output_dir, exist_ok=True)

        run_batch(
            batch_file=batch_file,
            optimization_calculator="mace_omol",
            energy_calculators=["mace_mp"],
            dipole_calculators=["espaloma"],
            skip_dft_baseline=True,
            output_dir=output_dir,
        )

        manifest_path = Path(output_dir) / "batch_manifest.json"
        assert manifest_path.exists()
        manifest = json.loads(manifest_path.read_text())
        assert "molecules" in manifest
        assert "water" in manifest["molecules"]
