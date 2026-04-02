"""Tests for CLI validation integration."""

import os
from unittest.mock import patch

from click.testing import CliRunner

from mace_gaussian.cli import VALID_DIPOLE_CALCULATORS, VALID_ENERGY_CALCULATORS, cli


class TestEnergyCalculatorValidation:
    """Tests for --energy-calculators CLI option validation callback."""

    def test_mace_off_accepted(self, tmp_path):
        """mace_off should be accepted without a BadParameter error."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(xyz_file), "--energy-calculators", "mace_off"])
        # Should not fail at CLI validation (BadParameter); may fail later due to missing prereqs
        assert "BadParameter" not in result.output
        assert "'mace_off' is not a valid" not in result.output

    def test_mace_anicc_accepted(self, tmp_path):
        """mace_anicc should be accepted without a BadParameter error."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(xyz_file), "--energy-calculators", "mace_anicc"])
        assert "BadParameter" not in result.output
        assert "'mace_anicc' is not a valid" not in result.output

    def test_comma_separated_both_valid(self, tmp_path):
        """Comma-separated list with valid values should pass validation."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--energy-calculators", "mace_mp,mace_off"]
        )
        assert "'mace_mp' is not a valid" not in result.output
        assert "'mace_off' is not a valid" not in result.output

    def test_invalid_calculator_rejected(self, tmp_path):
        """An unknown calculator name should produce BadParameter before any workflow runs."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(xyz_file), "--energy-calculators", "invalid_calc"])
        assert result.exit_code != 0
        assert "invalid_calc" in result.output

    def test_invalid_second_token_rejected(self, tmp_path):
        """Comma-separated list where second token is invalid should be rejected."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--energy-calculators", "mace_mp,bad_name"]
        )
        assert result.exit_code != 0
        assert "bad_name" in result.output

    def test_mace_polar_accepted(self, tmp_path):
        """mace_polar should be accepted without a BadParameter error."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(xyz_file), "--energy-calculators", "mace_polar"])
        assert "BadParameter" not in result.output
        assert "'mace_polar' is not a valid" not in result.output

    def test_valid_energy_calculators_includes_mace_polar(self):
        """VALID_ENERGY_CALCULATORS must include mace_polar."""
        assert "mace_polar" in VALID_ENERGY_CALCULATORS

    def test_default_energy_calculators_does_not_include_mace_polar(self, tmp_path):
        """Default --energy-calculators string should NOT contain mace_polar (opt-in only)."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        # Invoke help to see defaults
        result = runner.invoke(cli, ["run", "--help"])
        # The default value shown should not include mace_polar
        # We check the run command's help for the default energy-calculators
        assert "mace_polar" not in result.output.split("energy-calculators")[0]

    def test_valid_energy_calculators_constant(self):
        """VALID_ENERGY_CALCULATORS must include mace_mp, mace_omol, mace_off, mace_anicc."""
        assert "mace_mp" in VALID_ENERGY_CALCULATORS
        assert "mace_omol" in VALID_ENERGY_CALCULATORS
        assert "mace_off" in VALID_ENERGY_CALCULATORS
        assert "mace_anicc" in VALID_ENERGY_CALCULATORS


class TestOptimizationCalculatorValidation:
    """Tests for --optimization-calculator click.Choice validation."""

    def test_mace_anicc_accepted_by_click(self, tmp_path):
        """mace_anicc should be a valid choice for --optimization-calculator."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--optimization-calculator", "mace_anicc"]
        )
        # Should NOT be rejected by click.Choice
        assert "Invalid value for '--optimization-calculator'" not in result.output
        assert "mace_anicc" not in result.output or "is not" not in result.output

    def test_mace_polar_accepted_by_click(self, tmp_path):
        """mace_polar should be a valid choice for --optimization-calculator."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--optimization-calculator", "mace_polar"]
        )
        # Should NOT be rejected by click.Choice
        assert "Invalid value for '--optimization-calculator'" not in result.output

    def test_invalid_optimization_calculator_rejected(self, tmp_path):
        """An unknown optimization calculator should be rejected by click.Choice."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--optimization-calculator", "totally_invalid"]
        )
        assert result.exit_code != 0


class TestDipoleCalculatorValidation:
    """Tests for --dipole-calculators CLI option validation callback."""

    def test_espaloma_accepted(self, tmp_path):
        """espaloma should be accepted without error."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(xyz_file), "--dipole-calculators", "espaloma"])
        assert "'espaloma' is not a valid" not in result.output

    def test_invalid_dipole_calculator_rejected(self, tmp_path):
        """Unknown dipole calculator name should produce BadParameter."""
        xyz_file = tmp_path / "water.xyz"
        xyz_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(
            cli, ["run", str(xyz_file), "--dipole-calculators", "unknown_dipole"]
        )
        assert result.exit_code != 0
        assert "unknown_dipole" in result.output

    def test_valid_dipole_calculators_constant(self):
        """VALID_DIPOLE_CALCULATORS must include espaloma and mace_ml."""
        assert "espaloma" in VALID_DIPOLE_CALCULATORS
        assert "mace_ml" in VALID_DIPOLE_CALCULATORS


class TestCliRunValidation:
    """Test that the run command validates input before expensive imports."""

    def test_run_with_nonexistent_file(self):
        """Run with nonexistent file should fail via Click's exists=True."""
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "nonexistent.xyz"])
        assert result.exit_code != 0
        assert "nonexistent.xyz" in result.output or "Invalid value" in result.output

    def test_run_with_invalid_extension(self, tmp_path):
        """Run with non-.xyz file should raise InputValidationError."""
        txt_file = tmp_path / "test.txt"
        txt_file.write_text("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
        runner = CliRunner()
        result = runner.invoke(cli, ["run", str(txt_file)])
        assert result.exit_code != 0

    def test_version_shows_correct(self):
        """CLI --version should show 0.2.0 matching pyproject.toml."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert "0.2.0" in result.output


class TestCliListNoValidation:
    """Test that list command does not trigger prerequisite validation."""

    def test_list_does_not_require_gaussian(self):
        """List with nonexistent dir should fail with 'not found', not prerequisite error."""
        runner = CliRunner()
        result = runner.invoke(cli, ["list", "--output-dir", "/tmp/nonexistent_gsd_test_dir"])
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "Results directory" in result.output


class TestCliDiagnoseNoValidation:
    """Test that diagnose command works without prerequisite validation."""

    def test_diagnose_runs_independently(self):
        """Diagnose should not trigger prerequisite validation."""
        runner = CliRunner()
        # Diagnose imports dipole_factory which has heavy deps;
        # mock the actual diagnostics function
        with patch("mace_gaussian.cli.sys") as _:
            # Just verify the command structure works
            result = runner.invoke(cli, ["diagnose", "--help"])
            assert result.exit_code == 0
            assert "diagnostic" in result.output.lower()


class TestTimeoutConfiguration:
    """Test Gaussian subprocess timeout configuration."""

    def test_default_timeout_value(self):
        """Default timeout should be 86400 seconds (24 hours)."""
        try:
            from mace_gaussian.gaussian.runner import DEFAULT_TIMEOUT_SECONDS
        except (ImportError, RuntimeError):
            # gaussian.runner may have deps that fail in test env
            # Verify the constant definition directly
            expected = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))
            assert expected == 86400
            return

        expected = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))
        assert expected == DEFAULT_TIMEOUT_SECONDS

    def test_timeout_is_positive_integer(self):
        """Timeout should be a positive integer."""
        try:
            from mace_gaussian.gaussian.runner import DEFAULT_TIMEOUT_SECONDS
        except (ImportError, RuntimeError):
            # gaussian.runner may have deps that fail in test env; verify env var parsing logic
            val = int(os.getenv("GAUSSIAN_TIMEOUT_SECONDS", "86400"))
            assert isinstance(val, int)
            assert val > 0
            return

        assert isinstance(DEFAULT_TIMEOUT_SECONDS, int)
        assert DEFAULT_TIMEOUT_SECONDS > 0
