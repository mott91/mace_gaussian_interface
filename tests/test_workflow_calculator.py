"""Tests for mace_anicc branch in workflow.calculator() and element guard.

Tests are written to be GPU-free: the actual mace_anicc constructor is mocked
so no model download or CUDA device is required.

Covers:
- calculator("mace_anicc") returns a calculator without crashing (no TypeError
  from unexpected keyword arguments)
- _check_mace_anicc_elements passes for HCNO-only molecules
- _check_mace_anicc_elements raises ValueError for F, S, etc.
- Error message lists unsupported elements (sorted)
- Element guard is wired at run_frequency_calculation call site
- Element guard is wired at run_pipeline call site
- Other calculators (mace_mp, mace_off, mace_omol) are unaffected
"""

from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _atoms(symbols: list[str]) -> Atoms:
    """Build a minimal ASE Atoms object with the given element symbols."""
    n = len(symbols)
    positions = [[i * 1.5, 0.0, 0.0] for i in range(n)]
    return Atoms(symbols=symbols, positions=positions)


# ---------------------------------------------------------------------------
# Tests for _check_mace_anicc_elements
# ---------------------------------------------------------------------------


class TestCheckMaceAniccElements:
    """Unit tests for the element guard helper function."""

    def test_hcno_only_does_not_raise(self):
        """HCNO-only molecules pass the element guard without error."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        _check_mace_anicc_elements(_atoms(["H", "C", "N", "O"]))  # must not raise

    def test_single_hcno_does_not_raise(self):
        """Single supported element also passes."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        _check_mace_anicc_elements(_atoms(["H", "H", "O"]))  # water — must not raise

    def test_f_raises_value_error(self):
        """Molecule containing F raises ValueError."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        with pytest.raises(ValueError, match="mace_anicc only supports H/C/N/O"):
            _check_mace_anicc_elements(_atoms(["H", "C", "F"]))

    def test_f_listed_in_error_message(self):
        """The ValueError message explicitly lists the unsupported element F."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        with pytest.raises(ValueError, match=r"\['F'\]"):
            _check_mace_anicc_elements(_atoms(["H", "C", "F"]))

    def test_f_and_s_both_listed_sorted(self):
        """Multiple unsupported elements are listed sorted in the error message."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        with pytest.raises(ValueError, match=r"\['F', 'S'\]"):
            _check_mace_anicc_elements(_atoms(["H", "F", "S", "C"]))

    def test_s_only_raises(self):
        """A molecule with only S (and H) raises ValueError."""
        from mace_gaussian.workflow import _check_mace_anicc_elements

        with pytest.raises(ValueError, match="mace_anicc only supports H/C/N/O"):
            _check_mace_anicc_elements(_atoms(["H", "S"]))

    def test_supported_constant_exists(self):
        """The module-level constant _MACE_ANICC_SUPPORTED_ELEMENTS is importable."""
        from mace_gaussian.workflow import _MACE_ANICC_SUPPORTED_ELEMENTS

        assert _MACE_ANICC_SUPPORTED_ELEMENTS == frozenset({"H", "C", "N", "O"})


# ---------------------------------------------------------------------------
# Tests for calculator("mace_anicc")
# ---------------------------------------------------------------------------


class TestCalculatorMaceAnicc:
    """Tests for the mace_anicc branch in calculator()."""

    def test_mace_anicc_branch_exists_and_calls_correct_api(self):
        """calculator('mace_anicc') calls mace_anicc(device='cuda') only — no model= or default_dtype=."""
        mock_calc = MagicMock()
        mock_mace_anicc = MagicMock(return_value=mock_calc)

        with patch.dict("sys.modules", {"mace.calculators": MagicMock(mace_anicc=mock_mace_anicc)}):
            from importlib import reload

            import mace_gaussian.workflow as wf

            # Patch the import inside the function
            with patch("mace_gaussian.workflow.calculator") as mock_fn:
                # We want to test the real function, so we need to avoid the mock
                pass

        # Direct approach: patch the import inside calculator()
        with patch("builtins.__import__", side_effect=_make_import_interceptor(mock_mace_anicc)):
            from mace_gaussian.workflow import calculator

            result = calculator("mace_anicc")

        mock_mace_anicc.assert_called_once_with(device="cuda")
        assert result is mock_calc

    def test_mace_anicc_does_not_pass_model_kwarg(self):
        """mace_anicc is NOT called with model= kwarg (different API from other calcs)."""
        mock_calc = MagicMock()
        mock_mace_anicc = MagicMock(return_value=mock_calc)

        with patch("builtins.__import__", side_effect=_make_import_interceptor(mock_mace_anicc)):
            from mace_gaussian.workflow import calculator

            calculator("mace_anicc")

        call_kwargs = mock_mace_anicc.call_args.kwargs
        assert "model" not in call_kwargs
        assert "default_dtype" not in call_kwargs

    def test_mace_mp_still_works(self):
        """mace_mp branch is unaffected by the new mace_anicc branch."""
        mock_calc = MagicMock()
        mock_mace_mp = MagicMock(return_value=mock_calc)

        with patch("builtins.__import__", side_effect=_make_import_interceptor_mp(mock_mace_mp)):
            from mace_gaussian.workflow import calculator

            result = calculator("mace_mp")

        mock_mace_mp.assert_called_once()
        # mace_mp is called with model="large"
        assert mock_mace_mp.call_args.kwargs.get("model") == "large"

    def test_unknown_calculator_returns_none(self):
        """calculator() for an unknown name returns None (no branch matches)."""
        from mace_gaussian.workflow import calculator

        # Patch all known calculator imports so no network/GPU access occurs
        with patch.dict(
            "sys.modules",
            {
                "mace.calculators": MagicMock(
                    mace_mp=MagicMock(),
                    mace_off=MagicMock(),
                    mace_omol=MagicMock(),
                    mace_anicc=MagicMock(),
                ),
            },
        ):
            result = calculator("nonexistent_calculator")

        assert result is None


def _make_import_interceptor(mock_mace_anicc):
    """Return a side_effect for builtins.__import__ that intercepts mace.calculators."""
    original_import = __builtins__.__import__ if hasattr(__builtins__, "__import__") else __import__

    def _import(name, *args, **kwargs):
        if name == "mace.calculators":
            mock_module = MagicMock()
            mock_module.mace_anicc = mock_mace_anicc
            return mock_module
        return original_import(name, *args, **kwargs)

    return _import


def _make_import_interceptor_polar(mock_mace_polar):
    """Return a side_effect for builtins.__import__ that intercepts mace_polar."""
    original_import = __builtins__.__import__ if hasattr(__builtins__, "__import__") else __import__

    def _import(name, *args, **kwargs):
        if name == "mace.calculators":
            mock_module = MagicMock()
            mock_module.mace_polar = mock_mace_polar
            return mock_module
        return original_import(name, *args, **kwargs)

    return _import


def _make_import_interceptor_mp(mock_mace_mp):
    """Return a side_effect for builtins.__import__ that intercepts mace_mp."""
    original_import = __builtins__.__import__ if hasattr(__builtins__, "__import__") else __import__

    def _import(name, *args, **kwargs):
        if name == "mace.calculators":
            mock_module = MagicMock()
            mock_module.mace_mp = mock_mace_mp
            return mock_module
        return original_import(name, *args, **kwargs)

    return _import


# ---------------------------------------------------------------------------
# Tests for element guard at call sites
# ---------------------------------------------------------------------------


class TestElementGuardAtCallSites:
    """Tests verifying element guard fires at run_frequency_calculation and run_pipeline."""

    def test_run_frequency_calculation_raises_for_non_hcno_with_mace_anicc(self):
        """run_frequency_calculation raises ValueError before any model load for non-HCNO."""
        from mace_gaussian.workflow import run_frequency_calculation

        bad_atoms = _atoms(["H", "F"])  # F is not supported

        mock_results_mgr = MagicMock()

        # The guard should raise before calculator() is ever called
        with patch("mace_gaussian.workflow.calculator") as mock_calc_fn:
            with pytest.raises(ValueError, match="mace_anicc only supports H/C/N/O"):
                run_frequency_calculation(
                    bad_atoms,
                    molecule_name="test_mol",
                    energy_calculator_name="mace_anicc",
                    dipole_calculator_name="espaloma",
                    results_mgr=mock_results_mgr,
                )

            # calculator() must NOT have been called (guard fires first)
            mock_calc_fn.assert_not_called()

    def test_run_frequency_calculation_does_not_raise_for_hcno_with_mace_anicc(self):
        """run_frequency_calculation does NOT raise the element guard error for HCNO molecules."""
        from mace_gaussian.workflow import run_frequency_calculation

        good_atoms = _atoms(["H", "H", "O"])  # water — HCNO only

        mock_results_mgr = MagicMock()
        mock_calc = MagicMock()

        # We allow calculator() to be called and return a mock, then the function
        # will fail later (no real Gaussian), but the guard must not raise.
        with patch("mace_gaussian.workflow.calculator", return_value=mock_calc):
            # The function will fail somewhere after the guard — that's fine
            # We just need to confirm no ValueError from the element guard
            try:
                run_frequency_calculation(
                    good_atoms,
                    molecule_name="test_mol",
                    energy_calculator_name="mace_anicc",
                    dipole_calculator_name="espaloma",
                    results_mgr=mock_results_mgr,
                )
            except ValueError as e:
                if "mace_anicc only supports" in str(e):
                    pytest.fail(f"Element guard raised unexpectedly: {e}")
            except Exception:
                pass  # Other errors are expected (no real Gaussian available)

    def test_run_frequency_calculation_other_calculators_unaffected(self):
        """Element guard is NOT triggered for mace_mp (no restriction on element types)."""
        from mace_gaussian.workflow import run_frequency_calculation

        bad_for_anicc_atoms = _atoms(["H", "F"])  # F would fail anicc guard
        mock_results_mgr = MagicMock()
        mock_calc = MagicMock()

        with patch("mace_gaussian.workflow.calculator", return_value=mock_calc):
            # Should not raise ValueError from element guard (mace_mp has no restriction)
            try:
                run_frequency_calculation(
                    bad_for_anicc_atoms,
                    molecule_name="test_mol",
                    energy_calculator_name="mace_mp",
                    dipole_calculator_name="espaloma",
                    results_mgr=mock_results_mgr,
                )
            except ValueError as e:
                if "mace_anicc only supports" in str(e):
                    pytest.fail(f"Element guard triggered incorrectly for mace_mp: {e}")
            except Exception:
                pass  # Other errors expected

    def test_run_pipeline_raises_for_non_hcno_with_mace_anicc(self):
        """run_pipeline raises ValueError before model load when optimization_calculator is mace_anicc and molecule has non-HCNO elements."""
        from mace_gaussian.workflow import run_pipeline

        # Write a temp xyz file with F atom
        import tempfile
        import os

        xyz_content = "2\n\nH 0.0 0.0 0.0\nF 1.0 0.0 0.0\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as f:
            f.write(xyz_content)
            tmp_path = f.name

        try:
            mock_results_mgr_cls = MagicMock()
            mock_results_mgr = MagicMock()
            mock_results_mgr_cls.return_value = mock_results_mgr

            # opt_file must NOT exist so we reach the optimization branch
            mock_opt_dir = MagicMock()
            mock_opt_file = MagicMock()
            mock_opt_file.exists.return_value = False
            mock_opt_dir.__truediv__ = MagicMock(return_value=mock_opt_file)
            mock_results_mgr.create_optimization_directory.return_value = mock_opt_dir

            with patch("mace_gaussian.workflow.ResultsManager", mock_results_mgr_cls):
                with patch("mace_gaussian.workflow.calculator") as mock_calc_fn:
                    with pytest.raises(ValueError, match="mace_anicc only supports H/C/N/O"):
                        run_pipeline(
                            tmp_path,
                            optimization_calculator="mace_anicc",
                            energy_calculators=["mace_mp"],
                            dipole_calculators=["espaloma"],
                            include_dft_baselines=False,
                        )

                    # calculator() must NOT have been called (guard fires first)
                    mock_calc_fn.assert_not_called()
        finally:
            os.unlink(tmp_path)

    def test_run_pipeline_does_not_raise_guard_for_hcno_with_mace_anicc(self):
        """run_pipeline does NOT trigger element guard when molecule is HCNO-only."""
        from mace_gaussian.workflow import run_pipeline

        import tempfile
        import os

        xyz_content = "3\n\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as f:
            f.write(xyz_content)
            tmp_path = f.name

        try:
            mock_results_mgr_cls = MagicMock()
            mock_results_mgr = MagicMock()
            mock_results_mgr_cls.return_value = mock_results_mgr

            mock_opt_dir = MagicMock()
            mock_opt_file = MagicMock()
            mock_opt_file.exists.return_value = False
            mock_opt_dir.__truediv__ = MagicMock(return_value=mock_opt_file)
            mock_results_mgr.create_optimization_directory.return_value = mock_opt_dir

            mock_calc = MagicMock()

            with patch("mace_gaussian.workflow.ResultsManager", mock_results_mgr_cls):
                with patch("mace_gaussian.workflow.calculator", return_value=mock_calc):
                    try:
                        run_pipeline(
                            tmp_path,
                            optimization_calculator="mace_anicc",
                            energy_calculators=["mace_mp"],
                            dipole_calculators=["espaloma"],
                            include_dft_baselines=False,
                        )
                    except ValueError as e:
                        if "mace_anicc only supports" in str(e):
                            pytest.fail(f"Element guard raised unexpectedly: {e}")
                    except Exception:
                        pass  # Other errors expected (no real Gaussian)
        finally:
            os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Tests for calculator("mace_polar")
# ---------------------------------------------------------------------------


class TestCalculatorMacePolar:
    """Tests for the mace_polar branch in calculator()."""

    def test_mace_polar_branch_calls_correct_api(self):
        """calculator('mace_polar') calls mace_polar(model='polar-1-l', device='cuda', default_dtype='float64')."""
        mock_calc = MagicMock()
        mock_mace_polar = MagicMock(return_value=mock_calc)

        with patch(
            "builtins.__import__",
            side_effect=_make_import_interceptor_polar(mock_mace_polar),
        ):
            from mace_gaussian.workflow import calculator

            result = calculator("mace_polar")

        mock_mace_polar.assert_called_once_with(
            model="polar-1-l", device="cuda", default_dtype="float64"
        )
        assert result is mock_calc

    def test_mace_polar_does_not_pass_dispersion(self):
        """mace_polar is NOT called with dispersion kwarg (unlike mace_mp/mace_off/mace_omol)."""
        mock_calc = MagicMock()
        mock_mace_polar = MagicMock(return_value=mock_calc)

        with patch(
            "builtins.__import__",
            side_effect=_make_import_interceptor_polar(mock_mace_polar),
        ):
            from mace_gaussian.workflow import calculator

            calculator("mace_polar")

        call_kwargs = mock_mace_polar.call_args.kwargs
        assert "dispersion" not in call_kwargs

    def test_mace_polar_returns_calculator_instance(self):
        """calculator('mace_polar') returns the calculator object."""
        mock_calc = MagicMock()
        mock_mace_polar = MagicMock(return_value=mock_calc)

        with patch(
            "builtins.__import__",
            side_effect=_make_import_interceptor_polar(mock_mace_polar),
        ):
            from mace_gaussian.workflow import calculator

            result = calculator("mace_polar")

        assert result is mock_calc

    def test_mace_mp_still_works_after_polar(self):
        """mace_mp branch is unaffected by the new mace_polar branch."""
        mock_calc = MagicMock()
        mock_mace_mp = MagicMock(return_value=mock_calc)

        with patch(
            "builtins.__import__",
            side_effect=_make_import_interceptor_mp(mock_mace_mp),
        ):
            from mace_gaussian.workflow import calculator

            result = calculator("mace_mp")

        mock_mace_mp.assert_called_once()
        assert mock_mace_mp.call_args.kwargs.get("model") == "large"
