"""Tests for calculators/mace_loader.py -- safe MACE dipole loading mechanism.

Covers requirement: STRUCT-03 (MACE module monkey-patching replaced with safe loading).

All tests mock heavy dependencies (torch, mace_dipole_core) so they run
without GPU or ML packages installed.

Strategy: Pre-populate sys.modules with mocks for mace_dipole_core sub-packages
BEFORE importing calculators.mace_loader, then restore after import.
"""

from __future__ import annotations

import importlib
import inspect
import pickle
import sys
from io import BytesIO
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Pre-mock heavy dependencies before importing the module under test
# ---------------------------------------------------------------------------

_saved_modules: dict = {}
_heavy_deps = [
    "espaloma_charge",
    "rdkit",
    "rdkit.Chem",
    "xtb",
    "xtb.ase",
    "xtb.ase.calculator",
    "mace_dipole_core",
    "mace_dipole_core.modules",
    "mace_dipole_core.modules.models",
    "mace_dipole_core.calculators",
    "mace_dipole_core.calculators.mace",
]

for _dep in _heavy_deps:
    _saved_modules[_dep] = sys.modules.get(_dep)
    if _dep not in sys.modules:
        sys.modules[_dep] = MagicMock()

# Set up a mock class on the dipole models module so find_class can find it
_mock_dipole_models = sys.modules["mace_dipole_core.modules.models"]
_mock_dipole_models.AtomicDielectricMACE = type("AtomicDielectricMACE", (), {})

from mace_gaussian.calculators.mace_loader import (  # noqa: E402
    _DIPOLE_PKG_DIR,
    MACEDipoleCalculator,
    _DipoleModelUnpickler,
    _DipolePickleModule,
    _ensure_dipole_importable,
)

# Restore original module state
for _dep in _heavy_deps:
    if _saved_modules[_dep] is None:
        sys.modules.pop(_dep, None)
    else:
        sys.modules[_dep] = _saved_modules[_dep]


# ---------------------------------------------------------------------------
# TestDipoleModelUnpickler
# ---------------------------------------------------------------------------


class TestDipoleModelUnpickler:
    """Tests for _DipoleModelUnpickler class resolution remapping."""

    def test_remaps_mace_modules_models(self):
        """find_class remaps mace.modules.models to mace_dipole_core.modules.models."""
        sentinel_cls = type("AtomicDielectricMACE", (), {"_sentinel": True})
        mock_dipole_models = MagicMock()
        mock_dipole_models.AtomicDielectricMACE = sentinel_cls

        # Build a consistent mock package hierarchy so the import machinery
        # resolves the same object whether it uses sys.modules or parent attrs
        mock_core = MagicMock()
        mock_modules = MagicMock()
        mock_core.modules = mock_modules
        mock_modules.models = mock_dipole_models

        with patch.dict(
            sys.modules,
            {
                "mace_dipole_core": mock_core,
                "mace_dipole_core.modules": mock_modules,
                "mace_dipole_core.modules.models": mock_dipole_models,
            },
        ):
            unpickler = _DipoleModelUnpickler(BytesIO(b""))
            cls = unpickler.find_class("mace.modules.models", "AtomicDielectricMACE")

        assert cls is sentinel_cls

    def test_passes_through_other_modules(self):
        """find_class delegates to super() for non-mace modules."""
        unpickler = _DipoleModelUnpickler(BytesIO(b""))

        # torch.nn.Module should resolve via standard pickle resolution
        cls = unpickler.find_class("builtins", "dict")
        assert cls is dict

    def test_falls_through_if_attr_missing(self):
        """find_class falls through to super() when dipole module lacks the class."""
        mock_dipole_models = MagicMock(spec=[])  # Empty spec = no attributes
        mock_core = MagicMock()
        mock_modules = MagicMock()
        mock_core.modules = mock_modules
        mock_modules.models = mock_dipole_models

        with patch.dict(
            sys.modules,
            {
                "mace_dipole_core": mock_core,
                "mace_dipole_core.modules": mock_modules,
                "mace_dipole_core.modules.models": mock_dipole_models,
            },
        ):
            unpickler = _DipoleModelUnpickler(BytesIO(b""))

            # When the dipole module lacks the class, find_class falls through
            # to super().find_class() which resolves against the real
            # mace.modules.models (since mace-torch is installed).
            # Verify the fallthrough by checking we get a class from the real
            # module, not from the mock.
            cls = unpickler.find_class("mace.modules.models", "AtomicDielectricMACE")
            # The real mace.modules.models.AtomicDielectricMACE should be returned
            # (not from dipole mock, which has spec=[] so getattr returns None)
            import mace.modules.models as real_mace_models

            assert cls is real_mace_models.AtomicDielectricMACE


# ---------------------------------------------------------------------------
# TestDipolePickleModule
# ---------------------------------------------------------------------------


class TestDipolePickleModule:
    """Tests for _DipolePickleModule as a drop-in pickle replacement."""

    def test_has_unpickler_attribute(self):
        """_DipolePickleModule.Unpickler is _DipoleModelUnpickler."""
        pm = _DipolePickleModule()
        assert pm.Unpickler is _DipoleModelUnpickler

    def test_delegates_to_stdlib_pickle(self):
        """Other attributes are delegated to stdlib pickle."""
        pm = _DipolePickleModule()
        assert pm.loads is pickle.loads
        assert pm.dumps is pickle.dumps
        assert pm.HIGHEST_PROTOCOL is pickle.HIGHEST_PROTOCOL


# ---------------------------------------------------------------------------
# TestEnsureDipoleImportable
# ---------------------------------------------------------------------------


class TestEnsureDipoleImportable:
    """Tests for _ensure_dipole_importable sys.path setup."""

    def test_adds_dipole_pkg_to_sys_path(self):
        """_ensure_dipole_importable adds _DIPOLE_PKG_DIR to sys.path."""
        # Remove it first if present
        original_path = sys.path.copy()
        if _DIPOLE_PKG_DIR in sys.path:
            sys.path.remove(_DIPOLE_PKG_DIR)

        try:
            _ensure_dipole_importable()
            assert _DIPOLE_PKG_DIR in sys.path
        finally:
            sys.path[:] = original_path

    def test_idempotent(self):
        """Calling _ensure_dipole_importable twice adds the path only once."""
        original_path = sys.path.copy()
        if _DIPOLE_PKG_DIR in sys.path:
            sys.path.remove(_DIPOLE_PKG_DIR)

        try:
            _ensure_dipole_importable()
            _ensure_dipole_importable()
            assert sys.path.count(_DIPOLE_PKG_DIR) == 1
        finally:
            sys.path[:] = original_path


# ---------------------------------------------------------------------------
# TestMACEDipoleCalculatorWrapper
# ---------------------------------------------------------------------------


class TestMACEDipoleCalculatorWrapper:
    """Tests for MACEDipoleCalculator wrapper class."""

    def test_lazy_initialization(self):
        """Calculator attribute is None before first calculation."""
        calc = MACEDipoleCalculator(model_path="/fake/model", device="cpu")
        assert calc.calc is None

    def test_stores_config(self):
        """Constructor stores model_path and device."""
        calc = MACEDipoleCalculator(model_path="/fake/model.pt", device="cpu")
        assert calc.model_path == "/fake/model.pt"
        assert calc.device == "cpu"

    def test_no_cleanup_mace_modules_reference(self):
        """Module source must not contain cleanup_mace_modules."""
        source = inspect.getsource(importlib.import_module("mace_gaussian.calculators.mace_loader"))
        assert "cleanup_mace_modules" not in source

    def test_no_sys_modules_mutation(self):
        """Module source must not contain sys.modules["mace string."""
        source = inspect.getsource(importlib.import_module("mace_gaussian.calculators.mace_loader"))
        assert 'sys.modules["mace' not in source


# ---------------------------------------------------------------------------
# TestMACEDipoleCalculatorAutograd
# ---------------------------------------------------------------------------

import numpy as np


class TestMACEDipoleCalculatorAutograd:
    """Tests for MACEDipoleCalculator.calculate_dipole_derivatives() autograd override."""

    def _make_calculator_with_mock(self, use_polarizability=False):
        """Create a MACEDipoleCalculator with a pre-injected mock calc."""
        calc = MACEDipoleCalculator(model_path="/fake/model", device="cpu")
        mock_calc = MagicMock()
        mock_model = MagicMock(use_polarizability=use_polarizability)
        mock_calc.models = [mock_model]
        calc.calc = mock_calc  # bypass _ensure_calculator
        return calc, mock_calc

    def _make_mock_atoms(self, n_atoms=3):
        """Create a mock atoms object with n_atoms."""
        mock_atoms = MagicMock()
        mock_atoms.__len__ = lambda self: n_atoms
        return mock_atoms

    def test_calls_get_dielectric_derivatives(self):
        """Override calls get_dielectric_derivatives, not base class finite differences."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(3)
        dmu_dr_raw = np.random.randn(3, 3, 3)
        mock_calc.get_dielectric_derivatives.return_value = dmu_dr_raw

        calc.calculate_dipole_derivatives(mock_atoms)

        mock_calc.get_dielectric_derivatives.assert_called_once_with(mock_atoms)

    def test_no_polarizability_returns_correct_shape(self):
        """When use_polarizability=False, result is reshaped to (3*N, 3)."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(3)
        dmu_dr_raw = np.random.randn(3, 3, 3)  # (dipole_comp, atom, spatial_dir)
        mock_calc.get_dielectric_derivatives.return_value = dmu_dr_raw

        result = calc.calculate_dipole_derivatives(mock_atoms)

        assert result.shape == (9, 3)
        assert isinstance(result, np.ndarray)

    def test_polarizability_true_extracts_tuple(self):
        """When use_polarizability=True, dmu_dr is from index 0 of tuple."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=True)
        mock_atoms = self._make_mock_atoms(3)
        dmu_dr_raw = np.random.randn(3, 3, 3)
        dalpha_dr_raw = np.random.randn(3, 9, 3)
        mock_calc.get_dielectric_derivatives.return_value = (dmu_dr_raw, dalpha_dr_raw)

        result = calc.calculate_dipole_derivatives(mock_atoms)

        assert result.shape == (9, 3)
        assert calc._last_dalpha_dr is dalpha_dr_raw

    def test_no_polarizability_clears_dalpha(self):
        """When use_polarizability=False, _last_dalpha_dr is set to None."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(3)
        dmu_dr_raw = np.random.randn(3, 3, 3)
        mock_calc.get_dielectric_derivatives.return_value = dmu_dr_raw

        calc.calculate_dipole_derivatives(mock_atoms)

        assert calc._last_dalpha_dr is None

    def test_exception_propagates(self):
        """If get_dielectric_derivatives raises, it propagates (no try/except)."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(3)
        mock_calc.get_dielectric_derivatives.side_effect = RuntimeError("autograd failed")

        import pytest

        with pytest.raises(RuntimeError, match="autograd failed"):
            calc.calculate_dipole_derivatives(mock_atoms)

    def test_returned_dtype_is_float64(self):
        """Returned array dtype is float64-compatible."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(3)
        dmu_dr_raw = np.random.randn(3, 3, 3).astype(np.float64)
        mock_calc.get_dielectric_derivatives.return_value = dmu_dr_raw

        result = calc.calculate_dipole_derivatives(mock_atoms)

        assert result.dtype == np.float64

    def test_reshape_is_correct(self):
        """Verify the transpose+reshape produces correct layout."""
        calc, mock_calc = self._make_calculator_with_mock(use_polarizability=False)
        mock_atoms = self._make_mock_atoms(2)
        # (3, 2, 3) = (dipole_comp, atom, spatial_dir)
        dmu_dr_raw = np.arange(18).reshape(3, 2, 3).astype(np.float64)
        mock_calc.get_dielectric_derivatives.return_value = dmu_dr_raw

        result = calc.calculate_dipole_derivatives(mock_atoms)

        # transpose(1,2,0): (2, 3, 3) = (atom, spatial_dir, dipole_comp)
        # reshape(6, 3): row 3i+j = d(dipole)/d(atom_i, dir_j)
        expected = dmu_dr_raw.transpose(1, 2, 0).reshape(6, 3)
        np.testing.assert_array_equal(result, expected)

    def test_init_has_last_dalpha_dr(self):
        """__init__ initializes _last_dalpha_dr = None."""
        calc = MACEDipoleCalculator(model_path="/fake/model", device="cpu")
        assert hasattr(calc, "_last_dalpha_dr")
        assert calc._last_dalpha_dr is None
