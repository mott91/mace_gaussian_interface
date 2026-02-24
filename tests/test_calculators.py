"""Tests for the calculators/ package covering interface compliance,
factory pattern, and configuration constants.

Covers requirements: STRUCT-02 (calculator hierarchy extraction).
All tests mock heavy dependencies (espaloma, xtb, MACE models) so they
run without GPU or ML packages installed.

Strategy: We pre-populate sys.modules with mocks for all heavy dependencies
(espaloma_charge, rdkit, xtb, mace_dipole_core) BEFORE importing the
calculators package, so that _check_availability() calls don't trigger
real imports.
"""

import sys
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Pre-mock heavy dependencies before importing calculators
# ---------------------------------------------------------------------------

# Save any existing modules so we can restore them later
_saved_modules = {}
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

for dep in _heavy_deps:
    _saved_modules[dep] = sys.modules.get(dep)
    if dep not in sys.modules:
        sys.modules[dep] = MagicMock()

# Now import calculators - the _check_availability calls will hit our mocks
from mace_gaussian.calculators.base import DipoleCalculatorBase  # noqa: E402
from mace_gaussian.calculators.espaloma import EspalomaDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.mace_ml import MACEMLDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.xtb import XTBDipoleCalculator  # noqa: E402

# Restore original module state for heavy deps
for dep in _heavy_deps:
    if _saved_modules[dep] is None:
        sys.modules.pop(dep, None)
    else:
        sys.modules[dep] = _saved_modules[dep]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ALL_CALCULATOR_CLASSES = [EspalomaDipoleCalculator, MACEMLDipoleCalculator, XTBDipoleCalculator]


def _make_mock_calculator(name: str, available: bool) -> MagicMock:
    """Create a mock that quacks like a DipoleCalculatorBase instance."""
    mock = MagicMock(spec=DipoleCalculatorBase)
    mock.name = name
    mock.available = available
    return mock


# ---------------------------------------------------------------------------
# TestCalculatorInterface -- ABC / subclass / method checks
# ---------------------------------------------------------------------------


class TestCalculatorInterface:
    """Tests for calculator interface compliance."""

    @pytest.mark.parametrize("cls", ALL_CALCULATOR_CLASSES, ids=lambda c: c.__name__)
    def test_is_subclass(self, cls):
        """Every concrete calculator must be a subclass of DipoleCalculatorBase."""
        assert issubclass(cls, DipoleCalculatorBase)

    @pytest.mark.parametrize("cls", ALL_CALCULATOR_CLASSES, ids=lambda c: c.__name__)
    def test_has_required_methods(self, cls):
        """Every concrete calculator must implement the required methods."""
        for method_name in (
            "calculate_dipole",
            "calculate_dipole_derivatives",
            "_check_availability",
        ):
            assert hasattr(cls, method_name), f"{cls.__name__} missing {method_name}"
            assert callable(getattr(cls, method_name))

    def test_base_class_is_abstract(self):
        """DipoleCalculatorBase cannot be instantiated directly."""
        with pytest.raises(TypeError, match="abstract method"):
            DipoleCalculatorBase("test")


# ---------------------------------------------------------------------------
# TestDipoleCalculatorFactory -- factory pattern tests
# ---------------------------------------------------------------------------


class TestDipoleCalculatorFactory:
    """Tests for the DipoleCalculatorFactory class.

    Every test patches the three calculator constructors at their import
    location inside factory.py to avoid real heavy-dependency imports.
    """

    def _make_factory(self, espaloma_avail=True, xtb_avail=True, mace_ml_avail=True):
        """Create a factory with mocked calculators."""
        from mace_gaussian.calculators.factory import DipoleCalculatorFactory

        mock_esp = _make_mock_calculator("espaloma", espaloma_avail)
        mock_xtb = _make_mock_calculator("xtb", xtb_avail)
        mock_mace = _make_mock_calculator("mace_ml", mace_ml_avail)

        p1 = patch(
            "mace_gaussian.calculators.factory.EspalomaDipoleCalculator", return_value=mock_esp
        )
        p2 = patch("mace_gaussian.calculators.factory.XTBDipoleCalculator", return_value=mock_xtb)
        p3 = patch(
            "mace_gaussian.calculators.factory.MACEMLDipoleCalculator", return_value=mock_mace
        )
        with p1, p2, p3:
            factory = DipoleCalculatorFactory()

        return factory

    def test_get_known_calculator(self):
        """get_calculator returns the correct mock for a known name."""
        factory = self._make_factory()
        calc = factory.get_calculator("espaloma")
        assert calc.name == "espaloma"

    def test_get_unknown_calculator_raises(self):
        """get_calculator raises ValueError for an unregistered name."""
        factory = self._make_factory()
        with pytest.raises(ValueError, match="Unknown"):
            factory.get_calculator("nonexistent")

    def test_get_unavailable_calculator_raises(self):
        """get_calculator raises RuntimeError when the requested calc is unavailable."""
        factory = self._make_factory(espaloma_avail=False)
        with pytest.raises(RuntimeError, match="not available"):
            factory.get_calculator("espaloma")

    def test_auto_selects_first_available(self):
        """auto mode returns the first available calculator per preferred_order."""
        # preferred_order is [mace_ml, espaloma, xtb]
        # Make mace_ml unavailable so espaloma should be selected
        factory = self._make_factory(mace_ml_avail=False, espaloma_avail=True, xtb_avail=True)
        calc = factory.get_calculator("auto")
        assert calc.name == "espaloma"

    def test_auto_selects_mace_ml_when_available(self):
        """auto mode prefers mace_ml (first in preferred_order)."""
        factory = self._make_factory(mace_ml_avail=True, espaloma_avail=True, xtb_avail=True)
        calc = factory.get_calculator("auto")
        assert calc.name == "mace_ml"

    def test_auto_raises_when_none_available(self):
        """auto mode raises RuntimeError when no calculators are available."""
        factory = self._make_factory(espaloma_avail=False, xtb_avail=False, mace_ml_avail=False)
        with pytest.raises(RuntimeError, match="No dipole calculators available"):
            factory.get_calculator("auto")

    def test_list_available(self):
        """list_available returns a dict mapping names to availability booleans."""
        factory = self._make_factory(espaloma_avail=True, xtb_avail=False, mace_ml_avail=True)
        result = factory.list_available()

        assert isinstance(result, dict)
        assert result["espaloma"] is True
        assert result["xtb"] is False
        assert result["mace_ml"] is True
        assert len(result) == 3

    def test_preferred_order_is_set(self):
        """Factory has a preferred_order list with all three calculator names."""
        factory = self._make_factory()
        assert factory.preferred_order == ["mace_ml", "espaloma", "xtb"]


# ---------------------------------------------------------------------------
# TestDefaultMaceDipoleModel -- configuration constant tests
# ---------------------------------------------------------------------------


class TestDefaultMaceDipoleModel:
    """Tests for DEFAULT_MACE_DIPOLE_MODEL configuration constant."""

    def test_importable(self):
        """DEFAULT_MACE_DIPOLE_MODEL can be imported from mace_gaussian.calculators.mace_ml."""
        from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL

        assert isinstance(DEFAULT_MACE_DIPOLE_MODEL, str)

    def test_path_suffix(self):
        """Default model path ends with the expected suffix."""
        from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL

        assert DEFAULT_MACE_DIPOLE_MODEL.endswith("dipole_model/model_1.model")

    def test_env_var_override(self):
        """MACE_DIPOLE_MODEL_PATH env var overrides the default path."""
        import importlib

        import mace_gaussian.calculators.mace_ml as mace_ml_mod

        custom_path = "/tmp/custom_model.model"
        with patch.dict("os.environ", {"MACE_DIPOLE_MODEL_PATH": custom_path}):
            importlib.reload(mace_ml_mod)
            assert mace_ml_mod.DEFAULT_MACE_DIPOLE_MODEL == custom_path

        # Reload again to restore original value
        importlib.reload(mace_ml_mod)
