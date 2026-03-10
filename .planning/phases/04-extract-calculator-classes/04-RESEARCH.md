# Phase 4: Extract Calculator Classes - Research

**Researched:** 2026-02-17
**Domain:** Python class extraction, ABC pattern, factory pattern, module reorganization
**Confidence:** HIGH

## Summary

This phase extracts the dipole calculator class hierarchy (base class, 3 implementations, factory) from `gm_main.py` into a new `calculators/` package. The current code is self-contained in lines 86-317 of `gm_main.py`: the `DipoleCalculatorBase` ABC, `EspalomaDipoleCalculator`, `XTBDipoleCalculator`, `MACEMLDipoleCalculator`, `DipoleCalculatorFactory`, and the module-level `dipole_factory` singleton.

The extraction is straightforward because the calculator classes have minimal coupling to the rest of `gm_main.py`. The only dependency flowing inward is `MACEMLDipoleCalculator` importing `MACEDipoleCalculator` from `mace_calculators.py` and referencing the `DEFAULT_MACE_DIPOLE_MODEL` config constant. The only dependency flowing outward is that `gm_main.py` uses `dipole_factory.get_calculator()` in `run_next_calculation()` (line 635) and `dipole_factory.list_available()` in `print_diagnostics()` (line 1051). No other module in the codebase imports these classes directly.

**Primary recommendation:** Extract into `calculators/` package following the same flat layout pattern established by Phase 3's `utils/` package. Keep the `dipole_factory` singleton instantiation in `calculators/__init__.py` and re-export it, so `gm_main.py` changes from `dipole_factory = DipoleCalculatorFactory()` to `from calculators import dipole_factory`.

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions

All decisions are under Claude's Discretion for this phase. The following guidance from CONTEXT.md captures the agreed approach:

**Calculator boundary:**
- Extract: DipoleCalculatorBase (ABC), EspalomaCalculator, MaceMLCalculator, XtbCalculator, DipoleCalculatorFactory
- The `calculator()` function in gm_main.py that handles MACE standard model loading stays in gm_main.py -- it's deeply tied to the monkey-patching (Phase 5)
- Any helper that directly calls `cleanup_mace_modules()` stays in gm_main.py for now
- The calculators themselves import from mace/espaloma/xtb -- those imports move with them, but the sys.modules manipulation wrapper stays behind
- If a calculator's `__init__` or `calculate()` method calls cleanup, that call stays inline for now and Phase 5 will clean it up

**Factory design:**
- Keep the factory pattern simple -- a dict mapping or function that returns calculator instances
- No need to over-engineer (no plugin registry, no entry points) -- there are exactly 3 calculator types
- Factory should accept calculator name as string and return instantiated calculator
- Factory lives in `calculators/__init__.py` or `calculators/factory.py`

**Testing strategy:**
- Calculator tests should mock heavy dependencies (MACE models, espaloma, xtb)
- Test the factory pattern: correct class returned for each name, error on unknown name
- Test calculator interface compliance: all implementations have the required methods
- Integration tests (actually loading models) stay marked with @pytest.mark.gpu
- Use `unittest.mock.patch` for model loading, not actual model files

**Package layout:**
- `calculators/__init__.py` -- re-exports factory and base class
- `calculators/base.py` -- DipoleCalculatorBase ABC
- `calculators/espaloma.py` -- EspalomaCalculator
- `calculators/mace_ml.py` -- MaceMLCalculator
- `calculators/xtb.py` -- XtbCalculator (if it exists as a class)
- `calculators/factory.py` -- DipoleCalculatorFactory
- One class per file keeps things clean and matches the pattern established in Phase 3

### Claude's Discretion

User trusts Claude's expertise on all implementation decisions for this phase.

### Deferred Ideas (OUT OF SCOPE)

None -- discussion stayed within phase scope.

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-02 | DipoleCalculator classes extracted from gm_main.py into calculators/ package | Complete class inventory, dependency map, import graph, and extraction plan documented below |

</phase_requirements>

## Standard Stack

### Core

No new libraries needed. This phase uses only existing project infrastructure:

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Python ABC | stdlib | Abstract base class for calculator interface | Already used in gm_main.py |
| pytest | >=7.0.0 | Test extracted calculator modules | Already configured |
| unittest.mock | stdlib | Mock heavy ML dependencies in tests | Standard for testing classes with expensive dependencies |
| ruff | >=0.9.0 | Lint + format enforcement | Already configured in pyproject.toml |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | existing | Type hints and array operations in base class | Already a dependency |

### Alternatives Considered

None -- this is a pure extraction/refactoring operation with no technology choices to make.

## Architecture Patterns

### Recommended Package Structure

```
calculators/
    __init__.py          # Re-exports: DipoleCalculatorBase, DipoleCalculatorFactory,
                         #   EspalomaDipoleCalculator, XTBDipoleCalculator,
                         #   MACEMLDipoleCalculator, dipole_factory
    base.py              # DipoleCalculatorBase ABC
    espaloma.py          # EspalomaDipoleCalculator
    mace_ml.py           # MACEMLDipoleCalculator
    xtb.py               # XTBDipoleCalculator
    factory.py           # DipoleCalculatorFactory + dipole_factory singleton
```

### Pattern 1: ABC with Template Method (existing pattern, preserve)

**What:** `DipoleCalculatorBase` defines the interface with `calculate_dipole()` as abstract and `calculate_dipole_derivatives()` as a concrete template method using central differences.
**When to use:** Already the established pattern -- preserve as-is during extraction.
**Key detail:** The base class's `calculate_dipole_derivatives()` calls `self.calculate_dipole()` in a loop (6N calls for N atoms). This is the numerical differentiation workhorse.

```python
# calculators/base.py
from abc import ABC, abstractmethod
from typing import Optional, Tuple

import numpy as np

class DipoleCalculatorBase(ABC):
    """Abstract base class for dipole calculators"""

    def __init__(self, name: str):
        self.name = name
        self.available = None
        self._check_availability()

    @abstractmethod
    def _check_availability(self) -> bool:
        pass

    @abstractmethod
    def calculate_dipole(self, atoms, **kwargs) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        pass

    def calculate_dipole_derivatives(self, atoms, displacement=0.01, **kwargs) -> np.ndarray:
        # Existing implementation moves here unchanged
        ...
```

### Pattern 2: Factory with Simple Dict Registry (existing pattern, preserve)

**What:** `DipoleCalculatorFactory` eagerly instantiates all calculators, stores them in a dict, provides `get_calculator(method)` with "auto" selection by preferred order.
**When to use:** Already the established pattern -- preserve as-is.
**Key detail:** The factory eagerly instantiates ALL calculators at import time (via `_register_calculators()`). This means importing the factory triggers availability checks for espaloma, xtb, and MACE ML. This is intentional -- it front-loads failure detection.

```python
# calculators/factory.py
from calculators.base import DipoleCalculatorBase
from calculators.espaloma import EspalomaDipoleCalculator
from calculators.mace_ml import MACEMLDipoleCalculator
from calculators.xtb import XTBDipoleCalculator

class DipoleCalculatorFactory:
    def __init__(self):
        self.calculators = {}
        self.preferred_order = ["mace_ml", "espaloma", "xtb"]
        self._register_calculators()

    def _register_calculators(self):
        calculators = [
            EspalomaDipoleCalculator(),
            XTBDipoleCalculator(),
            MACEMLDipoleCalculator(),
        ]
        for calc in calculators:
            self.calculators[calc.name] = calc
    ...

# Module-level singleton
dipole_factory = DipoleCalculatorFactory()
```

### Pattern 3: Config Constants Migration

**What:** `DEFAULT_MACE_DIPOLE_MODEL` is currently defined at the top of `gm_main.py` and referenced by `MACEMLDipoleCalculator.__init__()`.
**Recommendation:** Move this constant into `calculators/mace_ml.py` since it is only used by `MACEMLDipoleCalculator`. It reads from `MACE_DIPOLE_MODEL_PATH` env var with a file path fallback. This keeps the calculator self-contained.

```python
# calculators/mace_ml.py
import os
from pathlib import Path

DEFAULT_MACE_DIPOLE_MODEL = os.getenv(
    "MACE_DIPOLE_MODEL_PATH",
    str(Path.home() / "mace_gaussian" / "dipole_model" / "model_1.model"),
)
```

### Anti-Patterns to Avoid

- **Circular imports:** `factory.py` imports all calculator implementations. Implementations must NOT import from `factory.py`. The `__init__.py` re-exports from all submodules but submodules should not import from `__init__.py`.
- **Lazy factory instantiation:** Don't defer the `dipole_factory` singleton creation -- the current eager initialization pattern is intentional for fail-fast behavior.
- **Breaking the monkey-patch boundary:** `MACEMLDipoleCalculator` uses `MACEDipoleCalculator` from `mace_calculators.py`, which internally calls `cleanup_mace_modules()`. This import chain must stay intact. Do NOT try to move `MACEDipoleCalculator` (the wrapper in `mace_calculators.py`) into the calculators package -- that belongs to Phase 5.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Calculator discovery | Plugin/entrypoint system | Simple dict in factory | Only 3 calculators, YAGNI |
| Config management | Config file parser | os.getenv with defaults | Single env var, already works |
| Availability detection | Separate probe system | Try-import in _check_availability | Existing pattern, proven reliable |

**Key insight:** The calculator hierarchy is simple (3 implementations, 1 factory). Keep it simple.

## Common Pitfalls

### Pitfall 1: Module-Level Side Effects from Eager Factory Instantiation

**What goes wrong:** When `calculators/__init__.py` imports `factory.py`, and `factory.py` creates the `dipole_factory` singleton at module level, all calculator `__init__` methods run, which trigger `_check_availability()`. This includes trying to import `espaloma_charge`, `xtb`, and loading the MACE dipole model.
**Why it happens:** Python executes module-level code on first import.
**How to avoid:** This is actually the DESIRED behavior (fail-fast). But tests must be aware: importing `from calculators import DipoleCalculatorFactory` will trigger these checks. Tests that don't need the factory should import individual classes directly: `from calculators.base import DipoleCalculatorBase`.
**Warning signs:** Tests failing because espaloma/xtb/MACE aren't available in the test environment.

### Pitfall 2: The `dipole_factory` Singleton Reference

**What goes wrong:** `gm_main.py` currently creates `dipole_factory` at module level (line 317). After extraction, `gm_main.py` must import it from `calculators`. If the import path is wrong, the singleton won't exist and `run_next_calculation()` (line 635) and `print_diagnostics()` (line 1051) will fail.
**How to avoid:** Update exactly two references in `gm_main.py`:
  1. Remove line 317 (`dipole_factory = DipoleCalculatorFactory()`)
  2. Add `from calculators import dipole_factory` to imports
  3. Verify `run_next_calculation` and `print_diagnostics` still work

### Pitfall 3: The `BOHR_TO_ANGSTROM` Import in Base Class

**What goes wrong:** `EspalomaDipoleCalculator.calculate_dipole()` uses `BOHR_TO_ANGSTROM` (imported in `gm_main.py` from `utils.units`). When extracted, the espaloma module needs its own import.
**How to avoid:** Add `from utils.units import BOHR_TO_ANGSTROM` to `calculators/espaloma.py`.

### Pitfall 4: Renaming Classes vs Preserving Names

**What goes wrong:** The CONTEXT.md mentions `EspalomaCalculator` but the actual class is `EspalomaDipoleCalculator`. If the class is renamed, the factory registration breaks because `calc.name` returns `"espaloma"` and is used as a dict key throughout the pipeline (CLI args like `--dipole-calculators espaloma`).
**How to avoid:** Preserve the exact class names from gm_main.py: `EspalomaDipoleCalculator`, `XTBDipoleCalculator`, `MACEMLDipoleCalculator`. The `.name` attribute (set in `__init__`) is what matters for the factory, but renaming classes adds unnecessary churn. If CONTEXT.md shortnames are desired, add them as aliases in `__init__.py`.

### Pitfall 5: `MACEDipoleCalculator` Import in `mace_ml.py`

**What goes wrong:** `MACEMLDipoleCalculator._check_availability()` calls `self.mace_calc = MACEDipoleCalculator(self.model_path)`. This import comes from `mace_calculators.py`. When the class moves to `calculators/mace_ml.py`, the import statement must be `from mace_calculators import MACEDipoleCalculator`.
**How to avoid:** Verify this import works from the new location. The import is at the module level in current `gm_main.py` (line 34), but in the extracted version it should remain in the class -- either at the top of `mace_ml.py` or lazily in `_check_availability()`. Lazy import is safer because `mace_calculators.py` has side effects (backs up `sys.modules["mace.modules.models"]`).

## Code Examples

### Example 1: `calculators/__init__.py` (re-export pattern matching utils/)

```python
"""Dipole calculator package for MACE-Gaussian interface."""

from calculators.base import DipoleCalculatorBase
from calculators.espaloma import EspalomaDipoleCalculator
from calculators.factory import DipoleCalculatorFactory, dipole_factory
from calculators.mace_ml import MACEMLDipoleCalculator
from calculators.xtb import XTBDipoleCalculator

__all__ = [
    "DipoleCalculatorBase",
    "DipoleCalculatorFactory",
    "EspalomaDipoleCalculator",
    "MACEMLDipoleCalculator",
    "XTBDipoleCalculator",
    "dipole_factory",
]
```

### Example 2: Updated imports in `gm_main.py`

```python
# REMOVE these lines from gm_main.py:
#   from abc import ABC, abstractmethod  (if only used by calculators)
#   from mace_calculators import MACEDipoleCalculator  (moved to mace_ml.py)
#   class DipoleCalculatorBase(ABC): ... (lines 86-139)
#   class EspalomaDipoleCalculator(...): ... (lines 142-205)
#   class XTBDipoleCalculator(...): ... (lines 208-238)
#   class MACEMLDipoleCalculator(...): ... (lines 241-267)
#   class DipoleCalculatorFactory: ... (lines 270-313)
#   dipole_factory = DipoleCalculatorFactory()  (line 317)

# ADD this import:
from calculators import dipole_factory
```

### Example 3: Test pattern for calculator interface compliance

```python
# tests/test_calculators.py
import pytest
from unittest.mock import patch, MagicMock
from calculators.base import DipoleCalculatorBase
from calculators.espaloma import EspalomaDipoleCalculator
from calculators.mace_ml import MACEMLDipoleCalculator
from calculators.xtb import XTBDipoleCalculator


class TestCalculatorInterface:
    """All implementations satisfy the DipoleCalculatorBase contract."""

    @pytest.mark.parametrize("cls", [
        EspalomaDipoleCalculator,
        MACEMLDipoleCalculator,
        XTBDipoleCalculator,
    ])
    def test_is_subclass(self, cls):
        assert issubclass(cls, DipoleCalculatorBase)

    @pytest.mark.parametrize("cls", [
        EspalomaDipoleCalculator,
        MACEMLDipoleCalculator,
        XTBDipoleCalculator,
    ])
    def test_has_required_methods(self, cls):
        assert hasattr(cls, 'calculate_dipole')
        assert hasattr(cls, 'calculate_dipole_derivatives')
        assert hasattr(cls, '_check_availability')


class TestFactory:
    """Factory correctly maps names to calculator instances."""

    @patch('calculators.factory.EspalomaDipoleCalculator')
    @patch('calculators.factory.XTBDipoleCalculator')
    @patch('calculators.factory.MACEMLDipoleCalculator')
    def test_get_known_calculator(self, mock_mace, mock_xtb, mock_esp):
        # Setup mocks to simulate available calculators
        for mock_cls in [mock_mace, mock_xtb, mock_esp]:
            instance = mock_cls.return_value
            instance.available = True
            instance.name = mock_cls.__name__.lower()
        ...

    def test_get_unknown_calculator_raises(self):
        from calculators.factory import DipoleCalculatorFactory
        factory = DipoleCalculatorFactory()
        with pytest.raises(ValueError, match="Unknown"):
            factory.get_calculator("nonexistent")
```

## Dependency Map (Current State)

```
gm_main.py
  |- DipoleCalculatorBase (ABC)          <- lines 86-139
  |    |- EspalomaDipoleCalculator       <- lines 142-205
  |    |    `- imports: espaloma_charge, rdkit, torch, utils.units.BOHR_TO_ANGSTROM
  |    |- XTBDipoleCalculator            <- lines 208-238
  |    |    `- imports: xtb.ase.calculator
  |    `- MACEMLDipoleCalculator         <- lines 241-267
  |         `- imports: mace_calculators.MACEDipoleCalculator, DEFAULT_MACE_DIPOLE_MODEL
  |- DipoleCalculatorFactory             <- lines 270-313
  |    `- imports: all 3 calculator classes above
  `- dipole_factory (singleton)          <- line 317
       `- used by: run_next_calculation (line 635), print_diagnostics (line 1051)

External consumers of these classes:
  - cli.py: imports run_workflow and print_diagnostics from gm_main (indirect)
  - tests/test_cli_validation.py: imports GAUSSIAN_TIMEOUT_SECONDS from gm_main (unrelated)
  - No module directly imports the calculator classes from gm_main
```

## Import Changes Required

| File | Current Import | New Import |
|------|---------------|------------|
| `gm_main.py` | `from mace_calculators import MACEDipoleCalculator` | Remove (moves to `calculators/mace_ml.py`) |
| `gm_main.py` | `from abc import ABC, abstractmethod` | Remove (if unused after extraction) |
| `gm_main.py` | (inline) `dipole_factory = DipoleCalculatorFactory()` | `from calculators import dipole_factory` |
| `calculators/espaloma.py` | (new file) | `from utils.units import BOHR_TO_ANGSTROM` |
| `calculators/mace_ml.py` | (new file) | `from mace_calculators import MACEDipoleCalculator` |
| `calculators/factory.py` | (new file) | `from calculators.base import DipoleCalculatorBase` + all 3 implementations |

## State of the Art

No technology changes needed. This is pure code organization.

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| All calculator code in gm_main.py (1285 lines) | Extract to calculators/ package | This phase | gm_main.py drops ~230 lines, gains modularity |

## Open Questions

1. **Should `dipole_factory` singleton stay in `calculators/__init__.py` or `calculators/factory.py`?**
   - What we know: Phase 3 put `ResultsManager` in `utils/results.py` and re-exported from `utils/__init__.py`. Same pattern works here: singleton in `factory.py`, re-exported from `__init__.py`.
   - Recommendation: Put singleton in `factory.py`, re-export from `__init__.py`. This avoids triggering all calculator imports when someone only needs `from calculators.base import DipoleCalculatorBase`.

2. **Should `ABC` import remain in `gm_main.py`?**
   - What we know: After removing the calculator classes, `gm_main.py` may still use `ABC` if any other abstract class exists there. Currently none do.
   - Recommendation: Remove the `ABC` import from `gm_main.py` if it becomes unused. Ruff will flag this automatically.

3. **Testing the `_check_availability` side effects**
   - What we know: `_check_availability()` is called in `__init__()`, which tries importing heavy deps (espaloma, xtb, MACE). Tests need to either mock these or accept they may not be available.
   - Recommendation: For unit tests, patch `_check_availability` to avoid import side effects. Test it separately with explicit mocks of the underlying libraries.

## Sources

### Primary (HIGH confidence)

- Direct code inspection of `gm_main.py` (lines 86-317) -- complete class hierarchy
- Direct code inspection of `mace_calculators.py` -- MACEDipoleCalculator wrapper and cleanup
- Direct code inspection of `cli.py` -- only consumer via `gm_main.run_workflow` and `gm_main.print_diagnostics`
- Grep across all `*.py` files -- verified no other module imports calculator classes directly
- Phase 3 `utils/` package structure -- established conventions for this project

### Secondary (MEDIUM confidence)

- None needed -- this is internal refactoring with no external library decisions

### Tertiary (LOW confidence)

- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, pure extraction
- Architecture: HIGH -- direct code inspection of all classes and their dependencies
- Pitfalls: HIGH -- every pitfall identified by tracing actual import chains in the codebase

**Research date:** 2026-02-17
**Valid until:** 2026-03-17 (stable -- no external dependencies changing)
