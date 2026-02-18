# Phase 5: Replace MACE Module Monkey-Patching - Research

**Researched:** 2026-02-18
**Domain:** Python module loading, PyTorch model deserialization, MACE ML potential packages
**Confidence:** HIGH

<user_constraints>

## User Constraints (from CONTEXT.md)

### Claude's Discretion

User trusts Claude's expertise on all implementation decisions for this phase. The following guidance captures the design intent and constraints for downstream agents:

**Loading strategy:**
- The core problem: standard MACE (`mace.modules.models`) and dipole MACE (`mace_dipole_core.modules.models`) collide in `sys.modules`. Current code hot-swaps the cache entry.
- Research should evaluate: (1) `importlib.util` isolated loading (load dipole variant without registering in `sys.modules`), (2) namespace/import tricks to avoid the collision entirely, (3) subprocess isolation as a fallback option.
- Priority order: simplest approach that works with CUDA first. `importlib.util` is the likely winner but needs verification that CUDA contexts survive isolated loading.
- Subprocess isolation is a last resort -- adds latency and GPU memory serialization complexity.
- The `fake_module_from_real()` function (shallow module copy) should be eliminated, not just replaced.

**Migration boundary:**
- `mace_calculators.py` should be absorbed into the `calculators/` package (likely as `calculators/mace_loader.py` or similar) since Phase 4 already extracted the calculator classes there.
- `MACEDipoleCalculator` wrapper class moves into calculators/ -- it's the natural home after Phase 4's extraction.
- `load_standard_mace_calculator()` and `load_dipole_mace_calculator()` get replaced with the new safe loading mechanism.
- `cleanup_mace_modules()` is deleted entirely -- the whole point is eliminating the need for cleanup.
- All `finally` blocks that call `cleanup_mace_modules()` are removed.
- The `calculator()` function in `gm_main.py` (lines 444-467) that loads standard MACE for energy calculations should also use the new safe loader.
- Update CLAUDE.md to remove the monkey-patching warning.

**Fallback behavior:**
- No fallback to monkey-patching. Clean break -- if the new loading doesn't work, it's a bug to fix, not a case to handle.
- Clear error messages if either MACE variant fails to load (which model, what went wrong, is it installed?).
- CUDA device placement must be explicit -- both models on the same GPU, with clear logging of device assignment.

**Testing approach:**
- Mock the actual model loading (torch, MACE internals) to test the loading mechanism in isolation.
- Verify: standard MACE and dipole MACE can be loaded in the same process without `sys.modules` conflicts.
- Verify: no `cleanup_mace_modules()` calls remain anywhere in the codebase.
- Verify: repeated dipole calculations don't accumulate state or leak memory.
- Integration tests (actually loading models) stay marked with `@pytest.mark.gpu`.
- The existing pipeline must still work with all 4 model combinations on water.

### Deferred Ideas (OUT OF SCOPE)

- Batching dipole derivative calculations (computing all 6N displacements in one GPU call) -- potential performance optimization but separate from fixing the loading mechanism. Could be a future phase or post-thesis work.

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-03 | MACE module monkey-patching replaced with safe loading pattern (lazy imports or process isolation) | Primary recommendation: `pickle_module` remapping in `torch.load`. Verified working -- see "Recommended Approach" section. No `sys.modules` mutation, no cleanup needed, standard MACE untouched. |

</phase_requirements>

## Summary

The monkey-patching in `mace_calculators.py` exists because the dipole model file (`dipole_model/model_1.model`) was saved with class reference `mace.modules.models.AtomicDielectricMACE`, but the actual class implementation needed at runtime is from the dipole fork (`mace_dipole_core.modules.models.AtomicDielectricMACE`). The standard mace-torch 0.3.14 package has its own `AtomicDielectricMACE` class, but its `forward()` method has different logic (requires polarizability, different tensor slicing, different helper functions) that is incompatible with models trained using the dipole fork.

The current code solves this by replacing `sys.modules["mace.modules.models"]` with a shallow copy of `mace_dipole_core.modules.models` before loading the dipole calculator, then cleaning up afterward. This is fragile because: (1) it mutates global state on every dipole calculation, (2) cleanup is mandatory or standard MACE breaks, (3) `finally` blocks are required everywhere, (4) `mace_dipole_core` must be on `sys.path` (currently requires `mace_dipole_pkg/` on the path).

I verified three replacement approaches experimentally. The recommended approach uses `torch.load`'s `pickle_module` parameter to intercept class resolution at deserialization time, remapping `mace.modules.models.*` to `mace_dipole_core.modules.models.*`. This was tested end-to-end: model loads correctly, dipole calculation produces correct results (`[0.21151979, -0.46528312, 0.0]` for water), and `sys.modules["mace.modules.models"]` remains completely untouched.

**Primary recommendation:** Use a custom `pickle_module` with `RemappingUnpickler` that redirects `mace.modules.models` class lookups to `mace_dipole_core.modules.models` during `torch.load`. This is the simplest, safest, and most surgical approach -- it requires no `sys.modules` manipulation, no cleanup, and no subprocess isolation.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| torch | 2.8.0+cu128 | Model loading via `torch.load` | Already installed; `pickle_module` param is the official extension point |
| mace-torch | 0.3.14 | Standard MACE energy calculators | Already installed in venv |
| mace_dipole_core | local fork | Dipole-enabled MACE (AtomicDielectricMACE) | Local package at `mace_dipole_pkg/mace_dipole_core/` |
| pickle (stdlib) | 3.12 | Custom `Unpickler.find_class` for class remapping | Standard Python; no extra deps |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| importlib (stdlib) | 3.12 | `importlib.import_module` for lazy loading of `mace_dipole_core` | Loading the dipole fork's module on first use |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `pickle_module` remapping | One-time class registration on `mace.modules.models` | Simpler but mutates the standard module's namespace permanently; `AtomicDielectricMACE` attribute on `mace.modules.models` would point to the wrong class |
| `pickle_module` remapping | `importlib.util.spec_from_file_location` isolated loading | More complex; the dipole module has deep dependency chains (`e3nn`, `torch`, `mace_dipole_core.tools.*`) that make truly isolated loading difficult |
| `pickle_module` remapping | Subprocess isolation | Last resort; adds latency, GPU memory duplication, IPC complexity |

## Architecture Patterns

### Recommended Project Structure

After this phase:
```
calculators/
    __init__.py          # Package exports
    base.py              # DipoleCalculatorBase (ABC)
    espaloma.py          # EspalomaDipoleCalculator
    mace_ml.py           # MACEMLDipoleCalculator (uses mace_loader)
    mace_loader.py       # NEW: Safe MACE loading (replaces mace_calculators.py)
    xtb.py               # XTBDipoleCalculator
    factory.py           # DipoleCalculatorFactory
mace_calculators.py      # DELETED
```

### Pattern 1: Custom pickle_module for torch.load Class Remapping

**What:** A custom pickle module class that overrides `Unpickler.find_class` to redirect class lookups from `mace.modules.models` to `mace_dipole_core.modules.models` during model deserialization.

**When to use:** When loading a PyTorch model whose saved class path differs from the runtime import path (common with forked packages).

**Verified working code:**
```python
import pickle
import sys
import torch

class _DipoleModelUnpickler(pickle.Unpickler):
    """Unpickler that remaps mace.modules.models classes to mace_dipole_core equivalents.

    The dipole model was saved with class path mace.modules.models.AtomicDielectricMACE,
    but the correct implementation lives in mace_dipole_core.modules.models. This unpickler
    intercepts class resolution during torch.load to use the dipole fork's classes.
    """
    def find_class(self, module: str, name: str):
        if module == "mace.modules.models":
            # Lazy import to avoid loading heavy deps until needed
            import mace_dipole_core.modules.models as dipole_models
            cls = getattr(dipole_models, name, None)
            if cls is not None:
                return cls
        return super().find_class(module, name)


class _DipolePickleModule:
    """Pickle module wrapper that uses _DipoleModelUnpickler for class resolution."""
    Unpickler = _DipoleModelUnpickler

    def __getattr__(self, name):
        return getattr(pickle, name)


def load_dipole_model(model_path: str, device: str = "cuda") -> torch.nn.Module:
    """Load a MACE dipole model using safe class remapping (no sys.modules mutation)."""
    return torch.load(
        f=model_path,
        map_location=device,
        weights_only=False,
        pickle_module=_DipolePickleModule(),
    )
```

**Source:** Verified experimentally on this codebase (2026-02-18). `torch.load` accepts `pickle_module` parameter since PyTorch 1.x. Confirmed working with torch 2.8.0.

### Pattern 2: MACEDipoleCalculator Absorption into calculators/mace_loader.py

**What:** Move `MACEDipoleCalculator` from `mace_calculators.py` into `calculators/mace_loader.py`, replacing internal `load_dipole_mace_calculator()` call with the new safe `load_dipole_model()` function. The calculator uses the dipole fork's `MACECalculator_dipole` for the ASE calculator interface, but loads the model via the safe remapping.

**Key insight:** The `MACECalculator_dipole` class from `mace_dipole_core.calculators.mace` calls `torch.load` internally in its constructor (line 139). To use our safe loading, we have two options:
1. **Patch `torch.load` on the dipole calculator module** temporarily during construction (tested, works)
2. **Pre-load the model and inject it** into the calculator (requires modifying the calculator constructor or using it differently)

Option 1 is simpler and was verified working:
```python
import mace_dipole_core.calculators.mace as dipole_calc_module

original_torch_load = torch.load

def _safe_torch_load(*args, **kwargs):
    kwargs["pickle_module"] = _DipolePickleModule()
    kwargs["weights_only"] = False
    return original_torch_load(*args, **kwargs)

# Temporarily patch during construction only
dipole_calc_module.torch.load = _safe_torch_load
calc = MACECalculator_dipole(model_paths=path, device=device, ...)
dipole_calc_module.torch.load = original_torch_load  # Restore immediately
```

This is a localized, scoped patch (not a global `sys.modules` mutation) and can be wrapped in a context manager.

### Pattern 3: Standard MACE Loading (No Changes Needed)

**What:** The `calculator()` function in `gm_main.py` (lines 444-467) loads standard MACE models for energy calculations using `from mace.calculators import mace_mp` etc. These already work correctly with no monkey-patching because they use the installed `mace-torch` package directly.

**When it matters:** The context document says this function should "use the new safe loader," but in practice it already uses standard imports and needs no changes. The safe loader is only needed for dipole models.

### Anti-Patterns to Avoid

- **Replacing the entire `sys.modules["mace.modules.models"]`:** This is what the current code does. Any replacement of the entire module object in sys.modules creates global state that requires cleanup and can silently corrupt other code.
- **`fake_module_from_real()` shallow copies:** Creates a types.ModuleType with attributes copied from another module. This breaks `isinstance` checks, module identity comparisons, and is fragile when the source module has lazy attributes or descriptors.
- **Global one-time class registration on standard module:** Replacing `mace.modules.models.AtomicDielectricMACE` permanently with the dipole fork version seems simple but creates version confusion -- the standard MACE's class would be silently shadowed, and any future code that expects the standard version would get the fork's version.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Class remapping during deserialization | Custom pickle stream parsing | `torch.load(pickle_module=...)` | PyTorch's official extension point; handles all edge cases (storage mapping, device placement, etc.) |
| Module isolation | Custom importlib machinery | Scoped `torch.load` patching on the dipole calc module | The dipole fork's deep dependency chain makes true module isolation impractical |
| CUDA device management | Custom device tracking | `torch.load(map_location=device)` + `model.to(device)` | Already handled by PyTorch |

**Key insight:** The problem looks like a module loading problem, but it's actually a **deserialization class resolution problem**. The fix belongs at the `torch.load` layer, not the `import` layer.

## Common Pitfalls

### Pitfall 1: mace_dipole_core Not on sys.path

**What goes wrong:** `import mace_dipole_core` fails with `ModuleNotFoundError` because the package lives at `mace_dipole_pkg/mace_dipole_core/` and `mace_dipole_pkg/` is not on `sys.path`.

**Why it happens:** The project root (`/home/mot/mace_gaussian`) is on `sys.path` via a `.pth` file, but `mace_dipole_pkg/` is a subdirectory not itself on the path. The package name is `mace_dipole_core`, not `mace_dipole_pkg.mace_dipole_core`.

**How to avoid:** The new `calculators/mace_loader.py` must add `mace_dipole_pkg/` to `sys.path` before importing `mace_dipole_core`. Do this once, lazily, on first dipole calculator use:
```python
import sys
from pathlib import Path

_DIPOLE_PKG_DIR = str(Path(__file__).resolve().parent.parent / "mace_dipole_pkg")

def _ensure_dipole_importable():
    if _DIPOLE_PKG_DIR not in sys.path:
        sys.path.insert(0, _DIPOLE_PKG_DIR)
```
**Warning signs:** `ModuleNotFoundError: No module named 'mace_dipole_core'` at runtime.

### Pitfall 2: Standard vs Dipole AtomicDielectricMACE Forward Incompatibility

**What goes wrong:** Model loads successfully but `forward()` raises `ValueError: Polarizability is not used in this model, but it is required for the AtomicDielectricMACE`.

**Why it happens:** Standard `mace.modules.models.AtomicDielectricMACE` (mace-torch 0.3.14) has a different `forward()` than the dipole fork's version. Key differences:
- Standard: requires `use_polarizability=True` (raises ValueError otherwise)
- Standard: uses `cutoff` parameter in interaction blocks
- Standard: indexes dipoles at `node_out[:, 2:5]`
- Fork: supports dipole-only mode (`use_polarizability=False, use_dipole=True`)
- Fork: no `cutoff` parameter in interaction blocks
- Fork: indexes dipoles at `node_out[:, 1:4]`
- Standard: calls `compute_fixed_charge_dipole_polar`
- Fork: calls `compute_fixed_charge_dipole`

**How to avoid:** The `pickle_module` remapping MUST be used when loading dipole models. Simply loading with standard `torch.load` will use the wrong `AtomicDielectricMACE` class.

**Warning signs:** Model loads without error but produces wrong results or crashes in `forward()`.

### Pitfall 3: TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD Environment Variable

**What goes wrong:** `torch.load` fails with `_pickle.UnpicklingError: Weights only load failed` when loading the dipole model.

**Why it happens:** PyTorch 2.6+ defaults `weights_only=True` for `torch.load`. The dipole model needs `weights_only=False` because it contains full model objects (not just state dicts). The current environment uses `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1` as a workaround.

**How to avoid:** Always pass `weights_only=False` explicitly in the safe loader function. Do NOT rely on the environment variable.

**Warning signs:** Unpickling errors mentioning "Weights only load failed" or "GLOBAL ... was not an allowed global."

### Pitfall 4: Forgetting to Remove All Cleanup References

**What goes wrong:** Dead `cleanup_mace_modules()` calls remain in the codebase, or the `mace_calculators` mock in tests becomes stale.

**Why it happens:** References to the old module exist in:
- `calculators/mace_ml.py` line 30: `from mace_calculators import MACEDipoleCalculator`
- `tests/test_calculators.py` line 32: `"mace_calculators"` in heavy deps mock list
- `CLAUDE.md`: monkey-patching warning
- `pyproject.toml` line 86: `mace_calculators` in coverage sources
- `pyproject.toml` line 65: `mace_calculators` in known-first-party

**How to avoid:** After migration, grep for `mace_calculators`, `cleanup_mace`, `fake_module_from_real`, `load_standard_mace`, `load_dipole_mace` across the entire codebase. Every hit should be addressed.

### Pitfall 5: Test Pre-Mock List Needs Updating

**What goes wrong:** `tests/test_calculators.py` pre-mocks `mace_calculators` in `sys.modules` (line 32). After this phase, that module no longer exists, but the import chain changes.

**How to avoid:** Update the test pre-mock list to mock `mace_dipole_core` (or whatever the new import target is in `calculators/mace_loader.py`) instead of `mace_calculators`.

## Code Examples

### Complete Safe Loader Module (calculators/mace_loader.py)

```python
"""Safe MACE model loading without sys.modules monkey-patching.

Replaces the old mace_calculators.py approach of swapping sys.modules entries.
Uses torch.load's pickle_module parameter for class remapping during deserialization.
"""

import logging
import pickle
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import torch

logger = logging.getLogger(__name__)

# Path to mace_dipole_pkg for lazy sys.path setup
_DIPOLE_PKG_DIR = str(Path(__file__).resolve().parent.parent / "mace_dipole_pkg")


def _ensure_dipole_importable():
    """Add mace_dipole_pkg to sys.path if not already present."""
    if _DIPOLE_PKG_DIR not in sys.path:
        sys.path.insert(0, _DIPOLE_PKG_DIR)


class _DipoleModelUnpickler(pickle.Unpickler):
    """Remaps mace.modules.models -> mace_dipole_core.modules.models during deserialization."""

    def find_class(self, module: str, name: str):
        if module == "mace.modules.models":
            _ensure_dipole_importable()
            import mace_dipole_core.modules.models as dipole_models
            cls = getattr(dipole_models, name, None)
            if cls is not None:
                return cls
        return super().find_class(module, name)


class _DipolePickleModule:
    """Drop-in pickle module replacement for torch.load."""
    Unpickler = _DipoleModelUnpickler

    def __getattr__(self, name):
        return getattr(pickle, name)


def load_dipole_calculator(model_path: str, device: str = "cuda",
                           model_type: str = "DipolePolarizabilityMACE",
                           default_dtype: str = "float64"):
    """Load a MACE dipole calculator with safe class remapping.

    Returns a MACECalculator_dipole instance ready for ASE calculations.
    No sys.modules cleanup needed.
    """
    _ensure_dipole_importable()
    from mace_dipole_core.calculators.mace import MACECalculator_dipole
    import mace_dipole_core.calculators.mace as dipole_calc_module

    # Temporarily redirect torch.load in the dipole calculator module
    # to use our class-remapping pickle module
    original_load = dipole_calc_module.torch.load

    def _safe_load(*args, **kwargs):
        kwargs["pickle_module"] = _DipolePickleModule()
        kwargs["weights_only"] = False
        return original_load(*args, **kwargs)

    try:
        dipole_calc_module.torch.load = _safe_load
        calc = MACECalculator_dipole(
            model_paths=model_path,
            device=device,
            model_type=model_type,
            default_dtype=default_dtype,
        )
    finally:
        dipole_calc_module.torch.load = original_load

    return calc
```

### Updated MACEDipoleCalculator (absorbed into calculators/mace_ml.py or mace_loader.py)

```python
class MACEDipoleCalculator:
    """Wrapper for MACE dipole calculator with safe loading."""

    def __init__(self, model_path, device="cuda"):
        self.model_path = model_path
        self.device = device
        self.calc = None

    def _ensure_calculator(self):
        if self.calc is None:
            from calculators.mace_loader import load_dipole_calculator
            self.calc = load_dipole_calculator(
                model_path=self.model_path,
                device=self.device,
            )

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole moment using MACE dipole model."""
        self._ensure_calculator()
        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc
        dipole_moment = atoms_copy.get_dipole_moment()

        # Handle different return formats (same as current code)
        if isinstance(dipole_moment, tuple):
            dipole_moment = dipole_moment[0]
        if isinstance(dipole_moment, torch.Tensor):
            dipole_moment = dipole_moment.detach().cpu().numpy()
        dipole_moment = np.atleast_1d(np.array(dipole_moment))
        if dipole_moment.shape == (1, 3):
            dipole_moment = dipole_moment.flatten()
        elif dipole_moment.shape != (3,):
            dipole_moment = dipole_moment.reshape(-1)[:3]

        return dipole_moment, None
        # NOTE: No cleanup_mace_modules() call needed!
```

## State of the Art

| Old Approach (Current) | New Approach (This Phase) | Impact |
|------------------------|--------------------------|--------|
| `sys.modules["mace.modules.models"] = fake_module(...)` on every call | `torch.load(pickle_module=_DipolePickleModule())` once at init | Eliminates global state mutation entirely |
| `cleanup_mace_modules()` in finally blocks | No cleanup needed | Removes fragile cleanup requirement |
| `fake_module_from_real()` shallow module copy | `_DipoleModelUnpickler.find_class()` override | Surgical class remapping vs wholesale module replacement |
| `mace_calculators.py` standalone module | `calculators/mace_loader.py` in package | Consistent with Phase 4 calculator package structure |
| Module reload on every load (`importlib.reload`) | No reloads needed | Simpler, faster |

## Open Questions

1. **Should `mace_dipole_pkg/` be installed as an editable package instead of sys.path manipulation?**
   - What we know: Currently not installed; relies on sys.path insertion. Installing via `pip install -e mace_dipole_pkg/` would make `mace_dipole_core` importable without path manipulation.
   - What's unclear: Whether editable install would conflict with the standard `mace-torch` package (both define similar module structures but under different names).
   - Recommendation: Keep sys.path approach for now (minimal disruption). Can revisit packaging in a future phase.

2. **What happens if the dipole model is re-saved with the standard mace-torch 0.3.14?**
   - What we know: The saved model uses `mace.modules.models.AtomicDielectricMACE` as the class path. If re-saved with standard MACE, it would use the standard class which has different forward() logic.
   - Recommendation: Document that dipole models must be trained/saved using the `mace_dipole_core` fork. The remapping in the loader handles the class path translation transparently.

## Sources

### Primary (HIGH confidence)
- Direct code inspection of `mace_calculators.py` -- full monkey-patching mechanism analyzed
- Direct code inspection of `mace_dipole_pkg/mace_dipole_core/modules/models.py` (lines 617-905) -- AtomicDielectricMACE fork version
- Direct code inspection of installed `mace.modules.models` (mace-torch 0.3.14) -- standard AtomicDielectricMACE
- **Experimental verification** (2026-02-18): Tested all three approaches in the project's Python environment:
  - `pickle_module` remapping: WORKS (full dipole calculation produces correct results)
  - One-time class registration: WORKS for loading but uses wrong forward() method
  - Direct torch.load without patching: Model loads but forward() raises ValueError
- `torch.load` documentation: `pickle_module` parameter confirmed since PyTorch 1.x, used for custom deserialization

### Secondary (MEDIUM confidence)
- `mace_dipole_pkg/mace_dipole_core/calculators/mace.py` (lines 128-142) -- torch.load call site in MACECalculator_dipole constructor
- `calculators/mace_ml.py` -- current import of MACEDipoleCalculator from mace_calculators
- `.pth` file at `.venv/lib/python3.12/site-packages/_mace_gaussian_interface.pth` -- confirms `/home/mot/mace_gaussian` on sys.path

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already in use; only the loading mechanism changes
- Architecture: HIGH -- experimentally verified the recommended approach produces correct results
- Pitfalls: HIGH -- discovered through hands-on testing (the forward() incompatibility was found by actually running the code)

**Research date:** 2026-02-18
**Valid until:** Indefinitely (core Python pickle/importlib and PyTorch torch.load APIs are stable)
