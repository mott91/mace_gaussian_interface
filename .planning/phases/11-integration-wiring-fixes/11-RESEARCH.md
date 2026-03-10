# Phase 11: Integration Wiring Fixes - Research

**Researched:** 2026-02-27
**Domain:** Python wiring gaps — warnings module, env var resolution, export plumbing, utils re-exports
**Confidence:** HIGH

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| REPR-02 | Calculation parameters captured alongside results | `workflow.py` calls `save_frequency_results()` without `calculation_parameters`; fix: assemble dict from `energy_calc_name` + `dipole_calc_name` + model constants and pass it in |
| ERR-03 | CUDA availability checked at startup with clear warning on CPU fallback | `detect_device()` uses `logger.warning` but never calls `warnings.warn(..., CUDANotAvailableWarning)`; fix: add the `warnings.warn` call alongside the logger call |
| ERR-04 | Env vars (MACE_DIPOLE_MODEL_PATH, MACE_HELPER_SCRIPT_PATH) validated at startup | `cli.py` calls `validate_all_prerequisites(dipole_model_path=None, helper_script_path=None)`; fix: resolve `os.getenv(...)` values and pass them |
</phase_requirements>

---

## Summary

Phase 11 closes three non-critical but concrete correctness gaps identified in the v1.0 audit. The gaps were not discovered during initial phases because the code was structured correctly at each layer in isolation — the failure is at the call-site integration layer, not in the implementation of individual functions. All three fixes are small, surgical, and self-contained. No new modules, no new dependencies.

Gap 1 (REPR-02): `workflow.py` already has all the information needed to populate `calculation_parameters` (energy calculator name, dipole calculator name, model config constants). The `save_frequency_results()` method already has a `calculation_parameters` parameter that accepts `dict | None`. The only missing step is assembling a dict and passing it at the call-site in `run_frequency_calculation()`.

Gap 2 (ERR-03): `CUDANotAvailableWarning` is fully defined in `utils/exceptions.py` and exported from `utils/__init__.py`. `detect_device()` only emits a `logger.warning`. The fix is one additional line: `warnings.warn("CUDA not available...", CUDANotAvailableWarning, stacklevel=2)`. The existing `logger.warning` line should stay for log file output; both channels serve different consumers.

Gap 3 (ERR-04 + bonus export gap): `cli.py` calls `validate_all_prerequisites(dipole_model_path=None, helper_script_path=None)`. The resolved env var paths are already available from module-level constants in `calculators/mace_ml.py` (`DEFAULT_MACE_DIPOLE_MODEL`) and `gaussian/io.py` (`DEFAULT_HELPER_SCRIPT`). The fix is to resolve the env vars at the cli.py call-site using `os.getenv()` with the same defaults. Additionally, `GaussianRunError` and `GaussianTimeoutError` are defined in `utils/exceptions.py` but not re-exported from `utils/__init__.py` — a one-line addition to the import list closes this.

**Primary recommendation:** One single plan (11-01) targeting all four fixes in `workflow.py`, `cli.py`, `utils/validation.py`, and `utils/__init__.py`. Each fix is 1-5 lines. Total implementation scope: approximately 15-20 lines changed across 4 files.

---

## Standard Stack

### Core (Python stdlib — no new dependencies)

| Module | Purpose | Notes |
|--------|---------|-------|
| `warnings` | Emit `CUDANotAvailableWarning` via `warnings.warn()` | Already imported in `workflow.py`; needs to be imported in `validation.py` |
| `os` | `os.getenv()` for env var resolution in `cli.py` | Already imported in `cli.py` implicitly; needs explicit import or use of existing import |

### Existing Project Modules (no new modules needed)

| Module | What to Use |
|--------|------------|
| `mace_gaussian.utils.exceptions` | `CUDANotAvailableWarning`, `GaussianRunError`, `GaussianTimeoutError` — all already defined |
| `mace_gaussian.utils.validation` | `detect_device()` — add `warnings.warn` call here |
| `mace_gaussian.utils.__init__` | Add `GaussianRunError`, `GaussianTimeoutError` to imports and `__all__` |
| `mace_gaussian.workflow` | `run_frequency_calculation()` — assemble and pass `calculation_parameters` dict |
| `mace_gaussian.cli` | `run()` command — resolve env var paths before calling `validate_all_prerequisites()` |

### No New Dependencies

All fixes use:
- Python stdlib: `warnings`, `os`
- Existing project constants: `DEFAULT_MACE_DIPOLE_MODEL` from `calculators/mace_ml.py`, `DEFAULT_HELPER_SCRIPT` from `gaussian/io.py`

**Installation:** None required.

---

## Architecture Patterns

### Current State (exactly as read from source)

```
cli.py run() command:
  validate_all_prerequisites(
      check_gaussian=True,
      check_formchk_tool=True,
      dipole_model_path=None,     # BUG: should be os.getenv("MACE_DIPOLE_MODEL_PATH", DEFAULT)
      helper_script_path=None,    # BUG: should be os.getenv("MACE_HELPER_SCRIPT_PATH", DEFAULT)
  )

validation.py detect_device():
  if not torch.cuda.is_available():
      logger.warning("CUDA not available, falling back to CPU")   # OK
      # MISSING: warnings.warn(..., CUDANotAvailableWarning)

workflow.py run_frequency_calculation():
  results_mgr.save_frequency_results(
      ...
      # MISSING: calculation_parameters={"energy_calculator": ..., "dipole_calculator": ...}
  )

utils/__init__.py:
  # MISSING: GaussianRunError, GaussianTimeoutError in imports and __all__
```

### Fixed State (target)

```
cli.py run() command:
  import os
  from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL
  from mace_gaussian.gaussian.io import DEFAULT_HELPER_SCRIPT
  validate_all_prerequisites(
      check_gaussian=True,
      check_formchk_tool=True,
      dipole_model_path=os.getenv("MACE_DIPOLE_MODEL_PATH", DEFAULT_MACE_DIPOLE_MODEL),
      helper_script_path=os.getenv("MACE_HELPER_SCRIPT_PATH", DEFAULT_HELPER_SCRIPT),
  )

validation.py detect_device():
  import warnings
  if not torch.cuda.is_available():
      logger.warning("CUDA not available, falling back to CPU")       # keep
      warnings.warn("CUDA not available, falling back to CPU",        # add
                    CUDANotAvailableWarning, stacklevel=2)

workflow.py run_frequency_calculation():
  calculation_parameters = {
      "energy_calculator": energy_calculator_name,
      "dipole_calculator": dipole_calculator_name,
  }
  results_mgr.save_frequency_results(
      ...
      calculation_parameters=calculation_parameters,
  )

utils/__init__.py:
  from .exceptions import (
      CUDANotAvailableWarning,
      GaussianParseError,
      GaussianRunError,          # ADD
      GaussianTimeoutError,      # ADD
      InputValidationError,
      MaceGaussianError,
      PrerequisiteError,
  )
  __all__ = [
      ...
      "GaussianRunError",        # ADD
      "GaussianTimeoutError",    # ADD
  ]
```

### Pattern: Minimal Surgery

Each fix targets exactly one location. No refactoring, no new abstractions. Read the exact line, change only the call-site.

### Anti-Patterns to Avoid

- **Do not duplicate env var resolution**: The env vars are already resolved via module-level constants in `calculators/mace_ml.py` and `gaussian/io.py`. Import and reuse those constants — don't hardcode paths in `cli.py`.
- **Do not remove `logger.warning`**: The `warnings.warn` call augments the existing logger call; it does not replace it. Both channels are used by different consumers (log files vs. `warnings.filterwarnings`).
- **Do not add model config details to `calculation_parameters` that require calculator instantiation**: Assembling `calculation_parameters` should use only the names/strings already available in `run_frequency_calculation()` scope — no new calculator calls.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Catchable GPU warning | Custom event/callback | `warnings.warn(..., CUDANotAvailableWarning)` | Python standard mechanism; compatible with `pytest.warns`, `warnings.filterwarnings` |
| Env var path resolution | Custom config parser | `os.getenv("KEY", default)` | Already used by the existing module-level constants in the codebase |
| Gaussian exception exports | Duplicate exception classes | Add to `utils/__init__.py` imports | Exceptions already exist in `utils/exceptions.py` |

---

## Common Pitfalls

### Pitfall 1: Importing Heavy Modules at CLI Top-Level

**What goes wrong:** Importing `DEFAULT_MACE_DIPOLE_MODEL` from `calculators/mace_ml.py` at the `cli.py` module level would trigger import of the calculators package, which loads heavy deps (DGL/espaloma). Previous phases decided to use lazy imports in the `run()` function body.

**Why it happens:** `mace_ml.py` has a module-level `os.getenv()` call that is safe, but the import itself may trigger `calculators/__init__.py` to load `EspalomaDipoleCalculator`.

**How to avoid:** Do not import from `calculators/mace_ml.py` or `gaussian/io.py` at cli.py module level. Instead, inline the `os.getenv()` calls with the same default paths directly inside the `run()` function body. This is cleaner and avoids the heavy-dep problem.

**Correct pattern in cli.py `run()` body:**
```python
import os
from pathlib import Path

dipole_model_path = os.getenv(
    "MACE_DIPOLE_MODEL_PATH",
    str(Path.home() / "mace_gaussian" / "mace4ir_models" / "pretrained_models" / "model_1_dipole.model")
)
helper_script_path = os.getenv(
    "MACE_HELPER_SCRIPT_PATH",
    str(Path(__file__).parent / "gm_helper.py")
)
```

### Pitfall 2: `warnings.warn` stacklevel

**What goes wrong:** Using `stacklevel=1` (the default) makes the warning appear to originate from inside `detect_device()`, not from the caller. Users cannot meaningfully filter it.

**How to avoid:** Use `stacklevel=2` so the warning appears to originate from the function calling `detect_device()`.

**Warning sign:** If `warnings.warn` is called with default `stacklevel`, the source file in warnings output will be `validation.py` rather than the module that called `detect_device()`.

### Pitfall 3: `validate_all_prerequisites` raises on `None`

**What goes wrong:** `check_dipole_model(None)` calls `Path(None)` which raises `TypeError` in Python 3.9+, not a clean `PrerequisiteError`. The current behavior is actually a deferred crash, not a clean early failure.

**The correct behavior:** The fix in ERR-04 is to pass the real resolved path. If the env var is unset and the default path doesn't exist, `check_dipole_model` will raise `PrerequisiteError` with a helpful message. This is the intended behavior.

**How to avoid:** Always pass the resolved path string (never `None`) to `validate_all_prerequisites`. The `if dipole_model_path is not None:` guard in `validate_all_prerequisites` is what allows skipping the check when `None` is passed — but this was the source of the bug.

### Pitfall 4: `calculation_parameters` dict keys

**What goes wrong:** Including keys whose values require additional computation or imports (e.g., model file path, model version, basis set) in `calculation_parameters` could add fragility or trigger additional imports.

**How to avoid:** Keep `calculation_parameters` minimal — just the string names of the calculators that are already in scope:
```python
calculation_parameters = {
    "energy_calculator": energy_calculator_name,
    "dipole_calculator": dipole_calculator_name,
}
```

### Pitfall 5: `from __future__ import annotations` compatibility

**What goes wrong:** `utils/results.py` has `dict[str, Any]` in a method signature without `from __future__ import annotations`. The test run showed this causes a `TypeError: 'type' object is not subscriptable` on Python 3.9 (pyproject.toml supports `>=3.9`, Python 3.9 does not support `dict[...]` as a type annotation at runtime without `__future__`).

**How to avoid:** When touching `utils/results.py`, add `from __future__ import annotations` at the top. This is an existing bug that should be fixed as part of the phase.

**Note:** This is a pre-existing bug discovered during test collection, not caused by Phase 11 changes. It blocks `test_validation.py` collection. Fix it as part of the wave 0 / first task.

---

## Code Examples

### Fix 1: REPR-02 — Populate calculation_parameters in workflow.py

```python
# Source: direct inspection of workflow.py run_frequency_calculation() and results.py
# In run_frequency_calculation(), before the save call:

calculation_parameters = {
    "energy_calculator": energy_calculator_name,
    "dipole_calculator": dipole_calculator_name,
}

results_mgr.save_frequency_results(
    molecule_name=molecule_name,
    energy_calculator=energy_calculator_name,
    dipole_calculator=dipole_calculator_name,
    calculator_type="ml",
    frequencies_data={...},
    energy=final_energy,
    dipole=parsed_data.get("dipole_moment"),
    runtime=runtime,
    gaussian_log=log_final if Path(log_final).exists() else None,
    gaussian_gjf=gjf_final if Path(gjf_final).exists() else None,
    timestamp=timestamp,
    calculation_parameters=calculation_parameters,   # ADD THIS
)
```

### Fix 2: ERR-03 — Emit CUDANotAvailableWarning in validation.py

```python
# Source: direct inspection of utils/validation.py detect_device()
# Change needed: add import and warnings.warn call

import warnings
from .exceptions import CUDANotAvailableWarning  # already in exceptions.py

def detect_device() -> str:
    try:
        import torch
        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            logger.info("CUDA available: %s", device_name)
            return "cuda"
        else:
            logger.warning("CUDA not available, falling back to CPU")  # keep
            warnings.warn(                                              # add
                "CUDA not available, falling back to CPU",
                CUDANotAvailableWarning,
                stacklevel=2,
            )
            return "cpu"
    except ImportError:
        logger.warning("PyTorch not installed, falling back to CPU")   # keep
        warnings.warn(                                                  # add
            "PyTorch not installed, falling back to CPU",
            CUDANotAvailableWarning,
            stacklevel=2,
        )
        return "cpu"
```

### Fix 3: ERR-04 — Resolve env var paths in cli.py

```python
# Source: direct inspection of cli.py run() and gaussian/io.py, calculators/mace_ml.py
# Add inside the run() function body, before validate_all_prerequisites call:

import os
from pathlib import Path

# Resolve env var paths using the same defaults as the calculator modules
dipole_model_path = os.getenv(
    "MACE_DIPOLE_MODEL_PATH",
    str(
        Path.home()
        / "mace_gaussian"
        / "mace4ir_models"
        / "pretrained_models"
        / "model_1_dipole.model"
    ),
)
helper_script_path = os.getenv(
    "MACE_HELPER_SCRIPT_PATH",
    str(Path(__file__).parent / "gm_helper.py"),
)

validate_all_prerequisites(
    check_gaussian=True,
    check_formchk_tool=True,
    dipole_model_path=dipole_model_path,
    helper_script_path=helper_script_path,
)
```

### Fix 4: Export GaussianRunError/GaussianTimeoutError from utils/__init__.py

```python
# Source: direct inspection of utils/__init__.py and utils/exceptions.py

from .exceptions import (
    CUDANotAvailableWarning,
    GaussianParseError,
    GaussianRunError,        # ADD
    GaussianTimeoutError,    # ADD
    InputValidationError,
    MaceGaussianError,
    PrerequisiteError,
)

__all__ = [
    # ... existing entries ...
    "GaussianRunError",      # ADD
    "GaussianTimeoutError",  # ADD
]
```

### Fix 5 (Pre-existing Bug): Add `from __future__ import annotations` to results.py

```python
# Source: test run error — TypeError: 'type' object is not subscriptable on Python 3.9
# dict[str, Any] annotation at line 170 of results.py requires __future__ on Python < 3.10

# Add at the top of utils/results.py:
from __future__ import annotations
```

---

## State of the Art

| Current Behavior | Target Behavior | Impact |
|-----------------|-----------------|--------|
| `detect_device()` logs to file only; GPU fallback not catchable | `warnings.warn(CUDANotAvailableWarning)` emitted | Users can use `warnings.filterwarnings` or `pytest.warns` to intercept |
| `calculation_parameters: {}` in every results JSON | Populated with `energy_calculator`, `dipole_calculator` | Results are self-describing; reproducibility improves |
| Env var paths skipped at CLI startup (passed as `None`) | Resolved paths validated at startup | Fails fast with helpful message if model/script missing |
| `GaussianRunError`/`GaussianTimeoutError` importable only from `mace_gaussian.utils.exceptions` | Importable from `mace_gaussian.utils` | Shorter import path; consistent with other exceptions |

---

## Open Questions

1. **Should `calculation_parameters` include more fields?**
   - What we know: Audit notes "parameter metadata (basis set, model config) is absent" but also notes "Calculator identity stored separately in `energy_calculator`/`dipole_calculator` fields"
   - What's unclear: Whether model file path or model version should also appear
   - Recommendation: Keep minimal for now — just the calculator names. The `version_info` field already captures model versions. Adding more fields risks triggering additional imports. This can be extended in v2.

2. **Should ERR-04 fix also validate dipole model / helper script on the `diagnose` command?**
   - What we know: Previous phases decided "Validation only on run command, not list/diagnose/compare/export" (STATE.md Phase 02-03 decision)
   - What's unclear: Nothing — this is a locked prior decision
   - Recommendation: ERR-04 fix applies only to the `run` command. Do not add validation to diagnose.

3. **Pre-existing `TypeError` in results.py on Python 3.9**
   - What we know: `dict[str, Any]` at line 170 causes `TypeError` without `from __future__ import annotations`; blocks `test_validation.py` collection
   - What's unclear: Whether this affects Python 3.10+ runtime (it doesn't — only an issue at 3.9)
   - Recommendation: Add `from __future__ import annotations` to `results.py` as part of Phase 11. The project supports `>=3.9` per pyproject.toml.

---

## Validation Architecture

`workflow.nyquist_validation` is not present in `.planning/config.json` — skipping this section as it is not enabled.

---

## Sources

### Primary (HIGH confidence)

- Direct code inspection: `mace_gaussian/utils/validation.py` — `detect_device()` implementation, confirmed uses `logger.warning` not `warnings.warn`
- Direct code inspection: `mace_gaussian/utils/exceptions.py` — `CUDANotAvailableWarning`, `GaussianRunError`, `GaussianTimeoutError` all defined
- Direct code inspection: `mace_gaussian/utils/__init__.py` — confirmed `GaussianRunError`/`GaussianTimeoutError` absent from imports and `__all__`
- Direct code inspection: `mace_gaussian/cli.py` lines 111-116 — `validate_all_prerequisites(dipole_model_path=None, helper_script_path=None)` confirmed
- Direct code inspection: `mace_gaussian/workflow.py` lines 442-459 — `save_frequency_results()` call confirmed missing `calculation_parameters` argument
- Direct code inspection: `mace_gaussian/utils/results.py` line 177 — `calculation_parameters: Optional[dict[str, Any]] = None` parameter exists and is wired to JSON at line 248
- Direct code inspection: `mace_gaussian/calculators/mace_ml.py` lines 12-21 — `DEFAULT_MACE_DIPOLE_MODEL` constant with `os.getenv()`
- Direct code inspection: `mace_gaussian/gaussian/io.py` lines 14-17 — `DEFAULT_HELPER_SCRIPT` constant with `os.getenv()`
- `.planning/v1.0-MILESTONE-AUDIT.md` — audit findings with detailed gap descriptions (REPR-02, ERR-03, ERR-04, GaussianRunError export gap)
- Test run: `pytest tests/test_validation.py tests/test_exceptions.py --collect-only` — confirmed pre-existing `TypeError` in `results.py` blocking test collection on Python 3.9

### Secondary (MEDIUM confidence)

- Python docs (well-known): `warnings.warn(msg, category, stacklevel)` — `stacklevel=2` makes warning appear at caller level; this is standard Python behavior verified by direct knowledge

---

## Metadata

**Confidence breakdown:**
- Gap identification: HIGH — audit document is authoritative and code inspection confirms each gap exactly as described
- Fix implementations: HIGH — all four fixes are 1-10 lines of straightforward Python; no new design decisions required
- Pitfalls: HIGH — confirmed by test collection failure (Python 3.9 `dict[...]` bug) and prior phase decisions in STATE.md

**Research date:** 2026-02-27
**Valid until:** 2026-03-29 (stable — pure Python wiring, no external dependencies)
