# Phase 3: Extract Utilities & Conventions - Research

**Researched:** 2026-02-17
**Domain:** Python module extraction, unit conversion consolidation, codebase conventions
**Confidence:** HIGH

## Summary

This phase is a pure refactoring operation: extract unit conversion constants and functions, relocate `validation.py` and `exceptions.py` into a `utils/` package, move `ResultsManager` into `utils/results.py`, and document code conventions. No new functionality is introduced.

The codebase currently defines unit conversion constants inline in multiple places across `gm_main.py` (lines 207, 419, 474-476, 563-568, 1006-1007), `dft_baseline.py` (line 354), `fchk_parser.py` (lines 188-189), and `charge_analysis.py` (lines 659-660). The values are slightly inconsistent (e.g., `0.52917721092` vs `0.529177210903` for Bohr-to-Angstrom). Centralizing these into `utils/units.py` eliminates duplication and inconsistency.

The `validation.py` and `exceptions.py` modules created in Phase 2 are imported by 6 files (`gm_main.py`, `cli.py`, `gaussian_parser.py`, `results_manager.py`, and 3 test files). Moving them into `utils/` requires updating all import statements. The `results_manager.py` module is already standalone and imported by 3 files (`gm_main.py`, `dft_baseline.py`, and tests).

**Primary recommendation:** Create `utils/` package with `units.py`, `exceptions.py`, `validation.py`, and `results.py`. Use compatibility re-exports at old locations temporarily if needed, but prefer clean import updates since the codebase is small (20 Python files).

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions

All decisions are under Claude's Discretion for this phase. The following guidance from CONTEXT.md captures the agreed approach:

**What gets extracted:**
- Unit conversion functions (e.g., hartree<->eV, bohr<->angstrom) -> `utils/units.py`
- Pure helper functions from gm_main.py that don't depend on Gaussian/MACE runtime state
- ResultsManager class -> `utils/results.py` (already somewhat standalone in `results_manager.py`)
- Boundary: if a function requires MACE models, Gaussian subprocess, or ZMQ -- it stays in place for now (Phases 4-6 handle those)

**Validation consolidation:**
- Phase 2 created `validation.py` and `exceptions.py` at the top level
- Move these into `utils/validation.py` and `utils/exceptions.py` respectively
- Update all imports across the codebase to use the new locations
- Keep the same API -- no behavioral changes, just relocation

**Import strategy:**
- `utils/__init__.py` should re-export key items for convenience (e.g., `from utils.units import hartree_to_ev`)
- Internal code uses direct imports: `from utils.units import ...`
- No deep nesting -- utils/ is flat (units.py, validation.py, exceptions.py, results.py)

**Convention scope:**
- Document conventions in a lightweight `docs/CONVENTIONS.md` or a section in existing `docs/DEVELOPMENT.md`
- Cover: naming conventions (snake_case functions, PascalCase classes), error handling pattern (use exceptions.py hierarchy), unit conventions (internal units: cm-1 for frequencies, Angstrom for distances), import ordering (stdlib -> third-party -> local)
- Keep it pragmatic -- conventions should reflect what the codebase already does well, not aspirational rules
- Ruff configuration in pyproject.toml enforces what can be automated

### Claude's Discretion

All implementation decisions are at Claude's discretion.

### Deferred Ideas (OUT OF SCOPE)

None -- discussion stayed within phase scope.

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-01 | Utility functions (unit conversions, validation) extracted from gm_main.py into separate modules | Unit conversion audit (6+ inline definitions), validation/exceptions relocation plan, ResultsManager extraction, import dependency map -- all documented below |

</phase_requirements>

## Standard Stack

### Core

No new libraries needed. This phase uses only existing project infrastructure:

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| ruff | >=0.9.0 | Lint + format + isort enforcement | Already configured in pyproject.toml |
| pytest | >=7.0.0 | Test extracted modules | Already configured with markers and test paths |
| Python stdlib | 3.9+ | pathlib, typing, logging | Already used throughout codebase |

### Supporting

No additional libraries needed. This is pure code reorganization.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Flat `utils/` package | Nested `utils/conversion/`, `utils/io/` | Over-engineering for ~4 modules; flat is correct here |
| Compatibility shims at old locations | Direct import updates only | Shims add complexity; codebase is small enough for direct updates |
| `pint` for unit conversions | Hand-written constants | Pint is overkill; these are 5 fixed physical constants, not a general unit system |

## Architecture Patterns

### Target Project Structure

```
mace_gaussian/
├── utils/
│   ├── __init__.py         # Re-exports key items
│   ├── units.py            # Physical constants + conversion functions
│   ├── exceptions.py       # Exception hierarchy (moved from top-level)
│   ├── validation.py       # Prerequisite checks, XYZ validation (moved from top-level)
│   └── results.py          # ResultsManager class (moved from results_manager.py)
├── gm_main.py              # Slimmed: imports from utils/ instead of inline constants
├── cli.py                  # Updated imports
├── gaussian_parser.py      # Updated imports
├── dft_baseline.py         # Updated imports
├── tests/
│   ├── test_units.py       # NEW: tests for unit conversion functions
│   ├── test_exceptions.py  # Updated imports
│   └── test_validation.py  # Updated imports
└── docs/
    └── CONVENTIONS.md       # NEW: documented code conventions
```

### Pattern 1: Centralized Physical Constants Module

**What:** Single source of truth for all physical constants and unit conversion functions.
**When to use:** Any time code needs to convert between unit systems (eV<->Hartree, Bohr<->Angstrom, etc.).

```python
# utils/units.py
"""Physical constants and unit conversion utilities.

All constants use CODATA 2018 recommended values.
Internal units: eV (energy), Angstrom (distance), cm-1 (frequency).
"""

# --- Physical Constants (CODATA 2018) ---
BOHR_TO_ANGSTROM: float = 0.529177210903
ANGSTROM_TO_BOHR: float = 1.0 / BOHR_TO_ANGSTROM

HARTREE_TO_EV: float = 27.211386245988
EV_TO_HARTREE: float = 1.0 / HARTREE_TO_EV

# Note: some code uses 27.211386246 (fewer decimals) -- standardize on CODATA 2018
# Note: charge_analysis.py uses ANGSTROM_TO_BOHR = 1.8897259886 which is the inverse
#       of 0.529177210903 = 1.8897259886. This is consistent.


# --- Convenience Functions ---
def hartree_to_ev(energy_hartree: float) -> float:
    """Convert energy from Hartree to eV."""
    return energy_hartree * HARTREE_TO_EV


def ev_to_hartree(energy_ev: float) -> float:
    """Convert energy from eV to Hartree."""
    return energy_ev * EV_TO_HARTREE


def bohr_to_angstrom(distance_bohr: float) -> float:
    """Convert distance from Bohr to Angstrom."""
    return distance_bohr * BOHR_TO_ANGSTROM


def angstrom_to_bohr(distance_angstrom: float) -> float:
    """Convert distance from Angstrom to Bohr."""
    return distance_angstrom * ANGSTROM_TO_BOHR
```

### Pattern 2: Package Re-exports via `__init__.py`

**What:** Convenience imports so callers can do `from utils import MaceGaussianError` or `from utils.exceptions import MaceGaussianError`.
**When to use:** For the most commonly used items (exceptions, key constants).

```python
# utils/__init__.py
"""Utility package for MACE-Gaussian interface."""

from utils.exceptions import (
    CUDANotAvailableWarning,
    GaussianParseError,
    InputValidationError,
    MaceGaussianError,
    PrerequisiteError,
)
from utils.units import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HARTREE,
    HARTREE_TO_EV,
    angstrom_to_bohr,
    bohr_to_angstrom,
    ev_to_hartree,
    hartree_to_ev,
)
```

### Pattern 3: Import Update Strategy

**What:** Systematic find-and-replace of all import statements.
**When to use:** When relocating modules.

Current imports that must change:

| File | Current Import | New Import |
|------|---------------|------------|
| `gm_main.py` | `from results_manager import ResultsManager` | `from utils.results import ResultsManager` |
| `gm_main.py` | `from validation import detect_device` | `from utils.validation import detect_device` |
| `cli.py` | `from exceptions import ...` | `from utils.exceptions import ...` |
| `cli.py` | `from validation import ...` | `from utils.validation import ...` |
| `gaussian_parser.py` | `from exceptions import GaussianParseError` | `from utils.exceptions import GaussianParseError` |
| `dft_baseline.py` | `from results_manager import ResultsManager` | `from utils.results import ResultsManager` |
| `results_manager.py` -> `utils/results.py` | `from validation import collect_version_metadata` | `from utils.validation import collect_version_metadata` |
| `validation.py` -> `utils/validation.py` | `from exceptions import ...` | `from utils.exceptions import ...` |
| `tests/test_validation.py` | `from exceptions import ...` | `from utils.exceptions import ...` |
| `tests/test_validation.py` | `from validation import ...` | `from utils.validation import ...` |
| `tests/test_exceptions.py` | `from exceptions import ...` | `from utils.exceptions import ...` |
| `tests/test_regression.py` | `from results_manager import ResultsManager` | `from utils.results import ResultsManager` |

### Anti-Patterns to Avoid

- **Leaving old modules as stubs:** Don't leave `exceptions.py` at top level with `from utils.exceptions import *`. Delete the old files after migration -- clean break. The codebase is small enough that lingering compatibility shims add confusion.
- **Circular imports in utils/:** `validation.py` imports from `exceptions.py`. Both will be in `utils/`. Make sure `exceptions.py` has zero local imports (it currently doesn't import anything local -- good).
- **Over-extracting from gm_main.py:** Functions like `parse_gaussian_input()`, `write_gaussian_output()`, `zmq_server()` touch Gaussian I/O and ZMQ -- these belong to Phases 4-6. Only extract unit constants and pure helpers.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Physical constants | Typed-out numbers inline | `utils/units.py` constants | Already have 3 slightly different Bohr values in the codebase |
| Import sorting | Manual ordering | `ruff check --select I --fix` | Already configured in pyproject.toml with known-first-party |
| Code formatting | Manual style | `ruff format` | Already configured, line-length 100 |

**Key insight:** The unit conversion problem looks trivial but is already causing inconsistency (two different Bohr-to-Angstrom values in the codebase). Centralizing it is the correct fix.

## Common Pitfalls

### Pitfall 1: Inconsistent Bohr-to-Angstrom Constants

**What goes wrong:** `gm_main.py` line 419 uses `0.52917721092` (10 decimals) while line 207 uses `0.529177210903` (12 decimals). These differ by ~1e-12 -- negligible for chemistry but indicates copy-paste drift.
**Why it happens:** Constants are defined inline wherever needed.
**How to avoid:** Single definition in `utils/units.py`. Use CODATA 2018 value: `0.529177210903`.
**Warning signs:** `grep -r "0.5291" *.py` returning multiple definitions.

### Pitfall 2: Breaking the `known-first-party` Ruff Configuration

**What goes wrong:** After moving modules into `utils/`, the ruff isort configuration in `pyproject.toml` still lists `"exceptions"`, `"validation"`, `"results_manager"` as known first-party. This will cause ruff to mis-sort imports.
**Why it happens:** Configuration not updated alongside code moves.
**How to avoid:** Update `[tool.ruff.lint.isort] known-first-party` to include `"utils"` and remove old module names.
**Warning signs:** `ruff check --select I` showing import ordering violations after migration.

### Pitfall 3: Breaking Lazy Imports in results_manager.py

**What goes wrong:** `results_manager.py` uses lazy imports (`from validation import collect_version_metadata` inside methods, wrapped in try/except ImportError). When moved to `utils/results.py`, this import path changes.
**Why it happens:** The lazy import pattern uses string-based module resolution.
**How to avoid:** When moving to `utils/results.py`, change to `from utils.validation import collect_version_metadata`. The try/except ImportError pattern should remain since it's a deliberate fallback.
**Warning signs:** `version_info` silently becoming `{}` in results JSON after migration.

### Pitfall 4: Forgetting to Update test_cli_validation.py

**What goes wrong:** Test file imports break silently or are missed in the sweep.
**Why it happens:** Test files are in a subdirectory and easy to overlook.
**How to avoid:** Run full test suite (`pytest tests/`) after every import update batch. Check all 7 test files.
**Warning signs:** `ImportError` in test collection.

### Pitfall 5: Missing `__init__.py` for utils Package

**What goes wrong:** Python cannot find `utils.units` if `utils/__init__.py` doesn't exist.
**Why it happens:** Forgetting to create the package init file.
**How to avoid:** Create `utils/__init__.py` as the first step before moving any modules.
**Warning signs:** `ModuleNotFoundError: No module named 'utils'`.

### Pitfall 6: Coverage Configuration Drift

**What goes wrong:** `[tool.coverage.run] source` in pyproject.toml lists `"results_manager"` but after the move it should be `"utils"` or `"utils/results"`.
**Why it happens:** Configuration not updated alongside code moves.
**How to avoid:** Update coverage source list to include `"utils"` and remove `"results_manager"`.

## Code Examples

### Unit Constant Replacement in gm_main.py

Before (scattered inline):
```python
# gm_main.py line 419
BOHR_TO_ANGSTROM = 0.52917721092

# gm_main.py line 474
ANGSTROM_TO_BOHR = 0.52917721092
EV_TO_HARTREE = 27.211386246

# gm_main.py line 563
EV_TO_HARTREE = 27.211386246
```

After (centralized import):
```python
# gm_main.py - top of file
from utils.units import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HARTREE,
)

# Then use constants directly in functions - no more local definitions
```

### Test for Unit Conversions

```python
# tests/test_units.py
import pytest
from utils.units import (
    BOHR_TO_ANGSTROM,
    HARTREE_TO_EV,
    bohr_to_angstrom,
    ev_to_hartree,
    hartree_to_ev,
)


class TestConstants:
    """Verify physical constant values against CODATA 2018."""

    def test_bohr_to_angstrom(self):
        assert BOHR_TO_ANGSTROM == pytest.approx(0.529177210903, rel=1e-10)

    def test_hartree_to_ev(self):
        assert HARTREE_TO_EV == pytest.approx(27.211386245988, rel=1e-10)

    def test_roundtrip_energy(self):
        """Converting eV -> Hartree -> eV should be identity."""
        original = 13.6
        assert hartree_to_ev(ev_to_hartree(original)) == pytest.approx(original, rel=1e-12)

    def test_roundtrip_distance(self):
        """Converting Bohr -> Angstrom -> Bohr should be identity."""
        from utils.units import angstrom_to_bohr
        original = 1.0
        assert angstrom_to_bohr(bohr_to_angstrom(original)) == pytest.approx(original, rel=1e-12)
```

### Convention Documentation Structure

```markdown
# docs/CONVENTIONS.md

## Naming
- Functions/variables: `snake_case`
- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Private functions: `_leading_underscore`

## Units (internal representation)
- Energy: eV (ASE default)
- Distance: Angstrom (ASE default)
- Frequency: cm-1
- Dipole: e*Bohr (Gaussian convention)
- All conversions via `utils.units` -- never define constants inline

## Error Handling
- Use exception hierarchy from `utils.exceptions`
- Catch `MaceGaussianError` for broad error handling
- Catch specific subclasses (`PrerequisiteError`, etc.) for targeted handling
- Never bare `except:` -- always specify exception type

## Imports (enforced by ruff)
1. Standard library
2. Third-party packages
3. Local imports (`utils`, `gm_main`, etc.)
- Separated by blank lines, alphabetized within groups

## File Organization
- One class per module when class is large (e.g., ResultsManager)
- Related pure functions grouped in themed modules (e.g., units.py)
- Tests mirror source structure: `utils/units.py` -> `tests/test_units.py`
```

## Detailed Extraction Inventory

### Functions/Constants to Extract from gm_main.py

| Item | Location | Target | Reason |
|------|----------|--------|--------|
| `BOHR_TO_ANGSTROM` constant | lines 419, 474 | `utils/units.py` | Duplicated, inconsistent values |
| `EV_TO_HARTREE` constant | lines 475, 563, 1006 | `utils/units.py` | Duplicated 3 times |
| `ANGSTROM_TO_BOHR` constant | lines 474, 567 | `utils/units.py` | Duplicated |
| `0.529177210903` (dipole conversion) | line 207 | `utils/units.py` as `BOHR_TO_ANGSTROM` | Inline magic number |

### Files to Move Wholesale

| Current Location | New Location | Imports to Update |
|-----------------|--------------|-------------------|
| `exceptions.py` | `utils/exceptions.py` | 6 files (see Pattern 3 table) |
| `validation.py` | `utils/validation.py` | 4 files |
| `results_manager.py` | `utils/results.py` | 4 files |

### Functions that Stay in gm_main.py (Phase 4-6 scope)

| Function | Why It Stays |
|----------|-------------|
| `zmq_server()` | ZMQ runtime dependency (Phase 6) |
| `is_calc_finished()` | ZMQ + subprocess dependency |
| `parse_gaussian_input()` | Gaussian I/O (Phase 6) |
| `write_gaussian_output()` | Gaussian I/O (Phase 6) |
| `run_next_calculation()` | Orchestrates Gaussian external interface |
| `DipoleCalculatorBase` and subclasses | Calculator hierarchy (Phase 4) |
| `DipoleCalculatorFactory` | Calculator hierarchy (Phase 4) |
| `calculator()` function | MACE imports (Phase 5) |
| `run_frequency_calculation()` | Gaussian subprocess + ZMQ |
| `run_geometry_optimization()` | ASE optimizer + MACE |
| `run_workflow()` | Top-level orchestrator |

### Configuration Updates Required

| File | Section | Change |
|------|---------|--------|
| `pyproject.toml` | `[tool.ruff.lint.isort] known-first-party` | Add `"utils"`, remove `"exceptions"`, `"validation"`, `"results_manager"` |
| `pyproject.toml` | `[tool.coverage.run] source` | Add `"utils"`, remove `"results_manager"` |
| `pyproject.toml` | `[tool.coverage.run] source` | Keep `"gm_main"`, `"cli"`, etc. |

## Open Questions

1. **Should old top-level files be deleted or kept as shims?**
   - What we know: Only ~20 Python files in the project, 12 import statements to update.
   - What's unclear: Whether any external scripts or notebooks import from the old locations.
   - Recommendation: Delete old files. The codebase is small. If anything breaks, the error message clearly says what to import from where. Add a one-line comment in git commit message noting the move.

2. **Which CODATA value to standardize on?**
   - What we know: Two slightly different Bohr values exist (`0.52917721092` vs `0.529177210903`). Both are from CODATA but different precision levels.
   - Recommendation: Use CODATA 2018 full precision: `0.529177210903` for Bohr, `27.211386245988` for Hartree-to-eV. This matches ASE's internal values.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis: `gm_main.py` (1311 lines), `validation.py`, `exceptions.py`, `results_manager.py`, `cli.py`, all test files
- `pyproject.toml` ruff and coverage configuration
- CODATA 2018 values verified against NIST (Bohr radius: 0.529177210903 Angstrom, Hartree: 27.211386245988 eV)

### Secondary (MEDIUM confidence)
- ASE documentation for internal unit conventions (eV, Angstrom)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - no new libraries, pure reorganization of existing code
- Architecture: HIGH - clear extraction targets identified from direct code analysis
- Pitfalls: HIGH - all identified from actual codebase inspection (import chains, config files, constant inconsistencies)

**Research date:** 2026-02-17
**Valid until:** 2026-03-17 (stable -- this is internal refactoring, no external API dependencies)
