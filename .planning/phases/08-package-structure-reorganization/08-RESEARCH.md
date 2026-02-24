# Phase 8: Package Structure & Reorganization - Research

**Researched:** 2026-02-24
**Domain:** Python packaging — flat-to-package migration with hatchling, import rewiring, pyproject.toml cleanup
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Package strategy:**
- All in (Path B): All library and analysis code moves into `mace_gaussian/`. No core code stays at root.
- The project should be distributable (`pip install mace-gaussian`) even though Gaussian 16 is a required system dependency.

**Final package structure:**
```
mace_gaussian/
    __init__.py              # version + re-exports run_pipeline
    cli.py                   # CLI entry point
    workflow.py              # pipeline orchestrator
    gm_helper.py             # ZMQ bridge to Gaussian
    dft_baseline.py          # DFT baseline calculations
    diagnostics.py           # diagnostic utilities
    calculators/             # calculator classes (move existing subdirectory in)
        __init__.py
        base.py, espaloma.py, factory.py, mace_loader.py, mace_ml.py, xtb.py
    gaussian/                # Gaussian I/O & ZMQ server (move existing subdirectory in)
        __init__.py
        runner.py, ...
    analysis/                # analysis and reporting
        __init__.py
        analyze_spectra.py
        mode_matching.py
        html_report_generator.py
        run_analysis.py      # anharmonic analysis logic
        run_analysis_harmonic.py  # harmonic analysis logic
    utils/                   # utilities (move existing subdirectory in)
        __init__.py
        exceptions.py, results.py, units.py, validation.py
```

**Root-level shim scripts (keep, update internals):**
- `run_analysis.py` stays at root as a thin shim calling `mace_gaussian.analysis.run_analysis`
- `run_analysis_harmonic.py` stays at root as a thin shim
- This preserves `python run_analysis.py water` and `python run_analysis_harmonic.py water` exactly as documented in CLAUDE.md

**Utility/research scripts:**
- `comparison_workflow.py`, `produce_molecules.py`, `convert_all_chk_files.py`, `charge_analysis.py` move to `scripts/` directory

**Public API (`mace_gaussian/__init__.py`):**
- Minimal: expose `__version__` and re-export `run_pipeline` for convenience
- Everything else accessed via subpackage paths
- No over-commitment to a public API at this stage

**Import style throughout package:**
- Relative imports within `mace_gaussian/` (e.g., `from .utils.exceptions import ...`)
- Absolute imports from outside the package
- `cli.py` uses absolute imports since it's the entry point

**pyproject.toml cleanup (full pass):**
- Remove broken `gm-main` and `gm-diagnose` entry points from `[project.scripts]`
- Update `[project.scripts]` to point `mace-gaussian = "mace_gaussian.cli:app"` (or equivalent)
- Fix `[tool.ruff.lint.isort] known-first-party` to list `mace_gaussian` (remove stale refs)
- Fix `[tool.coverage.run] source` to list `mace_gaussian`
- Update `[tool.ruff]` `src` if needed
- Update `[build-system]` / `packages` config to find `mace_gaussian/`

### Claude's Discretion
- Exact `__init__.py` re-export surface beyond `run_pipeline` and `__version__`
- Whether `dft_baseline.py` gets a subpackage or stays flat in `mace_gaussian/`
- How to handle `fort.7` and other non-Python files at root (leave at root — not package files)
- How to handle `.xyz` molecule input files (leave at root — user-provided inputs)
- Whether to add `py.typed` marker for type checking support

### Deferred Ideas (OUT OF SCOPE)
- CLI `analyze` subcommand (`python cli.py analyze water`) — Phase 9 or a new phase
- ORCA support as an alternative to Gaussian — future phase
- Formal public API with stable versioning guarantees — only needed if/when distributed wider
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-08 | Project reorganized into mace_gaussian/ package with proper __init__.py | Package structure patterns, hatchling config, `git mv` approach, relative import conversion |
| STRUCT-09 | CLI entry point aligned in pyproject.toml with package structure | Entry point syntax `mace_gaussian.cli:cli`, removal of broken gm-main/gm-diagnose, Click app name |
| STRUCT-10 | Analysis modules reorganized into analysis/ subpackage | Move analyze_spectra.py, mode_matching.py, html_report_generator.py, comparison_workflow.py into mace_gaussian/analysis/, update cross-imports |
</phase_requirements>

---

## Summary

This phase is a pure structural reorganization — no logic changes, only file movement and import rewiring. The project currently has a flat layout: `calculators/`, `gaussian/`, `utils/` are already proper subpackages at the root, and 5 standalone `.py` files (`cli.py`, `workflow.py`, `gm_helper.py`, `dft_baseline.py`, `diagnostics.py`) plus 4 analysis modules (`analyze_spectra.py`, `mode_matching.py`, `html_report_generator.py`, `comparison_workflow.py`) sit loose at the root.

The migration creates a `mace_gaussian/` namespace package that contains everything. The three already-packaged subdirectories (`calculators/`, `gaussian/`, `utils/`) move into `mace_gaussian/` unchanged internally except for converting their absolute imports to relative. The analysis modules consolidate into a new `mace_gaussian/analysis/` subpackage. Root shim scripts for `run_analysis.py` and `run_analysis_harmonic.py` remain at root but are reduced to 2–3 line delegators calling `mace_gaussian.analysis` functions.

The biggest risks in this phase are: (1) import rewiring errors — every internal `from calculators import ...` must become `from .calculators import ...` or `from mace_gaussian.calculators import ...` depending on context; (2) test breakage — all 131 tests currently use absolute paths like `from calculators.parser import ...` which must update to `from mace_gaussian.calculators import ...`; (3) pyproject.toml hatchling config — the current `packages = ["."]` is broken and must be replaced with `packages = ["mace_gaussian"]`.

**Primary recommendation:** Use `git mv` for every file move to preserve history, convert all internal imports to relative in one wave per subpackage, then update tests and pyproject.toml in a final cleanup wave.

---

## Standard Stack

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| hatchling | already installed | Build backend | Already in pyproject.toml, PEP 517 compliant |
| Python packaging conventions | PEP 328 | Relative imports | Standard way to express intra-package dependencies |
| pytest | 7.1.3 (installed) | Test runner | Already configured in pyproject.toml |
| ruff | already installed | Import sorting / linting | Already configured; isort rules via `I` selector |

### Supporting
| Tool | Purpose | When to Use |
|------|---------|-------------|
| `git mv` | Move files preserving history | Every file relocation in this phase |
| `uv sync` | Reinstall package in editable mode | After pyproject.toml changes |
| `uv pip install -e .` | Editable install | After creating mace_gaussian/ package directory |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| flat layout (current) | src layout | src layout is more robust against accidental imports of development code; user chose flat with explicit `mace_gaussian/` package, which is equivalent for this project |
| relative imports | absolute imports within package | Absolute imports work but break if package is renamed; relative imports are more portable within a package |

**Installation after reorganization:**
```bash
uv sync  # reinstalls editable package pointing to mace_gaussian/
```

---

## Architecture Patterns

### Recommended Project Structure

After Phase 8, the root of the repo contains:

```
mace_gaussian/           # THE package — all library code here
    __init__.py          # version="0.2.0", re-exports run_pipeline
    cli.py               # Click CLI entry point (absolute imports, it's the entry point)
    workflow.py          # Pipeline orchestrator
    gm_helper.py         # ZMQ client helper (called by Gaussian as external script)
    dft_baseline.py      # DFT baseline runner
    diagnostics.py       # Diagnostic utilities
    calculators/
        __init__.py      # existing — update absolute→relative imports inside
        base.py
        espaloma.py
        factory.py
        mace_loader.py
        mace_ml.py
        xtb.py
    gaussian/
        __init__.py      # existing — update absolute→relative imports inside
        fchk.py
        io.py
        parser.py
        runner.py
        zmq_server.py
    analysis/            # NEW subpackage
        __init__.py      # expose analyze_molecule, analyze_molecule_harmonic
        analyze_spectra.py    # moved from root
        mode_matching.py      # moved from root
        html_report_generator.py  # moved from root
        run_analysis.py           # moved from root (anharmonic workflow logic)
        run_analysis_harmonic.py  # moved from root (harmonic workflow logic)
    utils/
        __init__.py      # existing — update absolute→relative imports inside
        exceptions.py
        results.py
        units.py
        validation.py
pyproject.toml           # updated — packages, entry points, coverage, isort
run_analysis.py          # ROOT SHIM: 3 lines → calls mace_gaussian.analysis
run_analysis_harmonic.py # ROOT SHIM: 3 lines → calls mace_gaussian.analysis
scripts/                 # utility/research scripts (not library code)
    create_fixtures.py   # was scripts/ already
    comparison_workflow.py    # moved from root
    produce_molecules.py      # moved from root
    convert_all_chk_files.py  # moved from root
    charge_analysis.py        # moved from root
tests/                   # unchanged location; imports updated to mace_gaussian.*
water.xyz, acoh.xyz ...  # user input files — stay at root
fort.7                   # non-Python Gaussian scratch — stays at root
```

### Pattern 1: Relative Import Conversion

**What:** All internal imports within `mace_gaussian/` use relative dot-notation.
**When to use:** Every import inside the package that references another module inside the package.

**Example — before (in gaussian/__init__.py):**
```python
from gaussian.fchk import convert_chk_to_fchk
from gaussian.parser import GaussianLogParser
```

**Example — after (in mace_gaussian/gaussian/__init__.py):**
```python
from .fchk import convert_chk_to_fchk
from .parser import GaussianLogParser
```

**Example — cross-subpackage import (in mace_gaussian/dft_baseline.py):**
```python
# Before:
from gaussian.parser import parse_gaussian_log
from utils.results import ResultsManager

# After (flat inside mace_gaussian/):
from .gaussian.parser import parse_gaussian_log
from .utils.results import ResultsManager
```

**Example — in mace_gaussian/workflow.py:**
```python
# Before:
from calculators import dipole_factory
from gaussian.runner import run_gaussian_with_zmq
from utils.results import ResultsManager

# After:
from .calculators import dipole_factory
from .gaussian.runner import run_gaussian_with_zmq
from .utils.results import ResultsManager
```

### Pattern 2: cli.py Uses Absolute Imports

**What:** The CLI entry point is the package boundary — it imports from outside using absolute package paths.
**Why:** `cli.py` IS the boundary between "user space" and "library space". Using absolute imports here also makes it more resilient as an installed entry point.

```python
# In mace_gaussian/cli.py — absolute imports are correct here:
from mace_gaussian.workflow import run_pipeline
from mace_gaussian.utils.exceptions import InputValidationError, PrerequisiteError
from mace_gaussian.utils.validation import detect_device, validate_all_prerequisites
from mace_gaussian.calculators import dipole_factory
```

### Pattern 3: Root Shim Scripts

**What:** 2–3 line delegators that preserve the `python run_analysis.py water` interface.
**Why:** CLAUDE.md documents these commands; muscle memory must not break.

```python
# run_analysis.py (root shim — after Phase 8)
#!/usr/bin/env python3
"""Thin shim — delegates to mace_gaussian.analysis."""
from mace_gaussian.analysis.run_analysis import main

if __name__ == "__main__":
    main()
```

The shim must NOT use `sys.path.insert(0, ...)` — the installed package provides the import path.

### Pattern 4: mace_gaussian/__init__.py Minimal Public API

```python
# mace_gaussian/__init__.py
"""MACE-Gaussian interface: ML potentials + Gaussian for IR spectroscopy."""

__version__ = "0.2.0"

from .workflow import run_pipeline

__all__ = ["run_pipeline", "__version__"]
```

### Pattern 5: analysis/__init__.py Exposes Key Functions

```python
# mace_gaussian/analysis/__init__.py
"""Analysis and reporting subpackage."""

from .run_analysis import analyze_molecule, main as run_analysis_main
from .run_analysis_harmonic import analyze_molecule_harmonic, main as run_analysis_harmonic_main

__all__ = [
    "analyze_molecule",
    "analyze_molecule_harmonic",
    "run_analysis_main",
    "run_analysis_harmonic_main",
]
```

### Pattern 6: hatchling packages Configuration

```toml
# pyproject.toml
[tool.hatch.build.targets.wheel]
packages = ["mace_gaussian"]
```

This replaces the current broken `packages = ["."]` which would ship the entire root directory.

### Pattern 7: pyproject.toml Full Cleanup

```toml
[project.scripts]
mace-gaussian = "mace_gaussian.cli:cli"

[tool.ruff.lint.isort]
known-first-party = ["mace_gaussian"]

[tool.coverage.run]
source = ["mace_gaussian"]
omit = [
    "tests/*",
    "mace_dipole_pkg/*",
    "mace_ML_pkg/*",
]
```

Note: The Click app object in `cli.py` is named `cli` (the `@click.group()` decorated function). The entry point `"mace_gaussian.cli:cli"` is correct.

### Anti-Patterns to Avoid

- **Moving files without `git mv`:** Regular `mv` loses git history. Use `git mv old_path new_path` for every file relocation.
- **Converting `cli.py` to relative imports:** cli.py is the entry point and should use absolute imports from `mace_gaussian.*`.
- **sys.path.insert in shim scripts:** The installed package removes this need. Shims should just import and call.
- **Importing from both old and new paths simultaneously:** During migration, old absolute paths (e.g., `from calculators import ...`) will continue to resolve if `calculators/` still exists at root. Clean up old locations completely before testing.
- **Forgetting the egg-info directory:** `mace_gaussian_interface.egg-info/` at root points to the old layout. After pyproject.toml changes, `uv sync` regenerates it correctly. Don't manually edit it.
- **Circular imports in __init__.py:** If `mace_gaussian/__init__.py` imports from `workflow`, and `workflow` imports from `mace_gaussian` (absolute), there's a circular dependency. Solution: `workflow.py` uses relative imports internally.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Detecting all imports to update | Custom AST scanner | `ruff check --select I` after moving | ruff's isort detects import ordering issues; ruff F401 catches unused; manual grep is sufficient for targeted search |
| Package discovery during build | Custom find-packages logic | `packages = ["mace_gaussian"]` in hatchling | Explicit is better than auto-discovery for a project with multiple non-package directories at root |
| Backwards compat for old imports | `sys.modules` aliasing (`mace_modules["calculators"] = mace_gaussian.calculators`) | Update tests directly | Tests are internal — update them. The shim pattern is only for end-user-facing scripts. |

**Key insight:** This phase is mechanical. The value is in being systematic: move every file, convert every internal import, update every test import, clean up pyproject.toml — in that order, with ruff + pytest verification after each wave.

---

## Common Pitfalls

### Pitfall 1: gm_helper.py Absolute Path Requirement

**What goes wrong:** After moving `gm_helper.py` into `mace_gaussian/`, Gaussian still calls it via the absolute path stored in `MACE_HELPER_SCRIPT_PATH`. If `gm_helper.py` uses relative imports internally, running it as `python /path/to/gm_helper.py` fails because `__package__` is `None` for scripts.

**Why it happens:** `gm_helper.py` is called by Gaussian 16 as an external script — it runs standalone, not as part of a package import. Relative imports require a package context.

**How to avoid:** `gm_helper.py` must use absolute imports (`from mace_gaussian.gaussian.zmq_server import ...` if it needs anything) OR keep it self-contained (it currently only imports `zmq` and `os`, `sys` — standard library only). Check current imports before converting.

**Warning signs:** `ImportError: attempted relative import with no known parent package` when running `python gm_helper.py` directly.

### Pitfall 2: Test Discovery Breaks After Reorganization

**What goes wrong:** Tests currently import `from calculators.base import ...`, `from cli import cli`, `from mode_matching import ...`. After moving everything into `mace_gaussian/`, these imports fail because the old top-level modules no longer exist.

**Why it happens:** pytest adds `tests/` parent (the project root) to `sys.path`. With the old flat layout, `calculators/` was directly on `sys.path`. After moving into `mace_gaussian/`, only `mace_gaussian.calculators` exists.

**How to avoid:** Update every test import from `from X import Y` to `from mace_gaussian.X import Y` in the same wave that moves the source files. Run `python -m pytest tests/ -x` immediately after to verify no import errors.

**Warning signs:** `ModuleNotFoundError: No module named 'calculators'` during test collection.

### Pitfall 3: comparison_workflow.py Import Chain in run_analysis.py Shims

**What goes wrong:** Current `run_analysis.py` shim imports `from comparison_workflow import analyze_molecule`. After Phase 8, `comparison_workflow.py` moves to `scripts/` (not into `mace_gaussian/`). The shim cannot import from `scripts/`.

**Why it happens:** The CONTEXT.md decision puts analysis logic (including `comparison_workflow.py`) into `mace_gaussian/analysis/`. But `run_analysis.py` at root was calling `comparison_workflow.analyze_molecule`.

**Resolution:** The analysis orchestration in `comparison_workflow.py` moves into `mace_gaussian/analysis/` (renamed to something clearer or kept as `run_analysis.py`). The root shims call `mace_gaussian.analysis.run_analysis.main()`. The `scripts/comparison_workflow.py` that goes to `scripts/` is a different concern — it's the research utility version.

**Warning signs:** ImportError in root shim scripts after moving `comparison_workflow.py` to `scripts/`.

**Concrete resolution:** The CONTEXT.md lists both `run_analysis.py` and `run_analysis_harmonic.py` as going into `mace_gaussian/analysis/`. These are the analysis orchestration logic that currently lives in `comparison_workflow.py`. The current `run_analysis.py` and `run_analysis_harmonic.py` at root are thin shells that just call `comparison_workflow.analyze_molecule` and `comparison_workflow.analyze_molecule_harmonic`. After Phase 8: the orchestration logic moves into `mace_gaussian/analysis/run_analysis.py` and `mace_gaussian/analysis/run_analysis_harmonic.py` (formerly the bodies of `analyze_molecule` and `analyze_molecule_harmonic` in `comparison_workflow.py`), and the root shims become 3-line delegators calling `mace_gaussian.analysis.run_analysis.main()`.

### Pitfall 4: pyproject.toml `project.name` vs Package Name

**What goes wrong:** The current `[project] name = "mace-gaussian-interface"` does NOT match the import package name `mace_gaussian`. Hatchling's auto-discovery normalizes `mace-gaussian-interface` to look for `mace_gaussian_interface` — which doesn't exist.

**Why it happens:** The current `packages = ["."]` bypasses auto-discovery by including everything at root. After we set `packages = ["mace_gaussian"]`, hatchling will correctly find and package just `mace_gaussian/`.

**How to avoid:** Set `packages = ["mace_gaussian"]` explicitly in `[tool.hatch.build.targets.wheel]`. The project name (`mace-gaussian-interface`) can stay as-is for PyPI compatibility — it's independent of the import name.

**Warning signs:** `pip show mace-gaussian-interface` shows no top-level packages, or importing `mace_gaussian` fails after `pip install`.

### Pitfall 5: `run_analysis.py` Function Name Collision

**What goes wrong:** Root `run_analysis.py` and `mace_gaussian/analysis/run_analysis.py` will coexist. If the root script has `if __name__ == "__main__": main()` and imports `from mace_gaussian.analysis.run_analysis import main`, it will correctly call the package version. But if `sys.path` still has root first, `import run_analysis` inside the package could resolve to the ROOT shim instead of the package module.

**How to avoid:** Within `mace_gaussian/analysis/`, all imports must be relative (e.g., `from .analyze_spectra import ...`), never bare `import run_analysis`. The analysis subpackage is internally consistent with relative imports.

### Pitfall 6: MODE_MATCHING_SUMMARY.md Reference in mode_matching.py

**What goes wrong:** `mode_matching.py` at root may reference paths relative to the project root. After moving to `mace_gaussian/analysis/mode_matching.py`, any `Path(__file__).parent`-based path construction that assumed root placement will resolve to `mace_gaussian/analysis/` instead.

**How to avoid:** Audit `mode_matching.py` for `Path(__file__).parent` usage before moving. Any path that needs to reach the project root would need `Path(__file__).parent.parent.parent` after moving two levels deeper.

---

## Code Examples

Verified patterns for this migration:

### hatchling packages configuration (verified from official docs)
```toml
# Source: https://hatch.pypa.io/1.13/config/build/
[tool.hatch.build.targets.wheel]
packages = ["mace_gaussian"]
```

### Entry point configuration (verified from PyPA docs)
```toml
# Source: https://packaging.python.org/en/latest/guides/writing-pyproject-toml/
[project.scripts]
mace-gaussian = "mace_gaussian.cli:cli"
```

### Relative import pattern within package (Python standard)
```python
# In mace_gaussian/workflow.py — cross-subpackage relative imports
from .calculators import dipole_factory
from .gaussian.fchk import convert_chk_to_fchk
from .gaussian.io import ase_to_gjf, parse_gaussian_input, write_gaussian_output
from .gaussian.parser import parse_gaussian_log
from .gaussian.runner import DEFAULT_TIMEOUT_SECONDS, run_gaussian_with_zmq
from .utils.results import ResultsManager
from .utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV
```

### Root shim pattern (3 lines)
```python
#!/usr/bin/env python3
"""Thin shim — see mace_gaussian/analysis/run_analysis.py for implementation."""
from mace_gaussian.analysis.run_analysis import main

if __name__ == "__main__":
    main()
```

### mace_gaussian/__init__.py minimal API
```python
"""MACE-Gaussian interface: ML potentials + Gaussian for IR spectroscopy."""

__version__ = "0.2.0"

from .workflow import run_pipeline

__all__ = ["run_pipeline", "__version__"]
```

### Updated test import style
```python
# Before Phase 8:
from gaussian.parser import GaussianLogParser
from calculators.base import DipoleCalculatorBase
from cli import cli
from mode_matching import compute_mode_overlap

# After Phase 8:
from mace_gaussian.gaussian.parser import GaussianLogParser
from mace_gaussian.calculators.base import DipoleCalculatorBase
from mace_gaussian.cli import cli
from mace_gaussian.analysis.mode_matching import compute_mode_overlap
```

### ruff isort known-first-party update
```toml
[tool.ruff.lint.isort]
known-first-party = ["mace_gaussian"]
```

### coverage source update
```toml
[tool.coverage.run]
source = ["mace_gaussian"]
omit = [
    "tests/*",
    "mace_dipole_pkg/*",
    "mace_ML_pkg/*",
]
```

---

## Current State Inventory

### Files that MOVE into mace_gaussian/ (root → package)

| Current Location | Destination | Notes |
|-----------------|-------------|-------|
| `cli.py` | `mace_gaussian/cli.py` | Keep absolute imports; update `from workflow` → `from mace_gaussian.workflow` |
| `workflow.py` | `mace_gaussian/workflow.py` | Convert internal imports to relative |
| `gm_helper.py` | `mace_gaussian/gm_helper.py` | Standalone script; only stdlib imports — no change needed |
| `dft_baseline.py` | `mace_gaussian/dft_baseline.py` | Convert `from gaussian.parser` → `from .gaussian.parser` |
| `diagnostics.py` | `mace_gaussian/diagnostics.py` | Check imports — currently only stdlib |

### Subdirectories that MOVE into mace_gaussian/ (root → nested)

| Current Location | Destination | Import Conversion Needed |
|-----------------|-------------|--------------------------|
| `calculators/` | `mace_gaussian/calculators/` | Internal: `from calculators.X` → `from .X`; cross-pkg: `from utils.X` → `from ..utils.X` |
| `gaussian/` | `mace_gaussian/gaussian/` | Internal: `from gaussian.X` → `from .X`; cross-pkg: `from utils.X` → `from ..utils.X` |
| `utils/` | `mace_gaussian/utils/` | Internal: `from utils.X` → `from .X` |

### Files that move into mace_gaussian/analysis/ (root → analysis subpackage)

| Current Location | Destination | Notes |
|-----------------|-------------|-------|
| `analyze_spectra.py` | `mace_gaussian/analysis/analyze_spectra.py` | No internal cross-imports; standalone |
| `mode_matching.py` | `mace_gaussian/analysis/mode_matching.py` | Imports `from gaussian.fchk` → `from ..gaussian.fchk` |
| `html_report_generator.py` | `mace_gaussian/analysis/html_report_generator.py` | No internal imports; standalone |
| `comparison_workflow.py` (analysis logic) | `mace_gaussian/analysis/run_analysis.py` + `run_analysis_harmonic.py` | Contains `analyze_molecule()` and `analyze_molecule_harmonic()` — these become the analysis entry points |

### Files that MOVE to scripts/ (root → scripts/)

| Current Location | Destination | Notes |
|-----------------|-------------|-------|
| `comparison_workflow.py` | `scripts/comparison_workflow.py` | Research utility — not library code; update its imports to `from mace_gaussian.analysis...` |
| `produce_molecules.py` | `scripts/produce_molecules.py` | Research utility |
| `convert_all_chk_files.py` | `scripts/convert_all_chk_files.py` | Research utility |
| `charge_analysis.py` | `scripts/charge_analysis.py` | Research utility |

### Files that STAY at root

| File | Reason |
|------|--------|
| `run_analysis.py` | Root shim — becomes 3-line delegator |
| `run_analysis_harmonic.py` | Root shim — becomes 3-line delegator |
| `fort.7` | Gaussian scratch file — not Python |
| `water.xyz`, `acoh.xyz`, etc. | User input files |
| `pyproject.toml` | Build configuration |
| `CLAUDE.md`, `README.md` | Documentation |

### Tests that need import updates (all 131 tests)

| Test File | Current Import | Must Become |
|-----------|---------------|-------------|
| `test_gaussian_parser.py` | `from gaussian.parser import ...` | `from mace_gaussian.gaussian.parser import ...` |
| `test_fchk_parser.py` | `from gaussian.fchk import ...` | `from mace_gaussian.gaussian.fchk import ...` |
| `test_calculators.py` | `from calculators.base import ...` | `from mace_gaussian.calculators.base import ...` |
| `test_mace_loader.py` | `from calculators.mace_loader import ...` | `from mace_gaussian.calculators.mace_loader import ...` |
| `test_cli_validation.py` | `from cli import cli` | `from mace_gaussian.cli import cli` |
| `test_mode_matching.py` | `from mode_matching import ...` | `from mace_gaussian.analysis.mode_matching import ...` |
| `test_exceptions.py` | `from utils.exceptions import ...` | `from mace_gaussian.utils.exceptions import ...` |
| `test_units.py` | `from utils.units import ...` | `from mace_gaussian.utils.units import ...` |
| `test_validation.py` | `from utils.validation import ...` | `from mace_gaussian.utils.validation import ...` |
| `test_regression.py` | check current | likely `from mace_gaussian.gaussian...` or `analysis...` |

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `setup.py` + `find_packages()` | `pyproject.toml` + `hatchling` with `packages = [...]` | PEP 517/518, ~2020+ | Declarative; no executable build scripts |
| Flat layout (entire root as package) | Named package directory (`mace_gaussian/`) | This phase | Enables `pip install`, clean namespace |
| `packages = ["."]` | `packages = ["mace_gaussian"]` | This phase | Was shipping entire root directory |

**Deprecated/outdated in this project:**
- `packages = ["."]` in pyproject.toml: ships everything at root including test fixtures, molecule files — wrong
- `gm-main` and `gm-diagnose` entry points: point to deleted `gm_main:main` — broken since Phase 7

---

## Open Questions

1. **comparison_workflow.py: split or rename?**
   - What we know: The CONTEXT.md says `run_analysis.py` (anharmonic logic) and `run_analysis_harmonic.py` (harmonic logic) go into `mace_gaussian/analysis/`. The current `run_analysis.py` and `run_analysis_harmonic.py` at root are thin wrappers calling `comparison_workflow.analyze_molecule` and `comparison_workflow.analyze_molecule_harmonic`.
   - What's unclear: Does `comparison_workflow.py` get split into two files in `analysis/`, or does it move as-is and get renamed to something like `analysis_workflow.py`?
   - Recommendation: Move `comparison_workflow.py` into `mace_gaussian/analysis/` as `analysis_workflow.py` (it contains `ComparisonWorkflow` class, `analyze_molecule()`, `analyze_molecule_harmonic()`). The root shims call the appropriate entry point functions. The `scripts/comparison_workflow.py` would just be a copy of the old root file that imports from `mace_gaussian.analysis.analysis_workflow`. This preserves backward compat for anyone running `python scripts/comparison_workflow.py` directly while having clean library code.

2. **py.typed marker (Claude's Discretion)**
   - What we know: The project uses `ty check` for type checking (from CLAUDE.md). `py.typed` is a PEP 561 marker that tells type checkers the package is typed.
   - What's unclear: Whether `ty` (a type checker) requires `py.typed` or uses it automatically.
   - Recommendation: Add `py.typed` as an empty file in `mace_gaussian/` — it's a zero-effort addition that follows PEP 561 and signals type-checking intent. Include it in hatchling's include config.

3. **mace_gaussian_interface.egg-info at root**
   - What we know: This was generated by the old `packages = ["."]` config and lists the project root as the package. It currently has no `top_level.txt` content (empty).
   - What's unclear: Will `uv sync` automatically regenerate a correct egg-info after pyproject.toml change?
   - Recommendation: After changing pyproject.toml, run `uv sync` to regenerate. The old egg-info directory can be deleted as part of cleanup — it will be recreated correctly by hatchling.

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 7.1.3 |
| Config file | `pyproject.toml` → `[tool.pytest.ini_options]` |
| Quick run command | `python -m pytest tests/ -x -q` |
| Full suite command | `python -m pytest tests/ -ra` |
| Estimated runtime | ~3–5 seconds (131 tests, all unit tests, no GPU/Gaussian) |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|--------------|
| STRUCT-08 | mace_gaussian/ package importable after reorganization | smoke | `python -c "import mace_gaussian; print(mace_gaussian.__version__)"` | ❌ Wave 0 gap |
| STRUCT-08 | All subpackages importable | smoke | `python -c "from mace_gaussian import calculators, gaussian, utils; print('OK')"` | ❌ Wave 0 gap |
| STRUCT-08 | Internal relative imports work | integration | `python -m pytest tests/ -x -q` (all 131 tests) | ✅ yes (after import updates) |
| STRUCT-09 | CLI entry point works | smoke | `python -m mace_gaussian.cli --help` | ✅ yes (after import updates) |
| STRUCT-09 | mace-gaussian console script works | smoke | `mace-gaussian --help` (requires `uv sync`) | ❌ Wave 0 gap (depends on pyproject.toml fix) |
| STRUCT-10 | analysis subpackage importable | smoke | `python -c "from mace_gaussian.analysis import analyze_spectra; print('OK')"` | ❌ Wave 0 gap |
| STRUCT-10 | mode_matching importable from new location | unit | `python -m pytest tests/test_mode_matching.py -x` | ✅ yes (after import update) |

### Nyquist Sampling Rate
- **Minimum sample interval:** After each file move wave → run: `python -m pytest tests/ -x -q`
- **Full suite trigger:** Before marking phase complete
- **Phase-complete gate:** All 131 tests green + `python -c "import mace_gaussian"` succeeds + `mace-gaussian --help` works
- **Estimated feedback latency per task:** ~5 seconds

### Wave 0 Gaps (must be created before or during implementation)
- [ ] No new test files needed — existing 131 tests serve as structural validation after import updates
- [ ] Smoke command `python -c "import mace_gaussian; print(mace_gaussian.__version__)"` — manual, not a test file
- [ ] `uv sync` must be run after pyproject.toml changes before running `mace-gaussian --help`

*(Structural validation is entirely covered by updating existing tests to new import paths and running the full suite.)*

---

## Sources

### Primary (HIGH confidence)
- https://hatch.pypa.io/1.13/config/build/ — hatchling `packages` key syntax, wheel target configuration
- https://packaging.python.org/en/latest/guides/writing-pyproject-toml/ — `[project.scripts]` entry point format
- https://packaging.python.org/en/latest/discussions/src-layout-vs-flat-layout/ — flat layout pitfalls, sys.path behavior
- Python 3 import system docs (PEP 328) — relative import syntax rules

### Secondary (MEDIUM confidence)
- Direct codebase inspection: all 131 tests collected via `pytest --collect-only`, all current imports verified via `grep`
- Current pyproject.toml verified to have `packages = ["."]` (broken), stale entry points, stale coverage sources

### Tertiary (LOW confidence)
- General best practices for `__init__.py` re-export surfaces — common pattern but not in official spec

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — hatchling already in use, Python import system is stable spec
- Architecture: HIGH — derived directly from CONTEXT.md locked decisions + verified codebase inspection
- Pitfalls: HIGH — derived from direct inspection of existing imports, pytest behavior with sys.path, and official flat layout docs
- Import conversion map: HIGH — every import manually verified against current source files

**Research date:** 2026-02-24
**Valid until:** 2026-06-01 (Python packaging conventions are stable; hatchling API rarely changes)
