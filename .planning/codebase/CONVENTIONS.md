# Coding Conventions

**Analysis Date:** 2026-02-16

## Naming Patterns

**Files:**
- Python scripts use `snake_case` with `.py` extension
- Main entry points: `cli.py`, `gm_main.py` (module-level)
- Functional scripts: `run_analysis.py`, `run_analysis_harmonic.py`, `dft_baseline.py`
- Utility/library modules: `gaussian_parser.py`, `fchk_parser.py`, `mace_calculators.py`
- Analysis/workflow modules: `analyze_spectra.py`, `comparison_workflow.py`, `mode_matching.py`
- Manager modules: `results_manager.py`, `charge_analysis.py`

**Functions:**
- Use `snake_case` for function names
- Public functions: lowercase starting word (e.g., `run_workflow()`, `parse_gaussian_input()`)
- Private/internal functions: prefix with underscore (e.g., `_check_availability()`)
- Factory/builder patterns use `_<action>_<object>()` convention (e.g., `_register_calculators()`)
- Context managers use `@contextmanager` decorator (e.g., `zmq_server()`)

**Variables:**
- Use `snake_case` for local variables and parameters
- Constants: `UPPER_CASE` (e.g., `BOHR_TO_ANGSTROM`, `DEFAULT_HELPER_SCRIPT`)
- Boolean flags: prefixed with `is_`, `has_`, `should_`, `can_` (e.g., `is_calc_finished()`, `has_opt()`)
- Private/module-level constants: leading underscore (e.g., `_mace_params`)

**Types:**
- Class names: `PascalCase` (e.g., `DipoleCalculatorBase`, `EspalomaDipoleCalculator`)
- Abstract base classes: prefix `Base` (e.g., `DipoleCalculatorBase`)
- Data classes: use `@dataclass` decorator with descriptive names (e.g., `SpectrumData`, `ComparisonMetrics`)
- Enum members: `UPPER_CASE`
- Type hints: used throughout for function signatures and class attributes

## Code Style

**Formatting:**
- Line length: **100 characters** (enforced by ruff)
- Indentation: 4 spaces
- String quotes: double quotes preferred for consistency
- Trailing commas: used in multi-line function signatures for readability

**Linting:**
- Tool: `ruff` (version >= 0.9.0)
- Configuration: `pyproject.toml` with `[tool.ruff]` section
- Rules enabled: `["E", "F", "W", "I", "UP"]` (style, logical, warnings, imports, upgrades)
- Key settings:
  - `line-length = 100`
  - `target-version = "py39"`

**Type Checking:**
- Tool: `ty` (type checker, version `0.0.1a1`)
- Types used in function signatures and docstrings
- Optional types: `Optional[Type]` for nullable values
- Union types: documented in docstrings when needed

## Import Organization

**Order:**
1. Standard library imports (e.g., `import os`, `import sys`, `from pathlib import Path`)
2. Third-party imports (e.g., `import numpy as np`, `import torch`, `from ase.io import read`)
3. Local application imports (modules in same project)
4. Conditional imports wrapped in try/except for optional features

**Path Aliases:**
- Configured in `pyproject.toml` under `[tool.ruff.lint.isort]`
- Known first-party modules: `mace_calculators`, `gm_main`, `gaussian_parser`, `analyze_spectra`, `comparison_workflow`, `mode_matching`, `fchk_parser`
- No import aliases used; full module paths preferred for clarity

**Style Examples:**
```python
# Standard library first
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Third-party
import numpy as np
import torch
from ase.io import read
from ase.atoms import Atoms
import logging

# Local
from mace_calculators import MACEDipoleCalculator
from gaussian_parser import parse_gaussian_log
```

## Error Handling

**Patterns:**
- `try/except/finally` blocks with specific exception types
- Exceptions caught with narrowest type possible to avoid masking unrelated errors
- Generic exceptions caught at CLI boundaries (main entry points)
- Custom error messages include context (variable names, file paths)

**Common Patterns:**

```python
# Pattern 1: Graceful degradation with fallback
try:
    dipole, partial_charges = dipole_calc.calculate_dipole(atoms)
except Exception as e:
    logger.error(f"Dipole calculation failed: {e}")
    logger.warning("Falling back to zero dipole (IR intensities will be zero)")
    dipole = np.zeros(3)
    dipole_derivatives = np.zeros((3 * natoms, 3))
    return dipole, dipole_derivatives, None

# Pattern 2: Resource cleanup with context managers
with zmq_server(".ipc_file") as socket:
    while not is_calc_finished(proc, socket):
        msg = socket.recv_string()
        # Process message
        socket.send_string("ready")

# Pattern 3: File existence checks before operations
if not Path(log_file).exists():
    raise FileNotFoundError(f"Log file not found: {log_file}")

# Pattern 4: Validation with early return
if not calc:
    raise ValueError(f"Unknown dipole calculator: {method}")
if not calc.available:
    raise RuntimeError(f"Dipole calculator '{method}' not available")
```

## Logging

**Framework:** Python's standard `logging` module

**Setup pattern** (module level):
```python
import logging

logger = logging.getLogger(__name__)
```

**Configuration** (entry points like `cli.py`, `gm_main.py`):
```python
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
```

**When to Log:**
- `logger.info()`: Major workflow steps (e.g., "Geometry optimization started", "✓ Espaloma-charge dipole calculator available")
- `logger.debug()`: Detailed diagnostic info (e.g., "Espaloma dipole: [0.1, 0.2, 0.3]")
- `logger.warning()`: Recoverable issues (e.g., "Charge analysis module not available", "Could not convert .chk to .fchk")
- `logger.error()`: Failure states that block computation (e.g., "Energy/forces calculation failed: {e}")

**Examples:**
```python
logger.info(f"✓ Espaloma-charge dipole calculator available and tested")
logger.warning(f"⚠ xTB dipole calculator failed: {e}")
logger.error(f"Frequency calculation failed: {e}")
logger.debug(f"Espaloma dipole: {dipole}, charges sum: {np.sum(charges):.6f}")
```

## Comments

**When to Comment:**
- Complex algorithms explained with algorithm name (e.g., "Central difference derivative")
- Non-obvious unit conversions (e.g., "Convert from e*Angstrom to e*Bohr (Gaussian units)")
- Workarounds for known issues (e.g., "Temporarily set torch default dtype to float32")
- Important constraints or gotchas (marked with `# NOTE:` or `# WARNING:`)
- Deprecated patterns with migration guidance

**JSDoc/TSDoc:**
Python uses docstrings (not JSDoc). Convention:

```python
def calculate_dipole_derivatives(
    atoms, displacement=0.01, **kwargs
) -> np.ndarray:
    """
    Calculate dipole derivatives numerically.

    Uses central difference formula: (f(x+h) - f(x-h)) / (2h)

    Parameters
    ----------
    atoms : ase.Atoms
        Molecular structure
    displacement : float
        Finite difference step size in Angstroms (default: 0.01)
    **kwargs
        Additional arguments passed to calculator

    Returns
    -------
    np.ndarray
        Dipole derivatives, shape (3*natoms, 3)

    Raises
    ------
    Exception
        If calculation fails, returns zero array
    """
```

**Module-level docstrings:**
All non-trivial modules have docstrings at the top:

```python
"""
Gaussian log file parser for extracting frequencies and IR intensities.
"""
```

## Function Design

**Size:**
- Functions should be focused on a single responsibility
- Typical range: 10-50 lines
- Complex workflows broken into smaller functions with descriptive names
- Multi-step processes documented with numbered steps in comments

**Parameters:**
- Use type hints (e.g., `atoms: Atoms`, `natoms: int`)
- Default values for optional parameters (e.g., `force_optimization: bool = False`)
- Keyword-only arguments for clarity when function has many parameters
- Document parameter constraints in docstrings (e.g., "must be positive")

**Return Values:**
- Single return values preferred; tuples only for related values
- Document return types with `->` annotation (e.g., `-> Tuple[np.ndarray, Optional[np.ndarray]]`)
- Use dataclasses for multiple related return values (`@dataclass` pattern)

**Examples:**
```python
# Pattern 1: Tuple return for related multi-values
def parse_gaussian_input(infile: str) -> Tuple[int, int, int, int, np.ndarray, list]:
    """Returns (natoms, deriv, charge, spin, coordinates, atomnames)"""

# Pattern 2: Dataclass for structured returns
@dataclass
class ComparisonMetrics:
    mae_freq: float
    rmse_freq: float
    r2_freq: float
    num_peaks: int

# Pattern 3: Single return with type hint
def extract_mode_data_from_checkpoint(chk_or_fchk_file: str) -> Tuple[...]:
    """Returns modes array with detailed documentation"""
```

## Module Design

**Exports:**
- All public functions and classes at module level
- Private functions prefixed with underscore
- Modules act as cohesive units (e.g., `gaussian_parser.py` handles all Gaussian parsing)

**Barrel Files:**
- Not used in this codebase
- Each module is imported directly by name

**Module Responsibilities:**
- `mace_calculators.py`: MACE calculator wrapper and dipole factory pattern
- `gaussian_parser.py`: Parsing Gaussian output files
- `fchk_parser.py`: Parsing Gaussian checkpoint files
- `gm_main.py`: Main workflow orchestration
- `gm_helper.py`: Helper script for Gaussian external interface
- `analyze_spectra.py`: Spectral analysis and comparison
- `mode_matching.py`: Vibrational mode matching via eigenvector overlap
- `comparison_workflow.py`: End-to-end analysis pipeline
- `results_manager.py`: Output directory and metadata management
- `charge_analysis.py`: Partial charge analysis and visualization
- `dft_baseline.py`: DFT reference calculations
- `diagnostics.py`: Environment and dependency checking
- `cli.py`: Command-line interface

---

*Convention analysis: 2026-02-16*
