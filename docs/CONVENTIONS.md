# Code Conventions

Pragmatic conventions reflecting actual codebase patterns. Reference, not textbook.

## Naming

- **Functions/variables:** `snake_case`
- **Classes:** `PascalCase`
- **Constants:** `UPPER_SNAKE_CASE`
- **Private:** `_leading_underscore`
- **Modules:** `snake_case.py`, one class per module for large classes

## Units

Internal representation:

| Quantity  | Unit      |
|-----------|-----------|
| Energy    | eV        |
| Distance  | Angstrom  |
| Frequency | cm-1      |
| Dipole    | e*Bohr    |

All conversions go through `utils.units`. Never define constants inline.

```python
# Good
from utils.units import HARTREE_TO_EV
energy_ev = energy_hartree * HARTREE_TO_EV

# Bad
energy_ev = energy_hartree * 27.2114
```

Constants use CODATA 2018 recommended values.

## Error Handling

Use the exception hierarchy from `utils.exceptions`:

- `MaceGaussianError` -- base for all project exceptions
- `InputValidationError` -- bad user input (files, parameters)
- `GaussianParseError` -- Gaussian output parsing failures
- `PrerequisiteError` -- missing tools, models, environment
- `CUDANotAvailableWarning` -- GPU not found (warning, not error)

Catch `MaceGaussianError` for broad handling, specific subclasses for targeted handling.
Never use bare `except:`. Prefer `except SpecificError`.

## Imports

Order (enforced by ruff isort):

1. Standard library (`import os`, `from pathlib import Path`)
2. Third-party (`import numpy as np`, `from ase.io import read`)
3. Local (`from utils.units import HARTREE_TO_EV`, `from gaussian_parser import ...`)

Utils package imports use `from utils.X import Y`. Heavy dependencies (torch, mace, DGL)
use lazy imports inside functions to keep module loading fast.

## File Organization

- Related pure functions grouped thematically (e.g., `utils/units.py`)
- Large classes get their own module (e.g., `utils/results.py`)
- Tests mirror source structure: `utils/units.py` -> `tests/test_units.py`

## Code Quality

- **Linting/formatting:** ruff (line-length 100)
- **Type hints:** On public function signatures
- **Docstrings:** NumPy style on public functions
- **Tests:** pytest, class-based grouping, descriptive method names
