---
phase: 08-package-structure-reorganization
plan: 01
subsystem: infra
tags: [python-packaging, imports, refactor, git-mv]

# Dependency graph
requires:
  - phase: 07-extract-workflow-orchestrator
    provides: calculators/, gaussian/, utils/ subpackages at root level
provides:
  - mace_gaussian/ package namespace with __version__ = "0.2.0"
  - mace_gaussian/calculators/, mace_gaussian/gaussian/, mace_gaussian/utils/ subpackages
  - All internal imports converted to relative dot-notation
affects:
  - 08-02 (workflow.py move into mace_gaussian/)
  - 08-03 (cli.py and root-level callers updating imports to mace_gaussian.*)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Relative import pattern: intra-package uses from .module, cross-package uses from ..package.module"
    - "git mv for history-preserving directory relocation"

key-files:
  created:
    - mace_gaussian/__init__.py
  modified:
    - mace_gaussian/calculators/__init__.py
    - mace_gaussian/calculators/factory.py
    - mace_gaussian/calculators/espaloma.py
    - mace_gaussian/calculators/mace_ml.py
    - mace_gaussian/calculators/xtb.py
    - mace_gaussian/gaussian/__init__.py
    - mace_gaussian/gaussian/io.py
    - mace_gaussian/gaussian/parser.py
    - mace_gaussian/gaussian/runner.py
    - mace_gaussian/gaussian/fchk.py
    - mace_gaussian/utils/__init__.py
    - mace_gaussian/utils/validation.py

key-decisions:
  - "workflow.run_pipeline import commented out in mace_gaussian/__init__.py — workflow.py not yet inside package (Plan 02 activates it)"
  - "ruff --select I,F401 --fix applied to auto-sort 2 import blocks after manual relative import conversion"
  - "mace_loader.py lazy import inside _check_availability converted to from .mace_loader — same lazy semantics, now relative"

patterns-established:
  - "Cross-subpackage imports use two-dot relative: from ..utils.units import X"
  - "Intra-subpackage imports use one-dot relative: from .base import X"

requirements-completed: [STRUCT-08]

# Metrics
duration: 2min
completed: 2026-02-24
---

# Phase 8 Plan 01: Package Foundation Summary

**mace_gaussian/ package created via git mv with all three subpackages (calculators, gaussian, utils) relocated and internal imports converted to relative dot-notation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-24T23:09:23Z
- **Completed:** 2026-02-24T23:11:51Z
- **Tasks:** 2
- **Files modified:** 13 (1 created + 12 modified)

## Accomplishments
- Created `mace_gaussian/` package directory with `__init__.py` exposing `__version__ = "0.2.0"`
- Relocated all three subpackages via `git mv` (17 renamed files, git history preserved)
- Converted all 16 internal imports across 12 files from absolute to relative dot-notation
- Verified `from mace_gaussian.calculators`, `from mace_gaussian.gaussian`, and `from mace_gaussian.utils` all import successfully

## Task Commits

Each task was committed atomically:

1. **Task 1: Create mace_gaussian/ package and move subpackages via git mv** - `3ee4105` (chore)
2. **Task 2: Convert all internal imports to relative within the three moved subpackages** - `76ddb1e` (feat)

**Plan metadata:** (docs commit below)

## Files Created/Modified
- `mace_gaussian/__init__.py` - Package entry point, `__version__ = "0.2.0"`, run_pipeline import commented pending Plan 02
- `mace_gaussian/calculators/__init__.py` - from calculators.* → from .*
- `mace_gaussian/calculators/factory.py` - from calculators.* → from .*
- `mace_gaussian/calculators/espaloma.py` - from .base, from ..utils.units
- `mace_gaussian/calculators/mace_ml.py` - from .base, lazy from .mace_loader
- `mace_gaussian/calculators/xtb.py` - from .base
- `mace_gaussian/gaussian/__init__.py` - from gaussian.* → from .*
- `mace_gaussian/gaussian/io.py` - from ..utils.units
- `mace_gaussian/gaussian/parser.py` - from ..utils.exceptions
- `mace_gaussian/gaussian/runner.py` - from .zmq_server, from ..utils.exceptions
- `mace_gaussian/gaussian/fchk.py` - from ..utils.units
- `mace_gaussian/utils/__init__.py` - from utils.* → from .*
- `mace_gaussian/utils/validation.py` - from .exceptions

## Decisions Made
- `from .workflow import run_pipeline` left commented out in `mace_gaussian/__init__.py` — workflow.py hasn't moved into the package yet (Plan 02's job); activating it now would cause ImportError
- ruff auto-fix (`--select I,F401 --fix`) applied after manual import conversion to sort 2 import blocks in `espaloma.py` and `runner.py`
- `mace_ml.py` lazy import inside `_check_availability` (inside try/except) converted from `from calculators.mace_loader` to `from .mace_loader` — relative import semantics are identical in this context

## Deviations from Plan

None - plan executed exactly as written. The ruff import sorting was expected (the plan said to run ruff after changes); auto-fix was used to resolve the 2 flagged blocks rather than manually editing.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `mace_gaussian/` package namespace is established and all three subpackages import correctly
- Plan 02 can now move `workflow.py` into `mace_gaussian/` and activate `from .workflow import run_pipeline`
- Plan 03 will update all root-level callers (cli.py, gm_helper.py, etc.) from old absolute imports to `mace_gaussian.*` imports

---
*Phase: 08-package-structure-reorganization*
*Completed: 2026-02-24*
