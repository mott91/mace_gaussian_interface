# Phase 7: Extract Workflow Orchestrator - Research

**Researched:** 2026-02-22
**Domain:** Python module extraction / workflow orchestration pattern
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**workflow.py interface shape**
- Expose separate per-stage functions plus a top-level coordinator:
  - `run_geometry_optimization()` — stage 1, already exists in gm_main.py
  - `run_dft_baselines()` — stage 2, wraps `dft_baseline.run_all_dft_baselines()`
  - `run_ml_calculations()` — stage 3, loops over energy/dipole combinations
  - `run_pipeline()` — top-level coordinator, calls the three stages in order
- CLI calls `run_pipeline()` (or individual stage functions) directly — no subprocess
- `workflow.py` is a pure Python module, not a CLI script

**DFT / ML independence (cluster seam)**
- DFT baseline and ML calculations are independent stages — both depend only on the optimized geometry output, not on each other
- Each stage reads inputs from and writes results to disk (via ResultsManager) so stages can be run independently
- This creates a natural split point for future cluster execution: run DFT on cluster → copy results → continue from `run_ml_calculations()` or `run_pipeline()` with `include_dft_baselines=False`
- Full cluster automation (SLURM job submission, SSH, file sync) is deferred — design the seam now, implement cluster mechanics later

**gm_main.py removal**
- Hard delete on day one — no deprecation shim, no redirect, clean break
- Low-level helpers in gm_main.py (`update_molecule_geometry`, `calculate_energy_and_forces`, `calculate_hessian`, `calculate_dipole_properties`, `run_next_calculation`, `geometry_optimisation`, `calculator`, `analyze_molecular_charges`, `setup_output_directory`) either move to `workflow.py` or are inlined — researcher to determine best home
- `print_diagnostics()` stays with CLI or moves to `diagnostics.py` if it doesn't already exist there

**Stage control / skip flags**
- `run_pipeline()` accepts `include_dft_baselines: bool` (already in `run_workflow()`) to skip DFT
- Calculator selection (`energy_calculators`, `dipole_calculators`) remain as parameters
- `force_optimization: bool` remains to bypass cached geometry
- CLI flags map 1:1 to `run_pipeline()` parameters — no logic in CLI layer

### Claude's Discretion
- Where exactly each low-level helper ends up (workflow.py vs. another module)
- Whether `run_geometry_optimization` and `run_frequency_calculation` remain in workflow.py or become internal helpers
- Exact module-level docstrings and logging output format

### Deferred Ideas (OUT OF SCOPE)
- Cluster execution: Automatic SLURM job submission for DFT baseline with SSH file sync — future phase or milestone
- True parallelism: Running DFT and ML simultaneously in the same invocation — deferred pending cluster design
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| STRUCT-06 | Workflow orchestration extracted into workflow.py as thin coordinator | Analysis of gm_main.py `run_workflow()` function (lines 595-787) shows it is the direct 1:1 source for `run_pipeline()`. The three PHASE comment blocks map directly to the three stage functions. All dependencies (dft_baseline, ResultsManager, calculator, run_geometry_optimization, run_frequency_calculation) are already imported or local. |
</phase_requirements>

## Summary

This phase is a structural extraction with no behavioral changes. The entire job is to move `run_workflow()` from `gm_main.py` into a new `workflow.py` module, rename it `run_pipeline()`, split its three internal phases into named stage functions, update the CLI to import from `workflow` instead of `gm_main`, and hard-delete `gm_main.py`.

The codebase is already well-prepared for this. Phase 6 completed all gaussian/* extractions. All the dependencies `run_workflow()` needs (`dft_baseline`, `ResultsManager`, `gaussian.*`, `calculators`, `utils.*`) already exist as clean, standalone modules. The `gm_main.py` file is now a thin shell: its top-level imports all point to external packages, and `run_workflow()` itself already calls `dft_baseline.run_all_dft_baselines()` via an internal lazy import. The function is almost verbatim-portable.

The main decisions for the planner are: (1) which low-level helpers from gm_main.py belong in workflow.py vs. stay elsewhere, and (2) how to update `cli.py` to call workflow functions. The diagnostics question is already resolved: `diagnostics.py` already exists as a module but contains only Python environment and espaloma-specific checks — `print_diagnostics()` from gm_main.py is a calculator availability check that belongs in `cli.py` directly or as a thin wrapper there.

**Primary recommendation:** Create `workflow.py` by extracting and reorganizing content from `gm_main.py`. The stage functions map 1:1 to the PHASE 1/2/3 comment blocks already in `run_workflow()`. Update `cli.py` import from `gm_main` to `workflow`. Delete `gm_main.py`.

## Standard Stack

### Core (already in project — no new dependencies)
| Module | Version | Purpose | Why Standard |
|--------|---------|---------|--------------|
| `dft_baseline` | project | DFT baseline stage — `run_all_dft_baselines()` | Already exists, called by current `run_workflow()` |
| `gaussian.*` | project | Gaussian I/O, parsing, ZMQ, runner | Phase 6 complete — all callers already use this |
| `calculators` | project | Energy/dipole calculator factory and implementations | Phase 4 complete |
| `utils.results.ResultsManager` | project | Disk-based results persistence | Used by all stages |
| `utils.validation` | project | `detect_device()`, `validate_all_prerequisites()` | Already called in `run_workflow()` and CLI |
| `ase` | project | ASE Atoms, `read()`, LBFGS optimizer | Required for geometry optimization |
| `click` | project | CLI framework in cli.py | No changes to CLI library |

### No New Dependencies
This phase introduces zero new third-party libraries. It is a pure code reorganization within the existing project.

## Architecture Patterns

### Recommended Module Layout After Phase 7

```
mace_gaussian/
├── workflow.py          # NEW: thin orchestrator (was gm_main.py run_workflow)
├── cli.py               # UPDATED: imports from workflow, not gm_main
├── dft_baseline.py      # UNCHANGED: run_all_dft_baselines() stays here
├── gaussian/            # UNCHANGED: Phase 6 output
├── calculators/         # UNCHANGED: Phase 4 output
├── utils/               # UNCHANGED
└── gm_main.py           # DELETED: hard delete, no shim
```

### Pattern 1: Stage Function Decomposition

The existing `run_workflow()` body has three clearly-delineated phases. Each becomes a named stage function in workflow.py.

**Stage 1 — `run_geometry_optimization()` (already exists in gm_main.py, lines 338-395)**

This function already exists with the right signature. Move it verbatim to workflow.py. It depends only on `geometry_optimisation()` (the lower-level LBFGS runner, lines 226-239) and `ResultsManager`. Both can follow it into workflow.py.

```python
# workflow.py — Stage 1
def run_geometry_optimization(
    atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    calculator_name: str = "mace_omol",
) -> "ase.Atoms":
    """Run LBFGS geometry optimization and save results to disk."""
    ...
```

**Stage 2 — `run_dft_baselines()` (thin wrapper)**

Currently inside `run_workflow()` as an inline if-block calling `run_all_dft_baselines()`. Becomes a named function that just wraps the dft_baseline module call. Keep the signature aligned with `run_pipeline()` parameters.

```python
# workflow.py — Stage 2
def run_dft_baselines(
    optimized_atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    charge: int = 0,
    multiplicity: int = 1,
) -> dict[str, bool]:
    """Run DFT baseline calculations. Wraps dft_baseline.run_all_dft_baselines()."""
    from dft_baseline import run_all_dft_baselines
    return run_all_dft_baselines(
        optimized_atoms, molecule_name, results_mgr, charge, multiplicity, skip_if_exists=True
    )
```

**Stage 3 — `run_ml_calculations()` (extracted loop)**

Currently the nested for-loop at lines 741-753 in `run_workflow()`. Becomes a named function that loops over energy/dipole combinations and calls `run_frequency_calculation()`.

```python
# workflow.py — Stage 3
def run_ml_calculations(
    optimized_atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    energy_calculators: list[str],
    dipole_calculators: list[str],
    charge: int = 0,
    multiplicity: int = 1,
) -> list[dict]:
    """Run all ML energy+dipole calculator combinations via Gaussian."""
    ...
```

**Top-level — `run_pipeline()` (rename of `run_workflow()`)**

Calls the three stage functions in sequence. Accepts all the same parameters as current `run_workflow()`. CLI's call site changes only the import and function name.

```python
# workflow.py — Coordinator
def run_pipeline(
    input_file: str,
    optimization_calculator: str = "mace_omol",
    energy_calculators: list | None = None,
    dipole_calculators: list | None = None,
    force_optimization: bool = False,
    include_dft_baselines: bool = True,
    base_output_dir: str = "comparison_results",
) -> dict:
    """Coordinate all pipeline stages. Returns summary dict."""
    ...
```

### Pattern 2: Helper Function Disposition

gm_main.py contains several low-level helpers. Here is where each belongs:

| Function | Lines | Disposition | Rationale |
|----------|-------|-------------|-----------|
| `update_molecule_geometry()` | 70-83 | Move to workflow.py | Used only by `run_next_calculation()` which stays in workflow.py |
| `calculate_energy_and_forces()` | 85-98 | Move to workflow.py | Used only by `run_next_calculation()` |
| `calculate_hessian()` | 101-118 | Move to workflow.py | Used only by `run_next_calculation()` |
| `calculate_dipole_properties()` | 121-168 | Move to workflow.py | Used only by `run_next_calculation()` |
| `run_next_calculation()` | 173-223 | Move to workflow.py | ZMQ callback handler, called from `run_frequency_calculation()` |
| `geometry_optimisation()` | 226-239 | Move to workflow.py | Lower-level LBFGS runner, called by `run_geometry_optimization()` |
| `calculator()` | 242-266 | Move to workflow.py | Calculator factory for energy calcs, called in `run_frequency_calculation()` |
| `analyze_molecular_charges()` | 271-307 | Move to workflow.py OR inline | Only called within run_next_calculation flow; can stay in workflow.py or be inlined |
| `setup_output_directory()` | 310-335 | DROP | Not called anywhere in run_workflow(). Unused in the new design. |
| `run_geometry_optimization()` | 338-395 | Move to workflow.py | Stage 1 function, already has right signature |
| `run_frequency_calculation()` | 398-563 | Move to workflow.py | Called by run_ml_calculations() |
| `print_diagnostics()` | 578-587 | Move to cli.py inline | CLI-only function, does not belong in workflow logic |
| `run_workflow()` | 595-787 | Becomes `run_pipeline()` in workflow.py | Direct source |

**Key insight:** `setup_output_directory()` is defined but never called in the current run_workflow() code path. It can be safely dropped.

### Pattern 3: CLI Update

`cli.py` currently imports from `gm_main`:

```python
# BEFORE (cli.py lines 78-79, 366-367)
from gm_main import run_workflow
from gm_main import print_diagnostics
```

After phase 7:

```python
# AFTER (cli.py)
from workflow import run_pipeline
# print_diagnostics inlined into cli.py diagnose command, or imported from diagnostics.py
```

The `run` command call site changes from `run_workflow(...)` to `run_pipeline(...)` with identical arguments. The `diagnose` command absorbs `print_diagnostics()` logic directly since it's only 9 lines.

### Pattern 4: Module-Level Imports in workflow.py

Follow gm_main.py's import structure at the top of workflow.py. The lazy import of `dft_baseline` (currently inside `run_workflow()` body) should remain lazy to avoid heavy DGL/espaloma imports at module load time:

```python
# workflow.py top-level imports (always needed)
from __future__ import annotations

import json
import logging
import os
import shutil
import time
from pathlib import Path
from typing import Optional

import numpy as np
from ase.io import read
from ase.optimize import LBFGS

from calculators import dipole_factory
from gaussian.fchk import convert_chk_to_fchk
from gaussian.io import ase_to_gjf, parse_gaussian_input, write_gaussian_output
from gaussian.parser import parse_gaussian_log
from gaussian.runner import DEFAULT_TIMEOUT_SECONDS, run_gaussian_with_zmq
from utils.results import ResultsManager
from utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV

# Lazy import inside run_dft_baselines() to avoid heavy deps at import time:
# from dft_baseline import run_all_dft_baselines
```

### Anti-Patterns to Avoid

- **Leaving a stub in gm_main.py**: The decision is hard delete. No `from workflow import *`, no redirects. Anything referencing `gm_main` after this phase is a bug to fix.
- **Adding logic to CLI**: CLI flags must map 1:1 to `run_pipeline()` parameters. No branching, validation, or defaults in the CLI layer beyond what Click already provides.
- **Changing function signatures**: `run_pipeline()` must accept identical parameters to current `run_workflow()`. This is a rename + reorganization, not a redesign.
- **Moving dft_baseline.py contents**: `dft_baseline.run_all_dft_baselines()` stays in `dft_baseline.py`. Only the orchestration call moves to workflow.py.
- **Eager import of dft_baseline at module level**: This triggers heavy DGL/espaloma imports. Keep `from dft_baseline import run_all_dft_baselines` inside the stage function body.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Geometry optimizer | Custom optimization loop | `ase.optimize.LBFGS` (already used) | ASE handles convergence, line search, step limits correctly |
| Results persistence | Ad-hoc JSON writing | `ResultsManager` (already exists) | Handles directory structure, file naming, metadata |
| Calculator loading | Direct model instantiation in workflow | `calculators.dipole_factory` (already exists) | Handles lazy imports, availability checks |
| Gaussian subprocess | Custom process management in workflow | `gaussian.runner.run_gaussian_with_zmq` (already exists) | Handles ZMQ loop, SIGKILL timeout, error types |

**Key insight:** Every complex operation is already encapsulated in a dedicated module from previous phases. workflow.py is genuinely thin — it sequences calls to these modules, nothing more.

## Common Pitfalls

### Pitfall 1: Importing gm_main transitively
**What goes wrong:** Some file other than cli.py still imports from gm_main. After deletion, runtime ImportError appears only on that code path.
**Why it happens:** gm_main was the original monolith — it's plausible other scripts referenced it directly.
**How to avoid:** Run `grep -r "from gm_main\|import gm_main" .` before deleting. Verify zero hits besides cli.py.
**Warning signs:** Tests or scripts that hit gm_main import paths during CI.

### Pitfall 2: Breaking the diagnose command
**What goes wrong:** `cli.py diagnose` calls `from gm_main import print_diagnostics`. After deletion this ImportError is immediate on `python cli.py diagnose`.
**Why it happens:** `print_diagnostics()` is in gm_main.py, not in diagnostics.py (which exists but contains only Python environment and espaloma diagnostics, not calculator availability).
**How to avoid:** Either inline `print_diagnostics()` in the CLI `diagnose` command, or move it to `diagnostics.py`. The function is only 9 lines; inlining is simplest.
**Warning signs:** `python cli.py diagnose` fails with ImportError during smoke test.

### Pitfall 3: Losing the lazy dft_baseline import
**What goes wrong:** Moving `from dft_baseline import run_all_dft_baselines` to the top of workflow.py causes heavy DGL/espaloma imports at module load time, breaking tests that pre-mock sys.modules.
**Why it happens:** `dft_baseline.py` imports from `gaussian.parser` and `utils.results` which are safe, but espaloma calculators may be transitively imported if the import order changes. Current code puts this import inside the function body.
**How to avoid:** Keep `from dft_baseline import run_all_dft_baselines` inside `run_dft_baselines()` function body, matching the current pattern in `run_workflow()` (line 629).
**Warning signs:** `ImportError: No module named 'dgl'` appearing in tests that don't have GPU/espaloma mocks.

### Pitfall 4: Forgetting charge/multiplicity extraction after optimization
**What goes wrong:** `run_pipeline()` calls `run_geometry_optimization()`, then needs charge/spin for subsequent stages. If the extraction logic (lines 715-717 in gm_main.py) is lost during extraction, stages 2 and 3 get wrong charge/multiplicity.
**Why it happens:** The `int(optimized_mol.info.get("charge", 0))` lines are in `run_workflow()` body, not in a helper function. Easy to overlook when decomposing.
**How to avoid:** Trace the charge/multiplicity flow explicitly when writing `run_pipeline()`. After `run_geometry_optimization()` returns `optimized_atoms`, extract `charge = int(optimized_mol.info.get("charge", 0))` before passing to stages 2 and 3.
**Warning signs:** DFT or ML calculations run with charge=0 even for charged molecules.

### Pitfall 5: setup_output_directory() carried over unnecessarily
**What goes wrong:** `setup_output_directory()` is defined in gm_main.py but is NOT called anywhere in `run_workflow()`. If it is carried over to workflow.py, it adds dead code.
**Why it happens:** It exists in gm_main.py alongside the actual pipeline functions and looks related.
**How to avoid:** Verify it has no callers (`grep -r "setup_output_directory" .`). Drop it entirely.
**Warning signs:** If callers are found, they are likely legacy scripts that can be updated to use ResultsManager directly.

### Pitfall 6: rdkit/warnings suppression in gm_main.py header
**What goes wrong:** gm_main.py has top-of-file warning suppression (`warnings.filterwarnings`, `RDLogger.DisableLog`, `os.environ["PYTHONWARNINGS"]`). If not reproduced in workflow.py, rdkit warnings may appear during calculation runs.
**Why it happens:** The warning suppression is at module level in gm_main.py, not in any function.
**How to avoid:** Carry the warning suppression block over to workflow.py module level. Or move it to cli.py (more appropriate since it affects process-level output).
**Recommendation:** Move to cli.py since suppression is a presentation concern, not a workflow concern.

## Code Examples

Verified patterns from codebase analysis:

### run_pipeline() Skeleton (direct extraction from run_workflow())
```python
# Source: gm_main.py lines 595-787
def run_pipeline(
    input_file: str,
    optimization_calculator: str = "mace_omol",
    energy_calculators: list | None = None,
    dipole_calculators: list | None = None,
    force_optimization: bool = False,
    include_dft_baselines: bool = True,
    base_output_dir: str = "comparison_results",
) -> dict:
    """
    Run the complete MACE-Gaussian comparison pipeline.

    Sequences three independent stages:
      1. Geometry optimization (MACE ML)
      2. DFT baseline calculations (Gaussian, pure DFT)
      3. ML frequency calculations (Gaussian + MACE ZMQ interface)

    Stages 2 and 3 both depend only on the optimized geometry from stage 1
    and write results to disk independently via ResultsManager. This creates
    a natural seam for future cluster execution: run stage 2 on cluster,
    copy results, resume from stage 3 with include_dft_baselines=False.
    """
    from utils.validation import detect_device
    ...
    results_mgr = ResultsManager(base_output_dir=base_output_dir)

    # Stage 1
    optimized_atoms = run_geometry_optimization(...)

    # Stage 2 (independent of stage 3)
    charge = int(optimized_atoms.info.get("charge", 0))
    multiplicity = int(optimized_atoms.info.get("spin", 1))
    dft_results = {}
    if include_dft_baselines:
        dft_results = run_dft_baselines(optimized_atoms, molecule_name, results_mgr, charge, multiplicity)

    # Stage 3 (independent of stage 2)
    ml_results = run_ml_calculations(
        optimized_atoms, molecule_name, results_mgr,
        energy_calculators, dipole_calculators, charge, multiplicity,
    )

    return {"dft_baselines": dft_results, "ml_calculations": ml_results, "molecule_name": molecule_name}
```

### CLI Update Pattern
```python
# Source: cli.py lines 78-82, 366-370 — BEFORE
from gm_main import run_workflow
results = run_workflow(input_file=..., ...)

from gm_main import print_diagnostics
print_diagnostics()

# AFTER
from workflow import run_pipeline
results = run_pipeline(input_file=..., ...)

# print_diagnostics inlined in diagnose command:
from calculators import dipole_factory
for name, available in dipole_factory.list_available().items():
    status = "OK" if available else "UNAVAILABLE"
    click.echo(f"  {status}: {name}")
```

### Grep Commands for Safety Checks
```bash
# Before deleting gm_main.py — must return zero results
grep -r "from gm_main\|import gm_main" /home/mot/mace_gaussian/ --include="*.py"
# setup_output_directory callers — expect zero
grep -r "setup_output_directory" /home/mot/mace_gaussian/ --include="*.py"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| All pipeline logic in gm_main.py | Extracted modules (gaussian/, calculators/, utils/) + orchestrator in workflow.py | Phases 3-7 | gm_main.py becomes deletable |
| CLI imports from gm_main | CLI imports from workflow | Phase 7 | Single authoritative pipeline entry point |
| run_workflow() as monolith | run_pipeline() + named stage functions | Phase 7 | Stages callable independently for cluster split |

**Deprecated/outdated after this phase:**
- `gm_main.py`: Deleted. Any reference to it is a bug.
- `from gm_main import run_workflow` in cli.py: Replaced with `from workflow import run_pipeline`.
- `from gm_main import print_diagnostics` in cli.py: Replaced with inlined logic.

## Open Questions

1. **charge_analysis / ChargeAnalyzer import at top of gm_main.py**
   - What we know: `gm_main.py` has `from charge_analysis import ChargeAnalyzer` at module level (line 29) AND again inside a try/except (line 53). The `analyze_molecular_charges()` helper uses it but it is NOT called in `run_workflow()` path.
   - What's unclear: Whether charge analysis should be carried into workflow.py or dropped entirely for this phase.
   - Recommendation: Drop `analyze_molecular_charges()` and the `charge_analysis` import from workflow.py entirely. It is not called in the pipeline. If needed later, it can be added back. Avoids pulling in an optional dependency.

2. **Duplicate `from charge_analysis import ChargeAnalyzer` in gm_main.py**
   - What we know: Line 29 imports it unconditionally, line 53-59 imports it again inside try/except setting `CHARGE_ANALYSIS_AVAILABLE`. This is a bug in gm_main.py — the unconditional import at line 29 will fail if `charge_analysis` is missing, making the try/except irrelevant.
   - What's unclear: Does `charge_analysis` module exist in the project?
   - Recommendation: When writing workflow.py, omit the unconditional import entirely. The charge analysis path is dead code in the pipeline.

3. **mace_calculators.py/__pycache__ entries in git status**
   - What we know: Git status shows `D __pycache__/mace_calculators.cpython-310.pyc` and `.cpython-312.pyc` — deleted bytecode, not source.
   - What's unclear: Whether `mace_calculators.py` itself still exists as source.
   - Recommendation: Not relevant to this phase. The pycache deletions are from Phase 5's mace_calculators.py deletion.

## Sources

### Primary (HIGH confidence)
- Direct code analysis: `/home/mot/mace_gaussian/gm_main.py` — full read, line-by-line function mapping
- Direct code analysis: `/home/mot/mace_gaussian/cli.py` — import sites identified (lines 78-82, 366-370)
- Direct code analysis: `/home/mot/mace_gaussian/dft_baseline.py` — `run_all_dft_baselines()` signature confirmed
- Direct code analysis: `/home/mot/mace_gaussian/diagnostics.py` — confirmed it does NOT contain `print_diagnostics()` (calculator availability checker)
- Direct code analysis: `/home/mot/mace_gaussian/gaussian/__init__.py` — Phase 6 exports confirmed
- `.planning/phases/07-extract-workflow-orchestrator/07-CONTEXT.md` — locked decisions

### Secondary (MEDIUM confidence)
- N/A: This is an internal refactoring phase. No external library research was required.

### Tertiary (LOW confidence)
- N/A

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new dependencies, all existing modules confirmed by direct code read
- Architecture: HIGH — function-by-function mapping from actual source code
- Pitfalls: HIGH — derived from actual code analysis (e.g., duplicate imports, dead code), not speculation

**Research date:** 2026-02-22
**Valid until:** 2026-03-22 (stable internal refactoring; no external library changes affect this)
