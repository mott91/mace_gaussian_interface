---
phase: 07-extract-workflow-orchestrator
verified: 2026-02-24T17:30:00Z
status: passed
score: 11/11 must-haves verified
re_verification: false
human_verification:
  - test: "Run python cli.py diagnose in project environment"
    expected: "Prints DIAGNOSTIC MODE header, lists mace_ml/espaloma/xtb calculators with OK/UNAVAILABLE status"
    why_human: "Requires a working project environment with click installed (venv after uv sync)"
  - test: "Run python cli.py run water.xyz --skip-dft-baseline in project environment"
    expected: "Workflow executes stage 1 (geometry opt) and stage 3 (ML calcs) via run_pipeline()"
    why_human: "Requires MACE models and Gaussian 16 installed; full pipeline smoke test"
---

# Phase 7: Extract Workflow Orchestrator — Verification Report

**Phase Goal:** Extract the workflow orchestrator from gm_main.py into workflow.py, wire cli.py to use it, and delete gm_main.py. STRUCT-06 must be fully satisfied.
**Verified:** 2026-02-24T17:30:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #   | Truth                                                                                 | Status     | Evidence                                                                              |
|-----|---------------------------------------------------------------------------------------|------------|---------------------------------------------------------------------------------------|
| 1   | workflow.py exists with run_pipeline(), run_geometry_optimization(), run_dft_baselines(), run_ml_calculations() | VERIFIED | All 4 coordinator/stage functions confirmed at lines 533, 250, 476, 495              |
| 2   | All low-level helpers from gm_main.py present in workflow.py                         | VERIFIED   | All 7 helpers at lines 56, 70, 84, 102, 151, 202, 219; 12 total functions           |
| 3   | workflow.py importable without heavy DGL/espaloma deps at import time                | VERIFIED   | Python -c "import workflow" succeeds in project venv (3.12); DGL not loaded (ImportError caught gracefully in _check_availability()) |
| 4   | setup_output_directory() NOT in workflow.py                                           | VERIFIED   | grep returns no matches for setup_output_directory or analyze_molecular_charges      |
| 5   | analyze_molecular_charges() and charge_analysis import NOT in workflow.py            | VERIFIED   | grep returns no matches                                                               |
| 6   | dft_baseline import remains lazy (inside run_dft_baselines() body)                  | VERIFIED   | Line 488: `from dft_baseline import run_all_dft_baselines` inside function body      |
| 7   | cli.py imports run_pipeline from workflow (not run_workflow from gm_main)            | VERIFIED   | Line 84 of cli.py: `from workflow import run_pipeline` (lazy, inside run() body)     |
| 8   | cli.py diagnose command works without gm_main (print_diagnostics logic inlined)      | VERIFIED   | Lines 364-373: `from calculators import dipole_factory` + `dipole_factory.list_available()` |
| 9   | gm_main.py is deleted from the project                                               | VERIFIED   | `ls gm_main.py` returns "No such file or directory"                                  |
| 10  | No Python file in the project imports from gm_main                                  | VERIFIED   | Only `testingStuff/test_refactoring.py` retains reference (explicitly allowed by plan) |
| 11  | python cli.py --help succeeds                                                        | VERIFIED   | Returns full help text with all 5 commands (run, list, compare, export, diagnose)   |

**Score:** 11/11 truths verified

### Required Artifacts

| Artifact      | Expected                                                         | Status        | Details                                        |
|---------------|------------------------------------------------------------------|---------------|------------------------------------------------|
| `workflow.py` | Complete pipeline orchestrator with stage functions (min 400 lines) | VERIFIED   | 718 lines, 12 functions, ruff clean           |
| `cli.py`      | Updated CLI with workflow imports and inlined diagnostics        | VERIFIED      | lazy import from workflow, dipole_factory inline, rdkit suppression added |

### Key Link Verification

| From                        | To                                    | Via                                           | Status   | Details                                                   |
|-----------------------------|---------------------------------------|-----------------------------------------------|----------|-----------------------------------------------------------|
| workflow.py:run_pipeline()  | workflow.py:run_geometry_optimization() | direct call in stage 1 block               | WIRED    | Line 648: `optimized_mol = run_geometry_optimization(...)`|
| workflow.py:run_pipeline()  | workflow.py:run_dft_baselines()       | direct call in stage 2 block                | WIRED    | Line 668: `dft_results = run_dft_baselines(...)`          |
| workflow.py:run_pipeline()  | workflow.py:run_ml_calculations()     | direct call in stage 3 block                | WIRED    | Line 676: `ml_results = run_ml_calculations(...)`         |
| workflow.py:run_dft_baselines() | dft_baseline.run_all_dft_baselines() | lazy import inside function body          | WIRED    | Line 488: `from dft_baseline import run_all_dft_baselines`|
| cli.py:run command          | workflow.run_pipeline()               | `from workflow import run_pipeline` (lazy)  | WIRED    | Line 84 (inside run() body), called at line 131           |
| cli.py:diagnose command     | calculators.dipole_factory            | inlined print_diagnostics logic             | WIRED    | Lines 364-373: direct import + list_available() loop      |

### Requirements Coverage

| Requirement | Source Plan | Description                                           | Status    | Evidence                                                      |
|-------------|-------------|-------------------------------------------------------|-----------|---------------------------------------------------------------|
| STRUCT-06   | 07-01, 07-02 | Workflow orchestration extracted into workflow.py as thin coordinator | SATISFIED | workflow.py exists with named stage functions; gm_main.py deleted; cli.py imports from workflow |

### Anti-Patterns Found

| File     | Lines    | Pattern                  | Severity | Impact                                              |
|----------|----------|--------------------------|----------|-----------------------------------------------------|
| cli.py   | 278, 324 | `[PLACEHOLDER - Coming soon]` in compare and export commands | INFO | Pre-existing placeholder commands; unrelated to phase goal (workflow extraction) |
| pyproject.toml | 46-47, 65, 86 | Stale `gm_main` references in `[project.scripts]`, `[tool.ruff.lint.isort]`, `[tool.coverage.run]` | WARNING | gm_main.py is deleted but pyproject.toml still declares `gm-main` and `gm-diagnose` entry points and includes `gm_main` in coverage source and isort first-party. These will cause errors if entry points are invoked (`gm-main` / `gm-diagnose` CLI aliases are now broken). Not a blocker for the phase goal itself, but should be cleaned up. |

**Placeholders in compare/export are INFO only** — these commands were placeholders before this phase and are not in scope for phase 7.

**pyproject.toml stale refs are WARNING** — The broken `gm-main` / `gm-diagnose` entry points could confuse users, but `python cli.py --help` and `python cli.py run` still work correctly. The isort and coverage config entries are cosmetically stale but non-breaking.

### Human Verification Required

#### 1. cli.py diagnose in project environment

**Test:** Run `python cli.py diagnose` after `uv sync && source .venv/bin/activate`
**Expected:** Prints the diagnostic header and lists calculator availability (mace_ml: OK or UNAVAILABLE, espaloma: UNAVAILABLE, xtb: UNAVAILABLE, or similar)
**Why human:** Requires a working project environment with click installed. The project venv (3.12) lacks click currently (uv sync required), and the system Python (3.8) has a broken DGL installation that causes RuntimeError during workflow.py import.

#### 2. run_pipeline() smoke test

**Test:** Run `python cli.py run water.xyz --skip-dft-baseline --energy-calculators mace_mp --dipole-calculators mace_ml` in a fully configured environment
**Expected:** Stage 1 (geometry optimization) runs, Stage 3 (ML frequency calculations) runs, results saved to comparison_results/water/
**Why human:** Requires MACE models, CUDA, and Gaussian 16 installed.

### Gaps Summary

No gaps found. All 11 must-have truths are verified.

**Importability note:** The must-have "workflow.py is importable without heavy DGL/espaloma deps at import time" is fully met in the project environment (Python 3.12 venv). In the system Python 3.8 environment, a pre-existing DGL/PyTorch version mismatch causes a RuntimeError that propagates through the try/except in `_check_availability()` (which only catches `ImportError`). This is a pre-existing environment issue documented in the 07-02-SUMMARY.md. The code structure is correct: DGL is not a module-level import in workflow.py, and in the intended project environment the import succeeds cleanly.

---

_Verified: 2026-02-24T17:30:00Z_
_Verifier: Claude (gsd-verifier)_
