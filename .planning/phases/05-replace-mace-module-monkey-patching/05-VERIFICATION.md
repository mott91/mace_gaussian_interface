---
phase: 05-replace-mace-module-monkey-patching
verified: 2026-02-19T17:05:00Z
status: human_needed
score: 3/4 must-haves verified
re_verification: false
human_verification:
  - test: "Run full pipeline on water molecule end-to-end"
    expected: "Pipeline completes without calling cleanup_mace_modules() at any point; MACE energy and MACE dipole calculations run sequentially in same process"
    why_human: "Requires Gaussian 16 and GPU. Automated checks confirm no cleanup_mace_modules() calls exist in source, but E2E execution cannot be verified without the full environment."
---

# Phase 5: Replace MACE Module Monkey-Patching Verification Report

**Phase Goal:** Replace sys.modules manipulation with safe lazy import isolation pattern
**Verified:** 2026-02-19T17:05:00Z
**Status:** human_needed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | MACE standard and MACE dipole can be loaded independently without sys.modules cleanup | VERIFIED | `mace_loader.py` uses `pickle_module` remapping only; test `test_no_sys_modules_mutation` confirms no `sys.modules["mace` in source; test `test_no_cleanup_mace_modules_reference` confirms no cleanup calls |
| 2 | calculators/mace_loader.py provides safe loading mechanism tested in isolation | VERIFIED | File exists (209 lines), substantive implementation with `_DipoleModelUnpickler`, `_DipolePickleModule`, `load_dipole_calculator`, `MACEDipoleCalculator`; all 11 tests pass |
| 3 | CUDA device placement is handled correctly across MACE variants | VERIFIED | `load_dipole_calculator(device=...)` accepts explicit device; `MACEDipoleCalculator.__init__` stores and passes `device`; defaults to `"cuda"` matching project convention |
| 4 | Full pipeline runs on water without cleanup_mace_modules() calls | NEEDS HUMAN | Source verification passes (zero `cleanup_mace_modules` hits across all source files); E2E run requires Gaussian 16 + GPU |

**Score:** 3/4 truths verified (4th needs human confirmation)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `calculators/mace_loader.py` | Safe MACE dipole loading via pickle_module remapping + MACEDipoleCalculator wrapper | VERIFIED | 209 lines; contains `_DipoleModelUnpickler`, `_DipolePickleModule`, `load_dipole_calculator`, `MACEDipoleCalculator`; no stubs |
| `tests/test_mace_loader.py` | Unit tests for safe loading mechanism | VERIFIED | 222 lines; 11 tests across 4 test classes; all pass |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `calculators/mace_ml.py` | `calculators/mace_loader.py` | `from calculators.mace_loader import MACEDipoleCalculator` | WIRED | Line 30 of mace_ml.py imports from new location; import verified live |
| `calculators/mace_loader.py` | `mace_dipole_core.modules.models` | `_DipoleModelUnpickler.find_class` remapping | WIRED | Line 55: `import mace_dipole_core.modules.models as dipole_models` inside `find_class`; tested by `test_remaps_mace_modules_models` |
| `tests/test_calculators.py` | `mace_dipole_core` hierarchy | pre-mock list | WIRED | `_heavy_deps` contains all 5 `mace_dipole_core.*` entries |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| STRUCT-03 | 05-01-PLAN.md, 05-02-PLAN.md | MACE module monkey-patching replaced with safe loading pattern | SATISFIED | `mace_calculators.py` deleted; `calculators/mace_loader.py` implements pickle_module remapping; zero `cleanup_mace_modules` references in source; REQUIREMENTS.md marks it `[x]` Complete |

No orphaned requirements — only STRUCT-03 maps to Phase 5 in the traceability table.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | - | - | - | - |

No anti-patterns detected. Searched `calculators/mace_loader.py` and `calculators/mace_ml.py` for TODO/FIXME, placeholder returns, and console.log-only handlers. All implementations are substantive.

### Human Verification Required

#### 1. End-to-End Pipeline Run on Water

**Test:** Activate venv, run `python cli.py run water.xyz --energy-calculators mace_mp --dipole-calculators mace_ml`
**Expected:** Pipeline completes without error; no `cleanup_mace_modules()` reference appears in any traceback or log output; both MACE energy optimization (standard mace-torch) and MACE dipole calculation (mace_dipole_core) run in the same process without interference
**Why human:** Requires Gaussian 16 (`g16`) in PATH and a CUDA GPU. The automated checks confirm zero `cleanup_mace_modules` calls in source, but the actual runtime isolation between standard MACE and dipole MACE can only be confirmed by executing the pipeline.

### Gaps Summary

No blocking gaps. All source-verifiable aspects pass:

- `mace_calculators.py` is confirmed deleted (git rm'd)
- Zero references to `cleanup_mace_modules`, `fake_module_from_real`, `load_dipole_mace_calculator`, or `load_standard_mace` in any `.py` or `.toml` file outside `.planning/`
- Zero references to `mace_calculators` module in any `.py`, `.toml`, or `.md` file outside `.planning/`
- `calculators/mace_loader.py` is substantive: pickle_module class remapping fully implemented, scoped torch.load patching with finally-block restoration, MACEDipoleCalculator wrapper with lazy init
- `calculators/mace_ml.py` imports `MACEDipoleCalculator` from `calculators.mace_loader` (not `mace_calculators`)
- `CLAUDE.md` updated: monkey-patching gotcha replaced with pickle_module description
- `docs/ARCHITECTURE.md` updated: Dual MACE Package Architecture paragraph references `calculators/mace_loader.py`
- `pyproject.toml`: `known-first-party` and `coverage.run.source` both reference `calculators` (not `mace_calculators`)
- `tests/test_calculators.py`: `_heavy_deps` mocks `mace_dipole_core.*` hierarchy (not `mace_calculators`)
- All 11 new tests in `tests/test_mace_loader.py` pass
- 119 of 124 collected tests pass; 2 failures are pre-existing environment issues unrelated to Phase 5 (scipy `ObjSense` conflict on system Python, and a `test_optimization_results_contain_version_info` failure from Phase 2 REPR-01/02 work)

The only unverifiable item is the full pipeline E2E run which requires Gaussian 16 and GPU hardware.

---

_Verified: 2026-02-19T17:05:00Z_
_Verifier: Claude (gsd-verifier)_
