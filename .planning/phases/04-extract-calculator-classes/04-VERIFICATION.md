---
phase: 04-extract-calculator-classes
verified: 2026-02-17T15:30:00Z
status: passed
score: 11/11 must-haves verified
re_verification: false
---

# Phase 4: Extract Calculator Classes Verification Report

**Phase Goal:** Move calculator implementations from gm_main.py into dedicated calculators/ package
**Verified:** 2026-02-17T15:30:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #  | Truth                                                                                                          | Status     | Evidence                                                                 |
|----|----------------------------------------------------------------------------------------------------------------|------------|--------------------------------------------------------------------------|
| 1  | DipoleCalculatorBase ABC exists in calculators/base.py with calculate_dipole and calculate_dipole_derivatives | VERIFIED   | Both methods present and non-abstract base implements calculate_dipole_derivatives as template method |
| 2  | All 3 calculator implementations exist as separate modules in calculators/ package                            | VERIFIED   | calculators/espaloma.py, calculators/mace_ml.py, calculators/xtb.py all exist and are substantive |
| 3  | DipoleCalculatorFactory and dipole_factory singleton are importable from calculators package                  | VERIFIED   | calculators/factory.py defines both; calculators/__init__.py re-exports both |
| 4  | gm_main.py imports dipole_factory from calculators instead of defining it inline                              | VERIFIED   | Line 33: `from calculators import dipole_factory`; no class definitions remain |
| 5  | Calculator interface compliance is tested: all 3 are subclasses of DipoleCalculatorBase                       | VERIFIED   | test_calculators.py TestCalculatorInterface passes all 3 subclass assertions |
| 6  | All implementations have the required methods (calculate_dipole, calculate_dipole_derivatives, _check_availability) | VERIFIED | test_has_required_methods parametrized over all 3 classes — all 18 tests pass |
| 7  | Factory returns correct calculator type for each known name                                                    | VERIFIED   | test_get_known_calculator passes; factory dict keyed by .name attribute |
| 8  | Factory raises ValueError for unknown calculator name                                                         | VERIFIED   | test_get_unknown_calculator_raises asserts ValueError with "Unknown" |
| 9  | Factory auto-selection returns first available calculator in preferred order                                   | VERIFIED   | test_auto_selects_first_available and test_auto_selects_mace_ml_when_available pass |
| 10 | Factory list_available returns dict mapping names to availability booleans                                    | VERIFIED   | test_list_available asserts dict with correct True/False values per calculator |
| 11 | Existing tests still pass after extraction                                                                    | VERIFIED   | Full suite: 117 passed, 2 skipped, 1 xfailed (pre-existing) — 0 new failures |

**Score:** 11/11 truths verified

---

### Required Artifacts

#### Plan 04-01 Artifacts

| Artifact                   | Expected                                          | Status   | Details                                                             |
|----------------------------|---------------------------------------------------|----------|---------------------------------------------------------------------|
| `calculators/__init__.py`  | Package re-exports with `__all__`                 | VERIFIED | 17 lines; imports all 5 names + dipole_factory; has `__all__` list  |
| `calculators/base.py`      | DipoleCalculatorBase ABC                          | VERIFIED | 68 lines; ABC with @abstractmethod on _check_availability and calculate_dipole; calculate_dipole_derivatives implemented |
| `calculators/espaloma.py`  | EspalomaDipoleCalculator                          | VERIFIED | 77 lines; inherits base, lazy espaloma/rdkit imports inside methods  |
| `calculators/mace_ml.py`   | MACEMLDipoleCalculator + DEFAULT_MACE_DIPOLE_MODEL| VERIFIED | 47 lines; constant at module level, lazy MACEDipoleCalculator import inside _check_availability |
| `calculators/xtb.py`       | XTBDipoleCalculator                               | VERIFIED | 41 lines; inherits base, lazy xtb import inside methods              |
| `calculators/factory.py`   | DipoleCalculatorFactory + dipole_factory singleton| VERIFIED | 63 lines; class + `dipole_factory = DipoleCalculatorFactory()` at end |

#### Plan 04-02 Artifacts

| Artifact                    | Expected                              | Status   | Details                                             |
|-----------------------------|---------------------------------------|----------|-----------------------------------------------------|
| `tests/test_calculators.py` | Unit tests for calculator hierarchy   | VERIFIED | 214 lines (min 80 required); 18 tests across 3 classes; all pass |

---

### Key Link Verification

| From                       | To                                              | Via                               | Status   | Details                                                        |
|----------------------------|-------------------------------------------------|-----------------------------------|----------|----------------------------------------------------------------|
| calculators/factory.py     | calculators/espaloma.py, mace_ml.py, xtb.py   | import + instantiation in _register_calculators | VERIFIED | Lines 8-10: direct imports; instantiated in _register_calculators() |
| calculators/mace_ml.py     | mace_calculators.py                            | lazy import of MACEDipoleCalculator | VERIFIED | Line 30: `from mace_calculators import MACEDipoleCalculator` inside _check_availability — NOT at module level |
| gm_main.py                 | calculators/__init__.py                        | `from calculators import dipole_factory` | VERIFIED | Line 33 of gm_main.py; used at lines 392 and 808 |
| tests/test_calculators.py  | calculators/base.py                            | import DipoleCalculatorBase       | VERIFIED | Line 41: `from calculators.base import DipoleCalculatorBase` |
| tests/test_calculators.py  | calculators/factory.py                         | import DipoleCalculatorFactory    | VERIFIED | Line 113 inside _make_factory: `from calculators.factory import DipoleCalculatorFactory` |

---

### Requirements Coverage

| Requirement | Source Plan | Description                                                               | Status    | Evidence                                                                      |
|-------------|-------------|---------------------------------------------------------------------------|-----------|-------------------------------------------------------------------------------|
| STRUCT-02   | 04-01, 04-02 | DipoleCalculator classes extracted from gm_main.py into calculators/ package | SATISFIED | calculators/ package with 6 files exists; gm_main.py reduced by 246 lines (1285→1039); no class definitions remain in gm_main.py |

No orphaned requirements — only STRUCT-02 is mapped to Phase 4 in REQUIREMENTS.md.

---

### Anti-Patterns Found

| File                   | Line | Pattern            | Severity | Impact                                                              |
|------------------------|------|--------------------|----------|---------------------------------------------------------------------|
| calculators/base.py    | 28   | `Optional[np.ndarray]` (ruff UP045: should be `np.ndarray | None`) | Warning  | Style inconsistency; no functional impact; `from __future__ import annotations` makes it safe at runtime |

No blocker anti-patterns found. No TODOs, FIXMEs, placeholder implementations, or empty handlers in any new file.

---

### Human Verification Required

None. All observable truths can be verified programmatically via file content inspection and test execution.

---

### Gaps Summary

No gaps. All must-haves from both plans are fully satisfied.

The one minor finding is a non-blocking ruff UP045 warning in `calculators/base.py` line 28 (`Optional[np.ndarray]` should be `np.ndarray | None`). This is a style-only issue: the file already has `from __future__ import annotations` making the annotation safe at runtime, and it has no functional impact. It was not caught earlier because ruff format passes clean while ruff check catches it.

---

## Commit Verification

All commits documented in SUMMARYs verified to exist in git history:

- `89e9ab3` — feat(04-01): create calculators/ package with extracted dipole calculator classes
- `d47a0e6` — refactor(04-01): remove calculator classes from gm_main.py, import from calculators package
- `5aa277c` — test(04-02): add calculator interface, factory, and config unit tests

---

_Verified: 2026-02-17T15:30:00Z_
_Verifier: Claude (gsd-verifier)_
