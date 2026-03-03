---
phase: 13-calculator-expansion-acoh-bug-fix
plan: "02"
subsystem: calculators
tags: [mace_anicc, workflow, element-guard, ase, calculator]

requires:
  - phase: 12-test-infrastructure
    provides: pytest test infrastructure and patterns for mock-based testing

provides:
  - mace_anicc branch in workflow.calculator() using correct minimal API
  - _check_mace_anicc_elements() element guard for HCNO-only restriction
  - Guard wired at both run_frequency_calculation and run_pipeline call sites

affects:
  - Phase 14+ (any plan using mace_anicc as energy or optimization calculator)
  - Any benchmark plan combining mace_anicc with non-HCNO molecules (F, S, Cl, etc.)

tech-stack:
  added: []
  patterns:
    - "Element guard pattern: validate domain-specific constraints BEFORE expensive model load"
    - "Guard placed outside try/except so ValueError propagates (not silently caught)"
    - "mace_anicc API differs from other MACE calcs — no model= or default_dtype= args"

key-files:
  created:
    - tests/test_workflow_calculator.py
  modified:
    - mace_gaussian/workflow.py

key-decisions:
  - "mace_anicc uses mace_anicc(device='cuda') only — no model= or default_dtype= (TypeError if passed)"
  - "Element guard placed OUTSIDE try/except in run_frequency_calculation so ValueError propagates to caller"
  - "Guard checks elements at call sites (not inside calculator()) to fail fast before any model download"

patterns-established:
  - "Element guard pattern: if calc_name == 'restricted_calc': _check_elements(atoms) before calculator() call"
  - "TDD: write failing tests first, verify RED, implement minimal code, verify GREEN, lint"

requirements-completed:
  - CALC-02

duration: 3min
completed: 2026-03-03
---

# Phase 13 Plan 02: mace_anicc Calculator Branch Summary

**mace_anicc energy calculator added to workflow.py with HCNO element guard that raises ValueError before model download for unsupported elements**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-03-03T17:32:13Z
- **Completed:** 2026-03-03T17:34:47Z
- **Tasks:** 1 (TDD — RED + GREEN commits)
- **Files modified:** 2

## Accomplishments

- Added `_MACE_ANICC_SUPPORTED_ELEMENTS` frozenset constant at module level
- Added `_check_mace_anicc_elements()` helper raising ValueError with sorted unsupported element list
- Added `mace_anicc` branch to `calculator()` using `mace_anicc(device="cuda")` only (no `model=` or `default_dtype=`)
- Wired element guard at `run_frequency_calculation` call site (outside try/except so it propagates)
- Wired element guard at `run_pipeline` call site before `calculator()` call
- 16 new tests covering guard behavior, call sites, and non-regression of other calculators

## Task Commits

TDD task had two commits (RED then GREEN):

1. **Task 1 RED: Failing tests** - `6000896` (test)
2. **Task 1 GREEN: Implementation** - `c1af2f9` (feat)

## Files Created/Modified

- `tests/test_workflow_calculator.py` — 16 tests covering element guard, mace_anicc branch, and call site wiring (all mock-based, no GPU needed)
- `mace_gaussian/workflow.py` — Added constant, helper function, mace_anicc branch, guards at both call sites

## Decisions Made

- Element guard placed OUTSIDE the `try/except` block in `run_frequency_calculation` so `ValueError` propagates to the caller rather than being silently caught and returning `False`. This is the correct behavior: element mismatch is a programming error, not a recoverable runtime failure.
- `mace_anicc(device="cuda")` only — no `model=` or `default_dtype=` parameters, which differ from `mace_mp`/`mace_off`/`mace_omol` API and would raise `TypeError` if passed.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Element guard placement: moved outside try/except**
- **Found during:** Task 1 GREEN (test run after initial implementation)
- **Issue:** First placement put guard inside the `try:` block, causing `ValueError` to be caught by the `except Exception` handler and swallowed (function returned `False` instead of raising)
- **Fix:** Moved `_check_mace_anicc_elements(atoms)` to before the `try:` block so `ValueError` propagates naturally
- **Files modified:** `mace_gaussian/workflow.py`
- **Verification:** Test `test_run_frequency_calculation_raises_for_non_hcno_with_mace_anicc` passed after fix
- **Committed in:** c1af2f9 (Task 1 GREEN commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - bug fix)
**Impact on plan:** Fix was necessary for correctness — ValueError must propagate. No scope creep.

## Issues Encountered

- `pytest` is not installed in `mace4ir_v2` micromamba environment. Used `.venv/bin/pytest` (project's development environment) instead. This matches how existing tests are run per pyproject.toml.

## Next Phase Readiness

- `mace_anicc` is now available as an energy calculator and optimization calculator via CLI
- Any plan using `--energy-calculators mace_anicc` or `--optimization-calculator mace_anicc` will work for HCNO molecules and fail fast with a clear message for others
- Requirement CALC-02 complete

---
*Phase: 13-calculator-expansion-acoh-bug-fix*
*Completed: 2026-03-03*
