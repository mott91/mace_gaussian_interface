---
phase: 13-calculator-expansion-acoh-bug-fix
plan: "01"
subsystem: cli
tags: [click, cli, validation, mace_off, mace_anicc, calculator]

# Dependency graph
requires: []
provides:
  - CLI callback validation for --energy-calculators (mace_mp, mace_omol, mace_off, mace_anicc)
  - CLI callback validation for --dipole-calculators (espaloma, mace_ml)
  - mace_anicc added to --optimization-calculator click.Choice
  - VALID_ENERGY_CALCULATORS and VALID_DIPOLE_CALCULATORS constants exported from cli.py
affects:
  - Phase 13 plans that test or invoke calculator-specific workflows
  - Any downstream plan referencing CLI options for mace_off or mace_anicc

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Click callback pattern for comma-separated multi-value options (avoids click.Choice which atomically validates)"
    - "Module-level VALID_* constants for allowlist-based validation — importable for testing"

key-files:
  created:
    - "tests/test_cli_validation.py (11 new test methods added)"
  modified:
    - "mace_gaussian/cli.py — VALID_ENERGY_CALCULATORS, VALID_DIPOLE_CALCULATORS constants; _validate_energy_calculators and _validate_dipole_calculators callbacks; mace_anicc added to --optimization-calculator choice"

key-decisions:
  - "Use callback= not type=click.Choice for comma-separated options; click.Choice validates the whole string atomically and would reject 'mace_mp,mace_omol'"
  - "Export VALID_* constants at module level so tests can import them directly for whitebox assertions"

patterns-established:
  - "TDD flow confirmed: write tests first (ImportError on missing constants = RED), then implement (GREEN)"
  - "Click callback validation pattern: split on comma, check each token, raise BadParameter with list of valid choices"

requirements-completed: [CALC-01, CALC-02]

# Metrics
duration: 6min
completed: 2026-03-03
---

# Phase 13 Plan 01: CLI Calculator Validation Summary

**Click callback validation wired to --energy-calculators and --dipole-calculators, with mace_anicc added to --optimization-calculator choice, enabling fast-fail on invalid names before any workflow code runs**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-03T17:31:32Z
- **Completed:** 2026-03-03T17:37:41Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments

- Added `VALID_ENERGY_CALCULATORS = ["mace_mp", "mace_omol", "mace_off", "mace_anicc"]` and `VALID_DIPOLE_CALCULATORS = ["espaloma", "mace_ml"]` as module-level constants in `cli.py`
- Added `_validate_energy_calculators` and `_validate_dipole_calculators` click callbacks that split on comma and check each token individually, raising `click.BadParameter` with a clear message on failure
- Added `mace_anicc` to `--optimization-calculator` click.Choice list
- Wired callbacks to both options and updated help text to list valid choices
- 11 new TDD tests across 3 test classes; all 18 tests pass including 7 pre-existing ones

## Task Commits

Each task was committed atomically:

1. **Task 1: Add CLI validation callbacks for calculator options** - `4f5963e` (feat)

## Files Created/Modified

- `mace_gaussian/cli.py` - Added VALID_ENERGY_CALCULATORS, VALID_DIPOLE_CALCULATORS constants; _validate_energy_calculators and _validate_dipole_calculators callbacks; mace_anicc in --optimization-calculator click.Choice
- `tests/test_cli_validation.py` - Added TestEnergyCalculatorValidation (6 tests), TestOptimizationCalculatorValidation (2 tests), TestDipoleCalculatorValidation (3 tests)

## Decisions Made

- Used `callback=` instead of `type=click.Choice(...)` for the comma-separated options. click.Choice validates atomically and would reject "mace_mp,mace_omol" as a whole string. The callback approach splits on comma first, enabling per-token validation.
- Exported `VALID_*` constants at module level so tests can import them and make assertions like `assert "mace_anicc" in VALID_ENERGY_CALCULATORS` directly.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- `pytest` is not installed in the `mace4ir_v2` micromamba environment; tests run correctly using `.venv/bin/python -m pytest` instead. The plan's verification command used `micromamba run -n mace4ir_v2 python -m pytest` but `.venv` is the correct test environment per pyproject.toml. No fix needed — used the correct environment.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- CLI layer now validates `mace_off` and `mace_anicc` properly; Plan 13-02 can implement the actual calculator backends knowing the CLI will route them correctly
- `VALID_ENERGY_CALCULATORS` constant can be imported by workflow.py or tests to stay in sync with CLI validation

---
*Phase: 13-calculator-expansion-acoh-bug-fix*
*Completed: 2026-03-03*
