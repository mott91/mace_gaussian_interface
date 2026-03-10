---
phase: 13-calculator-expansion-acoh-bug-fix
plan: 03
subsystem: testing
tags: [parser, regex, gaussian, anharmonic, acoh, bug-fix]

# Dependency graph
requires:
  - phase: 13-calculator-expansion-acoh-bug-fix
    provides: "acoh fixture at tests/fixtures/acoh/ml_mace_mp_esp.log; xfail test documenting the bug"
provides:
  - "Dual-format anharmonic parser supporting both Format A (DFT, I(anharm)) and Format B (ML external, Status+E(anharm))"
  - "H/L-prefix-aware regex for Format B lines"
  - "Negative frequency handling in Format B regex"
  - "Promoted test_acoh_anharmonic_parsing from xfail to passing"
affects:
  - "Phase 13 analysis plans — acoh regression plots now have 18 anharmonic mode data points"
  - "Any future molecule run via ML external path produces Format B logs"

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Dual-format section detection: lookahead in 'Fundamental Bands' trigger checks for I(anharm)/DS(anharm) (Format A) vs E(anharm)+Status (Format B)"
    - "in_format_b flag drives separate data-line regex branches inside a single parser pass"
    - "Format B ir_intensity set to 0.0 — no intensity column in Vibrational Energies at Anharmonic Level output"

key-files:
  created: []
  modified:
    - mace_gaussian/gaussian/parser.py
    - tests/test_gaussian_parser.py

key-decisions:
  - "Format B entries return ir_intensity=0.0 — the ML external log format simply does not emit intensity; downstream consumers must handle this gracefully"
  - "Negative frequencies are valid in Format B (imaginary modes from unconverged geometry); regex updated to -?[\\d\\.]+"
  - "H/L prefix (high/low mode overlap indicator) is stripped; mode number is still captured correctly"

patterns-established:
  - "Parser deviation detection: when a section header is found, lookahead up to 10 lines to classify the format before processing data lines"

requirements-completed: [FIX-01]

# Metrics
duration: 2min
completed: 2026-03-03
---

# Phase 13 Plan 03: Acoh Anharmonic Parser Fix Summary

**Dual-format anharmonic frequency parser: adds Format B support for ML external Gaussian logs (Vibrational Energies at Anharmonic Level with H/L overlap indicators), fixing acoh and all ML-run molecules returning 0 anharmonic modes**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-03T17:32:14Z
- **Completed:** 2026-03-03T17:34:00Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments
- `parse_anharmonic_frequencies()` now handles both Format A (DFT logs with I(anharm) intensity column) and Format B (ML external logs with Status column, no intensity)
- H/L-prefixed lines (high/low overlap indicators) are correctly parsed with mode number extracted
- Negative frequencies (imaginary modes, e.g., mode 18 at -255.103 cm-1 in acoh) are handled
- `test_acoh_anharmonic_parsing` promoted from `@pytest.mark.xfail` to a normal passing test; suite goes from 20 passed + 1 xfailed to 21 passed

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix dual-format parser and remove xfail** - `6ccfd03` (fix)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `mace_gaussian/gaussian/parser.py` - Added `in_format_b` flag, dual lookahead trigger, Format B regex branch with H/L support and negative frequency handling
- `tests/test_gaussian_parser.py` - Removed `@pytest.mark.xfail` from `test_acoh_anharmonic_parsing`, updated docstring, added per-entry key assertions

## Decisions Made
- Format B entries return `ir_intensity=0.0` because the Gaussian ML external log format ("Vibrational Energies at Anharmonic Level") simply does not emit IR intensities — only rotational constants (Aa/Ba/Ca). Downstream plotting code that filters on `ir_intensity > 0` will correctly exclude these.
- Negative frequencies are valid and must be captured — mode 18 of acoh is an imaginary mode from the geometry not being perfectly converged; excluding it would cause the count to be 17 instead of 18.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Negative frequency in Format B regex**
- **Found during:** Task 1 (parser fix and test run)
- **Issue:** Initial Format B regex used `[\d\.]+` which doesn't match negative numbers; mode 18 (`-255.103  -250.532`) was silently skipped, giving 17 instead of 18 modes
- **Fix:** Updated regex to `-?[\d\.]+` for both E(harm) and E(anharm) capture groups
- **Files modified:** mace_gaussian/gaussian/parser.py
- **Verification:** 21 tests pass; acoh sanity check returns 18 modes
- **Committed in:** 6ccfd03 (same task commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug)
**Impact on plan:** Necessary correctness fix discovered during testing. No scope creep.

## Issues Encountered
- The acoh log contains an imaginary mode (negative frequency) for mode 18 due to a slightly unconverged geometry — this is expected physics, not a data error. The parser must capture it.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Acoh anharmonic parser is fixed; regression analysis plots for acoh ML runs will now have 18 data points instead of 0
- Format B support covers all future molecules run via the ML external calculation path
- No blockers for remaining Phase 13 plans

---
*Phase: 13-calculator-expansion-acoh-bug-fix*
*Completed: 2026-03-03*
