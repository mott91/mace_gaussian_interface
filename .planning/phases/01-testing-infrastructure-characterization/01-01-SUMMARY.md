---
phase: 01-testing-infrastructure-characterization
plan: 01
subsystem: testing
tags: [pytest, pytest-cov, fixtures, gaussian-parser, fchk-parser, conftest]

# Dependency graph
requires: []
provides:
  - "pytest configuration with gpu/gaussian/slow markers and strict mode"
  - "Test fixtures: water (DFT+ML log/fchk/json), CH4 (DFT+ML log/fchk/json), acoh (ML log)"
  - "conftest.py with 11 fixture functions for all test file paths"
  - "Coverage config targeting core modules (gaussian_parser, fchk_parser, mode_matching, etc.)"
  - "uv.lock tracked in git (REPR-04)"
affects: [01-02, 01-03, 01-04, 02-parser-hardening]

# Tech tracking
tech-stack:
  added: [pytest-cov>=4.0]
  patterns: [dependency-groups for dev deps, fixture trimming via scripts/create_fixtures.py]

key-files:
  created:
    - tests/conftest.py
    - tests/__init__.py
    - tests/fixtures/water/dft_b3lyp.log
    - tests/fixtures/water/dft_b3lyp.fchk
    - tests/fixtures/water/ml_mace_mp_esp.log
    - tests/fixtures/water/ml_mace_mp_esp.fchk
    - tests/fixtures/water/results.json
    - tests/fixtures/CH4_ase/dft_b3lyp.log
    - tests/fixtures/CH4_ase/dft_b3lyp.fchk
    - tests/fixtures/CH4_ase/ml_mace_mp_esp.log
    - tests/fixtures/CH4_ase/ml_mace_mp_esp.fchk
    - tests/fixtures/CH4_ase/results.json
    - tests/fixtures/acoh/ml_mace_mp_esp.log
    - scripts/create_fixtures.py
  modified:
    - .gitignore
    - pyproject.toml
    - uv.lock

key-decisions:
  - "Moved dev dependencies from orphaned [project] key to [dependency-groups] section for proper uv lockfile inclusion"
  - "Trimmed log files preserve exact section formatting for parser compatibility; verified output matches full originals"
  - "Acoh fixture preserves both bug root causes: missing Anharmonic IR section and H/L-prefixed lines"

patterns-established:
  - "Fixture naming: dft_b3lyp.{log,fchk} for DFT, ml_mace_mp_esp.{log,fchk} for ML"
  - "Fixture path access: FIXTURES_DIR constant + per-molecule fixture functions in conftest.py"
  - "Log trimming: extract only parser-relevant sections, verify round-trip with full originals"

# Metrics
duration: 9min
completed: 2026-02-16
---

# Phase 1 Plan 1: Testing Infrastructure Summary

**Pytest with 3 markers, pytest-cov, 13 committed test fixtures (water/CH4/acoh), and conftest.py with 11 fixture accessors**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-16T13:53:10Z
- **Completed:** 2026-02-16T14:02:12Z
- **Tasks:** 2
- **Files modified:** 17

## Accomplishments
- Configured pytest with gpu/gaussian/slow markers and strict-markers mode
- Added pytest-cov and coverage config targeting 10 core modules
- Created 13 test fixture files from comparison results (4 trimmed logs, 4 full fchk, 2 results.json, 1 acoh bug fixture)
- Built conftest.py with 11 fixture functions providing typed paths to all test files
- Fixed dev dependency configuration so uv.lock includes pytest/ruff/ty

## Task Commits

Each task was committed atomically:

1. **Task 1: Configure pytest, install pytest-cov, update .gitignore and pyproject.toml** - `cc96959` (chore)
2. **Task 2: Create test fixtures from comparison results and build conftest.py** - `43f259c` (feat)

## Files Created/Modified
- `.gitignore` - Removed uv.lock exclusion, added test fixture allowlist rules
- `pyproject.toml` - Added [dependency-groups], [tool.pytest.ini_options], [tool.coverage.run], [tool.coverage.report]
- `uv.lock` - Regenerated with dev dependencies (pytest, pytest-cov, ruff, ty)
- `tests/conftest.py` - 11 fixture functions for water/CH4/acoh file paths
- `tests/__init__.py` - Empty package marker
- `tests/fixtures/water/dft_b3lyp.log` - Trimmed water DFT log (88 lines, 3 harmonic + 3 anharmonic + 3 overtones + 6 combo bands)
- `tests/fixtures/water/dft_b3lyp.fchk` - Full water DFT formatted checkpoint (1701 lines, 3 atoms, 3 modes)
- `tests/fixtures/water/ml_mace_mp_esp.log` - Trimmed water ML log (89 lines)
- `tests/fixtures/water/ml_mace_mp_esp.fchk` - Full water ML formatted checkpoint (1289 lines)
- `tests/fixtures/water/results.json` - Water ML reference results
- `tests/fixtures/CH4_ase/dft_b3lyp.log` - Trimmed CH4 DFT log (230 lines, 4 harmonic + 9 anharmonic + 9 overtones + 72 combo)
- `tests/fixtures/CH4_ase/dft_b3lyp.fchk` - Full CH4 DFT formatted checkpoint (3114 lines, 5 atoms, 9 modes)
- `tests/fixtures/CH4_ase/ml_mace_mp_esp.log` - Trimmed CH4 ML log (231 lines)
- `tests/fixtures/CH4_ase/ml_mace_mp_esp.fchk` - Full CH4 ML formatted checkpoint (2235 lines)
- `tests/fixtures/CH4_ase/results.json` - CH4 ML reference results
- `tests/fixtures/acoh/ml_mace_mp_esp.log` - Trimmed acoh ML log (304 lines, demonstrates parser bug)
- `scripts/create_fixtures.py` - Reproducible fixture generation script

## Decisions Made
- **Moved dev deps to [dependency-groups]:** The existing `dev = [...]` under `[project]` was an orphaned TOML key -- uv never resolved pytest/ruff/ty into the lockfile. Moving to `[dependency-groups]` is the correct PEP 735 / uv approach and ensures `uv sync` installs dev dependencies.
- **Trimmed log strategy:** Extract header, frequency blocks with displacement data, thermochemistry, and full Anharmonic IR sections. Verified every trimmed fixture produces identical parser output to the full original.
- **Acoh bug preservation:** The acoh fixture contains both root causes -- missing "Anharmonic Infrared Spectroscopy" section (no I(anharm)/DS(anharm) columns) and H/L-prefixed lines in "Vibrational Energies at Anharmonic Level" section. Parser returns 0 anharmonic/overtone/combo results from this file.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed dev dependency configuration**
- **Found during:** Task 1 (pytest-cov installation)
- **Issue:** `dev = [...]` was an orphaned key under `[project]` table, not recognized by uv. Dev dependencies (pytest, ruff, ty) were never included in uv.lock or installed by `uv sync`.
- **Fix:** Moved dev dependency list to `[dependency-groups]` section (PEP 735 standard, supported by uv). Regenerated uv.lock which now includes pytest, pytest-cov, ruff, ty and their transitive deps.
- **Files modified:** pyproject.toml, uv.lock
- **Verification:** `uv tree --group dev` shows pytest and pytest-cov in dependency tree
- **Committed in:** cc96959 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Essential fix -- without it, dev dependencies would not be installed. No scope creep.

## Issues Encountered
- `dgl==2.2.1` (transitive dependency via espaloma-charge) has no Linux wheel. Worked around with `uv sync --no-install-package dgl`. This is a pre-existing issue, not introduced by this plan.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Testing infrastructure is fully configured and ready for test writing
- All fixture files committed and verified against parsers
- Plans 01-02 through 01-04 can now write tests that use these fixtures via conftest.py fixtures
- Coverage reporting ready via `pytest --cov`

## Self-Check: PASSED

All 17 created/modified files verified present on disk. Both task commits (cc96959, 43f259c) verified in git log.

---
*Phase: 01-testing-infrastructure-characterization*
*Completed: 2026-02-16*
