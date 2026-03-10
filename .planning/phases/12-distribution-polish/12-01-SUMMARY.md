---
phase: 12-distribution-polish
plan: "01"
subsystem: docs
tags: [pyproject, ci, readme, quickstart, entry-point, distribution]

# Dependency graph
requires:
  - phase: 08-distribution-packaging
    provides: mace-gaussian entry point wired in pyproject.toml scripts
  - phase: 10-documentation
    provides: quickstart.md and README.md created
  - phase: 09-ci-testing
    provides: .github/workflows/ci.yml created
provides:
  - Real GitHub repository URLs in pyproject.toml [project.urls]
  - CI ruff version aligned with pyproject.toml dev dep floor (>=0.9.0)
  - All user-facing docs (quickstart.md, README.md) show mace-gaussian entry point
  - Developer docs (CLAUDE.md) Commands section shows mace-gaussian entry point
affects:
  - external researchers installing from repository
  - CI reproducibility across ruff releases

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Entry point consistency: all doc references use installed mace-gaussian command, not python cli.py"
    - "CI version floor matching: CI ruff pin matches pyproject.toml dev dep floor"

key-files:
  created: []
  modified:
    - pyproject.toml
    - .github/workflows/ci.yml
    - docs/quickstart.md
    - README.md
    - CLAUDE.md

key-decisions:
  - "Replace yourusername placeholder in [project.urls] with mott91/mace_gaussian_interface (real GitHub coordinates)"
  - "Use ruff>=0.9.0 floor in CI instead of ==0.15.1 exact pin — removes maintenance burden for research project (audit: LOW risk)"
  - "Replace all python cli.py with mace-gaussian in user-facing and developer docs, preserving python run_analysis.py and python run_analysis_harmonic.py unchanged"

patterns-established:
  - "Distribution readiness: no placeholder URLs, no CI version skew, consistent entry-point docs"

requirements-completed: []

# Metrics
duration: 2min
completed: 2026-02-27
---

# Phase 12 Plan 01: Distribution Polish — Placeholder URLs, CI Ruff Pin, Entry-Point Docs Summary

**Removed all four distribution-readiness gaps: replaced yourusername placeholder URLs with real mott91/mace_gaussian_interface coordinates, aligned CI ruff pin from ==0.15.1 to >=0.9.0, and replaced 18 python cli.py occurrences across quickstart.md, README.md, and CLAUDE.md with the installed mace-gaussian entry point.**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-27T10:56:31Z
- **Completed:** 2026-02-27T10:58:37Z
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments
- pyproject.toml [project.urls] now points to https://github.com/mott91/mace_gaussian_interface (Homepage, Repository, Issues — 3 lines)
- .github/workflows/ci.yml ruff pin changed from exact ==0.15.1 to floor >=0.9.0, matching pyproject.toml dev dep declaration
- docs/quickstart.md updated: 5 occurrences of `python cli.py` replaced with `mace-gaussian`
- README.md updated: 8 occurrences of `python cli.py` replaced with `mace-gaussian`, run_analysis.py references preserved
- CLAUDE.md Commands section updated: 5 occurrences of `python cli.py` replaced with `mace-gaussian`, run_analysis scripts left unchanged

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix pyproject.toml placeholder URLs and align CI ruff pin** - `8132213` (chore)
2. **Task 2: Update user-facing docs to mace-gaussian entry point** - `a11f62c` (docs)
3. **Task 3: Update CLAUDE.md Commands section to mace-gaussian entry point** - `c65fe00` (docs)

## Files Created/Modified
- `pyproject.toml` - Three [project.urls] entries replaced with real GitHub repository coordinates
- `.github/workflows/ci.yml` - ruff install changed from exact pin to floor version
- `docs/quickstart.md` - 5 python cli.py occurrences replaced with mace-gaussian
- `README.md` - 8 python cli.py occurrences replaced with mace-gaussian
- `CLAUDE.md` - 5 python cli.py occurrences replaced with mace-gaussian

## Decisions Made
- Use ruff>=0.9.0 in CI (floor match) rather than an exact pin — removes the maintenance burden of keeping CI and pyproject.toml in sync for a research project where the audit classified this LOW risk
- Preserve python run_analysis.py and python run_analysis_harmonic.py in all docs — these are standalone scripts with no entry-point equivalent and must remain as-is

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 12 Plan 01 complete — four distribution-readiness gaps from the v1.0 audit closed
- Package now presents correctly to external researchers: real URLs, consistent entry point, no version skew
- Phase 12 further plans (if any) can proceed; project is distribution-ready

## Self-Check: PASSED

All modified files present: pyproject.toml, .github/workflows/ci.yml, docs/quickstart.md, README.md, CLAUDE.md, 12-01-SUMMARY.md.
All task commits verified: 8132213 (Task 1), a11f62c (Task 2), c65fe00 (Task 3).

---
*Phase: 12-distribution-polish*
*Completed: 2026-02-27*
