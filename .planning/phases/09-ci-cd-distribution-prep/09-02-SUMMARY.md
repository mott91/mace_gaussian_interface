---
phase: 09-ci-cd-distribution-prep
plan: 02
subsystem: infra
tags: [github-actions, ci, ruff, pytest, coverage]

requires:
  - phase: 09-01
    provides: clean-ruff-baseline with zero violations across mace_gaussian/ and tests/

provides:
  - GitHub Actions CI workflow (.github/workflows/ci.yml)
  - Automated lint checks on every push (ruff format + ruff check)
  - Automated unit tests on every push (pytest with coverage)

affects: [all-future-branches, pr-workflow]

tech-stack:
  added: [github-actions]
  patterns: [parallel-ci-jobs, lightweight-dep-install-for-ci, marker-based-test-exclusion]

key-files:
  created:
    - .github/workflows/ci.yml
  modified: []

key-decisions:
  - "on: push with no branch filter — CI runs on every push to every branch (locked decision)"
  - "pip install -e . --no-deps avoids dgl==2.2.1 Windows-only wheels that block uv sync on ubuntu-latest"
  - "ruff==0.15.1 pinned to match uv.lock version for reproducibility"
  - "lint and test jobs run in parallel (no needs: dependency) for fastest CI"
  - "coverage is informational only — no --cov-fail-under threshold (locked decision)"

patterns-established:
  - "CI pattern: lightweight install (--no-deps + explicit test deps) to avoid heavy ML dependency wheels"
  - "Marker exclusion: not gpu and not gaussian and not slow keeps CI fast without GPU/Gaussian"

requirements-completed: [CI-01, CI-02]

duration: 5min
completed: 2026-02-26
---

# Phase 09 Plan 02: GitHub Actions CI Workflow Summary

**Single-file GitHub Actions CI with parallel lint (ruff) and test (pytest) jobs that run on every push without dgl/espaloma/mace_torch wheel failures.**

## Performance

- **Duration:** ~5 min
- **Started:** 2026-02-26T11:01:25Z
- **Completed:** 2026-02-26T11:06:00Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments

- Created `.github/workflows/ci.yml` with two parallel jobs: lint and test
- Lint job installs pinned `ruff==0.15.1` and runs `ruff format --check` then `ruff check` on `mace_gaussian/` and `tests/`
- Test job installs package without heavy deps (`--no-deps`) plus lightweight test dependencies, then runs pytest excluding gpu/gaussian/slow markers with terminal coverage report

## Task Commits

Each task was committed atomically:

1. **Task 1: Create .github/workflows/ci.yml** - `b8251f1` (feat)

**Plan metadata:** (final docs commit follows)

## Files Created/Modified

- `.github/workflows/ci.yml` - GitHub Actions CI workflow with lint and test parallel jobs

## Decisions Made

- `pip install -e . --no-deps` chosen over `uv sync --locked` because dgl==2.2.1 has Windows-only wheels, causing `uv sync` to fail on `ubuntu-latest`. This is the key pitfall identified in 09-RESEARCH.md.
- `ruff==0.15.1` pinned (matches uv.lock) for reproducibility.
- No `--cov-fail-under` threshold — coverage is informational only, locked decision from research phase.
- Parallel jobs (no `needs:` field) for fastest possible CI feedback.
- `python-version: "3.10"` matches the mace4ir_v2 production environment.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - straightforward file creation. YAML validated successfully on first attempt.

## User Setup Required

None - no external service configuration required. The workflow will activate automatically when pushed to a GitHub repository.

## Next Phase Readiness

- CI workflow is in place and ready to enforce lint/test quality on all future pushes
- Plans 09-03 through 09-05 can proceed (pre-commit hooks, pyproject.toml metadata, release notes)

---
*Phase: 09-ci-cd-distribution-prep*
*Completed: 2026-02-26*
