---
phase: 09-ci-cd-distribution-prep
plan: "03"
subsystem: docs
tags: [readme, installation, documentation, reproducibility]

# Dependency graph
requires:
  - phase: 09-01
    provides: pyproject.toml wired with mace-gaussian entry point and package structure
provides:
  - README.md Installation section with explicit step-by-step instructions (5 commands visible)
  - CI-04 documented install procedure satisfied
affects: [future contributors, thesis reviewers, CI documentation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Show manual steps first, convenience scripts second — never hide required commands inside a script reference"

key-files:
  created: []
  modified:
    - README.md

key-decisions:
  - "Manual pip install steps shown explicitly before mentioning install_mace_packages.sh convenience script"
  - "Prerequisites note Python 3.10 requirement with mace4ir_v2 environment name for clarity"
  - "install_mace_packages.sh retained as convenience shortcut in Environment Files section, not removed"

patterns-established:
  - "Installation docs: numbered bold headings (1-5), one code block per step, brief explanatory prose after each block"

requirements-completed:
  - CI-04

# Metrics
duration: 2min
completed: 2026-02-26
---

# Phase 9 Plan 03: README Installation Section Summary

**README Installation section rewritten with five explicit pip install steps visible — researchers can reproduce the environment without reading the install script**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-26T11:06:42Z
- **Completed:** 2026-02-26T11:07:45Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments

- Replaced "Quick Start" (which hid pip install steps inside install_mace_packages.sh) with numbered "Step-by-Step Installation" showing each command explicitly
- `pip install -e mace_ML_pkg` and `pip install -e mace_dipole_pkg` are now visible separate steps
- `environment.yml` referenced by name so researchers know exactly which file to use
- `python cli.py diagnose` is the final verification step
- `install_mace_packages.sh` kept but demoted to convenience note in Environment Files section

## Task Commits

1. **Task 1: Rewrite README Installation section with explicit manual steps** - `a5f2a54` (feat)

**Plan metadata:** (docs commit below)

## Files Created/Modified

- `README.md` — Installation section replaced (lines 14-62 replaced with explicit 5-step procedure)

## Decisions Made

- Manual pip install steps shown first before mentioning the convenience script — avoids the "Pitfall 5" anti-pattern of hiding required commands in a script reference
- Prerequisites bullet updated to note `mace4ir_v2` environment name explicitly (matches what users will activate)
- install_mace_packages.sh retained as convenience shortcut after the manual steps (not removed — still useful for repetitive installs)

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- CI-04 documented install procedure is satisfied
- README Installation section is now suitable for thesis appendix and GitHub front page
- Ready for plan 09-04 (remaining CI/CD distribution prep plans)

---
*Phase: 09-ci-cd-distribution-prep*
*Completed: 2026-02-26*
