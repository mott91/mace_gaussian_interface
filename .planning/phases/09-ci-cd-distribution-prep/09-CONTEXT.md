# Phase 9: CI/CD & Distribution Prep - Context

**Gathered:** 2026-02-25
**Status:** Ready for planning

<domain>
## Phase Boundary

Set up GitHub Actions CI that runs automatically on every push, and provide a clear
install procedure for a general researcher wanting to use the project on their own machine.
This does not include PyPI publishing, containers, or hosted docs — those are v2.

</domain>

<decisions>
## Implementation Decisions

### CI triggers
- Run on every push to **any branch** (not just main)
- Single workflow file is sufficient — no separate PR workflow needed

### Test scope in CI
- Run only tests that do **not** require custom local packages (`mace_ML_pkg`, `mace_dipole_pkg`),
  GPU, or Gaussian 16
- The existing markers already cover this: skip `gpu`, `gaussian`, and `slow` tests
- Unit tests for parsers, mode matching, utils, and validation run fine without the ML packages
- Do NOT attempt to install mace_ML_pkg or mace_dipole_pkg in CI — they are not on PyPI
  and the CI install would be fragile

### Ruff rule expansion (CI-03)
- Add rules B, SIM, PTH, RUF to the existing ruff config
- **Fix all existing violations now** with `ruff check --fix` before CI enforces them
- CI should then enforce: no new violations allowed (non-zero exit = fail)
- Both `ruff check` and `ruff format --check` run in CI

### Coverage reporting
- Add `pytest-cov` to CI run — report coverage for `mace_gaussian/` module
- **Do not set a failure threshold** — coverage is informational only
- Output as terminal summary (not uploaded to external service)

### Distribution / install procedure (CI-04)
- Target: a general researcher on their own machine (not the same HPC cluster)
- Deliverable: clear README section with step-by-step install instructions
- Must cover:
  1. Cloning the repo
  2. Creating a conda/micromamba environment (reference `environment.yml`)
  3. Installing the local custom packages (`mace_ML_pkg`, `mace_dipole_pkg`) with `pip install -e`
  4. Installing the main package with `pip install -e .`
  5. Verifying with `python cli.py diagnose`
- No bash install script needed — documented steps in README are sufficient

### Claude's Discretion
- Exact Python version(s) to test against in CI (pick one stable version, e.g. 3.10)
- Whether to cache pip/conda dependencies in the CI workflow
- Exact coverage report format (term-missing or just term)

</decisions>

<specifics>
## Specific Ideas

- The project already has `environment.yml` — the install docs should reference it
- `python cli.py diagnose` is the natural "did it work?" check for a new install
- CI should be fast — skipping GPU/Gaussian/slow tests is correct, not a compromise

</specifics>

<deferred>
## Deferred Ideas

- None came up — discussion stayed within phase scope

</deferred>

---

*Phase: 09-ci-cd-distribution-prep*
*Context gathered: 2026-02-25*
