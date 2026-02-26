---
phase: 09-ci-cd-distribution-prep
verified: 2026-02-26T12:30:00Z
status: passed
score: 10/10 must-haves verified
re_verification: false
---

# Phase 09: CI/CD and Distribution Prep Verification Report

**Phase Goal:** Automate testing, linting, and prepare infrastructure for distribution to other researchers
**Verified:** 2026-02-26T12:30:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | ruff check with B, SIM, PTH, RUF rules exits 0 on mace_gaussian/ and tests/ | VERIFIED | `ruff check mace_gaussian/ tests/` → "All checks passed!" |
| 2 | ruff format --check exits 0 on mace_gaussian/ and tests/ | VERIFIED | `ruff format --check mace_gaussian/ tests/` → "41 files already formatted" |
| 3 | All 128 tests pass (test_path_suffix stale test fixed) | VERIFIED | `pytest tests/ -m "not gpu and not gaussian and not slow" -q` → "128 passed, 2 deselected, 1 xfailed in 4.03s" |
| 4 | pyproject.toml ruff select includes B, SIM, PTH, RUF | VERIFIED | Line 62: `select = ["E", "F", "W", "I", "UP", "B", "SIM", "PTH", "RUF"]` |
| 5 | .github/workflows/ci.yml exists and is valid YAML | VERIFIED | File exists at 48 lines; `python3 -c "import yaml; yaml.safe_load(...)"` → "YAML valid" |
| 6 | Workflow triggers on every push to any branch (no branch filter) | VERIFIED | `on:\n  push:` with no `branches:` key — confirmed by grep |
| 7 | lint job runs ruff format --check and ruff check | VERIFIED | Lines 21-24 of ci.yml: `ruff format --check mace_gaussian/ tests/` then `ruff check mace_gaussian/ tests/` |
| 8 | test job excludes gpu, gaussian, slow markers; runs coverage to terminal only | VERIFIED | Lines 46-48: `-m "not gpu and not gaussian and not slow" --cov=mace_gaussian --cov-report=term-missing`; no `--cov-fail-under` |
| 9 | README Installation section shows all 5 explicit steps | VERIFIED | All 5 strings present: `pip install -e mace_ML_pkg`, `pip install -e mace_dipole_pkg`, `environment.yml`, `python cli.py diagnose`, `pip install -e .` |
| 10 | README references environment.yml by name and install_mace_packages.sh as convenience only | VERIFIED | `environment.yml` at line 34 and 82; `install_mace_packages.sh` demoted to Environment Files section (line 83) |

**Score:** 10/10 truths verified

---

## Required Artifacts

### Plan 09-01 Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `pyproject.toml` | Expanded ruff rule set | VERIFIED | `select = ["E", "F", "W", "I", "UP", "B", "SIM", "PTH", "RUF"]` at line 62 |
| `tests/test_calculators.py` | Fixed stale test_path_suffix assertion | VERIFIED | Contains `mace4ir_models/pretrained_models/model_1_dipole.model` at line 209 |

### Plan 09-02 Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `.github/workflows/ci.yml` | GitHub Actions CI workflow | VERIFIED | 48-line file exists; valid YAML; created in commit b8251f1 |

### Plan 09-03 Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `README.md` | Updated Installation section with manual step-by-step | VERIFIED | Numbered 5-step section present; `pip install -e mace_ML_pkg` at line 53; commit a5f2a54 |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `pyproject.toml [tool.ruff.lint]` | `ruff check mace_gaussian/ tests/` | ruff reads pyproject.toml for rule config | WIRED | ruff 0.15.1 reads select from pyproject.toml; live run confirmed 0 violations |
| `.github/workflows/ci.yml (lint job)` | `ruff check mace_gaussian/ tests/` | GitHub Actions step run | WIRED | Line 24: `run: ruff check mace_gaussian/ tests/` present verbatim |
| `.github/workflows/ci.yml (test job)` | `pytest tests/ -m "not gpu and not gaussian and not slow"` | GitHub Actions step run | WIRED | Lines 46-47 contain exact marker exclusion string |
| `README.md Installation section` | `environment.yml` | `micromamba env create -f environment.yml` reference | WIRED | Line 41: `micromamba env create -f environment.yml` |

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| CI-01 | 09-02 | GitHub Actions workflow runs ruff check and ruff format --check on push | SATISFIED | `.github/workflows/ci.yml` lint job runs both checks on every push |
| CI-02 | 09-02 | GitHub Actions workflow runs pytest (unit tests only, no GPU/Gaussian) on push | SATISFIED | `.github/workflows/ci.yml` test job runs pytest with `-m "not gpu and not gaussian and not slow"` |
| CI-03 | 09-01 | Ruff rules expanded to include B, SIM, PTH, RUF | SATISFIED | `pyproject.toml` select includes all four rule sets; live ruff check exits 0 |
| CI-04 | 09-03 | Install script or documented procedure for custom MACE packages | SATISFIED | README Installation section shows explicit `pip install -e mace_ML_pkg` and `pip install -e mace_dipole_pkg` steps |

**Orphaned requirements:** None. All CI-01 through CI-04 appear in plan frontmatter and REQUIREMENTS.md maps them to Phase 9.

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `mace_gaussian/cli.py` | 278 | `[PLACEHOLDER - Coming soon]` in `compare` command docstring | Info | Pre-existing stub for unimplemented `compare` CLI command; unrelated to phase 09 goals; ruff passes; no test coverage gap for this phase |
| `mace_gaussian/cli.py` | 324 | `[PLACEHOLDER - Coming soon]` in `export` command docstring | Info | Pre-existing stub for unimplemented `export` CLI command; same as above |

**Assessment:** Both placeholders are in docstrings of CLI commands that were not targets of phase 09. They were present before this phase and are not blockers for the CI/CD or distribution goals. Ruff does not flag docstring content. Neither command is tested, which is consistent with their stub nature.

---

## Commit Verification

All commits documented in SUMMARY files were verified to exist in git history:

| Commit | Plan | Description |
|--------|------|-------------|
| `10dae2a` | 09-01 Task 1 | chore(09-01): expand ruff rules to B,SIM,PTH,RUF and apply auto-fixes |
| `26a98ff` | 09-01 Task 2 | fix(09-01): manually fix all remaining ruff violations (PTH, SIM, B, RUF, UP) |
| `b8251f1` | 09-02 Task 1 | feat(09-02): add GitHub Actions CI workflow |
| `a5f2a54` | 09-03 Task 1 | feat(09-03): rewrite README Installation section with explicit manual steps |

---

## Human Verification Required

### 1. CI Workflow Execution on GitHub

**Test:** Push the current branch to GitHub and observe the Actions tab.
**Expected:** Both "Lint" and "Test" jobs run in parallel; Lint job completes with green checkmark; Test job shows "128 passed" in logs.
**Why human:** Cannot trigger GitHub Actions locally — requires a push to a GitHub-hosted repository.

### 2. Installation Procedure Clarity for External Researchers

**Test:** Ask a researcher unfamiliar with the project to follow the README Installation section on a fresh machine.
**Expected:** They can clone, create the environment, install all packages, and run `python cli.py diagnose` without needing to read any other file.
**Why human:** Clarity and completeness of prose cannot be verified programmatically — requires a human reader.

---

## Summary

Phase 09 achieved its goal. All three plans executed cleanly:

- **Plan 09-01 (CI-03):** Ruff rule set expanded to B/SIM/PTH/RUF; 55 violations fixed across 23 files; 128 tests passing including the previously-stale `test_path_suffix`. Live `ruff check` and `ruff format --check` both exit 0.

- **Plan 09-02 (CI-01, CI-02):** `.github/workflows/ci.yml` created with two parallel jobs. Lint job pins ruff==0.15.1 and checks format then violations on `mace_gaussian/` and `tests/`. Test job installs without heavy ML deps (`--no-deps`) to avoid dgl Windows-wheel failure on ubuntu-latest, then runs pytest with correct marker exclusions and terminal coverage (no failure threshold). YAML validates cleanly.

- **Plan 09-03 (CI-04):** README Installation section rewritten from "Quick Start (run the script)" to a numbered 5-step procedure showing each pip install command explicitly. `environment.yml` referenced by name. `install_mace_packages.sh` demoted to convenience note.

The codebase is ready for CI enforcement: any push will trigger automated lint and test checks without requiring GPU, Gaussian 16, or heavy ML dependency wheels.

---

_Verified: 2026-02-26T12:30:00Z_
_Verifier: Claude (gsd-verifier)_
