---
phase: 12-distribution-polish
verified: 2026-02-28T00:00:00Z
status: passed
score: 5/5 must-haves verified
re_verification: false
---

# Phase 12: Distribution Polish Verification Report

**Phase Goal:** Complete distribution-readiness: real GitHub URLs, aligned CI versions, correct docs entry point, and clean planning metadata
**Verified:** 2026-02-28
**Status:** PASSED
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | pyproject.toml [project.urls] references `mott91/mace_gaussian_interface` (not yourusername) | VERIFIED | Lines 50-52: all three URL fields contain `mott91/mace_gaussian_interface`; grep for `yourusername` returns no matches |
| 2 | CI workflow installs `ruff>=0.9.0` (not `==0.15.1`), matching pyproject.toml dev dep floor | VERIFIED | `.github/workflows/ci.yml` line 18: `pip install "ruff>=0.9.0"`; grep for `ruff==` returns no matches |
| 3 | docs/quickstart.md shows `mace-gaussian <subcommand>` everywhere `cli.py` was referenced | VERIFIED | All CLI commands in quickstart.md use `mace-gaussian run`, `mace-gaussian list`, `mace-gaussian diagnose`; grep for `python cli.py` returns no matches in this file |
| 4 | README.md shows `mace-gaussian <subcommand>` everywhere `cli.py` was referenced (run_analysis.py references unchanged) | VERIFIED | All CLI commands use `mace-gaussian`; `python run_analysis.py water` appears at lines 99 and 151 (preserved); grep for `python cli.py` returns no matches |
| 5 | CLAUDE.md Commands section shows `mace-gaussian` for all cli.py invocations | VERIFIED | Lines 13-25 of CLAUDE.md use `mace-gaussian diagnose`, `mace-gaussian run`, `mace-gaussian list`; `python run_analysis.py` and `python run_analysis_harmonic.py` intentionally preserved |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `pyproject.toml` | Real GitHub URLs in [project.urls] | VERIFIED | Lines 50-52 contain Homepage, Repository, Issues all pointing to `https://github.com/mott91/mace_gaussian_interface` |
| `.github/workflows/ci.yml` | Aligned ruff version constraint | VERIFIED | Line 18: `pip install "ruff>=0.9.0"` — matches floor in `[dependency-groups] dev` |
| `docs/quickstart.md` | Correct entry-point in user-facing quickstart | VERIFIED | Contains `mace-gaussian run` at lines 23, 35, 43; `mace-gaussian list water` at line 43; `mace-gaussian diagnose` at lines 13, 111 |
| `README.md` | Correct entry-point in user-facing README | VERIFIED | Contains `mace-gaussian run`, `mace-gaussian list`, `mace-gaussian diagnose` throughout; no `python cli.py` occurrences |
| `CLAUDE.md` | Correct entry-point in developer docs | VERIFIED | Contains `mace-gaussian diagnose` at line 13; all 5 former cli.py invocations replaced |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| pyproject.toml [project.urls] | https://github.com/mott91/mace_gaussian_interface | string replacement of yourusername placeholder | WIRED | Three URL lines verified: Homepage (line 50), Repository (line 51), Issues (line 52) |
| docs/quickstart.md | mace-gaussian entry point | replace all python cli.py occurrences | WIRED | Pattern `mace-gaussian \w` found at multiple lines; no `python cli.py` residue |

### Requirements Coverage

No explicit requirement IDs declared in PLAN frontmatter `requirements` field (empty list). Phase was driven by success criteria derived from the v1.0 milestone audit. All five audit findings are resolved.

### Anti-Patterns Found

Scanned the five modified files for TODO/FIXME/placeholder patterns.

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `pyproject.toml` | 11 | `email = "your.email@domain.com"` | Info | Placeholder email in authors field — not in scope for this phase (URLs were the audit finding); no impact on distribution functionality |

No blockers or warnings found in the five target files.

Note: `python cli.py` references remain in out-of-scope files (`docs/DEVELOPMENT.md`, `docs/examples/water/README.md`, `.planning/` history files). These were not in the plan's `files_modified` list and do not affect distribution readiness for the identified gaps.

### Human Verification Required

None. All five must-haves are unambiguously verifiable via file content inspection.

### Gaps Summary

No gaps. All five truths verified, all five artifacts substantive and in place, both key links wired.

The one minor observation (placeholder email in `pyproject.toml` author field) was not in scope for this phase and does not block any of the stated must-haves. It may be worth a separate one-line fix before a public release, but it is not a distribution-polish gap per the phase definition.

---

_Verified: 2026-02-28_
_Verifier: Claude (gsd-verifier)_
