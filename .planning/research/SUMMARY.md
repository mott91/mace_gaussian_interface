# Project Research Summary

**Project:** MACE-Gaussian Refactoring & Distribution
**Domain:** Scientific Python Package for Computational Chemistry
**Researched:** 2026-02-16
**Confidence:** HIGH

## Executive Summary

The MACE-Gaussian project is a functional research tool that bridges ML potentials (MACE) with Gaussian 16 for molecular IR spectroscopy. The codebase was built with AI assistance and is currently thesis-ready but not distribution-ready. The primary structural issue is a monolithic main module (gm_main.py, ~1000 lines) that mixes workflow orchestration, I/O, ZMQ server management, and calculator implementations. The highest-risk component is module monkey-patching used to swap MACE implementations at runtime.

The recommended approach is **incremental refactoring protected by tests**, not a rewrite. The codebase has good foundations (uv, pyproject.toml, ruff) but lacks the defensive infrastructure (tests, error handling, documentation) needed for distribution. The critical success factor is establishing a comprehensive test suite BEFORE refactoring anything — this protects against silent behavior changes in numerical code and fragile parsers.

Key risks are (1) breaking the working pipeline during refactoring, (2) losing reproducibility of thesis results, and (3) MACE dependency issues when replacing monkey-patching. Mitigation: commit reference outputs before refactoring, maintain backward compatibility throughout, and keep the tool functional at every step.

## Key Findings

### Recommended Stack

The project already uses modern Python packaging tools (uv, pyproject.toml, hatchling) which are the current standard. The main additions needed are testing infrastructure and CI/CD for maintainability.

**Core technologies:**
- **pytest + pytest-cov**: Testing framework with coverage reporting — industry standard for scientific Python, supports fixtures for reference data
- **pytest-mock + unittest.mock**: Mocking external dependencies (Gaussian, CUDA) — allows testing parsers without installing licensed software
- **GitHub Actions CI**: Automated linting and unit tests — follows ASE/cclib pattern of running parser tests in CI while excluding integration tests
- **ruff (expand rules)**: Add B, SIM, PTH, RUF rules for better code quality — fast enough to run manually without pre-commit hooks
- **MkDocs + mkdocs-material**: User documentation — simpler than Sphinx, easier for non-programmers to maintain
- **Test fixtures**: Committed Gaussian .log/.fchk snippets — standard pattern for tools with external dependencies (ASE, MDAnalysis, cclib all do this)

**Type checking note:** Current ty (alpha) is acceptable for thesis scope; switch to mypy only if ty causes issues. Don't over-invest in type annotations — focus on critical interfaces (parsers, calculators).

### Expected Features

The tool already implements the core differentiators (automated pipeline, mode matching, multi-model comparison, HTML reports). The main gaps are in table-stakes features that users expect from any distributed scientific tool.

**Must have (table stakes):**
- **Unit tests with reference data** — users need confidence that parsers and mode matching produce correct results
- **Explicit failure over silent corruption** — if a parser returns empty data, raise an error with context; never silently return empty results
- **Prerequisite checking** — validate Gaussian, formchk, CUDA, dipole model exist before starting multi-hour calculations
- **Documentation** — README with installation + quickstart, worked example, method documentation for papers
- **CLI experience** — consistent exit codes, progress indication, verbose/quiet modes
- **Reproducibility metadata** — version, Python version, model versions, config in results.json

**Should have (competitive):**
- **Test markers** — `@pytest.mark.gpu`, `@pytest.mark.gaussian` to separate integration tests from unit tests
- **Config file support** — for HPC job scripts where CLI args are awkward (ConfigArgParse already in dependencies)
- **Coverage reporting** — pytest-cov to identify untested code paths

**Defer (v2+):**
- Web interface, database backend, multi-user support, automatic model training, cloud deployment, plugin system, real-time visualization, support for non-Gaussian QM packages — all explicitly anti-features for thesis scope.

### Architecture Approach

The refactoring goal is to decompose the monolithic gm_main.py into focused modules without breaking the working pipeline. The rest of the codebase (parsers, analysis, reports) is already reasonably modular and just needs reorganization.

**Major components:**
1. **calculators/** — Extract DipoleCalculator classes (base, espaloma, mace_ml, xtb, factory, mace_loader) from gm_main.py. Replace monkey-patching with lazy import isolation or process isolation.
2. **gaussian/** — Extract I/O (parse_gaussian_input, write outputs), runner (subprocess management), parser (move GaussianLogParser), fchk (move FCHK parsing), zmq_server (context manager for message loop) from gm_main.py.
3. **analysis/** — Move existing modules (spectra, mode_matching, comparison, report) into package structure. Already clean, just reorganize.
4. **utils/** — Extract pure functions (units, validation, results) from gm_main.py.
5. **workflow.py** — Thin orchestrator that sequences phases by calling modular components. Extract from gm_main.py last.

**Critical constraint:** Each refactoring step must keep the tool functional. Never have a broken state. Pattern: write tests → refactor → verify tests pass → commit.

### Critical Pitfalls

Research identified 7 pitfalls with varying severity. The top 5 that directly impact roadmap planning:

1. **Refactoring without tests (CRITICAL)** — Silent behavior changes in numerical code and parsers. Mitigation: Add characterization tests with water/CH4 reference outputs BEFORE any refactoring. Use numpy.testing for float comparisons. Run full pipeline after every change.

2. **Breaking the working pipeline (HIGH)** — ZMQ + Gaussian subprocess integration is complex and hard to debug. Mitigation: Extract from gm_main.py, don't rewrite it. Each extraction: move code → update imports → verify pipeline runs on water. Don't change helper script mechanism until everything else is stable.

3. **Losing reproducibility (HIGH)** — Thesis results must be reproducible but refactoring can change floating-point accumulation or import order. Mitigation: Commit reference outputs before refactoring. After refactoring, verify outputs match within tolerance. Pin dependencies (uv.lock helps). Record versions in output metadata.

4. **MACE dependency hell (HIGH)** — Custom local packages have implicit dependencies on import order and CUDA initialization. Monkey-patching exists because of these complexities. Mitigation: Document exact import sequence before changing. Test MACE loading in isolation. Start with lazy imports (less disruptive than process isolation). Keep rollback plan.

5. **Over-engineering (MEDIUM but HIGH likelihood)** — Risk of turning simple file moves into architecture redesigns. Mitigation: Only refactor what's actually broken. If existing code works and isn't fragile, move it but don't redesign it. Three concrete duplications before abstracting. Ask "does this make the thesis easier?" for every change.

**Additional pitfalls:** AI code inconsistencies (different conventions across sessions, inconsistent units/error handling), acetic acid parsing bug (documented but unfixed, needs test case).

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Testing Infrastructure & Characterization
**Rationale:** Tests are the safety net for all subsequent refactoring. Without tests, refactoring risks silent behavior changes in numerical code. This is the highest priority from pitfalls research.
**Delivers:** pytest setup, test fixtures (committed Gaussian .log/.fchk snippets), characterization tests for water/CH4 that capture exact current behavior, test markers for gpu/gaussian-dependent tests.
**Addresses:** Table stakes features (tests), critical pitfall #1 (refactoring without tests), pitfall #5 (losing reproducibility — establishes baseline).
**Avoids:** Silent breakage during refactoring.
**Research flag:** Standard pattern — pytest with fixtures is well-documented. No phase-specific research needed.

### Phase 2: Error Handling & Input Validation
**Rationale:** Can be done before structural refactoring. Adds defensive programming that makes refactoring safer. Addresses table stakes features without touching the monolithic structure yet.
**Delivers:** Explicit failures with context in parsers, prerequisite checking (Gaussian/CUDA/dipole model exist), improved error messages, environment diagnostics expansion.
**Addresses:** Table stakes features (error handling, prerequisite checking), pitfall #2 (breaking pipeline — better errors make debugging easier).
**Uses:** Pure additions to existing code, no structural changes yet.
**Research flag:** Standard patterns — no additional research needed.

### Phase 3: Extract Utilities & Conventions
**Rationale:** Pure functions (units, validation) have zero coupling risk. Extracting them first builds confidence in the refactoring process. Also addresses AI code inconsistencies early.
**Delivers:** utils/units.py (unit conversions), utils/validation.py (input checks), utils/results.py (move ResultsManager), documented conventions for naming/error handling/units.
**Addresses:** Pitfall #4 (AI code inconsistencies — establishes conventions), architecture goal (decompose gm_main.py).
**Avoids:** Starting with high-risk refactoring (calculators, MACE loading).
**Research flag:** No research needed — straightforward extraction.

### Phase 4: Extract Calculator Classes
**Rationale:** Calculators are self-contained with a clear factory pattern already in place. Extract before tackling MACE loading complexity.
**Delivers:** calculators/ package (base, espaloma, mace_ml, xtb, factory), updated imports, tests for calculator factory.
**Addresses:** Architecture goal (decompose gm_main.py), prepares for next phase (MACE loading).
**Uses:** Factory pattern already exists, just moving code.
**Research flag:** No research needed — extraction only, not redesign.

### Phase 5: Replace MACE Module Monkey-Patching
**Rationale:** Highest-risk refactoring task. Must be done after tests are in place and calculators are extracted. Start with lazy import isolation (less disruptive than process isolation).
**Delivers:** calculators/mace_loader.py with lazy import isolation, cleanup of sys.modules manipulation, tests for MACE loading in isolation.
**Addresses:** Pitfall #7 (MACE dependency hell — the root cause of monkey-patching), architecture fragility (highest-risk component).
**Avoids:** Breaking MACE imports by testing in isolation first.
**Research flag:** NEEDS RESEARCH — importlib.util patterns for isolated module loading, CUDA initialization state management, rollback strategies if lazy imports fail.

### Phase 6: Extract Gaussian I/O & ZMQ Server
**Rationale:** Tightly coupled components (I/O, subprocess runner, ZMQ server) should be extracted together. Do this after MACE loading is stable.
**Delivers:** gaussian/ package (io, runner, parser, fchk, zmq_server), moved GaussianLogParser and FCHK parser into package.
**Addresses:** Architecture goal (decompose gm_main.py), modularizes the Gaussian integration.
**Uses:** Moves existing modules + extracts I/O from gm_main.py.
**Research flag:** Standard patterns — no additional research needed.

### Phase 7: Extract Workflow Orchestrator
**Rationale:** Do this last — workflow.py becomes a thin orchestrator calling modular components. By this point, all components are extracted and tested.
**Delivers:** workflow.py (phase sequencing), updated CLI to use workflow, gm_main.py deprecated or removed.
**Addresses:** Architecture goal (complete decomposition of gm_main.py).
**Avoids:** Pitfall #2 (breaking pipeline — workflow is last to change so pipeline stays functional throughout).
**Research flag:** No research needed — thin wrapper around extracted components.

### Phase 8: Package Structure & Reorganization
**Rationale:** Final cleanup phase — organize into mace_gaussian/ package, update all imports, reorganize analysis modules.
**Delivers:** mace_gaussian/ package layout, updated imports everywhere, analysis/ subpackage organized, CLI entry point aligned in pyproject.toml.
**Addresses:** Architecture goal (final package structure), stack recommendation (entry points).
**Uses:** Existing modular code, just reorganizing.
**Research flag:** No research needed — standard package layout.

### Phase 9: CI/CD & Distribution Prep
**Rationale:** After refactoring is complete, set up automation and distribution infrastructure.
**Delivers:** GitHub Actions CI (lint + unit tests), expanded ruff rules (B, SIM, PTH, RUF), pytest-cov in CI, install script for custom MACE packages.
**Addresses:** Stack recommendations (CI/CD), table stakes features (test coverage).
**Uses:** GitHub Actions (well-documented), ruff (already in use).
**Research flag:** Standard patterns — no additional research needed.

### Phase 10: Documentation
**Rationale:** Document after refactoring is stable. Tests provide verified examples for documentation.
**Delivers:** Updated README with installation + quickstart, worked example (water end-to-end), method documentation, MkDocs setup if desired, improved CLI help text.
**Addresses:** Table stakes features (documentation), enables distribution.
**Uses:** MkDocs + mkdocs-material (recommended), Google-style docstrings.
**Research flag:** No research needed — documentation of working code.

### Phase Ordering Rationale

- **Tests first (Phase 1)** protects everything else — non-negotiable based on pitfalls research.
- **Error handling (Phase 2)** can be done early without structural changes — makes refactoring safer.
- **Low-risk first (Phases 3-4)** builds confidence — utilities and calculator extraction have minimal coupling.
- **High-risk isolated (Phase 5)** — MACE monkey-patching is the most fragile component, do it after safety net is in place and with dedicated research.
- **Tightly coupled together (Phase 6)** — Gaussian I/O, runner, ZMQ server are interdependent, extract as a unit.
- **Orchestrator last (Phase 7)** — workflow touches everything, so do it after components are stable.
- **Structure after function (Phase 8)** — reorganize into package layout after modular structure is proven.
- **Automation last (Phases 9-10)** — CI/CD and docs document the refactored codebase, not the work-in-progress.

This ordering avoids pitfall #2 (breaking pipeline) by keeping the tool functional throughout, and addresses pitfall #1 (refactoring without tests) by front-loading testing infrastructure.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 5 (MACE loading):** Complex integration with custom packages, importlib.util patterns for isolated loading, CUDA device placement, rollback strategies. This is the highest-risk refactoring task.

Phases with standard patterns (skip research-phase):
- **Phases 1, 2, 3, 4, 6, 7, 8, 9, 10:** All use well-documented patterns (pytest with fixtures, error handling, package organization, CI with GitHub Actions, MkDocs). No phase-specific research needed.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | pytest + uv + pyproject.toml are industry standards for scientific Python. CI patterns verified with ASE, cclib examples. |
| Features | HIGH | Table stakes features identified from distributed package norms. Differentiators already implemented. Anti-features clearly scoped. |
| Architecture | HIGH | Refactoring pattern (extract from monolith incrementally) is well-established. Component boundaries are clear from existing code structure. |
| Pitfalls | HIGH | Pitfalls are grounded in existing fragility (CLAUDE.md documents monkey-patching issues) and standard scientific code risks (numerical reproducibility). |

**Overall confidence:** HIGH

Research is based on established patterns in scientific Python community (ASE, MDAnalysis, cclib, RDKit all follow similar testing/packaging approaches) and analysis of existing codebase fragility points documented in CLAUDE.md.

### Gaps to Address

- **MACE loading implementation details** — Phase 5 needs dedicated research on importlib.util patterns for isolated module loading and CUDA initialization state management. This is the only area where standard patterns don't directly apply due to custom package dependencies.

- **Test coverage targets** — Research recommends pytest-cov but doesn't specify coverage targets. During Phase 1, establish pragmatic coverage goals (focus on parsers, mode matching, calculators — not necessarily 100%).

- **Type checking strategy** — Research suggests keeping ty but switching to mypy if issues arise. Decision should be made during Phase 2-3 based on actual experience with ty during refactoring.

- **Acetic acid bug resolution** — Documented in CLAUDE.md (commit a4384c4) but root cause unclear. Should be addressed during Phase 1 (testing) with a dedicated test case, but may need parser-specific research if the fix is non-obvious.

## Sources

### Primary (HIGH confidence)
- **Existing codebase analysis** — gm_main.py structure, CLAUDE.md documented fragility (monkey-patching, acetic acid bug), pyproject.toml current setup
- **Scientific Python packaging norms** — ASE, MDAnalysis, cclib, RDKit all use pytest + fixtures + CI patterns documented in research
- **pytest documentation** — testing strategy with fixtures for external dependencies (Gaussian .log/.fchk files)
- **uv documentation** — modern Python packaging with pyproject.toml + lockfile
- **ruff documentation** — rule expansion for code quality

### Secondary (MEDIUM confidence)
- **MkDocs vs Sphinx trade-offs** — MkDocs is gaining ground but Sphinx still dominates scientific Python; either works for thesis scope
- **Type checker selection** — ty is alpha but may be sufficient; mypy is more battle-tested but slower

### Tertiary (LOW confidence)
- **MACE import isolation patterns** — will need Phase 5-specific research; existing research identifies problem but not solution details
- **Test coverage targets** — no specific research on what coverage percentage is appropriate for scientific code; community practice varies

---
*Research completed: 2026-02-16*
*Ready for roadmap: yes*
