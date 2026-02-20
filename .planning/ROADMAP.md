# Roadmap: MACE-Gaussian Refactoring & Distribution

## Overview

Transform MACE-Gaussian from a functional thesis tool into a distribution-ready scientific Python package. The journey starts with establishing a test suite to protect against regressions, adds defensive error handling, incrementally extracts components from the monolithic main module (with special care around MACE module loading), reorganizes into proper package structure, and culminates in CI/CD automation and comprehensive documentation. Each phase keeps the pipeline functional—no broken states—while building toward a package that other computational chemistry researchers can clone, install, and trust.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Testing Infrastructure & Characterization** - Establish pytest suite with fixtures and reference outputs
- [ ] **Phase 2: Error Handling & Input Validation** - Add defensive programming and prerequisite checking
- [ ] **Phase 3: Extract Utilities & Conventions** - Move pure functions and establish code conventions
- [ ] **Phase 4: Extract Calculator Classes** - Modularize calculator implementations into dedicated package
- [x] **Phase 5: Replace MACE Module Monkey-Patching** - Fix highest-risk fragility with safe loading pattern (completed 2026-02-19)
- [ ] **Phase 6: Extract Gaussian I/O & ZMQ Server** - Modularize Gaussian integration components
- [ ] **Phase 7: Extract Workflow Orchestrator** - Create thin coordinator for pipeline phases
- [ ] **Phase 8: Package Structure & Reorganization** - Reorganize into proper mace_gaussian/ package
- [ ] **Phase 9: CI/CD & Distribution Prep** - Automate testing and prepare for distribution
- [ ] **Phase 10: Documentation** - Complete user docs and distribution materials

## Phase Details

### Phase 1: Testing Infrastructure & Characterization
**Goal**: Establish comprehensive test suite that protects against regressions during refactoring
**Depends on**: Nothing (first phase)
**Requirements**: TEST-01, TEST-02, TEST-03, TEST-04, TEST-05, TEST-06, TEST-07, TEST-08, TEST-09, REPR-04
**Success Criteria** (what must be TRUE):
  1. User can run pytest and see passing tests for Gaussian parsers using committed fixtures
  2. Water and CH4 reference outputs are committed and can be regenerated to verify reproducibility
  3. Test markers allow separating unit tests (run in CI) from integration tests (require GPU/Gaussian)
  4. Acetic acid parser bug is captured in a test that documents expected vs actual behavior
  5. Coverage report shows which code paths are tested vs untested
**Plans**: 4 plans

Plans:
- [ ] 01-01-PLAN.md — Test infrastructure: pytest config, fixtures, conftest.py, .gitignore updates, uv.lock
- [ ] 01-02-PLAN.md — Gaussian log parser tests (harmonic, anharmonic, overtones, combination bands, acoh bug)
- [ ] 01-03-PLAN.md — FCHK parser and mode matching tests
- [ ] 01-04-PLAN.md — Reference output regression tests and coverage verification

### Phase 2: Error Handling & Input Validation
**Goal**: Add defensive programming that makes failures explicit and catches problems before multi-hour calculations
**Depends on**: Phase 1
**Requirements**: ERR-01, ERR-02, ERR-03, ERR-04, ERR-05, ERR-06, REPR-01, REPR-02, REPR-03
**Success Criteria** (what must be TRUE):
  1. Parser failures raise exceptions with context instead of returning empty data silently
  2. User gets clear error at startup if Gaussian, formchk, or dipole model are missing
  3. CUDA availability is checked and logged with warning if falling back to CPU
  4. Results JSON includes version metadata (tool version, Python version, model versions, calculation parameters)
  5. Optimization step count is tracked correctly and appears in results
**Plans**: 3 plans

Plans:
- [ ] 02-01-PLAN.md — Exception hierarchy and validation module with tests
- [ ] 02-02-PLAN.md — Parser error hardening, optimization step tracking, and results metadata
- [ ] 02-03-PLAN.md — CLI validation integration, subprocess timeout, and version fix

### Phase 3: Extract Utilities & Conventions
**Goal**: Extract pure functions from gm_main.py and establish documented code conventions
**Depends on**: Phase 2
**Requirements**: STRUCT-01
**Success Criteria** (what must be TRUE):
  1. Unit conversion functions exist in utils/units.py and are tested
  2. Input validation functions exist in utils/validation.py and are reused consistently
  3. ResultsManager is extracted to utils/results.py
  4. Code conventions are documented (naming, error handling, units) and followed
**Plans**: 2 plans

Plans:
- [ ] 03-01-PLAN.md — Create utils/ package with units.py, move exceptions.py and validation.py, update all imports
- [ ] 03-02-PLAN.md — Move ResultsManager to utils/results.py, add unit tests, document conventions

### Phase 4: Extract Calculator Classes
**Goal**: Move calculator implementations from gm_main.py into dedicated calculators/ package
**Depends on**: Phase 3
**Requirements**: STRUCT-02
**Success Criteria** (what must be TRUE):
  1. DipoleCalculator base class and implementations (espaloma, mace_ml, xtb) exist in calculators/ package
  2. Calculator factory pattern is cleanly separated and tested
  3. Full pipeline still runs on water molecule using all 4 model combinations
**Plans**: 2 plans

Plans:
- [ ] 04-01-PLAN.md — Extract calculator hierarchy into calculators/ package, update gm_main.py imports
- [ ] 04-02-PLAN.md — Unit tests for calculator interface, factory pattern, and config constants

### Phase 5: Replace MACE Module Monkey-Patching
**Goal**: Replace sys.modules manipulation with safe lazy import isolation pattern
**Depends on**: Phase 4
**Requirements**: STRUCT-03
**Success Criteria** (what must be TRUE):
  1. MACE standard and MACE dipole can be loaded independently without sys.modules cleanup
  2. calculators/mace_loader.py provides safe loading mechanism tested in isolation
  3. CUDA device placement is handled correctly across MACE variants
  4. Full pipeline runs on water without cleanup_mace_modules() calls
**Plans**: 2 plans

Plans:
- [ ] 05-01-PLAN.md — Create safe mace_loader.py with pickle_module remapping, update mace_ml.py import, add tests
- [ ] 05-02-PLAN.md — Delete mace_calculators.py, clean up all references (CLAUDE.md, pyproject.toml, test mocks)

### Phase 6: Extract Gaussian I/O & ZMQ Server
**Goal**: Modularize Gaussian integration into focused components in gaussian/ package
**Depends on**: Phase 5
**Requirements**: STRUCT-04, STRUCT-05, STRUCT-07
**Success Criteria** (what must be TRUE):
  1. Gaussian I/O functions (parse_gaussian_input, write outputs) are in gaussian/io.py
  2. GaussianLogParser and FCHK parser are in gaussian/parser.py and gaussian/fchk.py
  3. ZMQ server is a context manager in gaussian/zmq_server.py with proper socket cleanup
  4. Gaussian subprocess runner in gaussian/runner.py handles timeouts
  5. Full pipeline still runs with modular Gaussian integration
**Plans**: 5 plans

Plans:
- [ ] 06-01-PLAN.md — Add GaussianRunError/GaussianTimeoutError exceptions, create gaussian/io.py
- [ ] 06-02-PLAN.md — Create gaussian/parser.py and gaussian/fchk.py (move top-level parser files)
- [ ] 06-03-PLAN.md — Create gaussian/zmq_server.py with GaussianZMQServer class (LINGER=0 fix)
- [ ] 06-04-PLAN.md — Create gaussian/runner.py with run_gaussian_with_zmq (SIGKILL timeout)
- [ ] 06-05-PLAN.md — Wire everything: gaussian/__init__.py, update all callers, delete old files

### Phase 7: Extract Workflow Orchestrator
**Goal**: Create thin workflow.py that sequences pipeline phases by calling modular components
**Depends on**: Phase 6
**Requirements**: STRUCT-06
**Success Criteria** (what must be TRUE):
  1. workflow.py exists as the main entry point coordinating all pipeline phases
  2. gm_main.py is deprecated or removed
  3. CLI delegates to workflow.py functions
  4. Full pipeline runs with new orchestrator on water and CH4
**Plans**: TBD

Plans:
- [ ] 07-01: TBD

### Phase 8: Package Structure & Reorganization
**Goal**: Reorganize entire codebase into proper mace_gaussian/ package with clean imports
**Depends on**: Phase 7
**Requirements**: STRUCT-08, STRUCT-09, STRUCT-10
**Success Criteria** (what must be TRUE):
  1. Code is organized as mace_gaussian/ package with proper __init__.py files
  2. Analysis modules are in mace_gaussian/analysis/ subpackage
  3. CLI entry point in pyproject.toml aligns with package structure
  4. All imports are relative within package, absolute from outside
  5. Full pipeline runs from installed package location
**Plans**: TBD

Plans:
- [ ] 08-01: TBD
- [ ] 08-02: TBD

### Phase 9: CI/CD & Distribution Prep
**Goal**: Automate testing, linting, and prepare infrastructure for distribution to other researchers
**Depends on**: Phase 8
**Requirements**: CI-01, CI-02, CI-03, CI-04
**Success Criteria** (what must be TRUE):
  1. GitHub Actions workflow runs ruff check and ruff format --check on every push
  2. GitHub Actions runs pytest unit tests (no GPU/Gaussian) on every push
  3. Ruff rules include B, SIM, PTH, RUF for better code quality
  4. Install script or documented procedure exists for custom MACE packages
  5. pytest-cov reports coverage in CI for core modules
**Plans**: TBD

Plans:
- [ ] 09-01: TBD
- [ ] 09-02: TBD

### Phase 10: Documentation
**Goal**: Provide complete documentation for installation, usage, and method description
**Depends on**: Phase 9
**Requirements**: DOC-01, DOC-02, DOC-03, DOC-04, DOC-05
**Success Criteria** (what must be TRUE):
  1. README includes step-by-step installation (including custom MACE packages)
  2. Quickstart guide (clone -> install -> run water -> view results) is tested and works
  3. Worked example with expected output is committed to repo
  4. CLI help text is complete and accurate for all commands
  5. Method description suitable for thesis methods section is available
**Plans**: TBD

Plans:
- [ ] 10-01: TBD
- [ ] 10-02: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 -> 9 -> 10

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Testing Infrastructure & Characterization | 0/4 | Planning complete | - |
| 2. Error Handling & Input Validation | 0/3 | Planning complete | - |
| 3. Extract Utilities & Conventions | 0/2 | Planning complete | - |
| 4. Extract Calculator Classes | 0/2 | Planning complete | - |
| 5. Replace MACE Module Monkey-Patching | 0/2 | Complete    | 2026-02-19 |
| 6. Extract Gaussian I/O & ZMQ Server | 3/5 | In Progress|  |
| 7. Extract Workflow Orchestrator | 0/TBD | Not started | - |
| 8. Package Structure & Reorganization | 0/TBD | Not started | - |
| 9. CI/CD & Distribution Prep | 0/TBD | Not started | - |
| 10. Documentation | 0/TBD | Not started | - |
