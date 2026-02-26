# Requirements: MACE-Gaussian Refactoring

**Defined:** 2026-02-16
**Core Value:** Reliable, reproducible IR spectral predictions using ML potentials that can be validated against DFT reference data.

## v1 Requirements

Requirements for the refactoring milestone. Each maps to roadmap phases.

### Testing

- [ ] **TEST-01**: Parser unit tests exist for harmonic frequency extraction with committed .log fixtures
- [ ] **TEST-02**: Parser unit tests exist for anharmonic frequency, overtone, and combination band extraction
- [ ] **TEST-03**: FCHK parser unit tests exist with committed .fchk fixtures
- [ ] **TEST-04**: Mode matching unit tests verify eigenvector overlap and mode pairing logic
- [ ] **TEST-05**: Reference outputs for water molecule committed as regression baseline
- [ ] **TEST-06**: Reference outputs for CH4 molecule committed as regression baseline
- [ ] **TEST-07**: Test markers (@pytest.mark.gpu, @pytest.mark.gaussian) separate unit from integration tests
- [ ] **TEST-08**: Acetic acid parser edge case has a dedicated test (reproduces known bug)
- [ ] **TEST-09**: pytest-cov configured and reports coverage for core modules

### Error Handling

- [ ] **ERR-01**: Gaussian parser raises explicit error with context when expected section is missing (not empty return)
- [ ] **ERR-02**: Prerequisite validation checks Gaussian (g16), formchk, and dipole model file exist before starting calculations
- [ ] **ERR-03**: CUDA availability checked at startup with clear warning on CPU fallback
- [ ] **ERR-04**: Environment variables (MACE_DIPOLE_MODEL_PATH, MACE_HELPER_SCRIPT_PATH) validated at startup
- [ ] **ERR-05**: Input .xyz file validated for basic sanity (exists, parseable, reasonable atom count)
- [x] **ERR-06**: Subprocess timeout added to Gaussian process to prevent indefinite hangs

### Code Structure

- [ ] **STRUCT-01**: Utility functions (unit conversions, validation) extracted from gm_main.py into separate modules
- [ ] **STRUCT-02**: DipoleCalculator classes extracted from gm_main.py into calculators/ package
- [x] **STRUCT-03**: MACE module monkey-patching replaced with safe loading pattern (lazy imports or process isolation)
- [x] **STRUCT-04**: Gaussian I/O functions extracted from gm_main.py into gaussian/ package
- [x] **STRUCT-05**: ZMQ server context manager extracted into dedicated module
- [x] **STRUCT-06**: Workflow orchestration extracted into workflow.py as thin coordinator
- [x] **STRUCT-07**: ZMQ socket cleanup race condition fixed (proper LINGER settings)
- [x] **STRUCT-08**: Project reorganized into mace_gaussian/ package with proper __init__.py
- [x] **STRUCT-09**: CLI entry point aligned in pyproject.toml with package structure
- [x] **STRUCT-10**: Analysis modules reorganized into analysis/ subpackage

### Reproducibility

- [ ] **REPR-01**: Results JSON includes tool version, Python version, and MACE model version
- [ ] **REPR-02**: Calculation parameters captured alongside results
- [ ] **REPR-03**: Optimization step count tracked properly (not hardcoded to 0)
- [ ] **REPR-04**: uv.lock committed and documented as reproducibility mechanism

### CI/CD

- [x] **CI-01**: GitHub Actions workflow runs ruff check and ruff format --check on push
- [x] **CI-02**: GitHub Actions workflow runs pytest (unit tests only, no GPU/Gaussian) on push
- [x] **CI-03**: Ruff rules expanded to include B, SIM, PTH, RUF
- [x] **CI-04**: Install script or documented procedure for custom MACE packages

### Documentation

- [ ] **DOC-01**: README updated with complete installation steps (including custom MACE packages)
- [ ] **DOC-02**: Quickstart guide: clone → install → run water → view results
- [ ] **DOC-03**: Worked example with expected output committed to repo
- [ ] **DOC-04**: CLI help text complete for all commands and options
- [ ] **DOC-05**: Method description suitable for citing in thesis methods section

## v2 Requirements

Deferred to future milestone. Tracked but not in current roadmap.

### Distribution

- **DIST-01**: pip-installable package on PyPI or conda-forge
- **DIST-02**: Docker/Singularity container for reproducible environment
- **DIST-03**: MkDocs documentation site hosted on GitHub Pages

### Features

- **FEAT-01**: Config file support for HPC job scripts (YAML/TOML)
- **FEAT-02**: Automatic retry logic for transient Gaussian failures
- **FEAT-03**: Checkpoint/restart mechanism for interrupted calculations
- **FEAT-04**: Verbose/quiet CLI modes (-v/-q flags)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Web interface / GUI | CLI is appropriate for research audience |
| Database backend | JSON files sufficient for thesis-scale data |
| Support for non-Gaussian QM packages | Single-target simplifies validation |
| New ML model backends | Focus on cleaning existing 4 combinations |
| Automatic model training | Separate research problem |
| Cloud deployment | Gaussian licensing prevents this |
| Plugin/extension system | Premature abstraction for thesis scope |
| Performance optimization (>100 atoms) | Thesis molecules are smaller |
| Real-time visualization | Static plots in reports are sufficient |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| TEST-01 | Phase 1 | Pending |
| TEST-02 | Phase 1 | Pending |
| TEST-03 | Phase 1 | Pending |
| TEST-04 | Phase 1 | Pending |
| TEST-05 | Phase 1 | Pending |
| TEST-06 | Phase 1 | Pending |
| TEST-07 | Phase 1 | Pending |
| TEST-08 | Phase 1 | Pending |
| TEST-09 | Phase 1 | Pending |
| ERR-01 | Phase 2 | Pending |
| ERR-02 | Phase 2 | Pending |
| ERR-03 | Phase 2 | Pending |
| ERR-04 | Phase 2 | Pending |
| ERR-05 | Phase 2 | Pending |
| ERR-06 | Phase 2 | Complete |
| STRUCT-01 | Phase 3 | Pending |
| STRUCT-02 | Phase 4 | Pending |
| STRUCT-03 | Phase 5 | Complete |
| STRUCT-04 | Phase 6 | Complete |
| STRUCT-05 | Phase 6 | Complete |
| STRUCT-06 | Phase 7 | Complete |
| STRUCT-07 | Phase 6 | Complete |
| STRUCT-08 | Phase 8 | Complete |
| STRUCT-09 | Phase 8 | Complete |
| STRUCT-10 | Phase 8 | Complete |
| REPR-01 | Phase 2 | Pending |
| REPR-02 | Phase 2 | Pending |
| REPR-03 | Phase 2 | Pending |
| REPR-04 | Phase 1 | Pending |
| CI-01 | Phase 9 | Complete |
| CI-02 | Phase 9 | Complete |
| CI-03 | Phase 9 | Complete |
| CI-04 | Phase 9 | Complete |
| DOC-01 | Phase 10 | Pending |
| DOC-02 | Phase 10 | Pending |
| DOC-03 | Phase 10 | Pending |
| DOC-04 | Phase 10 | Pending |
| DOC-05 | Phase 10 | Pending |

**Coverage:**
- v1 requirements: 38 total
- Mapped to phases: 38
- Unmapped: 0

---
*Requirements defined: 2026-02-16*
*Last updated: 2026-02-16 after roadmap creation (traceability confirmed)*
