# Features Research

**Research Date:** 2026-02-16
**Domain:** Distributable scientific Python packages for computational chemistry

## Table Stakes

Features researchers expect from any distributed scientific tool — without these, they won't trust or cite it.

### Testing & Validation

- **Unit tests for core logic** — parsers, mode matching, metric calculations must have tests with known-good reference data. Researchers need confidence that the code does what the paper claims.
- **Regression tests with reference outputs** — run a small molecule (water) end-to-end and compare outputs against committed reference data. Detects silent breakage.
- **Test data fixtures** — committed sample Gaussian .log and .fchk files so tests don't need Gaussian installed.
- **Complexity:** Medium — need to extract test fixtures from existing runs and write pytest tests.
- **Dependencies:** None — can be done first.

### Error Handling & Messages

- **Explicit failure over silent corruption** — if a parser returns empty data, raise an error with context (which file, what was expected, what was found). Never silently return empty results.
- **Prerequisite checking** — validate that Gaussian, formchk, CUDA, dipole model file all exist before starting a multi-hour calculation. Fail fast with clear messages.
- **Environment diagnostics** — `cli.py diagnose` already exists but should cover all failure modes.
- **Complexity:** Low-medium — mostly adding validation and better error messages to existing code.
- **Dependencies:** None.

### Documentation

- **README with installation and quickstart** — step-by-step from clone to first result. Include expected output.
- **Method documentation** — what the tool does scientifically (for the paper's methods section and for users to understand what they're running).
- **CLI help text** — every command and option documented via Click's built-in help.
- **Example workflow** — a complete worked example with a small molecule showing input → output.
- **Complexity:** Medium — writing, not coding. But critical for adoption.
- **Dependencies:** Tests should pass first (so documented examples are verified).

### CLI Experience

- **Clear progress indication** — which phase is running, what molecule, estimated progress.
- **Consistent exit codes** — 0 for success, non-zero for failure. Critical for HPC job scripts.
- **Verbose/quiet modes** — `-v` for debug output, `-q` for minimal output.
- **Complexity:** Low — Click supports all of this natively.
- **Dependencies:** None.

### Reproducibility

- **Pinned dependency versions** — uv.lock already exists, good. Ensure it's committed.
- **Random seed management** — if any ML inference has stochastic elements, seeds must be settable.
- **Version in output** — results.json should include tool version, Python version, MACE model versions.
- **Configuration capture** — save the exact parameters used for each run alongside results.
- **Complexity:** Low — mostly adding metadata to existing output.
- **Dependencies:** None.

### Configuration

- **Sensible defaults** — tool should work out of the box for the common case (water, B3LYP/6-31G(d,p), mace_omol/mace_ml).
- **Config file support** — for HPC environments where CLI args are awkward in job scripts.
- **Environment variable documentation** — MACE_DIPOLE_MODEL_PATH, MACE_HELPER_SCRIPT_PATH, etc.
- **Complexity:** Low — ConfigArgParse already in dependencies.
- **Dependencies:** None.

## Differentiators

Features that make this tool stand out vs doing the ML+Gaussian workflow manually:

- **Automated pipeline** — single command from .xyz to full spectral analysis with plots and reports. Manual workflow requires 10+ separate steps.
- **Mode matching** — eigenvector overlap matching between DFT and ML modes is novel and not available in standard tools. This is the thesis contribution.
- **Multi-model comparison** — run 4 ML combinations against DFT in one workflow. Manual comparison would be extremely tedious.
- **HTML reports** — interactive, self-contained reports that can be shared and viewed in a browser. Better than raw data files.
- **Harmonic + anharmonic support** — most ML spectroscopy tools only do harmonic. Anharmonic via Gaussian external interface is unique.

## Anti-Features

Things to deliberately NOT build for a thesis-scope project:

- **Web interface / GUI** — CLI is appropriate for the audience. Web adds massive complexity for no research value.
- **Database backend** — JSON files in directories are fine for thesis-scale data (<100 molecules).
- **Multi-user support** — single-researcher tool. No auth, no permissions, no collaboration features.
- **Automatic model training** — use pre-trained MACE models only. Training is a separate research problem.
- **Cloud deployment** — Gaussian 16 licensing makes cloud deployment impractical anyway.
- **Plugin system** — hardcoded calculator backends are fine. Extensibility can come later if the tool is adopted.
- **Real-time visualization** — static plots in reports are sufficient. No need for interactive dashboards.
- **Support for non-Gaussian QM packages** — Psi4, ORCA, etc. are future work if the approach proves valuable.

## Feature Dependencies

```
Testing & Validation ──→ Documentation (examples need verified outputs)
Error Handling ──→ CLI Experience (errors must surface cleanly)
Reproducibility ──→ Documentation (reproducibility claims need docs)
Configuration ──→ CLI Experience (config and CLI must be consistent)
```

Testing has no dependencies and should be done first — it protects all subsequent changes.

## Recommendations

**Priority order for refactoring:**

1. **Testing infrastructure** — add fixtures, unit tests for parsers and mode matching. This protects everything else.
2. **Error handling** — fail fast with clear messages. Replace silent empty returns.
3. **Input validation** — check prerequisites before long calculations.
4. **Code modularity** — split gm_main.py into focused modules (refactoring protected by tests).
5. **Reproducibility metadata** — version and config capture in output.
6. **Documentation** — README, quickstart, worked example.
7. **CLI polish** — progress, verbose/quiet, exit codes.
8. **Configuration** — config file support, defaults documentation.

---
*Features research: 2026-02-16*
