# Phase 8: Package Structure & Reorganization - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Reorganize the flat collection of top-level `.py` files and existing subdirectories into a proper `mace_gaussian/` installable Python package. The goal is a clean importable structure suitable for eventual `pip install` distribution. CLI expansion (adding `analyze` subcommands) is deferred to a later phase. This phase is structural only — no new features.

</domain>

<decisions>
## Implementation Decisions

### Package strategy
- **All in (Path B):** All library and analysis code moves into `mace_gaussian/`. No core code stays at root.
- The project should be distributable (`pip install mace-gaussian`) even though Gaussian 16 is a required system dependency — this is standard in computational chemistry tools.

### Final package structure
```
mace_gaussian/
    __init__.py              # version + re-exports run_pipeline
    cli.py                   # CLI entry point
    workflow.py              # pipeline orchestrator
    gm_helper.py             # ZMQ bridge to Gaussian
    dft_baseline.py          # DFT baseline calculations
    diagnostics.py           # diagnostic utilities
    calculators/             # calculator classes (move existing subdirectory in)
        __init__.py
        base.py, espaloma.py, factory.py, mace_loader.py, mace_ml.py, xtb.py
    gaussian/                # Gaussian I/O & ZMQ server (move existing subdirectory in)
        __init__.py
        runner.py, ...
    analysis/                # analysis and reporting
        __init__.py
        analyze_spectra.py
        mode_matching.py
        html_report_generator.py
        run_analysis.py      # anharmonic analysis logic
        run_analysis_harmonic.py  # harmonic analysis logic
    utils/                   # utilities (move existing subdirectory in)
        __init__.py
        exceptions.py, results.py, units.py, validation.py
```

### Root-level shim scripts (keep, update internals)
- `run_analysis.py` stays at root as a thin shim calling `mace_gaussian.analysis.run_analysis`
- `run_analysis_harmonic.py` stays at root as a thin shim
- This preserves `python run_analysis.py water` and `python run_analysis_harmonic.py water` exactly as documented in CLAUDE.md — no muscle memory breakage

### Utility/research scripts
- `comparison_workflow.py`, `produce_molecules.py`, `convert_all_chk_files.py`, `charge_analysis.py` move to `scripts/` directory
- These are research utilities, not library code — `scripts/` is the standard home

### Public API (`mace_gaussian/__init__.py`)
- Minimal: expose `__version__` and re-export `run_pipeline` for convenience
- Everything else accessed via subpackage paths (e.g., `from mace_gaussian.workflow import run_pipeline`)
- No over-commitment to a public API at this stage

### Import style throughout package
- Relative imports within `mace_gaussian/` (e.g., `from .utils.exceptions import ...`)
- Absolute imports from outside the package
- `cli.py` uses absolute imports since it's the entry point

### pyproject.toml cleanup (full pass)
- Remove broken `gm-main` and `gm-diagnose` entry points from `[project.scripts]`
- Update `[project.scripts]` to point `mace-gaussian = "mace_gaussian.cli:app"` (or equivalent)
- Fix `[tool.ruff.lint.isort] known-first-party` to list `mace_gaussian` (remove stale `gm_main`, `gaussian_parser`, `fchk_parser` refs)
- Fix `[tool.coverage.run] source` to list `mace_gaussian`
- Update `[tool.ruff]` `src` if needed to reflect new package location
- Update `[build-system]` / `packages` config to find `mace_gaussian/`

### Claude's Discretion
- Exact `__init__.py` re-export surface beyond `run_pipeline` and `__version__`
- Whether `dft_baseline.py` gets a subpackage or stays flat in `mace_gaussian/`
- How to handle `fort.7` and other non-Python files at root (leave at root — not package files)
- How to handle `.xyz` molecule input files (leave at root — user-provided inputs)
- Whether to add `py.typed` marker for type checking support

</decisions>

<specifics>
## Specific Ideas

- Distribution motivation: even though Gaussian 16 is a closed-source system dependency, the Python package should be `pip install`-able (like cclib, ase, or psi4 wrappers). Thesis collaborators with Gaussian licenses should be able to `pip install mace-gaussian`.
- The `scripts/` directory already exists in the repo — the utility scripts move there rather than creating a new directory.
- `run_analysis.py` shims: should be literally 2-3 lines — import and call `main()` from the moved module. Not wrapper logic, just delegation.

</specifics>

<deferred>
## Deferred Ideas

- CLI `analyze` subcommand (`python cli.py analyze water`) — Phase 9 or a new phase between 8 and 9
- ORCA support as an alternative to Gaussian — future phase when ORCA path is confirmed
- Formal public API with stable versioning guarantees — only needed if/when distributed to a wider audience

</deferred>

---

*Phase: 08-package-structure-reorganization*
*Context gathered: 2026-02-24*
