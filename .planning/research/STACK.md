# Stack Research

**Research Date:** 2026-02-16
**Domain:** Tooling stack for distributable scientific Python packages

## Current State

The project already has a reasonable foundation:
- **Package manager:** uv with uv.lock (good — modern, fast, reproducible)
- **Build backend:** hatchling (good — simple, standard)
- **Project metadata:** pyproject.toml (good — PEP 621 compliant)
- **Linting:** ruff with E, F, W, I, UP rules at line-length 100
- **Testing:** pytest (minimal — few or no tests currently)
- **Type checking:** ty (alpha — should evaluate stability)
- **CLI:** Click (good — standard for scientific tools)

## Testing

### Framework: pytest (keep)
Already in dev dependencies. The standard for scientific Python.

### Key additions needed:

**pytest-cov** — Coverage reporting. Essential to know which code paths are tested.
```
pip install pytest-cov
pytest --cov=. --cov-report=html
```

**Test fixtures from real data** — Commit sanitized Gaussian .log and .fchk file snippets as test data. This is the standard pattern for scientific packages with external dependencies: test against committed reference data, not against live external tools.

**Mocking strategy for external deps:**
- `unittest.mock` (stdlib) for mocking Gaussian subprocess calls
- `pytest.mark.skipif` for tests requiring CUDA or Gaussian
- `pytest-mock` for cleaner mock syntax
- Separate test markers: `@pytest.mark.gpu`, `@pytest.mark.gaussian` for integration tests that need the real tools

**Test organization:**
```
tests/
  conftest.py          # shared fixtures
  test_data/           # committed reference .log, .fchk snippets
  test_gaussian_parser.py
  test_fchk_parser.py
  test_mode_matching.py
  test_mace_calculators.py
  test_analyze_spectra.py
  test_integration.py  # end-to-end, marked as slow/requires-gaussian
```

**Confidence:** HIGH — pytest + committed test data is the universal pattern for scientific Python packages (ASE, MDAnalysis, cclib all use this approach).

## Packaging & Distribution

### Current setup (keep, with improvements):
- **uv + pyproject.toml + hatchling** — this is the modern standard. No changes needed to the build system.
- **uv.lock** — must be committed for reproducibility.

### Improvements needed:

**Version management** — currently hardcoded `version = "0.2.0"` in pyproject.toml. Fine for thesis scope. If distributing later, consider `hatch-vcs` for git-tag-based versioning.

**Entry points** — currently defines `gm-main` and `gm-diagnose` scripts. The CLI entry point via `cli.py` uses Click but isn't registered as a console script. Should align:
```toml
[project.scripts]
mace-gaussian = "cli:cli"
```

**Custom MACE packages** — `mace_torch` and `mace_dipole_core` are local packages. For distribution:
- Document exact installation steps in README
- Consider adding an install script or Makefile target
- Long-term: publish to PyPI or provide as git submodules

**Confidence:** HIGH — pyproject.toml + uv is the established standard.

## Documentation

### Recommended: MkDocs with mkdocs-material

**Why MkDocs over Sphinx:**
- Markdown-based (easier for non-programmers to maintain)
- mkdocs-material theme is polished and widely used
- Simpler configuration than Sphinx
- GitHub Pages deployment is trivial

**For thesis scope, minimum docs:**
1. README.md with installation + quickstart (already exists, needs updating)
2. docs/ARCHITECTURE.md (already exists)
3. docs/DEVELOPMENT.md (already exists)
4. A worked example (water molecule, end-to-end)

**Docstrings:** Use Google-style docstrings. ruff can enforce with `D` rules (pydocstyle).

**Confidence:** MEDIUM — MkDocs is gaining ground but Sphinx still dominates in scientific Python. Either works. MkDocs is easier for a non-programmer.

## CI/CD

### GitHub Actions (recommended if using GitHub)

**What CAN run in CI (no Gaussian/GPU):**
- Linting: `ruff check && ruff format --check`
- Type checking: `ty check` (or mypy)
- Unit tests: `pytest -m "not gpu and not gaussian"` — parser tests, mode matching with fixtures, metric calculations
- Build check: `uv build`

**What CANNOT run in CI:**
- Integration tests requiring Gaussian 16 (licensed software)
- GPU tests requiring CUDA
- Full end-to-end workflow

**Minimal CI workflow:**
```yaml
# .github/workflows/ci.yml
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: astral-sh/setup-uv@v4
      - run: uv sync --dev
      - run: uv run ruff check
      - run: uv run ruff format --check
      - run: uv run pytest -m "not gpu and not gaussian" --cov
```

**Confidence:** HIGH — this is exactly what ASE, cclib, and similar packages do.

## Code Quality

### Keep: ruff (expand rules)

Current rules: `["E", "F", "W", "I", "UP"]`

**Recommended additions:**
- `"B"` — bugbear (common Python pitfalls)
- `"SIM"` — simplification suggestions
- `"PTH"` — prefer pathlib over os.path (consistent with existing code)
- `"RUF"` — ruff-specific rules

Don't add: `"D"` (docstrings) until documentation phase — too noisy for refactoring.

### Type checking: evaluate ty vs mypy

`ty` is currently in alpha (0.0.1a1). Options:
- **Keep ty** if it works well enough for the project — it's fast and developing rapidly
- **mypy** is the battle-tested standard but slower and more verbose
- **pyright** is faster than mypy but more strict

For thesis scope: keep ty if it's working, switch to mypy if ty causes issues. Don't invest heavily in type annotations — focus on the critical paths (parser return types, calculator interfaces).

**Confidence:** HIGH for ruff expansion, MEDIUM for type checker choice.

## Recommendations Summary

### Adopt now (high priority):
1. **pytest-cov** + test fixtures from real Gaussian output
2. **Expand ruff rules** to include B, SIM, PTH, RUF
3. **Test markers** for gpu/gaussian-dependent tests
4. **Align CLI entry points** in pyproject.toml

### Adopt during refactoring:
5. **Google-style docstrings** on public functions (as you touch them)
6. **Type hints** on function signatures (as you touch them)

### Adopt for distribution:
7. **GitHub Actions CI** with lint + unit test jobs
8. **MkDocs** for user-facing documentation
9. **Install script** for custom MACE packages

### Don't adopt (thesis scope):
- Full Sphinx documentation build
- 100% test coverage target
- Pre-commit hooks (ruff is fast enough to run manually)
- Containerization (Docker/Singularity) — users have their own HPC setups

---
*Stack research: 2026-02-16*
