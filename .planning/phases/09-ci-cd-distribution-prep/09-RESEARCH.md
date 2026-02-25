# Phase 9: CI/CD & Distribution Prep - Research

**Researched:** 2026-02-25
**Domain:** GitHub Actions CI, Ruff linting, pytest-cov, Python install procedures
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**CI triggers**
- Run on every push to **any branch** (not just main)
- Single workflow file is sufficient — no separate PR workflow needed

**Test scope in CI**
- Run only tests that do **not** require custom local packages (`mace_ML_pkg`, `mace_dipole_pkg`),
  GPU, or Gaussian 16
- The existing markers already cover this: skip `gpu`, `gaussian`, and `slow` tests
- Unit tests for parsers, mode matching, utils, and validation run fine without the ML packages
- Do NOT attempt to install mace_ML_pkg or mace_dipole_pkg in CI — they are not on PyPI
  and the CI install would be fragile

**Ruff rule expansion (CI-03)**
- Add rules B, SIM, PTH, RUF to the existing ruff config
- **Fix all existing violations now** with `ruff check --fix` before CI enforces them
- CI should then enforce: no new violations allowed (non-zero exit = fail)
- Both `ruff check` and `ruff format --check` run in CI

**Coverage reporting**
- Add `pytest-cov` to CI run — report coverage for `mace_gaussian/` module
- **Do not set a failure threshold** — coverage is informational only
- Output as terminal summary (not uploaded to external service)

**Distribution / install procedure (CI-04)**
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

### Deferred Ideas (OUT OF SCOPE)
- None came up — discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| CI-01 | GitHub Actions workflow runs ruff check and ruff format --check on push | GitHub Actions YAML patterns, ruff CI flags documented below |
| CI-02 | GitHub Actions workflow runs pytest unit tests (no GPU/Gaussian) on push | pytest -m "not gpu and not gaussian and not slow" confirmed working; pip install strategy documented |
| CI-03 | Ruff rules expanded to include B, SIM, PTH, RUF | 65 existing violations inventoried; extend-select pattern documented; all are fixable before enforcement |
| CI-04 | Install script or documented procedure for custom MACE packages | install_mace_packages.sh already exists; README needs updated step-by-step section |
</phase_requirements>

## Summary

Phase 9 is a pure infrastructure phase: create one GitHub Actions YAML, expand ruff rules, fix existing violations, add pytest-cov to CI, and update the README installation section. No new Python code is written — only configuration files and documentation.

The critical CI constraint is that `uv sync --locked` **cannot** be used in CI because `dgl==2.2.1` (a transitive dependency of `espaloma_charge`) has Windows-only wheels pinned in `uv.lock`. The correct approach is `pip install -e . --no-deps` to install the package without pulling in the heavy ML dependencies, then `pip install pytest pytest-cov ruff`. This is sufficient because all tests that run in CI (parsers, mode matching, utils, validation) use only lazy-imported heavy deps, which gracefully fail to load at runtime without affecting test collection or execution.

There are exactly 65 ruff violations under the new B/SIM/PTH/RUF rules. 8 are auto-fixable with `ruff check --fix`; the remaining 57 require manual edits (mostly replacing `open()` with `Path.open()` and `os.path.*` with `Path.*` methods). These must be fixed before CI is activated or CI will fail on first push.

**Primary recommendation:** Use `pip install -e . --no-deps && pip install pytest pytest-cov ruff` for CI, run pytest with `-m "not gpu and not gaussian and not slow"`, and fix all 65 ruff violations before pushing the workflow file.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| GitHub Actions | N/A | CI/CD automation | Native GitHub integration, free for public repos |
| `astral-sh/setup-uv` | v7 | Install uv in CI | Official action from uv maintainers, handles cache |
| `actions/checkout` | v4 | Check out repo | Standard GitHub Actions checkout |
| `actions/setup-python` | v5 | Set Python version | Official GitHub action for Python setup |
| ruff | 0.15.1 (pinned in uv.lock) | Lint + format check | Already in dev deps; single tool replaces flake8+black |
| pytest | 7.1.3 (pinned) | Test runner | Already configured in pyproject.toml |
| pytest-cov | 7.0.0 (pinned) | Coverage reporting | Already in dev deps; integrates with pytest |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `actions/cache` | v4 | Pip dependency caching | Reduces CI time by caching pip downloads |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `pip install --no-deps` | `uv sync --locked` | uv sync fails because dgl Windows-only wheel is in lock file; pip --no-deps is correct here |
| `pip install -e . --no-deps` | `pip install pytest pytest-cov ruff` only | The package must be installed so `mace_gaussian.*` imports resolve in tests |
| single combined job | separate lint + test jobs | Combined is faster and simpler for a single-dev research project |

**Installation (CI workflow — not local):**
```bash
# CI install strategy (NOT uv sync -- see pitfall below)
pip install -e . --no-deps
pip install pytest pytest-cov ruff
```

## Architecture Patterns

### Recommended Workflow File Location
```
.github/
└── workflows/
    └── ci.yml        # Single workflow for all CI checks
```

### Pattern 1: Single Workflow File, Two Jobs
**What:** One YAML file with a `lint` job and a `test` job running in parallel.
**When to use:** For a single-developer research project where simplicity matters over parallelism.

```yaml
# Source: https://docs.astral.sh/uv/guides/integration/github/
# Source: https://docs.github.com/en/actions/tutorials/build-and-test-code/python
name: CI

on:
  push:

jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"

      - name: Install ruff
        run: pip install ruff==0.15.1

      - name: Check formatting
        run: ruff format --check mace_gaussian/ tests/

      - name: Lint
        run: ruff check mace_gaussian/ tests/

  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"
          cache: "pip"

      - name: Install package (no heavy deps)
        run: pip install -e . --no-deps

      - name: Install test dependencies
        run: pip install pytest pytest-cov ruff

      - name: Run tests
        run: |
          pytest tests/ \
            -m "not gpu and not gaussian and not slow" \
            --cov=mace_gaussian \
            --cov-report=term-missing
```

**Key decisions baked in:**
- `on: push` with no branch filter = runs on every push to any branch (locked decision)
- `python-version: "3.10"` = matches production environment (Claude's discretion)
- `pip install -e . --no-deps` = avoids dgl Windows-only wheel problem (see Pitfall 1)
- `--cov-report=term-missing` = shows uncovered line numbers in terminal (Claude's discretion)
- No failure threshold = `--cov-fail-under` is absent (locked decision)

### Pattern 2: Ruff Rule Configuration in pyproject.toml
**What:** Use `extend-select` to add new rules alongside existing ones.
**When to use:** When a base `select` is already defined (current situation).

```toml
# Source: https://docs.astral.sh/ruff/linter/
[tool.ruff.lint]
select = ["E", "F", "W", "I", "UP", "B", "SIM", "PTH", "RUF"]
```

Note: `extend-select` is also valid but `select` with full list makes the complete rule set explicit and is preferred by the ruff docs.

### Pattern 3: Coverage Configuration in pyproject.toml
**What:** `[tool.coverage.run]` and `[tool.coverage.report]` sections already exist.
**When to use:** Always configure coverage source and omit patterns in pyproject.toml, not via CLI flags.

```toml
# Already configured in pyproject.toml — no changes needed
[tool.coverage.run]
source = ["mace_gaussian"]
omit = ["tests/*", "mace_dipole_pkg/*", "mace_ML_pkg/*"]

[tool.coverage.report]
show_missing = true
skip_empty = true
```

### Anti-Patterns to Avoid
- **`uv sync --locked` in CI**: Fails because `dgl==2.2.1` has Windows-only wheels in `uv.lock`. Never use `uv sync` in CI for this project.
- **Installing mace_ML_pkg or mace_dipole_pkg in CI**: These are not on PyPI and require large GPU-specific binaries. Explicitly forbidden in locked decisions.
- **`--cov-fail-under` threshold**: Explicitly forbidden in locked decisions — coverage is informational only.
- **Branch-filtered triggers**: `on: push` with no branches filter is required (every branch, locked decision).
- **Uploading coverage to external service**: Terminal output only (locked decision).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Caching pip downloads | Custom cache logic | `actions/setup-python@v5` with `cache: "pip"` | Built-in to setup-python action |
| Python version matrix | Multiple workflow files | `strategy.matrix.python-version` (if needed) | GitHub Actions native feature |
| Coverage threshold enforcement | Custom script | `--cov-fail-under=N` flag (NOT USED — locked decision) | Locked out of scope |
| Ruff auto-fix in CI | Custom pre-commit | `ruff check --fix` locally before push | CI enforces, dev fixes locally |

**Key insight:** This phase is configuration-only. Every tool already exists in `dev` dependencies. The challenge is the correct install strategy in CI, not building anything new.

## Common Pitfalls

### Pitfall 1: `uv sync --locked` Fails Due to dgl Windows-Only Wheel
**What goes wrong:** `uv sync --locked` exits with error: `Distribution dgl==2.2.1 can't be installed because it doesn't have a source distribution or wheel for the current platform (linux)`.
**Why it happens:** `espaloma_charge` depends on `dgl`, and `dgl==2.2.1` only has Windows wheels. The lock file was created on a Windows machine (or the wheel metadata was resolved for Windows). CI runs on `ubuntu-latest`.
**How to avoid:** Use `pip install -e . --no-deps` instead. This installs the `mace_gaussian` package for importability without pulling any transitive dependencies.
**Warning signs:** Any attempt to `pip install espaloma-charge` or `uv sync` will fail on Linux CI.

### Pitfall 2: 65 Ruff Violations Must Be Fixed BEFORE Pushing the Workflow
**What goes wrong:** If `.github/workflows/ci.yml` is pushed before ruff violations are fixed, CI fails immediately on first push.
**Why it happens:** The 65 violations under B/SIM/PTH/RUF rules all exist now. CI enforces "no violations" on push.
**How to avoid:** Task ordering is critical — fix ruff violations (Task 1) before creating the workflow (Task 2). Use `ruff check --fix` for the 8 auto-fixable ones; manually fix the remaining 57.
**Warning signs:** `ruff check --select B,SIM,PTH,RUF mace_gaussian/ tests/` shows 65 errors on the current codebase.

### Pitfall 3: test_path_suffix Test Fails — Stale Test
**What goes wrong:** `tests/test_calculators.py::TestDefaultMaceDipoleModel::test_path_suffix` fails because it asserts the default model path ends with `"dipole_model/model_1.model"` but the actual default is `mace4ir_models/pretrained_models/model_1_dipole.model`.
**Why it happens:** The `DEFAULT_MACE_DIPOLE_MODEL` constant in `mace_gaussian/calculators/mace_ml.py` was updated to point to the new MACE4IR model, but the test was not updated.
**How to avoid:** Fix this test before CI runs. Update the assertion to match the actual default path suffix.
**Warning signs:** Currently 1 failing test in the test suite: `FAILED tests/test_calculators.py::TestDefaultMaceDipoleModel::test_path_suffix`.

### Pitfall 4: Heavy Dep Import Warnings at Test Collection Time
**What goes wrong:** When pytest collects tests, `mace_gaussian/__init__.py` triggers lazy calculator checks that emit `✗ Espaloma-charge dipole calculator failed: No module named 'dgl'` to stdout/stderr.
**Why it happens:** The calculators package `__init__.py` imports `EspalomaDipoleCalculator` which calls `_check_availability()` at module load. Without `dgl` installed, this logs warnings.
**How to avoid:** This is cosmetic only — tests still collect and run correctly. The warnings are not test failures. However they may clutter CI output. Acceptable per locked decision (do not install heavy deps).
**Warning signs:** `✗ Espaloma-charge dipole calculator failed` lines in CI output — not a failure, just noise.

### Pitfall 5: README Install Section Already Partially Exists
**What goes wrong:** The README already has an installation section that references `install_mace_packages.sh`. The locked decision says "documented steps in README are sufficient — no bash install script needed". But the script already exists and is referenced.
**Why it happens:** The script `install_mace_packages.sh` was created earlier. The CONTEXT.md decision is about what the README documents, not whether the script exists.
**How to avoid:** Keep the script (it's useful), but rewrite the README section to show the explicit manual steps. The script can be mentioned as an alternative. The README must cover all 5 steps explicitly (clone, conda env, `pip install -e mace_ML_pkg`, `pip install -e mace_dipole_pkg`, `pip install -e .`, diagnose).
**Warning signs:** README install section currently says "run `install_mace_packages.sh`" without showing the individual pip commands.

## Code Examples

Verified patterns from official sources:

### Complete GitHub Actions Workflow (annotated)
```yaml
# Source: https://docs.astral.sh/uv/guides/integration/github/
# Source: https://docs.github.com/en/actions/tutorials/build-and-test-code/python
name: CI

on:
  push:            # Every branch, every push (no branch filter = locked decision)

jobs:
  lint:
    name: Lint
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"

      - name: Install ruff
        run: pip install "ruff==0.15.1"

      - name: Check code formatting
        run: ruff format --check mace_gaussian/ tests/

      - name: Lint code
        run: ruff check mace_gaussian/ tests/

  test:
    name: Test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.10"
          cache: "pip"            # Caches pip downloads (Claude's discretion: yes, cache)

      - name: Install package without heavy deps
        run: pip install -e . --no-deps

      - name: Install test tools
        run: pip install "pytest>=7.0.0" "pytest-cov>=4.0" "numpy>=1.20.0" "pyzmq>=22.0.0" "ase>=3.22.0"

      - name: Run unit tests
        run: |
          pytest tests/ \
            -m "not gpu and not gaussian and not slow" \
            --cov=mace_gaussian \
            --cov-report=term-missing
```

Note on deps: `pip install -e . --no-deps` skips all declared dependencies. Tests for parsers, mode matching, units, and validation only need `numpy`, `ase`, `pyzmq`. Install these explicitly.

### Ruff pyproject.toml Configuration (expanded)
```toml
# Source: https://docs.astral.sh/ruff/linter/
[tool.ruff.lint]
select = ["E", "F", "W", "I", "UP", "B", "SIM", "PTH", "RUF"]

[tool.ruff.lint.isort]
known-first-party = ["mace_gaussian"]
```

Change from current `select = ["E", "F", "W", "I", "UP"]` to the above.

### Ruff Violation Fix Strategy
```bash
# Step 1: Auto-fix the 8 fixable violations
ruff check --select B,SIM,PTH,RUF --fix mace_gaussian/

# Step 2: Check what remains
ruff check --select B,SIM,PTH,RUF mace_gaussian/
# Remaining violations by type:
# - PTH123 (open() → Path.open()): ~20 occurrences in analysis/, gaussian/, cli.py, utils/
# - PTH110 (os.path.exists() → Path.exists()): ~8 occurrences in workflow.py, utils/results.py
# - PTH100 (os.path.abspath() → Path.resolve()): ~5 occurrences
# - PTH107 (os.remove() → Path.unlink()): 2 occurrences in zmq_server.py
# - SIM108 (if/else → ternary): 4 occurrences in analyze_spectra.py
# - SIM105 (try/except/pass → contextlib.suppress): 1 occurrence in zmq_server.py
# - SIM117 (nested with → single with): 1 occurrence in gm_helper.py
# - B904 (raise from err): 3 occurrences in gaussian/fchk.py
# - RUF013 (implicit Optional): 5 occurrences in analyze_spectra.py
# - RUF059 (unused unpacked var): 1 occurrence in mode_matching.py
# - B007 (unused loop variable): 2 occurrences
# - RUF002 (ambiguous char in docstring): 1 occurrence
```

### pytest-cov Command (no threshold)
```bash
# Source: https://pytest-cov.readthedocs.io/en/latest/reporting.html
pytest tests/ \
  -m "not gpu and not gaussian and not slow" \
  --cov=mace_gaussian \
  --cov-report=term-missing
# NO --cov-fail-under flag (locked decision: coverage is informational only)
```

### Fix for test_path_suffix
```python
# tests/test_calculators.py — update the stale assertion
def test_path_suffix(self):
    """Default model path ends with the expected suffix."""
    from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL

    # Updated: points to MACE4IR model (old path was dipole_model/model_1.model)
    assert DEFAULT_MACE_DIPOLE_MODEL.endswith("mace4ir_models/pretrained_models/model_1_dipole.model")
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| flake8 + black + isort separately | ruff as unified tool | 2023 | Single tool, 10-100x faster, same rules |
| `setup-python` with pip requirements | `astral-sh/setup-uv` with uv.lock | 2024 | Reproducible, faster installs |
| codecov.io upload for coverage | Terminal output only | N/A (project choice) | Simpler, no external service |
| `requirements.txt` for CI deps | `pyproject.toml` dependency groups | 2024 | Single source of truth |

**Deprecated/outdated:**
- `actions/setup-python` with `cache: "pip"` is still valid but `astral-sh/setup-uv` is preferred for uv projects. However for this project, `setup-python` + direct pip is correct because uv sync fails (dgl issue).

## Open Questions

1. **What Python version should CI use?**
   - What we know: Production uses Python 3.10 (`mace4ir_v2` env). The venv uses 3.12. pyproject.toml declares `requires-python = ">=3.9"`.
   - What's unclear: Whether tests pass on Python 3.10 vs 3.12 vs other versions.
   - Recommendation: Use Python 3.10 (matches production `mace4ir_v2` env). Single version, not a matrix (Claude's discretion per CONTEXT.md).

2. **Should numpy/ase/pyzmq be installed in CI?**
   - What we know: `pip install -e . --no-deps` skips all declared deps. Some test imports need numpy (gaussian parser uses it), ase (mode_matching), pyzmq (zmq_server tests).
   - What's unclear: Exact minimum set of lightweight deps needed for the 127 passing tests.
   - Recommendation: Install `numpy ase pyzmq` alongside pytest/pytest-cov/ruff. These have no transitive Linux issues.

3. **Does ruff format need `--check` or `--diff` in CI?**
   - What we know: `ruff format --check` exits non-zero if any file would be reformatted (fail the CI). `ruff format --diff` shows the diff but always exits 0 (informational only). Currently 12 files would be reformatted.
   - What's unclear: Whether to auto-format before CI enforcement or just enforce.
   - Recommendation: Use `--check` (fail CI on unformatted code). But run `ruff format mace_gaussian/ tests/` locally first to fix the 12 files before pushing the workflow.

## Sources

### Primary (HIGH confidence)
- https://docs.astral.sh/uv/guides/integration/github/ — uv GitHub Actions integration guide
- https://docs.astral.sh/ruff/linter/ — ruff rule configuration and extend-select syntax
- https://github.com/astral-sh/setup-uv — setup-uv@v7 parameters and usage
- https://docs.github.com/en/actions/tutorials/build-and-test-code/python — GitHub official Python CI docs
- https://pytest-cov.readthedocs.io/en/latest/reporting.html — pytest-cov report options

### Secondary (MEDIUM confidence)
- https://ber2.github.io/posts/2025_github_actions_python/ — 2025 GitHub Actions Python best practices
- Local codebase analysis: `ruff check --select B,SIM,PTH,RUF mace_gaussian/` confirmed 65 violations
- Local test run: `pytest tests/ -m "not gpu and not gaussian and not slow"` confirmed 127 passing, 1 failing (stale test)

### Tertiary (LOW confidence)
- None — all findings verified against official sources or direct codebase analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — verified against official docs and local codebase
- Architecture: HIGH — workflow YAML patterns from official uv/GitHub docs; dgl issue confirmed locally
- Pitfalls: HIGH — all pitfalls discovered through direct local testing (not speculation)

**Research date:** 2026-02-25
**Valid until:** 2026-03-25 (GitHub Actions syntax is stable; ruff is fast-moving but breaking changes rare)
