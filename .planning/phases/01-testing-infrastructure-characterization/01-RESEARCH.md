# Phase 1: Testing Infrastructure & Characterization - Research

**Researched:** 2026-02-16
**Domain:** Python testing (pytest) for scientific Gaussian output parsing and vibrational mode analysis
**Confidence:** HIGH

## Summary

This phase establishes the test suite that protects the codebase during subsequent refactoring. The project currently has **zero tests** -- no `tests/` directory, no `conftest.py`, no pytest configuration. The core modules to test are pure Python parsers (`gaussian_parser.py`, `fchk_parser.py`) and numerical algorithms (`mode_matching.py`, `analyze_spectra.py`) that operate on text/numerical data, making them highly amenable to unit testing with committed fixtures.

The key challenge is fixture management: Gaussian `.log` and `.fchk` files are currently `.gitignore`d (lines 29-33 of `.gitignore`), so test fixtures need explicit allowlisting via `!tests/fixtures/**`. File sizes are manageable -- water fixtures total ~430KB (log+fchk for DFT and ML), CH4 fixtures total ~940KB. These can be committed without concern.

A critical discovery is the **acetic acid parsing bug**: the parser's `parse_anharmonic_frequencies()` method requires an `I(anharm)` or `DS(anharm)` header to enter the Fundamental Bands section, but the acoh ML external calculation log files lack the "Anharmonic Infrared Spectroscopy" section entirely. Additionally, there is a "H" prefix on some lines (indicating high overlap with assigned state after variational correction, e.g., `H      4(1)      active`) that breaks the regex pattern `^\s*(\d+)\(1\)`. A secondary bug: `parse_combination_bands()` in the water DFT log captures combination bands twice (once from `I(anharm)` section, once from `DS(anharm)` section), yielding 6 entries instead of the expected 3.

**Primary recommendation:** Use pytest 7.x (already installed) with pytest-cov, create a `tests/` directory with `conftest.py` providing fixture loading, and commit trimmed Gaussian log/fchk excerpts as test fixtures. Register custom markers `gpu`, `gaussian`, and `slow` in `pyproject.toml`.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pytest | >=7.0 (7.1.3 installed) | Test framework | Already in dev deps, standard for Python |
| pytest-cov | >=4.0 | Coverage reporting | Standard pytest coverage plugin |
| numpy | (installed) | Numerical assertions | `np.testing.assert_allclose` for float comparisons |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest-xdist | >=3.0 | Parallel test execution | Optional, for future CI speedup |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| pytest-cov | coverage.py directly | pytest-cov integrates more cleanly with pytest CLI |
| Hand-crafted fixture strings | Committed real files | Real files ensure parser handles actual Gaussian output edge cases |

**Installation:**
```bash
uv add --dev pytest-cov
```

## Architecture Patterns

### Recommended Project Structure
```
tests/
├── conftest.py              # Shared fixtures, markers, fixture path helpers
├── fixtures/
│   ├── water/
│   │   ├── dft_b3lyp.log        # Trimmed water DFT log (anharmonic section)
│   │   ├── dft_b3lyp.fchk       # Full water DFT fchk (1701 lines, 128KB)
│   │   ├── ml_mace_mp_esp.log   # Trimmed water ML log
│   │   ├── ml_mace_mp_esp.fchk  # Water ML fchk
│   │   └── results.json         # Reference parsed output
│   ├── CH4_ase/
│   │   ├── dft_b3lyp.log        # Trimmed CH4 DFT log
│   │   ├── dft_b3lyp.fchk       # CH4 DFT fchk (3114 lines, 240KB)
│   │   ├── ml_mace_mp_esp.log   # Trimmed CH4 ML log
│   │   ├── ml_mace_mp_esp.fchk  # CH4 ML fchk
│   │   └── results.json         # Reference parsed output
│   └── acoh/
│       └── ml_mace_mp_esp.log   # Trimmed acoh ML log (demonstrates parsing bug)
├── test_gaussian_parser.py  # TEST-01, TEST-02, TEST-08
├── test_fchk_parser.py      # TEST-03
├── test_mode_matching.py    # TEST-04
└── test_regression.py       # TEST-05, TEST-06 (reference output validation)
```

### Pattern 1: Fixture Path Helper
**What:** Central fixture path resolution in conftest.py
**When to use:** All test modules that load fixture files
**Example:**
```python
# conftest.py
import pytest
from pathlib import Path

FIXTURES_DIR = Path(__file__).parent / "fixtures"

@pytest.fixture
def fixtures_dir():
    return FIXTURES_DIR

@pytest.fixture
def water_dft_log():
    return str(FIXTURES_DIR / "water" / "dft_b3lyp.log")

@pytest.fixture
def water_dft_fchk():
    return str(FIXTURES_DIR / "water" / "dft_b3lyp.fchk")
```

### Pattern 2: Known Expected Values as Test Constants
**What:** Hard-coded expected parse results derived from verified runs
**When to use:** Parser regression tests
**Example:**
```python
# test_gaussian_parser.py
WATER_DFT_HARMONIC_EXPECTED = [
    {"freq_cm": 1665.3001, "ir_intensity": 70.3304},
    {"freq_cm": 3799.2203, "ir_intensity": 1.6374},
    {"freq_cm": 3912.4174, "ir_intensity": 20.2278},
]

def test_parse_harmonic_frequencies(water_dft_log):
    parser = GaussianLogParser(water_dft_log)
    result = parser.parse_harmonic_frequencies()
    assert len(result) == 3
    for actual, expected in zip(result, WATER_DFT_HARMONIC_EXPECTED):
        assert actual["freq_cm"] == pytest.approx(expected["freq_cm"], abs=0.001)
        assert actual["ir_intensity"] == pytest.approx(expected["ir_intensity"], abs=0.001)
```

### Pattern 3: Parametrized Tests for Multiple Molecules
**What:** Same test logic applied to water and CH4 fixtures
**When to use:** Cross-molecule regression checks
**Example:**
```python
@pytest.mark.parametrize("molecule,expected_n_modes", [
    ("water", 3),
    ("CH4_ase", 9),  # 5 atoms -> 9 modes (3N-6)
])
def test_harmonic_mode_count(molecule, expected_n_modes, fixtures_dir):
    log_path = fixtures_dir / molecule / "dft_b3lyp.log"
    parser = GaussianLogParser(str(log_path))
    result = parser.parse_harmonic_frequencies()
    assert len(result) == expected_n_modes
```

### Pattern 4: xfail for Known Bugs
**What:** Mark acoh parser test as expected failure to document the bug
**When to use:** TEST-08 (acetic acid bug)
**Example:**
```python
@pytest.mark.xfail(reason="Acoh ML log lacks Anharmonic IR section; parser regex also misses H-prefixed lines")
def test_acoh_anharmonic_parsing(fixtures_dir):
    log_path = fixtures_dir / "acoh" / "ml_mace_mp_esp.log"
    parser = GaussianLogParser(str(log_path))
    result = parser.parse_anharmonic_frequencies()
    assert len(result) == 18  # Acoh has 18 modes (8 atoms -> 3*8-6=18)
```

### Anti-Patterns to Avoid
- **Do NOT use full untrimmed log files as fixtures:** The full water DFT log is 2710 lines with ~2000 lines of Gaussian boilerplate. Trim to the relevant parsing sections (frequencies, anharmonic, archive) to keep fixtures focused and fast to read.
- **Do NOT hard-code absolute paths in tests:** Use the `fixtures_dir` fixture for portability.
- **Do NOT test plotting/visualization in unit tests:** The `plot_*` functions in `mode_matching.py` and `analyze_spectra.py` produce matplotlib figures. These are integration-level. Focus unit tests on data parsing and numerical computations.
- **Do NOT import modules that trigger MACE/torch in unit tests:** `mode_matching.py` imports `fchk_parser` which is safe, but avoid importing `mace_calculators.py` or `gm_main.py` in unit tests as they trigger CUDA/torch initialization.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Float comparison | Custom epsilon logic | `pytest.approx()` or `np.testing.assert_allclose()` | Handles relative/absolute tolerance correctly |
| Test discovery | Manual test registration | pytest auto-discovery (`test_*.py`) | Standard, zero config |
| Fixture cleanup | Manual tmp file management | `tmp_path` pytest fixture | Built-in, auto-cleaned |
| Coverage reporting | Custom instrumentation | pytest-cov with `--cov` flag | HTML + terminal reports |
| Test selection | Custom skip logic | `pytest.mark` decorators + `-m` flag | Standard, composable |

**Key insight:** pytest's fixture system and marker system already solve every organizational challenge this phase faces. The only custom work is writing the actual test assertions and preparing the fixture files.

## Common Pitfalls

### Pitfall 1: .gitignore Blocks Test Fixtures
**What goes wrong:** `.gitignore` contains `*.log`, `*.fchk`, `*.chk`, `*.gjf` rules which will prevent test fixtures from being committed.
**Why it happens:** The project gitignore was written before a test suite existed.
**How to avoid:** Add explicit allowlist rules at the bottom of `.gitignore`:
```
# Allow test fixtures
!tests/fixtures/**/*.log
!tests/fixtures/**/*.fchk
!tests/fixtures/**/*.json
```
**Warning signs:** `git add tests/fixtures/` shows no files staged.

### Pitfall 2: Harmonic Frequencies Parsed Twice
**What goes wrong:** `parse_harmonic_frequencies()` returns duplicate entries because Gaussian log files repeat the frequency section -- once in the initial frequency calculation and again in the anharmonic pre-analysis.
**Why it happens:** The parser uses a `seen_freqs` set with `(round(freq, 4), round(intensity, 4))` keys to deduplicate, which works for exact duplicates. But the test should verify that this deduplication is actually working.
**How to avoid:** Tests should assert exact expected counts: water has 3 harmonic modes, CH4 has 9 (including degenerate modes).
**Warning signs:** Harmonic count is double expected.

### Pitfall 3: Combination Bands Parsed Twice from Different Sections
**What goes wrong:** `parse_combination_bands()` (and `parse_overtones()`) parse from both the `I(anharm)` section and the `DS(anharm)` section, because the parser enters the section again when it finds another "Combination Bands" header.
**Why it happens:** The Gaussian log repeats Fundamental Bands / Overtones / Combination Bands twice: once with `I(anharm)` intensities and once with `DS(anharm)` dipole strengths. The parser's break condition for combination bands (`Electric dipole :` or `Rotational Constants` or `==`) exits the first block, but then the same pattern appears again.
**How to avoid:** The test should document this as a known issue. Water DFT log produces 6 combination bands instead of 3.
**Warning signs:** Combination band count is 2x expected.

### Pitfall 4: fchk_parser Calls formchk Subprocess
**What goes wrong:** `convert_chk_to_fchk()` runs `subprocess.run(['formchk', ...])` which requires Gaussian 16 installed. `get_fchk_from_chk()` also calls this.
**Why it happens:** The function is designed for environments with Gaussian.
**How to avoid:** In unit tests, always use `.fchk` files directly (already converted). Never pass `.chk` files to test functions. Mark any test that requires Gaussian conversion with `@pytest.mark.gaussian`.
**Warning signs:** `FileNotFoundError: formchk command not found`.

### Pitfall 5: FCHK Section Parsing Sensitivity to Line Format
**What goes wrong:** `parse_fchk_section()` uses `re.match(r'^[A-Z]', line)` to detect new section boundaries. This can fail if data values happen to start at column 0 or if there are unexpected header formats.
**Why it happens:** The .fchk format is somewhat rigid but edge cases exist.
**How to avoid:** Test with actual committed .fchk files rather than synthetic strings.
**Warning signs:** `ValueError: Expected N values for 'Section', got M`.

### Pitfall 6: Acetic Acid Bug Has Two Distinct Root Causes
**What goes wrong:** The acoh anharmonic parsing fails completely (returns empty list).
**Why it happens:** TWO separate issues:
1. The ML external calculation log for acoh does NOT contain the "Anharmonic Infrared Spectroscopy" section (with `I(anharm)` / `DS(anharm)` headers). The parser looks for these headers to set `in_fundamental_section = True`, so it never enters the parsing loop.
2. Even if it did enter the section, the "Vibrational Energies at Anharmonic Level" section has lines prefixed with "H" or "L" overlap indicators (e.g., `H      4(1)      active`) that don't match the regex `^\s*(\d+)\(1\)`.
**How to avoid:** Test should be `pytest.mark.xfail` and document both root causes clearly in the test docstring.
**Warning signs:** `parse_anharmonic_frequencies()` returns empty list for acoh.

## Code Examples

### Verified Expected Values for Water DFT (from actual parser run)

```python
# Water B3LYP/6-31G(d,p) - verified 2026-02-16
# Source: comparison_results/water/b3lyp_6-31Gdp/gaussian_freq.log

WATER_HARMONIC = [
    {"freq_cm": 1665.3001, "ir_intensity": 70.3304},
    {"freq_cm": 3799.2203, "ir_intensity": 1.6374},
    {"freq_cm": 3912.4174, "ir_intensity": 20.2278},
]

WATER_ANHARMONIC = [
    {"mode": 1, "freq_cm": 3624.001, "ir_intensity": 0.62829225, "freq_harmonic": 3799.22},
    {"mode": 2, "freq_cm": 1615.153, "ir_intensity": 69.29467704, "freq_harmonic": 1665.3},
    {"mode": 3, "freq_cm": 3722.273, "ir_intensity": 16.09760189, "freq_harmonic": 3912.417},
]

WATER_OVERTONES = [
    {"mode": 1, "overtone_level": 2, "freq_harmonic": 7598.441,
     "freq_anharmonic": 7161.57, "ir_intensity": 0.77569767},
    {"mode": 2, "overtone_level": 2, "freq_harmonic": 3330.6,
     "freq_anharmonic": 3194.746, "ir_intensity": 0.80602582},
    {"mode": 3, "overtone_level": 2, "freq_harmonic": 7824.835,
     "freq_anharmonic": 7346.675, "ir_intensity": 0.04914658},
]

# NOTE: Parser currently returns 6 combination bands (duplicates from I(anharm) + DS(anharm))
# True unique count should be 3
WATER_COMBINATION_BANDS_RAW = 6  # Current parser behavior
WATER_COMBINATION_BANDS_UNIQUE = 3  # Expected after bug fix

WATER_ENERGY = 0.003705  # Hartrees (approximate, from "Energy=" pattern)
```

### Verified Expected Values for CH4 ML (from results.json)

```python
# CH4 MACE-MP/Espaloma - verified 2026-02-16
# Source: comparison_results/CH4_ase/mace_mp_espaloma/results.json

CH4_HARMONIC_COUNT = 4  # Only unique harmonic freqs (degenerate modes collapsed)
CH4_ANHARMONIC_COUNT = 9  # All 9 modes (3*5-6) including degenerates
CH4_OVERTONES_COUNT = 9
# CH4 has many combination bands (36 pairs for 9 modes)
```

### Mode Matching Unit Test Pattern

```python
import numpy as np
from mode_matching import compute_mode_overlap, match_modes

def test_identical_modes_perfect_overlap():
    """Two identical mode vectors should have overlap = 1.0"""
    mode = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    overlap = compute_mode_overlap(mode, mode)
    assert overlap == pytest.approx(1.0, abs=1e-10)

def test_opposite_modes_perfect_overlap():
    """Sign-flipped modes should still have overlap = 1.0 (abs value)"""
    mode = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    overlap = compute_mode_overlap(mode, -mode)
    assert overlap == pytest.approx(1.0, abs=1e-10)

def test_orthogonal_modes_zero_overlap():
    """Orthogonal modes should have overlap = 0.0"""
    mode1 = np.array([[1.0, 0.0, 0.0]])
    mode2 = np.array([[0.0, 1.0, 0.0]])
    overlap = compute_mode_overlap(mode1, mode2)
    assert overlap == pytest.approx(0.0, abs=1e-10)

def test_match_modes_identity():
    """Matching modes against themselves should give 1-to-1 mapping"""
    modes = np.random.randn(5, 3, 3)  # 5 modes, 3 atoms, 3 coords
    matches = match_modes(modes, modes, threshold=0.5)
    for calc_idx, (ref_idx, overlap) in matches.items():
        assert calc_idx == ref_idx
        assert overlap == pytest.approx(1.0, abs=1e-10)
```

### pyproject.toml Configuration

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "gpu: marks tests that require CUDA GPU",
    "gaussian: marks tests that require Gaussian 16 (g16/formchk)",
    "slow: marks tests that take >30 seconds",
]
addopts = [
    "--strict-markers",
    "-ra",
]

[tool.coverage.run]
source = [
    "gaussian_parser",
    "fchk_parser",
    "mode_matching",
    "analyze_spectra",
    "comparison_workflow",
    "results_manager",
    "mace_calculators",
    "gm_main",
    "cli",
    "dft_baseline",
]
omit = [
    "tests/*",
    "mace_dipole_pkg/*",
    "mace_ML_pkg/*",
]

[tool.coverage.report]
show_missing = true
skip_empty = true
```

### Fixture Trimming Strategy

```python
# Script to create trimmed fixtures from full logs
# Keep only: header (5 lines), frequency sections, anharmonic sections, archive section
# Goal: <100 lines per fixture for parseable content

# For harmonic: Keep lines containing "Frequencies --" through "IR Inten --" + context
# For anharmonic: Keep "Fundamental Bands" through "Combination Bands" end
# For energy: Keep "Energy=" lines
# For dipole: Keep "Dipole=" archive line
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| No tests | Need tests before refactoring | Now (Phase 1) | Enables safe refactoring in Phases 2-5 |
| `uv.lock` in `.gitignore` | `uv.lock` committed (REPR-04) | Phase 1 | Reproducible dev environment |
| pytest only in dev deps | pytest + pytest-cov in dev deps | Phase 1 | Coverage visibility |

**Deprecated/outdated:**
- The `.gitignore` currently excludes `uv.lock` (line 57) -- this needs to be removed per REPR-04

## Open Questions

1. **Fixture trimming vs full files**
   - What we know: Full log files are 72KB-528KB, feasible to commit but noisy
   - What's unclear: Whether trimmed excerpts might miss edge cases in section boundary detection
   - Recommendation: Commit FULL `.fchk` files (needed intact for numerical parsing), but TRIM `.log` files to relevant sections. Include enough context lines around section boundaries to test parser's section detection logic.

2. **Combination band duplication bug -- document or fix?**
   - What we know: `parse_combination_bands()` returns duplicates from I(anharm) and DS(anharm) sections
   - What's unclear: Whether downstream code depends on this behavior (the results.json files show duplicate combination bands in water)
   - Recommendation: Write the test to document current behavior, add a comment noting it's a bug. Fixing it is a refactoring task (Phase 2+), not Phase 1.

3. **Scope of "reference outputs" (TEST-05, TEST-06)**
   - What we know: The `results.json` files in `comparison_results/` contain the full parsed output
   - What's unclear: Whether "reference outputs" means committed JSON snapshots or just hard-coded expected values in tests
   - Recommendation: Commit the `results.json` files as fixtures AND have hard-coded expected values in tests. The JSON provides the full regression baseline; the hard-coded values make test intent explicit.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `gaussian_parser.py`, `fchk_parser.py`, `mode_matching.py`, `analyze_spectra.py` -- direct source code reading
- Actual parser execution on water DFT log -- verified expected values
- Actual parser execution on acoh ML log -- confirmed bug (0 anharmonic frequencies parsed)
- Git commit `a4384c4` -- confirmed acoh bug documented in commit message
- `.gitignore` analysis -- confirmed fixture blocking rules
- `pyproject.toml` -- confirmed pytest >=7.0 in dev deps, no test config exists
- [pytest documentation: How to use fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html)
- [pytest documentation: How to mark test functions](https://docs.pytest.org/en/stable/how-to/mark.html)
- [pytest-cov configuration docs](https://pytest-cov.readthedocs.io/en/latest/config.html)

### Secondary (MEDIUM confidence)
- [pytest documentation: Configuration](https://docs.pytest.org/en/stable/reference/customize.html) -- pyproject.toml `[tool.pytest.ini_options]`
- [pytest-cov PyPI](https://pypi.org/project/pytest-cov/) -- version 7.0.0 available

### Tertiary (LOW confidence)
- None. All findings verified from primary sources.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- pytest 7.x already installed, pytest-cov is the universal choice
- Architecture: HIGH -- test directory structure is standard Python convention
- Pitfalls: HIGH -- all pitfalls verified by direct code analysis and test execution
- Acoh bug analysis: HIGH -- reproduced and root-caused from actual log file inspection

**Research date:** 2026-02-16
**Valid until:** 2026-06-16 (stable domain, pytest/coverage tooling changes slowly)
