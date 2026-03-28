# Technology Stack: v1.2 Analysis Quality Overhaul

**Project:** MACE-Gaussian Interface
**Researched:** 2026-03-28
**Scope:** Stack additions for new analysis quality features only

## Existing Stack (NOT changing)

Python 3.10, ASE, torch, numpy, scipy, matplotlib, seaborn, ZMQ. Already uses `time.perf_counter()` for per-call timing in `workflow.py`.

---

## Recommended Stack Additions

### 1. Lorentzian Broadening: NO new dependency

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| numpy + scipy (existing) | -- | Lorentzian broadening | Trivial math: `(gamma/pi) / ((x - x0)**2 + gamma**2)`. 10-15 lines of code. No library needed. |

**Rationale:** Lorentzian broadening is a direct analytical formula applied to stick spectra. The function `L(x; x0, gamma) = (1/pi) * gamma / ((x - x0)^2 + gamma^2)` is summed over all modes weighted by intensity. numpy vectorization handles this in a single broadcast. scipy is only needed if you want to use `scipy.signal.fftconvolve` for performance on extremely dense grids (unlikely for IR with <100 modes).

**Implementation pattern:**
```python
def lorentzian_broaden(frequencies, intensities, x_range, gamma=4.0):
    """Apply Lorentzian broadening to stick spectrum.

    Args:
        frequencies: Mode frequencies in cm^-1
        intensities: IR intensities in km/mol
        x_range: Wavenumber grid (e.g., np.linspace(400, 4000, 3601))
        gamma: Half-width at half-maximum in cm^-1 (default 4.0)

    Returns:
        Broadened spectrum as 1D array matching x_range
    """
    spectrum = np.zeros_like(x_range)
    for freq, inten in zip(frequencies, intensities):
        spectrum += inten * (gamma / np.pi) / ((x_range - freq)**2 + gamma**2)
    return spectrum
```

Typical HWHM values: 4 cm^-1 for gas-phase comparison, 8-15 cm^-1 for solution/solid-state.

**Alternatives considered:**
| Package | Why Not |
|---------|---------|
| [galore](https://github.com/SMTG-Bham/galore) | Overkill for this; designed for XPS/DOS broadening, pulls in extra deps |
| [RADIS](https://radis.readthedocs.io/) | Atmospheric line-by-line code; massive dependency for a 10-line function |

### 2. Experimental Spectra Comparison: nistchempy + jcamp

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| nistchempy | 1.0.5 | Fetch IR spectra from NIST Chemistry WebBook | Only Python package that wraps NIST WebBook. Extracts IR spectra by compound name/CAS/InChI. MIT license. Active (updated Mar 2026). |
| jcamp | 1.3.0 | Parse JCAMP-DX spectral files | Standard JCAMP-DX parser. Returns dict with `x` (wavenumber) and `y` (absorbance/transmittance) arrays. MIT license. Handles NIST's JCAMP format directly. |

**How they work together:**
1. `nistchempy.run_search(name="water")` returns compound objects
2. Compound has `.ir_specs` property with IR spectrum URLs
3. Download spectrum data (JCAMP-DX format from NIST)
4. Parse with `jcamp.jcamp_reader()` to get `x`, `y` arrays
5. Compare broadened ML spectrum against experimental

**NIST access details (IMPORTANT):**
- NIST Chemistry WebBook has NO official API. NistChemPy is a web scraper.
- Rate limiting is required: add delays between requests (1-2 seconds).
- Not all compounds have IR spectra. Coverage is ~16,000 compounds.
- Spectra are in transmittance (%) or absorbance. Must normalize for comparison.
- Gas-phase spectra are most comparable to computed (no solvent effects).

**SDBS (Japan AIST) alternative:**
- Contains ~34,000 FT-IR spectra (more than NIST)
- **NO API whatsoever.** No programmatic access. Web-only interface.
- Terms of use explicitly prohibit automated access.
- **Verdict: Do not use.** Stick with NIST via nistchempy.

**Confidence:** MEDIUM. NistChemPy works but is unofficial scraping -- could break if NIST changes their HTML. The `jcamp` library is stable and well-tested.

### 3. VPT2 Alternatives: Deferred / Research-only

| Technology | Version | Purpose | Verdict |
|------------|---------|---------|---------|
| PyVPT2 | latest (conda) | Standalone VPT2 computations | **DO NOT ADD** -- requires Psi4 1.8+ as hard dependency |
| Psience/VPT2 | latest (git) | McCoy Group standalone VPT2 | **INVESTIGATE ONLY** -- may accept external force constants |

**PyVPT2 analysis (HIGH confidence):**
- Published Jan 2025 in J. Chem. Phys. by Nelson & Sherrill (Georgia Tech)
- BSD-3-Clause license, conda-forge package
- **Hard dependency on Psi4 1.8+** -- imports `psi4` at module level, uses `psi4.core.Wavefunction` throughout
- Cannot consume pre-computed force constants from Gaussian or MACE
- Designed to compute its own cubic/quartic force constants via finite differences through Psi4/QCEngine
- Would require installing Psi4 (~2GB conda environment) alongside existing MACE stack
- **Verdict: Incompatible with MACE-Gaussian architecture.** The project already gets anharmonic frequencies from Gaussian's VPT2 implementation. PyVPT2 would be a replacement for Gaussian, not a complement, and it can't use ML potentials.

**Psience VPT2 analysis (LOW confidence -- needs deeper investigation):**
- McCoy Group (U. Washington), MIT license, actively maintained (Mar 2026)
- Pure Python VPT2 implementation with `McUtils` as only non-standard dependency
- Has `DegeneracySpecs.py` -- handles Fermi resonances and degeneracies natively
- Accepts `Molecule` objects with potential energy surfaces
- **Key question:** Can it consume pre-computed Hessian + cubic/quartic force constants from Gaussian .fchk files? The `Molecule` class may accept these, but this needs hands-on testing.
- **If yes:** Could provide independent VPT2 validation against Gaussian's implementation using the same force constants -- valuable for thesis.
- **If no:** Same problem as PyVPT2 -- needs its own QC backend.
- **Verdict: Flag for phase-specific research.** Worth a 2-hour spike to test if Psience can ingest Gaussian force constants.

**Prior ORCA research (from v1.1-orca-vpt2.md):**
- ORCA lacks external hooks for VPT2 -- technically blocked.
- Already marked out of scope in PROJECT.md.

### 4. Degenerate Mode Handling: NO new dependency

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| scipy.spatial.distance (existing) | -- | Frequency clustering for degeneracy detection | Already available. Use threshold-based grouping (e.g., modes within 1 cm^-1). |

**Rationale:** Degenerate mode handling is algorithmic, not a library problem. The codebase already has degenerate mode awareness in the Gaussian parser (see `test_ch4_harmonic_degenerate_deduplication`). What's needed:

1. **Detection:** Group modes whose frequencies differ by < threshold (e.g., 1 cm^-1 for exact degeneracies, 5-10 cm^-1 for near-degeneracies)
2. **Intensity aggregation:** Sum intensities of degenerate modes for broadened spectrum (individual modes may split slightly in anharmonic treatment)
3. **Mode matching:** Hungarian algorithm already handles this; ensure degenerate modes are properly represented in the overlap matrix rather than being collapsed

No new dependencies needed. This is a logic fix in the parser and mode matching code.

### 5. Wall-Clock Timing: NO new dependency

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| time.perf_counter (stdlib) | -- | Wall-clock timing per pipeline stage | Already used in workflow.py. Extend pattern to full pipeline. |
| contextlib.contextmanager (stdlib) | -- | Timer context manager for clean instrumentation | Standard pattern for wrapping code blocks. |

**Rationale:** The project already uses `time.perf_counter()` for per-Gaussian-call timing in `workflow.py`. What's needed is:

1. A `Timer` context manager or decorator that collects named durations
2. Wrap each pipeline stage: optimization, harmonic freq, anharmonic freq, DFT baseline
3. Write timing data to results JSON alongside existing calculation_parameters
4. Add timing summary to HTML report

**Implementation pattern:**
```python
import time
from contextlib import contextmanager

class PipelineTimer:
    def __init__(self):
        self.timings = {}

    @contextmanager
    def stage(self, name):
        t0 = time.perf_counter()
        yield
        self.timings[name] = time.perf_counter() - t0

    def to_dict(self):
        return {**self.timings, "total": sum(self.timings.values())}
```

**Alternatives considered:**
| Package | Why Not |
|---------|---------|
| `line_profiler` | For dev profiling, not production instrumentation |
| `scalene` | Same -- profiling tool, not timing data collection |
| `timeit` | For microbenchmarks, not pipeline stage timing |

stdlib `time.perf_counter()` is the right tool. No dependency needed.

---

## Summary: What to Install

### Required new packages (2 packages)

```bash
pip install nistchempy==1.0.5 jcamp==1.3.0
```

Or in `pyproject.toml`:
```toml
dependencies = [
    # ... existing ...
    "nistchempy>=1.0.5",
    "jcamp>=1.3.0",
]
```

### Optional (investigate in spike)

```bash
pip install McUtils Psience  # McCoy Group VPT2 -- for spike only
```

### NOT installing

| Package | Reason |
|---------|--------|
| PyVPT2 / psi4 | Hard Psi4 dependency, can't use ML potentials |
| galore | Overkill for Lorentzian broadening |
| RADIS | Overkill, massive deps |
| Any SDBS scraper | No API, prohibited by ToS |

---

## Alternatives Considered (Full Table)

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Broadening | numpy (existing) | galore, RADIS | 10-line formula vs. pulling in packages |
| Exp. spectra source | NIST via nistchempy | SDBS | SDBS prohibits automated access |
| JCAMP parsing | jcamp 1.3.0 | pyjdx | jcamp is more mature, MIT, better documented |
| VPT2 | Gaussian (existing) | PyVPT2, Psience | PyVPT2 needs Psi4; Psience needs investigation |
| Timing | time.perf_counter (stdlib) | line_profiler, scalene | Profiling tools, not production instrumentation |
| Degeneracy | scipy (existing) | -- | Algorithmic fix, not library problem |

---

## Integration Points

### nistchempy integration
- New module: `mace_gaussian/analysis/experimental_comparison.py`
- Fetch spectra during analysis phase (after ML/DFT calculations complete)
- Cache fetched spectra locally (JCAMP files in molecule result directory)
- Rate-limit requests: minimum 2 seconds between NIST hits
- Graceful fallback: if no experimental spectrum found, skip comparison (don't fail)

### jcamp integration
- Parse cached JCAMP-DX files to numpy arrays
- Convert transmittance to absorbance if needed: `A = -log10(T/100)`
- Interpolate onto common wavenumber grid for comparison with broadened ML spectrum
- Use `scipy.interpolate.interp1d` (already available) for grid alignment

### Lorentzian broadening integration
- New function in `mace_gaussian/analysis/spectral_comparison.py`
- Called from analysis workflow after frequency parsing
- Takes stick spectrum (frequencies + intensities) -> continuous spectrum
- Used for: (1) ML vs DFT visual comparison, (2) ML vs experimental comparison
- Default gamma=4.0 cm^-1 (gas-phase), configurable via CLI/config

### Timing integration
- Wrap existing `time.perf_counter()` calls in `workflow.py` with PipelineTimer
- Add timer to `run_pipeline()` in cli.py for full end-to-end timing
- Store in results.json under `"timings"` key
- Display in HTML report as timing breakdown table

---

## Sources

- [PyVPT2 GitHub](https://github.com/philipmnel/pyvpt2) -- BSD-3, requires Psi4 1.8+ (HIGH confidence)
- [PyVPT2 paper](https://pubmed.ncbi.nlm.nih.gov/39820337/) -- J. Chem. Phys. 162, 032501 (2025)
- [Psience/VPT2](https://github.com/McCoyGroup/Psience) -- McCoy Group, MIT, pure Python (LOW confidence on integration)
- [NistChemPy PyPI](https://pypi.org/project/nistchempy/) -- v1.0.5, MIT
- [NistChemPy GitHub](https://github.com/IvanChernyshov/NistChemPy) -- unofficial NIST WebBook API
- [NistChemPy Docs](https://ivanchernyshov.github.io/NistChemPy/) -- v1.0.3 docs (latest published)
- [jcamp PyPI](https://pypi.org/project/jcamp/) -- v1.3.0, MIT
- [jcamp GitHub](https://github.com/nzhagen/jcamp) -- JCAMP-DX parser
- [SDBS](https://sdbs.db.aist.go.jp) -- no API, web-only (HIGH confidence: not usable)
- [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/) -- no official API (HIGH confidence)
- [Galore](https://github.com/SMTG-Bham/galore) -- broadening package, overkill for IR

---
*Stack research for v1.2 milestone: 2026-03-28. Supersedes v1.0 stack research from 2026-02-16 for new features only; original stack decisions remain valid.*
