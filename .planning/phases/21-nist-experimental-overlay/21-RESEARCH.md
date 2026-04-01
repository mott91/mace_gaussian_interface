# Phase 21: NIST Experimental Overlay - Research

**Researched:** 2026-04-01
**Domain:** NIST WebBook IR data fetching, JCAMP-DX parsing, spectral overlay plotting
**Confidence:** HIGH

## Summary

This phase adds experimental IR spectrum overlay to the existing analysis pipeline. Two Python libraries -- `nistchempy` (v1.0.5) and `jcamp` (v1.2.2) -- provide the complete data pipeline: search NIST WebBook by molecule name, download JDX text, and parse into wavenumber/intensity arrays. Both libraries were installed and tested successfully in the `mace4ir_v2` environment.

The existing `SpectrumAnalyzer` plotting methods use a stacked offset layout (DFT on bottom, ML on top) with Lorentzian-broadened spectra. The experimental trace will be added as a third trace (black dashed line) on the existing plots. NIST data comes in wavenumber (cm-1) on x-axis, but y-units vary: ABSORBANCE, TRANSMITTANCE, or cross-section. Transmittance must be converted to absorbance before overlay.

**Primary recommendation:** Create a standalone `nist_fetcher.py` module in `mace_gaussian/analysis/` that handles search, download, caching, and parsing. Hook it into `analysis_workflow.py` at the start of `run_full_analysis()`. Pass parsed experimental data to plotting methods as an optional parameter.

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Use `nistchempy` library to fetch experimental spectra from NIST WebBook
- **D-02:** Parse JCAMP-DX format using `jcamp` library
- **D-03:** Gas-phase IR data only. If no gas-phase IR spectrum, skip silently
- **D-04:** Best-effort fetching -- NIST fetch failure never blocks analysis
- **D-05:** Cache in `comparison_results/{molecule}/experimental/`
- **D-06:** If cached data exists, skip re-download
- **D-07:** Experimental spectrum as third trace on existing DFT-above/ML-below offset plot (not separate panel)
- **D-08:** Both experimental and computed spectra normalized to [0, 1] for visual comparison
- **D-09:** Experimental trace: black dashed line
- **D-10:** No quantitative peak matching (NIST-03 deferred)
- **D-11:** Automatic fetch, no CLI flag. Graceful degradation on failure
- **D-12:** Both `run_analysis.py` and `run_analysis_harmonic.py` support overlay
- **D-13:** Batch report includes per-molecule experimental overlay when available

### Claude's Discretion
- JCAMP-DX parsing details (handling malformed files, unit detection)
- Legend placement and labeling for the experimental trace
- How to note "no experimental data available" in the report (if at all)
- Frequency range alignment between experimental and computed grids
- Cache file format (raw JCAMP-DX, parsed numpy arrays, or both)

### Deferred Ideas (OUT OF SCOPE)
- NIST-03 quantitative peak comparison with MAE/RMSE metrics
- SDBS as secondary spectral source
- Condensed-phase fallback with shift warning

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| NIST-01 | Fetch experimental IR spectrum from NIST WebBook by molecule name (cached locally) | nistchempy `run_search()` + `get_ir_spectra()` API verified working; jcamp parsing tested; caching via JDX file save |
| NIST-02 | Analysis report overlays experimental spectrum on computed spectra plot when available | Existing `plot_spectra_comparison()` and `plot_combined_spectra()` use offset layout; experimental trace adds at bottom with negative offset or overlaid on DFT panel |
| NIST-03 | Quantitative peak position comparison (experimental vs computed) with error metrics | **DEFERRED** per D-10 in CONTEXT.md |

</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| nistchempy | 1.0.5 | Search NIST WebBook, download IR spectra as JDX text | Only Python API for NIST Chemistry WebBook; locked by D-01 |
| jcamp | 1.2.2 | Parse JCAMP-DX format into numpy arrays | Standard JCAMP-DX parser; locked by D-02; returns `x`, `y` arrays with unit metadata |

### Supporting (already in project)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | 1.26.4 | Array operations for spectral data | Interpolation, normalization |
| matplotlib | (installed) | Plotting experimental overlay | Adding third trace to existing figures |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| nistchempy | Direct HTTP scraping | nistchempy already wraps NIST WebBook; scraping is fragile |
| jcamp | Manual JCAMP-DX parsing | JCAMP-DX format has edge cases (DIF/SQZ compression); library handles them |

**Installation:**
```bash
micromamba run -n mace4ir_v2 pip install nistchempy jcamp
```

**Version verification:** Both verified installed and working on 2026-04-01. nistchempy 1.0.5 (latest on PyPI), jcamp 1.2.2 (pip resolves to this from 1.3.0 due to Python 3.10 compatibility).

## Architecture Patterns

### Recommended Module Structure
```
mace_gaussian/analysis/
  nist_fetcher.py          # NEW: fetch, cache, parse experimental spectra
  analyze_spectra.py       # MODIFY: add experimental trace to plot methods
  analysis_workflow.py     # MODIFY: call nist_fetcher before plotting
  html_report_generator.py # MODIFY: add experimental overlay section
  batch_report.py          # MODIFY: per-molecule experimental overlay
```

### Pattern 1: NIST Fetcher Module
**What:** Standalone module with a single public function that returns parsed experimental data or None
**When to use:** Called from analysis_workflow.py before plotting

```python
# mace_gaussian/analysis/nist_fetcher.py

@dataclass
class ExperimentalSpectrum:
    """Parsed experimental IR spectrum from NIST."""
    wavenumbers: np.ndarray   # cm-1
    absorbance: np.ndarray    # normalized to [0, 1]
    source: str               # e.g., "NIST WebBook (Sadtler)"
    molecule_name: str
    cas_number: str

def fetch_experimental_spectrum(
    molecule_name: str,
    cache_dir: Path,
) -> ExperimentalSpectrum | None:
    """Fetch and cache experimental IR spectrum from NIST.

    Returns None if molecule not found or no gas-phase IR data.
    Never raises exceptions -- all errors are logged and return None.
    """
```

### Pattern 2: Cache-First with JDX File
**What:** Cache the raw JDX text file in `comparison_results/{molecule}/experimental/`. On subsequent runs, load from cached JDX file instead of re-fetching from NIST.
**Why both raw + parsed:** Raw JDX preserves provenance; parsing is fast (~1ms) so no need to cache parsed arrays separately.

```python
cache_path = cache_dir / "experimental" / "nist_ir.jdx"
if cache_path.exists():
    return _parse_jdx_file(cache_path)
else:
    # Fetch from NIST, save JDX, parse
    spectrum = _fetch_from_nist(molecule_name)
    if spectrum:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(spectrum.jdx_text)
        return _parse_jdx_file(cache_path)
    return None
```

### Pattern 3: Experimental Trace on Existing Plots
**What:** Add experimental spectrum as optional parameter to existing plot methods
**How:** The current plot layout uses fixed offset=1.5 (DFT at y=0, ML at y=1.5). Experimental data goes on the DFT panel (y=0) as an overlaid black dashed line, since we want to compare experimental vs DFT and experimental vs ML.

```python
def plot_spectra_comparison(
    self,
    ml_spectrum: SpectrumData,
    dft_spectrum: SpectrumData,
    ml_name: str,
    molecule_name: str | None = None,
    save_path: str | None = None,
    experimental: ExperimentalSpectrum | None = None,  # NEW
) -> plt.Figure:
    # ... existing DFT + ML plotting ...
    if experimental is not None:
        # Interpolate experimental data onto our frequency grid
        exp_interp = np.interp(self.freq_grid, experimental.wavenumbers, experimental.absorbance)
        # Normalize to [0, 1]
        exp_max = np.max(exp_interp)
        if exp_max > 0:
            exp_norm = exp_interp / exp_max
        else:
            exp_norm = exp_interp
        # Plot on both panels (DFT at y=0, ML at y=offset)
        ax1.plot(self.freq_grid, exp_norm, '--', color='black', linewidth=1.5,
                 label='Experimental (NIST)', alpha=0.7, zorder=5)
        ax1.plot(self.freq_grid, exp_norm + offset, '--', color='black', linewidth=1.5,
                 alpha=0.7, zorder=5)
```

### Anti-Patterns to Avoid
- **Fetching NIST data inside plot methods:** Keep data fetching separate from plotting. Fetch once in workflow, pass to all plot methods.
- **Blocking analysis on NIST failure:** Every NIST call must be wrapped in try/except. Return None, log warning, continue.
- **Caching parsed numpy arrays:** The JDX file is tiny (~5-40KB) and parses in milliseconds. Caching parsed data adds complexity without benefit.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NIST WebBook search | HTTP scraping + HTML parsing | `nistchempy.run_search()` | NIST page structure changes; library handles it |
| JCAMP-DX parsing | Regex-based parser | `jcamp.jcamp_readfile()` | JCAMP-DX has DIF/SQZ compressed formats; library handles decompression |
| Transmittance-to-absorbance | N/A | Simple formula: `A = -log10(T/100)` if T in %, or `A = -log10(T)` if T in [0,1] | Straightforward but must check units |
| Spectral interpolation | Custom interpolation | `np.interp()` | Standard linear interpolation; experimental data has 880-56000 points, much denser than our 0.5 cm-1 grid |

**Key insight:** The hard part is not the math -- it's handling the variety of NIST data formats (different y-units, different x-ranges, different providers) robustly.

## Common Pitfalls

### Pitfall 1: Y-Unit Mismatch
**What goes wrong:** NIST spectra come in different y-units: ABSORBANCE, TRANSMITTANCE, or cross-section `(micromol/mol)-1m-1 (base 10)`. Plotting transmittance data directly against absorbance-based computed spectra produces inverted or misscaled overlays.
**Why it happens:** Different NIST data providers use different units. Even for the same molecule (e.g., methanol), spectra 0-1 use TRANSMITTANCE and spectra 2-4 use cross-section.
**How to avoid:** Check `data['yunits']` after jcamp parsing. Convert transmittance to absorbance. For cross-section units, normalize to [0,1] directly (shape is correct, only scale differs).
**Warning signs:** Experimental trace appears inverted (peaks go down instead of up) or flat.

### Pitfall 2: Multiple Spectra Per Molecule
**What goes wrong:** NIST often provides 2-5 IR spectra per molecule from different sources (e.g., Sadtler, Coblentz, PNNL). They have different resolution, range, and quality.
**Why it happens:** NIST aggregates data from multiple providers.
**How to avoid:** Prefer the first gas-phase spectrum with ABSORBANCE or TRANSMITTANCE units. The high-resolution PNNL cross-section data (14K-56K points) is excellent but uses unusual units. Pick the simplest one for visual overlay.
**Warning signs:** Unexpected data range or density.

### Pitfall 3: Molecule Name Ambiguity
**What goes wrong:** `nistchempy.run_search('cocaine', 'name')` returns 2 compounds. Common names may match multiple molecules or no molecules.
**Why it happens:** NIST uses common names and synonyms; some names are ambiguous.
**How to avoid:** Take the first compound if `num_compounds >= 1`. If 0 compounds found, return None gracefully. Log warnings for ambiguous matches.
**Warning signs:** `search.num_compounds > 1` or `search.num_compounds == 0`.

### Pitfall 4: NIST Rate Limiting / Network Failures
**What goes wrong:** NIST WebBook may be slow, temporarily unavailable, or rate-limit rapid requests.
**Why it happens:** It's a government service, not a high-availability API.
**How to avoid:** Wrap all network calls in try/except with timeout. Cache aggressively. For batch runs, add a small delay between molecules. Never let network failure block analysis.
**Warning signs:** `requests.exceptions.Timeout` or HTTP 503.

### Pitfall 5: X-Range Mismatch
**What goes wrong:** Experimental spectrum may cover 450-3966 cm-1 while computed spectrum covers 400-4000 cm-1. `np.interp()` will extrapolate to 0 outside the experimental range.
**Why it happens:** Different instruments have different spectral ranges.
**How to avoid:** `np.interp()` naturally handles this -- values outside the experimental range become 0 (if `left=0, right=0` is passed). This is correct behavior: we simply don't have experimental data in that region.
**Warning signs:** Sharp cutoffs at edges of experimental trace.

### Pitfall 6: jcamp API Name
**What goes wrong:** Training data and some documentation references `jcamp_reader()` but the actual function in v1.2.2 is `jcamp_readfile()`.
**Why it happens:** API was renamed between versions.
**How to avoid:** Use `jcamp.jcamp_readfile(filepath)` -- verified working in v1.2.2.
**Warning signs:** `AttributeError: module 'jcamp' has no attribute 'jcamp_reader'`.

## Code Examples

### Verified: Complete NIST Fetch + Parse Pipeline
```python
# Verified working 2026-04-01 in mace4ir_v2 environment
import nistchempy as nist
import jcamp
import numpy as np
import tempfile
from pathlib import Path

def fetch_nist_ir(molecule_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Fetch gas-phase IR spectrum from NIST. Returns (wavenumbers, absorbance) or None."""
    search = nist.run_search(molecule_name, 'name')
    if not search.success or search.num_compounds == 0:
        return None

    compound = search.compounds[0]
    compound.get_ir_spectra()

    if not compound.ir_specs:
        return None

    # Find first gas-phase spectrum with absorbance or transmittance
    for spec in compound.ir_specs:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.jdx', delete=False) as f:
            f.write(spec.jdx_text)
            tmp_path = f.name

        data = jcamp.jcamp_readfile(tmp_path)
        Path(tmp_path).unlink()

        state = data.get('state', '').lower()
        if 'gas' not in state:
            continue

        x = data['x']  # wavenumbers (cm-1)
        y = data['y']
        yunits = data.get('yunits', '').upper()

        if 'TRANSMITTANCE' in yunits:
            # Convert T -> A: A = 2 - log10(T) if T in percent
            # or A = -log10(T) if T in [0, 1]
            if np.max(y) > 2:  # likely percent
                y = np.clip(y, 0.01, 100)
                y = 2.0 - np.log10(y)
            else:
                y = np.clip(y, 1e-6, 1.0)
                y = -np.log10(y)
        elif 'ABSORBANCE' in yunits:
            pass  # already absorbance
        else:
            # Cross-section or other: just use shape, normalize later
            pass

        # Ensure sorted by wavenumber (ascending)
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]

        # Normalize to [0, 1]
        y_max = np.max(y)
        if y_max > 0:
            y = y / y_max

        return x, y

    return None  # No gas-phase spectrum found
```

### Verified: jcamp Parsed Output Structure
```python
# Water spectrum (Sadtler source):
# data['x']:      array shape=(880,), range=[450.0, 3966.0]  (wavenumbers, cm-1)
# data['y']:      array shape=(880,), range=[0.0, 0.63]       (absorbance)
# data['xunits']: '1/CM'
# data['yunits']: 'ABSORBANCE'
# data['state']:  'gas'
# data['title']:  'Water'
# data['cas registry no']: '7732-18-5'

# Methanol spectrum (spec 0):
# data['yunits']: 'TRANSMITTANCE'
# data['state']:  'GAS (VAPOR)'
# npoints: 2741
```

### Verified: nistchempy Search API
```python
import nistchempy as nist

# Basic search
search = nist.run_search('water', 'name')
# search.success: True
# search.num_compounds: 1
# search.compounds: [NistCompound(ID=C7732185)]
# search.compound_ids: ['C7732185']

# Compound properties
compound = search.compounds[0]
# compound.ID: 'C7732185'
# compound.name: 'Water'
# compound.synonyms: ['Water vapor', 'Distilled water', 'Ice']

# IR spectra
compound.get_ir_spectra()
# compound.ir_specs: [Spectrum(...), Spectrum(...)]
# Each spec has: spec_type='IR', spec_idx='0', jdx_text='##TITLE=...'

# Save JDX to file
compound.ir_specs[0].save(name='water_ir.jdx', path_dir='/some/dir')
```

### Integration Point: Adding Experimental to plot_spectra_comparison
```python
# In analyze_spectra.py, plot_spectra_comparison() at line 563
# Current layout: DFT at y=0, ML at y=offset (1.5)
# Add experimental on BOTH panels for easy visual comparison

if experimental is not None:
    exp_interp = np.interp(
        self.freq_grid, experimental.wavenumbers, experimental.absorbance,
        left=0, right=0  # zero outside experimental range
    )
    # Plot on DFT panel
    ax1.plot(self.freq_grid, exp_interp, '--', color='black',
             linewidth=1.5, label='Experimental (NIST)', alpha=0.7, zorder=5)
    # Plot on ML panel
    ax1.plot(self.freq_grid, exp_interp + offset, '--', color='black',
             linewidth=1.5, alpha=0.7, zorder=5)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual NIST WebBook browsing | nistchempy automated search | nistchempy 1.0.0 (2023) | Enables programmatic access |
| jcamp_reader() | jcamp_readfile() | jcamp 1.2.x | API rename; old name no longer works |
| KDE broadening for IR | Lorentzian broadening | Phase 18 (2026-03) | Experimental overlay uses same Lorentzian infrastructure |

**Deprecated/outdated:**
- `jcamp.jcamp_reader()`: Renamed to `jcamp.jcamp_readfile()` in current version
- `jcamp.jcamp_read()`: Also exists but reads from file object, not path

## Open Questions

1. **Which NIST spectrum to prefer when multiple exist?**
   - What we know: Molecules often have 2-5 spectra from different providers. ABSORBANCE spectra are simpler to handle. PNNL cross-section data is highest resolution but unusual units.
   - What's unclear: Whether to always prefer ABSORBANCE over TRANSMITTANCE, or just take the first gas-phase one.
   - Recommendation: Take the first gas-phase spectrum. If it's TRANSMITTANCE, convert. This is simplest and matches the "best-effort" philosophy. If the first one produces bad results, the user can delete the cache and we'll try the next one in a future enhancement.

2. **Experimental trace on combined plot (all ML vs DFT)?**
   - What we know: `plot_combined_spectra()` stacks DFT + N ML traces with 1.5 offset each. Adding experimental to every panel would be visually cluttered.
   - What's unclear: Whether to show experimental only on the DFT (bottom) panel, or on all panels.
   - Recommendation: Show experimental only on the DFT (bottom) panel in the combined plot. Individual per-calculator plots can show it on both panels.

## Project Constraints (from CLAUDE.md)

- **Code quality:** `ruff check --fix && ruff format` (line length 100)
- **Type checking:** `ty check`
- **Plot conventions:** PNG, DPI=300, seaborn "colorblind" palette, Arial/Helvetica
- **Standard freq range:** 400-4000 cm-1
- **Environment:** Use `micromamba activate mace4ir_v2` (not `.venv`)
- **Dependencies:** nistchempy and jcamp are new deps, NOT yet in pyproject.toml -- must be added

## Sources

### Primary (HIGH confidence)
- nistchempy 1.0.5 -- installed and tested in mace4ir_v2; search, compound, spectrum APIs all verified working
- jcamp 1.2.2 -- installed and tested; `jcamp_readfile()` verified returning correct x/y arrays with metadata
- Existing codebase: `analyze_spectra.py`, `analysis_workflow.py` -- read and analyzed for integration points

### Secondary (MEDIUM confidence)
- [NistChemPy PyPI](https://pypi.org/project/nistchempy/) -- version and dependency info
- [NistChemPy GitHub](https://github.com/IvanChernyshov/NistChemPy) -- API structure and Spectrum class source
- [jcamp PyPI](https://pypi.org/project/jcamp/) -- version info
- [jcamp GitHub](https://github.com/nzhagen/jcamp) -- API and returned data structure

### Tertiary (LOW confidence)
- None -- all findings verified by direct testing

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- both libraries installed and tested end-to-end
- Architecture: HIGH -- integration points clearly identified in existing code
- Pitfalls: HIGH -- discovered by testing real molecules (water, methanol, cocaine)
- Y-unit handling: MEDIUM -- conversion logic is straightforward but edge cases in unusual NIST data possible

**Research date:** 2026-04-01
**Valid until:** 2026-05-01 (stable libraries, unlikely to change)
