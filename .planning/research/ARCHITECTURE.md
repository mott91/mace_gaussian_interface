# Architecture Patterns: v1.2 Analysis Quality Overhaul

**Domain:** Integration of new analysis features into existing MACE-Gaussian IR spectroscopy pipeline
**Researched:** 2026-03-28
**Confidence:** HIGH (based on direct codebase analysis, not external sources)

## Current Architecture Summary

The existing pipeline has three layers:

1. **Pipeline layer** -- `workflow.py` (single-molecule orchestrator), `batch.py` (multi-molecule with manifest), `slurm.py` (remote DFT offload)
2. **Analysis layer** -- `analysis/analysis_workflow.py` (ComparisonWorkflow orchestrator), `analysis/analyze_spectra.py` (SpectrumAnalyzer with broadening/regression/plotting), `analysis/mode_matching.py` (Hungarian algorithm eigenvector matching), `analysis/html_report_generator.py`, `analysis/batch_report.py`
3. **Data layer** -- `results.json` per calculation, `.fchk` checkpoint files, `batch_manifest.json`

Key data types: `SpectrumData` (frequencies, intensities, labels, mode_ids), `ComparisonMetrics` (MAE, RMSE, R2, match stats).

## Recommended Architecture for New Features

### Component Map: New vs Modified

| Component | Status | Location | Rationale |
|-----------|--------|----------|-----------|
| Lorentzian broadening | **MODIFY** `analyze_spectra.py` | `SpectrumAnalyzer` | Parallel to existing `broaden_spectrum()` Gaussian method |
| NIST fetcher | **NEW** `analysis/nist_spectra.py` | New module | Separate concern: network I/O, caching, JCAMP-DX parsing |
| NIST cache store | **NEW** `~/.cache/mace_gaussian/nist/` | User-level cache | Avoid polluting project results directories |
| VPT2 engine | **NEW** `vpt2/` subpackage | Top-level subpackage | Parallel pathway to `gaussian/`, not a sub-feature of analysis |
| Degenerate mode handling | **MODIFY** `analysis/mode_matching.py` | Extend `match_modes()` | Add subspace detection + overlap, keep API compatible |
| Zero-intensity filtering | **MODIFY** `analysis/analyze_spectra.py` | `calculate_metrics()` | Filter before intensity regression, not frequency regression |
| Timing instrumentation | **NEW** `utils/timing.py` | Utility module | Context manager + decorator pattern for reuse |
| Timing in results | **MODIFY** `workflow.py`, `batch.py` | Extend results dict | Add per-stage timing breakdown to existing runtime field |
| Early SLURM | **MODIFY** `batch.py` | Restructure main loop | Submit DFT after geom-opt, before ML combos (not after all) |
| Experimental overlay in plots | **MODIFY** `analyze_spectra.py` | `plot_spectra_comparison()` | Add optional experimental trace |
| Report: timing section | **MODIFY** `html_report_generator.py`, `batch_report.py` | New report section | Cost-accuracy tradeoff table/plot |
| Report: experimental section | **MODIFY** `html_report_generator.py` | New report section | Experimental overlay + peak comparison |

### Component Boundaries

| Component | Responsibility | Communicates With |
|-----------|---------------|-------------------|
| `analysis/nist_spectra.py` | Fetch, cache, parse JCAMP-DX experimental spectra | `analysis_workflow.py` (provides SpectrumData), `html_report_generator.py` (provides plot data) |
| `utils/timing.py` | `@timed` decorator and `TimingContext` context manager, timing aggregation | `workflow.py`, `batch.py` (instruments stages), `analysis_workflow.py` (reads timing from results) |
| `vpt2/__init__.py` | PyVPT2 engine wrapper: takes ASE atoms + calculator, returns frequencies/intensities | `workflow.py` (alternative to Gaussian anharmonic stage), `cli.py` (--vpt2-engine flag) |
| `analysis/mode_matching.py` (modified) | Degenerate subspace detection + group overlap scoring | `analysis_workflow.py` (receives overlap scores), `html_report_generator.py` (renders group confidence) |

### Data Flow: How New Features Integrate

```
Existing flow (unchanged):
  Input (.xyz) -> Geom Opt -> Freq Calc (Gaussian+ZMQ) -> Anharmonic -> Parse -> results.json

New timing flow (wraps existing):
  Each stage wrapped with TimingContext -> timing dict added to results.json
  results.json gains: {"timing": {"geom_opt_s": X, "freq_calc_s": Y, "gaussian_anharmonic_s": Z}}

New VPT2 flow (parallel alternative):
  Geom Opt -> [if --vpt2-engine=pyvpt2] -> vpt2.run(atoms, calculator) -> frequencies/intensities
  Output: same results.json schema (frequencies.anharmonic, frequencies.overtones)
  Skip: Gaussian ZMQ bridge, .gjf/.log parsing entirely

New NIST flow (analysis-time addon):
  analysis_workflow.py -> nist_spectra.fetch(molecule_name) -> SpectrumData
  Cache: ~/.cache/mace_gaussian/nist/{cas_number}.jdx
  Overlay: passed to plot_spectra_comparison() as optional experimental_spectrum param

New degenerate mode flow (enhanced existing):
  match_modes() -> detect_degenerate_groups(freqs, threshold=5.0) -> list[list[int]]
  For each group: compute_subspace_overlap(modes_calc[group], modes_ref[group]) -> float
  Return: matches dict unchanged + degenerate_groups metadata dict

Early SLURM flow (restructured batch loop):
  For each molecule: geom_opt -> submit_dft_to_slurm (non-blocking)
  Then: run all ML combos locally (while DFT runs on cluster)
  After all molecules: poll_and_retrieve_dft (blocking)
  Before: geom_opt ALL -> ML combos ALL -> submit_dft ALL -> poll
  After:  geom_opt -> submit_dft -> ML combos (per molecule, pipelined)
```

## Detailed Integration Points

### 1. Lorentzian Broadening -- MODIFY `analyze_spectra.py`

**Where:** Add `broaden_spectrum_lorentzian()` method to `SpectrumAnalyzer` class, parallel to existing `broaden_spectrum()` (which is Gaussian broadening).

**Why not a new module:** The broadening logic is 10 lines. It operates on the same `SpectrumData` and `freq_grid` already owned by `SpectrumAnalyzer`. Extracting it would create unnecessary indirection.

**Implementation:**
```python
def broaden_spectrum_lorentzian(self, spectrum: SpectrumData, hwhm: float = 10.0) -> np.ndarray:
    """Apply Lorentzian broadening: L(v) = (gamma/pi) / ((v - v0)^2 + gamma^2)"""
    broadened = np.zeros_like(self.freq_grid)
    for freq, intensity in zip(spectrum.frequencies, spectrum.intensities):
        broadened += intensity * (hwhm / np.pi) / ((self.freq_grid - freq)**2 + hwhm**2)
    return broadened
```

**API change:** Add `broadening_type: str = "gaussian"` parameter to `SpectrumAnalyzer.__init__()` (or to `plot_spectra_comparison()`). Accepted values: `"gaussian"`, `"lorentzian"`, `"voigt"`. Default stays `"gaussian"` for backward compatibility.

**Plot integration:** `plot_spectra_comparison()` already calls `self.broaden_spectrum()`. Add a dispatch: if `self.broadening_type == "lorentzian"`, call `broaden_spectrum_lorentzian()` instead. No changes to callers.

**Voigt option (optional, low priority):** Voigt = convolution of Gaussian and Lorentzian. Use `scipy.special.voigt_profile`. Adds one dependency already in the environment. Most realistic for gas-phase IR but Lorentzian alone is standard for computational comparison.

### 2. VPT2 Alternative -- NEW `vpt2/` Subpackage

**Where:** New top-level subpackage `mace_gaussian/vpt2/` -- NOT inside `analysis/` and NOT as a calculator type.

**Why a separate subpackage:** VPT2 is a complete alternative anharmonic pathway. It replaces the Gaussian ZMQ bridge + anharmonic calculation stage entirely. It has its own dependencies (PyVPT2/QCEngine), its own execution flow, and its own failure modes. Stuffing it into `workflow.py` or `calculators/` would blur the boundary between "how we talk to Gaussian" and "how we compute anharmonic frequencies."

**Structure:**
```
mace_gaussian/vpt2/
    __init__.py          # run_vpt2(atoms, energy_calc, dipole_calc) -> results dict
    pyvpt2_engine.py     # PyVPT2 + QCEngine adapter
    qcengine_harness.py  # QCEngine ASE harness (maps MACE calculator to QCEngine protocol)
```

**Integration with workflow.py:** Add a branch in `run_frequency_calculation()`:
```python
if vpt2_engine == "pyvpt2":
    from .vpt2 import run_vpt2
    results = run_vpt2(atoms, energy_calc, dipole_calc)
else:
    # Existing Gaussian ZMQ path
    results = run_gaussian_with_zmq(...)
```

**Output schema:** VPT2 results MUST produce the same `results.json` schema (frequencies.harmonic, frequencies.anharmonic, frequencies.overtones, frequencies.combination_bands) so that `analysis_workflow.py` works unchanged.

**CLI flag:** `--vpt2-engine gaussian|pyvpt2` on `mace-gaussian run` and `mace-gaussian batch`. Default: `gaussian`.

**Important:** PyVPT2 is research software. Expect integration friction. Plan a prototype phase that validates it works on water before wiring into the full pipeline. This is the highest-risk feature in v1.2.

### 3. NIST Spectra -- NEW `analysis/nist_spectra.py`

**Where:** New module `mace_gaussian/analysis/nist_spectra.py`.

**Why a new module:** This is a distinct concern (network I/O, caching, format parsing) that should not pollute `analyze_spectra.py`. But it lives in `analysis/` because it produces `SpectrumData` consumed by the analysis workflow.

**Structure:**
```python
# nist_spectra.py
class NISTFetcher:
    def __init__(self, cache_dir: Path = None):
        """Default cache: ~/.cache/mace_gaussian/nist/"""

    def fetch(self, molecule_name: str, cas_number: str = None) -> SpectrumData | None:
        """Fetch gas-phase IR spectrum from NIST WebBook. Returns None if not available."""

    def _resolve_cas(self, molecule_name: str) -> str | None:
        """Look up CAS number from molecule name via NIST search."""

    def _download_jcamp(self, cas_number: str) -> Path | None:
        """Download JCAMP-DX file, cache locally. Return cached path."""

    def _parse_jcamp(self, jdx_path: Path) -> SpectrumData:
        """Parse JCAMP-DX into SpectrumData (transmittance -> absorbance conversion)."""
```

**Cache strategy:** `~/.cache/mace_gaussian/nist/{cas_number}.jdx` for raw JCAMP-DX files. Separate from project `comparison_results/` because experimental data is molecule-global, not run-specific. Check cache before network request.

**Integration with analysis_workflow.py:**
```python
# In ComparisonWorkflow.run_full_comparison():
from .nist_spectra import NISTFetcher
fetcher = NISTFetcher()
experimental = fetcher.fetch(self.molecule_name)
if experimental:
    # Pass to plot functions as optional overlay
    self.analyzer.plot_spectra_comparison(..., experimental_spectrum=experimental)
```

**Integration with reports:** Add an "Experimental Comparison" section to `html_report_generator.py`. Show: (a) overlaid plot, (b) matched peak positions table, (c) note on data source (gas-phase vs condensed-phase).

**JCAMP-DX parsing:** The format is well-documented (IUPAC standard). Key fields: `##XUNITS=`, `##YUNITS=`, `##XYDATA=`. NIST uses transmittance (%T); convert to absorbance via `A = -log10(T/100)` for overlay with computed absorption spectra. Use no external JCAMP library -- the format is simple enough to parse with ~50 lines of Python. Avoid adding a dependency for this.

**Failure modes:** NIST may not have gas-phase IR for all molecules. The fetcher MUST return `None` gracefully so the analysis pipeline continues without experimental data. Never block analysis on NIST availability.

### 4. Degenerate Mode Handling -- MODIFY `analysis/mode_matching.py`

**Where:** Add `detect_degenerate_groups()` and `compute_subspace_overlap()` functions to `mode_matching.py`. Modify `match_modes()` return type to optionally include degeneracy metadata.

**Why in mode_matching.py:** Degeneracy is a property of the normal modes and their overlap matrix. The detection and subspace overlap computation belong with the existing eigenvector matching logic.

**Implementation approach:**

```python
def detect_degenerate_groups(
    frequencies: np.ndarray, threshold: float = 5.0
) -> list[list[int]]:
    """Group mode indices where frequencies are within threshold cm-1.

    Returns list of groups, each group is a list of mode indices.
    Singletons (non-degenerate) are included as 1-element lists.
    """

def compute_subspace_overlap(
    modes_a: np.ndarray, modes_b: np.ndarray
) -> float:
    """Compute subspace overlap for a degenerate group.

    Given modes_a shape (k, n_atoms, 3) and modes_b shape (k, n_atoms, 3),
    compute trace(M^T M) / k where M[i,j] = dot(a_i, b_j) / (||a_i|| ||b_j||).
    Returns value in [0, 1]. 1.0 = subspaces are identical.
    """
```

**API compatibility:** `match_modes()` return type stays `dict[int, tuple[int | None, float]]`. Add a new function `match_modes_with_degeneracy()` that returns an extended result including group metadata. This keeps backward compatibility with all existing callers.

**Where zero-intensity filtering goes:** In `analyze_spectra.py`'s `calculate_metrics()` method. After mode matching (which is intensity-agnostic), filter pairs where `dft_intensity < min_threshold` before computing intensity R2/RMSE. Frequency metrics use ALL matched pairs. Add `min_dft_intensity: float = 0.1` parameter.

### 5. Timing Instrumentation -- NEW `utils/timing.py`

**Where:** New module `mace_gaussian/utils/timing.py`.

**Why a utility module:** Timing is cross-cutting. Used in `workflow.py` (per-stage), `batch.py` (per-molecule), and read by `analysis/` (for reports). A shared utility avoids duplicating `time.time()` patterns.

**Pattern: context manager + decorator.**

```python
import time
from contextlib import contextmanager
from functools import wraps

@contextmanager
def timed(name: str, timings: dict):
    """Context manager that records wall-clock time into timings[name]."""
    start = time.monotonic()
    yield
    timings[name] = round(time.monotonic() - start, 2)

def timed_decorator(name: str):
    """Decorator version. Stores timing in return value's 'timing' key if dict."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time.monotonic()
            result = func(*args, **kwargs)
            elapsed = round(time.monotonic() - start, 2)
            if isinstance(result, dict):
                result.setdefault("timing", {})[name] = elapsed
            return result
        return wrapper
    return decorator
```

**Instrumentation points in workflow.py:**
```python
timings = {}
with timed("geom_opt", timings):
    optimized_atoms = run_geometry_optimization(...)
with timed("freq_calc", timings):
    run_gaussian_with_zmq(...)
# Store in results.json:
results["timing"] = timings
```

**Existing timing:** `workflow.py` already records `runtime_s` via `time.time()` (lines 412-434, 479-615). The new system should extend this, not replace it. Add per-stage breakdown alongside the existing total.

**Batch report integration:** `batch_report.py` already aggregates per-molecule data. Add a timing column to the leaderboard table and a cost-accuracy scatter plot (RMSE vs wall-clock time).

### 6. Early SLURM Submission -- MODIFY `batch.py`

**Current structure (sequential, DFT submitted AFTER all molecules):**
```
for molecule in molecules:
    geom_opt(molecule)
    if not dft_on_cluster:
        run_dft_local(molecule)   # blocking
    run_ml_combos(molecule)       # blocking

# After ALL molecules done:
if dft_on_cluster:
    submit_all_dft_to_slurm()     # <-- too late!
    poll_until_done()
    retrieve_results()
```

**Problem:** DFT jobs (hours each) are only submitted after ALL local ML calculations complete. If the batch has 25 molecules with 4 ML combos each, DFT waits for ~100 ML calculations before even being submitted.

**Proposed structure (pipelined, DFT submitted per-molecule after geom-opt):**
```
pending_slurm_jobs = {}

for molecule in molecules:
    geom_opt(molecule)

    if dft_on_cluster:
        job_id = submit_single_dft(molecule)   # non-blocking, ~1 second
        pending_slurm_jobs[molecule] = job_id

    run_ml_combos(molecule)   # blocking, but DFT already running on cluster

# After all molecules:
if pending_slurm_jobs:
    poll_and_retrieve(pending_slurm_jobs)
```

**Key changes to batch.py:**
1. Extract `submit_single_dft()` from the bulk submission block (lines 322-416). The existing `submit_dft_jobs()` in `slurm.py` already handles a list -- call it with a single-element list per molecule.
2. Move the `if dft_on_cluster: pass` block (line 220-222) to actually submit instead of deferring.
3. Keep the poll/retrieve block at the end (lines 389-415) but have it operate on the already-submitted jobs.
4. Manifest tracking per-molecule SLURM job_id stays the same.

**Risk:** The existing batch flow defers DFT to guarantee all geom-opts succeed first. With per-molecule submission, a failed geom-opt after DFT submission wastes a cluster job. Mitigation: only submit DFT if geom-opt succeeded (which is already the flow -- DFT is inside the try block after geom-opt).

**slurm.py changes:** Minimal. `submit_dft_jobs()` already accepts a list of molecules. Calling it per-molecule (list of 1) works. No API change needed.

## Patterns to Follow

### Pattern 1: Optional Enhancement (for NIST, timing, Lorentzian)
**What:** New features are always optional and additive. Existing behavior is the default.
**When:** Adding features that may fail (network), are not yet validated (VPT2), or change output format (Lorentzian).
**Rule:** Every new parameter has a default that preserves current behavior. `experimental_spectrum=None`, `broadening_type="gaussian"`, `min_dft_intensity=0.0` (no filtering).

### Pattern 2: Same Schema, Different Engine (for VPT2)
**What:** Alternative computation engines produce identical output schemas so downstream analysis is unchanged.
**When:** Replacing one computational method with another (Gaussian VPT2 -> PyVPT2).
**Rule:** The `results.json` schema is the contract. New engines must produce compliant output. Test with schema validation, not just "it runs."

### Pattern 3: Lazy Import (existing pattern in `analysis/__init__.py`)
**What:** Heavy dependencies (matplotlib, PyVPT2, network fetchers) imported only when called.
**When:** Modules that add optional dependencies or significant import time.
**Rule:** Follow existing `__init__.py` lazy wrapper pattern for new analysis modules.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Mixing Network I/O with Analysis Logic
**What:** Putting NIST HTTP requests inside `analyze_spectra.py`.
**Why bad:** Makes analysis depend on network availability. Breaks tests. Couples concerns.
**Instead:** `nist_spectra.py` handles all I/O. Returns `SpectrumData` or `None`. Analysis code never knows about HTTP.

### Anti-Pattern 2: Modifying match_modes() Return Type
**What:** Changing the return dict type to include degeneracy info, breaking all callers.
**Why bad:** `match_modes()` is called in `analysis_workflow.py` (line 361), `mode_matching.py` main block, and potentially tests. Changing its return breaks them all.
**Instead:** Add `match_modes_with_degeneracy()` as a new function. Or add degeneracy info as an optional second return value via a `return_degeneracy=False` parameter.

### Anti-Pattern 3: Tight VPT2-Gaussian Coupling
**What:** Threading PyVPT2 configuration through Gaussian-specific code paths.
**Why bad:** VPT2 via PyVPT2 has nothing to do with Gaussian. It should not import from `gaussian/` or know about ZMQ.
**Instead:** `vpt2/` subpackage is self-contained. `workflow.py` dispatches to either `gaussian/` or `vpt2/` based on `--vpt2-engine`.

## Build Order (Dependency-Aware)

The features have minimal cross-dependencies. Recommended order based on risk and value:

1. **Timing instrumentation** (`utils/timing.py` + `workflow.py` changes) -- Zero risk, immediate value, enables cost-accuracy analysis. No new dependencies.

2. **Lorentzian broadening** (`analyze_spectra.py`) -- 10-line implementation, improves all existing spectra plots immediately. No new dependencies.

3. **Zero-intensity filtering** (`analyze_spectra.py`) -- Small change to `calculate_metrics()`, fixes methane R2=0.34 intensity issue immediately.

4. **Degenerate mode handling** (`mode_matching.py`) -- Medium complexity, well-scoped. Requires testing with methane (T2 3-fold) and other symmetric molecules.

5. **NIST spectra** (`analysis/nist_spectra.py`) -- Independent module, can be developed in parallel. Needs network access for testing. No dependencies on 1-4.

6. **Early SLURM submission** (`batch.py`) -- Restructures the batch loop. Low risk (no new logic, just reordering existing calls). Biggest value for the benchmark campaign.

7. **VPT2 integration** (`vpt2/`) -- Highest risk, new external dependency (PyVPT2), research-grade software. Build last so other features are stable. May slip to v1.3 if PyVPT2 integration is harder than expected.

Items 1-3 can be done in a single phase. Items 4-6 can be parallelized. Item 7 is a standalone phase.

## Scalability Considerations

| Concern | Current (5 molecules) | At 25 molecules (benchmark) | At 100+ molecules |
|---------|----------------------|----------------------------|-------------------|
| NIST fetching | Negligible | ~25 HTTP requests, cache helps | Rate limit risk; add delay between requests |
| Analysis runtime | Seconds | Minutes (fine) | May need parallel analysis |
| Batch manifest size | Trivial | ~500 entries, fine | Still fine (JSON) |
| SLURM job submission | N/A or serial | 25 simultaneous jobs acceptable | May need job array batching |
| Timing data in results | ~3 fields | Same | Same |
| Report generation | Seconds | batch_report already handles | Already handles |

## Sources

- Direct codebase analysis: `mace_gaussian/analysis/analyze_spectra.py`, `analysis_workflow.py`, `mode_matching.py`, `batch.py`, `slurm.py`, `workflow.py`
- Existing todo files: `.planning/todos/pending/2026-03-26-*.md`, `.planning/todos/pending/2026-03-10-*.md`
- JCAMP-DX format: IUPAC standard (training data, HIGH confidence for format simplicity)
- Lorentzian broadening: standard computational spectroscopy (HIGH confidence)
- PyVPT2: as documented in the VPT2 todo file (MEDIUM confidence for integration feasibility -- research software)
