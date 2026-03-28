# Phase 18: Spectral Quality Foundations - Research

**Researched:** 2026-03-28
**Domain:** IR spectral broadening (Lorentzian line shapes), intensity filtering for regression metrics
**Confidence:** HIGH

## Summary

This phase makes two surgical changes to the analysis layer: (1) replace the Gaussian broadening kernel with a Lorentzian kernel in `broaden_spectrum()`, and (2) filter zero-intensity modes from intensity regression metrics. Both changes are well-contained -- the broadening swap touches one method and its initialization, while intensity filtering touches one method (`compare_spectra`) plus the intensity regression plot.

The Lorentzian line shape formula is standard physics (Cauchy distribution), requires no new dependencies, and the existing code structure (loop over frequency/intensity pairs, accumulate on grid) is directly reusable. The FWHM parameter is already accepted by `SpectrumAnalyzer.__init__()` but hardcoded at the `ComparisonWorkflow` level -- wiring it to CLI requires adding `--fwhm` to both entry points and threading it through `analyze_molecule()` / `analyze_molecule_harmonic()`.

**Primary recommendation:** Implement Lorentzian kernel swap first, then intensity filtering, then CLI wiring. Each is independently testable.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Show Lorentzian broadened spectrum only -- no stick spectrum for now (SPEC-03 stick display deferred)
- **D-02:** Keep existing DFT-above / ML-below offset layout, just swap Gaussian broadening for Lorentzian in `broaden_spectrum()`
- **D-03:** Filter mode pairs from intensity regression metrics when DFT reference intensity < 0.1 km/mol
- **D-04:** Frequency metrics (R-squared, RMSE for frequencies) still include all modes regardless of intensity
- **D-05:** Lorentzian spectrum plots still show all modes -- filtering is for statistics only, not visual display
- **D-06:** Default FWHM changed from 8 cm-1 to 10 cm-1 (standard in computational spectroscopy literature)
- **D-07:** FWHM configurable via `--fwhm` CLI argument on analysis scripts (run_analysis.py, run_analysis_harmonic.py). Already has `bandwidth_fwhm` parameter in SpectralAnalyzer -- just wire to CLI.

### Claude's Discretion
- Lorentzian implementation details (normalization, numerical cutoff)
- How to label/annotate filtered modes in report text
- Whether to mention filtering count in report methodology section

### Deferred Ideas (OUT OF SCOPE)
- Stick spectrum display (SPEC-03 partial) -- add later, possibly Phase 23
- PIPE-02: mace_polar dipole calculator investigation -- separate todo, not Phase 18
- Interactive Plotly spectrum viewer -- future phase
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SPEC-01 | Simulated IR spectra use Lorentzian line shapes instead of Gaussian KDE, with configurable FWHM (default 10 cm-1) | Lorentzian formula documented below; swap kernel in `broaden_spectrum()`, change default in `__init__`, wire `--fwhm` to CLI |
| SPEC-02 | Modes with IR intensity below 0.1 km/mol are filtered from intensity regression metrics (but retained in frequency metrics) | Filter in `compare_spectra()` before intensity `linregress`; DFT intensity column already available as `dft_int` array |
| SPEC-03 | Analysis reports display both stick spectrum and broadened Lorentzian spectrum | **DEFERRED by user decision D-01** -- Lorentzian-only display for now |
| PIPE-02 | Mace_polar dipole calculator failure mode investigated and either fixed or documented with skip logic | **DEFERRED by user request** -- becomes standalone todo |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- **Code quality:** `ruff check --fix && ruff format` (line length 100), `ty check` for type checking
- **Plots:** PNG format, DPI=300, seaborn "colorblind" palette, Arial/Helvetica font. Always include R-squared and RMSE in regression plots. Standard range 400-4000 cm-1.
- **Environment:** Use `micromamba activate mace4ir_v2` (not `.venv`)
- **Testing:** pytest with markers `gpu`, `gaussian`, `slow`; strict markers enabled

## Standard Stack

No new dependencies required. All changes use numpy (already imported in analyze_spectra.py).

### Core (already in project)
| Library | Purpose | Role in This Phase |
|---------|---------|-------------------|
| numpy | Numerical arrays | Lorentzian kernel computation on frequency grid |
| scipy.stats.linregress | Linear regression | Already used for intensity R-squared; filtering applied before call |
| matplotlib | Plotting | Already used; no changes to plot structure needed |
| argparse | CLI parsing | Will be added to entry points (currently use raw `sys.argv`) |

## Architecture Patterns

### Files to Modify

```
mace_gaussian/analysis/
  analyze_spectra.py          # broaden_spectrum(), __init__(), compare_spectra()
  analysis_workflow.py        # ComparisonWorkflow.__init__(), analyze_molecule(),
                              #   analyze_molecule_harmonic(), run_analysis_main(),
                              #   run_analysis_harmonic_main()
  html_report_generator.py    # Methodology text (line 485: "8 cm^-1" -> dynamic)
```

### Pattern 1: Lorentzian Kernel Swap

**What:** Replace Gaussian broadening `exp(-B * (v - v0)^2)` with Lorentzian `(gamma^2) / ((v - v0)^2 + gamma^2)` where `gamma = FWHM / 2`.

**Current code** (analyze_spectra.py line 324-348):
```python
def broaden_spectrum(self, spectrum: SpectrumData) -> np.ndarray:
    broadened = np.zeros_like(self.freq_grid)
    for freq, intensity in zip(spectrum.frequencies, spectrum.intensities):
        argument = -self.broad_param * (self.freq_grid - freq) ** 2
        mask = argument > -50
        broadened[mask] += intensity * np.exp(argument[mask])
    return broadened
```

**New code pattern:**
```python
def broaden_spectrum(self, spectrum: SpectrumData) -> np.ndarray:
    broadened = np.zeros_like(self.freq_grid)
    gamma = self.bandwidth_fwhm / 2.0
    gamma_sq = gamma ** 2
    for freq, intensity in zip(spectrum.frequencies, spectrum.intensities):
        lorentzian = gamma_sq / ((self.freq_grid - freq) ** 2 + gamma_sq)
        broadened += intensity * lorentzian
    return broadened
```

**Key details:**
- The Lorentzian naturally decays as `1/x^2` at large distances (much slower than Gaussian), so the `mask` optimization for underflow is no longer needed -- Lorentzian values at 2000 cm-1 away from center are ~6e-6 * gamma^2, negligible but not underflow.
- Optional: add a cutoff at ~100*FWHM for performance, but with typical grid sizes (7200 points for 400-4000 at 0.5 step) the full computation is fast enough.
- The `self.broad_param` precomputed in `__init__` is Gaussian-specific. Replace it with storing `self.bandwidth_fwhm` directly (already stored at line 87).
- Remove or simplify the `sigma`/`broad_param` computation in `__init__` (lines 96-100) since Lorentzian only needs `gamma = FWHM / 2`.

**Normalization choice (Claude's discretion):** Use unnormalized Lorentzian (peak = intensity at center). This matches the current Gaussian behavior where peak height equals intensity. Area-normalized Lorentzian would divide by `pi * gamma`, but that changes the visual scale and is unnecessary for comparison plots where both DFT and ML use the same broadening.

### Pattern 2: Intensity Filtering in compare_spectra()

**What:** Before computing intensity regression metrics, filter out mode pairs where DFT intensity < 0.1 km/mol.

**Current code** (analyze_spectra.py line 526-535):
```python
# Intensity metrics
int_errors = ml_int - dft_int
mae_intensity = np.mean(np.abs(int_errors))

if len(dft_int) > 0 and np.any(dft_int > 0):
    _, _, r_value_int, _, _ = linregress(dft_int, ml_int)
    r2_intensity = r_value_int**2
```

**New pattern:**
```python
# Intensity metrics -- filter out near-zero DFT modes (SPEC-02)
int_mask = dft_int >= 0.1  # km/mol threshold
dft_int_filtered = dft_int[int_mask]
ml_int_filtered = ml_int[int_mask]

if len(dft_int_filtered) > 1:
    int_errors = ml_int_filtered - dft_int_filtered
    mae_intensity = np.mean(np.abs(int_errors))
    _, _, r_value_int, _, _ = linregress(dft_int_filtered, ml_int_filtered)
    r2_intensity = r_value_int**2
else:
    mae_intensity = 0.0
    r2_intensity = 0.0
```

**Important:** `linregress` requires at least 2 data points. The current code checks `np.any(dft_int > 0)` but does not check length >= 2. The new code should guard against this.

**Where filtering does NOT apply (D-04, D-05):**
- Frequency metrics (lines 517-524): use all modes, no filter
- Spectrum plots: show all modes visually
- `num_peaks`, `num_matched`, etc.: count all modes

**New field suggestion:** Add `num_intensity_filtered: int` to `ComparisonMetrics` so the report can state "N modes excluded from intensity regression."

### Pattern 3: CLI Argument Wiring

**Current state:** `run_analysis_main()` and `run_analysis_harmonic_main()` use raw `sys.argv[1]` with no argparse. The `ComparisonWorkflow.__init__()` hardcodes `bandwidth_fwhm=8.0`.

**New pattern:** Convert both entry points to use `argparse.ArgumentParser` with:
- Positional `molecule_name`
- `--fwhm` with default 10.0, type float
- Thread FWHM through `analyze_molecule(fwhm=...)` -> `ComparisonWorkflow(bandwidth_fwhm=...)` -> `SpectrumAnalyzer(bandwidth_fwhm=...)`

This requires:
1. Add `bandwidth_fwhm` parameter to `ComparisonWorkflow.__init__()`
2. Add `fwhm` parameter to `analyze_molecule()` and `analyze_molecule_harmonic()`
3. Replace `sys.argv` parsing with argparse in both `run_*_main()` functions

### Pattern 4: HTML Report Methodology Text Update

**Current:** Line 485 in html_report_generator.py says "Gaussian convolution with 8 cm^-1 FWHM"

**New:** Change to "Lorentzian broadening with {fwhm} cm^-1 FWHM" -- the FWHM value should be passed from the workflow to the report generator.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Lorentzian function | Custom approximation | Direct formula `gamma^2 / (dx^2 + gamma^2)` | One-liner, no library needed, numerically stable |
| CLI argument parsing | Manual sys.argv slicing | `argparse.ArgumentParser` | Already in stdlib, handles help text and type validation |

## Common Pitfalls

### Pitfall 1: Lorentzian Long Tails
**What goes wrong:** Lorentzian has much heavier tails than Gaussian. A peak at 3000 cm-1 with FWHM=10 still contributes visibly at 400 cm-1 (about 0.001% of peak height). With many modes, these tails accumulate and create a nonzero baseline.
**Why it happens:** Lorentzian decays as 1/x^2 vs Gaussian exp(-x^2).
**How to avoid:** This is physically correct behavior for natural line shapes. Do NOT try to fix it -- the baseline shift is real and expected. If it becomes visually distracting, a future phase could subtract it, but for now it matches the physics.
**Warning signs:** Spectrum baseline is not zero between peaks. This is expected.

### Pitfall 2: Intensity Filter Edge Case with Few Modes
**What goes wrong:** If a molecule has very few modes (e.g., water with 3 modes) and most have low DFT intensity, filtering could leave 0 or 1 modes, making `linregress` fail.
**Why it happens:** `linregress` needs >= 2 points.
**How to avoid:** Guard with `len(dft_int_filtered) > 1` check before calling linregress. Return 0.0 for R-squared if insufficient data, and log a warning.
**Warning signs:** `ValueError` from scipy linregress.

### Pitfall 3: Forgetting to Update Both Entry Points
**What goes wrong:** Adding `--fwhm` to `run_analysis_main()` but forgetting `run_analysis_harmonic_main()`, or updating one `analyze_molecule` function but not the other.
**Why it happens:** There are two parallel paths (anharmonic and harmonic) with nearly identical boilerplate.
**How to avoid:** Change list: both `run_*_main()`, both `analyze_molecule*()`, `ComparisonWorkflow.__init__()`, and `SpectrumAnalyzer.__init__()` default value.

### Pitfall 4: HTML Report Hardcoded FWHM String
**What goes wrong:** The methodology text in the HTML report says "8 cm^-1 FWHM" even after default changes to 10.
**Why it happens:** String literal in `html_report_generator.py` line 485 is not parameterized.
**How to avoid:** Pass FWHM value to the report generator and interpolate into the methodology text.

## Code Examples

### Lorentzian Broadening (verified formula)
```python
# Lorentzian (Cauchy) line shape: L(v) = I * gamma^2 / ((v - v0)^2 + gamma^2)
# At v = v0: L = I (peak equals intensity)
# At v = v0 +/- gamma: L = I/2 (half-maximum, confirming FWHM = 2*gamma)
gamma = fwhm / 2.0
gamma_sq = gamma ** 2
lorentzian = gamma_sq / ((freq_grid - freq_center) ** 2 + gamma_sq)
broadened += intensity * lorentzian
```

### Intensity Threshold Filter
```python
# Filter near-zero DFT intensity modes from regression (SPEC-02)
INTENSITY_THRESHOLD = 0.1  # km/mol
int_mask = dft_int >= INTENSITY_THRESHOLD
n_filtered = (~int_mask).sum()
if n_filtered > 0:
    logger.info(f"Filtered {n_filtered} modes with DFT intensity < {INTENSITY_THRESHOLD} km/mol from intensity regression")
```

### Argparse for Analysis Entry Points
```python
import argparse

def run_analysis_main() -> None:
    parser = argparse.ArgumentParser(description="IR Spectral Analysis (anharmonic)")
    parser.add_argument("molecule_name", help="Name of molecule to analyze")
    parser.add_argument("--fwhm", type=float, default=10.0,
                        help="Full width at half maximum for Lorentzian broadening (cm-1)")
    args = parser.parse_args()
    results = analyze_molecule(args.molecule_name, fwhm=args.fwhm)
```

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest >= 7.0.0 |
| Config file | `pyproject.toml` [tool.pytest.ini_options] |
| Quick run command | `pytest tests/ -x -q --no-header` |
| Full suite command | `pytest tests/ -ra` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SPEC-01 | Lorentzian broadening produces correct peak shape (FWHM check) | unit | `pytest tests/test_spectral_broadening.py::test_lorentzian_fwhm -x` | Wave 0 |
| SPEC-01 | Lorentzian peak height equals input intensity | unit | `pytest tests/test_spectral_broadening.py::test_lorentzian_peak_height -x` | Wave 0 |
| SPEC-01 | Default FWHM is 10 cm-1 | unit | `pytest tests/test_spectral_broadening.py::test_default_fwhm -x` | Wave 0 |
| SPEC-02 | Intensity regression excludes modes with DFT intensity < 0.1 | unit | `pytest tests/test_spectral_broadening.py::test_intensity_filter_threshold -x` | Wave 0 |
| SPEC-02 | Frequency metrics include all modes (no filtering) | unit | `pytest tests/test_spectral_broadening.py::test_frequency_metrics_unfiltered -x` | Wave 0 |
| SPEC-02 | Edge case: all modes below threshold returns r2=0 | unit | `pytest tests/test_spectral_broadening.py::test_intensity_filter_all_below -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest tests/test_spectral_broadening.py -x -q`
- **Per wave merge:** `pytest tests/ -ra`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/test_spectral_broadening.py` -- covers SPEC-01, SPEC-02 (new file)
- [ ] No framework install needed -- pytest already configured

## Open Questions

1. **Report generator FWHM parameterization**
   - What we know: Line 485 of html_report_generator.py has a hardcoded string "8 cm^-1 FWHM"
   - What's unclear: How much of the report generator receives context from the workflow (need to check if a metadata dict is passed)
   - Recommendation: Thread FWHM through workflow -> report generator; if too invasive, at minimum update the hardcoded string to 10

## Sources

### Primary (HIGH confidence)
- Direct code reading of `mace_gaussian/analysis/analyze_spectra.py` -- current Gaussian broadening implementation
- Direct code reading of `mace_gaussian/analysis/analysis_workflow.py` -- FWHM hardcoding at line 73, CLI entry points
- Direct code reading of `mace_gaussian/analysis/html_report_generator.py` -- methodology text at line 485
- Lorentzian line shape formula -- standard physics (Cauchy distribution), verified against the project's existing Gaussian formula structure

### Secondary (MEDIUM confidence)
- FWHM default of 10 cm-1 as "standard in computational spectroscopy" -- per user decision D-06, standard practice in published DFT/IR work

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, purely internal code changes
- Architecture: HIGH -- all target files read and understood, change points identified
- Pitfalls: HIGH -- edge cases identified from direct code reading (linregress guard, dual entry points)

**Research date:** 2026-03-28
**Valid until:** 2026-04-28 (stable -- no external dependency changes)
