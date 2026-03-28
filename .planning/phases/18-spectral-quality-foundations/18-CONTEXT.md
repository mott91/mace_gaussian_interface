# Phase 18: Spectral Quality Foundations - Context

**Gathered:** 2026-03-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace Gaussian KDE broadening with physically correct Lorentzian line shapes in IR spectrum plots, filter zero-intensity modes from intensity regression metrics, and update default FWHM. The mace_polar dipole investigation (PIPE-02) is deferred out of this phase.

</domain>

<decisions>
## Implementation Decisions

### Spectrum display
- **D-01:** Show Lorentzian broadened spectrum only — no stick spectrum for now (SPEC-03 stick display deferred)
- **D-02:** Keep existing DFT-above / ML-below offset layout, just swap Gaussian broadening for Lorentzian in `broaden_spectrum()`

### Zero-intensity filtering
- **D-03:** Filter mode pairs from intensity regression metrics when DFT reference intensity < 0.1 km/mol
- **D-04:** Frequency metrics (R², RMSE for frequencies) still include all modes regardless of intensity
- **D-05:** Lorentzian spectrum plots still show all modes — filtering is for statistics only, not visual display

### FWHM & broadening
- **D-06:** Default FWHM changed from 8 cm⁻¹ to 10 cm⁻¹ (standard in computational spectroscopy literature)
- **D-07:** FWHM configurable via `--fwhm` CLI argument on analysis scripts (run_analysis.py, run_analysis_harmonic.py). Already has `bandwidth_fwhm` parameter in SpectralAnalyzer — just wire to CLI.

### Claude's Discretion
- Lorentzian implementation details (normalization, numerical cutoff)
- How to label/annotate filtered modes in report text
- Whether to mention filtering count in report methodology section

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Analysis code
- `mace_gaussian/analysis/analyze_spectra.py` — Core analysis: `SpectralAnalyzer` class with `broaden_spectrum()` (line 324), FWHM parameter (line 88), spectrum plotting (line 575)
- `mace_gaussian/analysis/html_report_generator.py` — HTML report generation, mentions "8 cm⁻¹ FWHM" in methodology text (line 485)
- `mace_gaussian/analysis/analysis_workflow.py` — Wires intensity data as `DFT_Intensity_km_mol` / `ML_Intensity_km_mol` (line 291)

### Analysis entry points
- `run_analysis.py` — Anharmonic analysis entry point
- `run_analysis_harmonic.py` — Harmonic analysis entry point

### Requirements
- `.planning/REQUIREMENTS.md` §v1.2 — SPEC-01 (Lorentzian), SPEC-02 (zero-intensity filter), SPEC-03 (stick+broadened, deferred)

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `SpectralAnalyzer.broaden_spectrum()`: Already implements broadening loop over frequency grid with FWHM — swap Gaussian kernel for Lorentzian kernel
- `SpectralAnalyzer.__init__()`: Already accepts `bandwidth_fwhm` parameter and creates frequency grid
- Intensity data columns (`DFT_Intensity_km_mol`, `ML_Intensity_km_mol`) already available in analysis data pipeline

### Established Patterns
- Broadening: iterate over (freq, intensity) pairs, accumulate on `self.freq_grid` array
- Plot layout: DFT above with offset, ML below, normalized to global max
- HTML report: methodology text describes broadening approach

### Integration Points
- `broaden_spectrum()` is called in spectrum plotting methods — changing the kernel affects all plots automatically
- Intensity regression is computed in `match_by_mode()` and related methods — add filter there
- CLI arguments in `run_analysis.py` / `run_analysis_harmonic.py` need `--fwhm` flag

</code_context>

<specifics>
## Specific Ideas

- User explicitly does not want stick spectrum displayed — Lorentzian only, matching current single-panel layout
- PIPE-02 (mace_polar dipole) deferred by user request — should become a standalone todo, not part of this phase

</specifics>

<deferred>
## Deferred Ideas

- Stick spectrum display (SPEC-03 partial) — add later when needed, possibly in Phase 23 report overhaul
- PIPE-02: mace_polar dipole calculator investigation — user wants this as a separate todo, not in Phase 18
- Interactive Plotly spectrum viewer (from todo backlog) — future phase

</deferred>

---

*Phase: 18-spectral-quality-foundations*
*Context gathered: 2026-03-28*
