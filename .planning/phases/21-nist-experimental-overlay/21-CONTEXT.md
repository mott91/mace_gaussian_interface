# Phase 21: NIST Experimental Overlay - Context

**Gathered:** 2026-04-01
**Status:** Ready for planning

<domain>
## Phase Boundary

Fetch experimental IR spectra from NIST WebBook, cache them locally per molecule, and overlay them on computed DFT/ML spectrum plots in the analysis reports. This phase covers NIST-01 (fetch + cache) and NIST-02 (visual overlay) only. NIST-03 (quantitative peak position comparison with error metrics) is deferred to a future phase.

</domain>

<decisions>
## Implementation Decisions

### Data sourcing
- **D-01:** Use `nistchempy` library to fetch experimental spectra from NIST WebBook. Already identified as a new dependency in STATE.md.
- **D-02:** Parse JCAMP-DX format using `jcamp` library (second new dependency from STATE.md).
- **D-03:** Gas-phase IR data only. If a molecule has no gas-phase IR spectrum in NIST, skip silently — no condensed-phase fallback.
- **D-04:** Best-effort fetching — NIST fetch failure never blocks or errors the analysis pipeline.

### Caching
- **D-05:** Cache experimental spectra in `comparison_results/{molecule}/experimental/` alongside existing DFT and ML result directories.
- **D-06:** If cached data exists, skip re-download on subsequent runs.

### Spectral overlay visuals
- **D-07:** Experimental spectrum appears as a third trace on the existing DFT-above/ML-below offset plot layout. Not a separate panel or figure.
- **D-08:** Both experimental and computed spectra normalized to [0, 1] range for visual comparison. Axis labeled "Normalized intensity".
- **D-09:** Experimental trace rendered as black dashed line — distinct from seaborn colorblind palette colors (blue DFT, orange ML). Convention in spectroscopy papers.

### Peak matching
- **D-10:** No quantitative peak matching in this phase. NIST-03 deferred. Comparison is visual overlay only.

### Pipeline integration
- **D-11:** Automatic — analysis scripts auto-fetch from NIST if spectrum not cached. No CLI flag required. If fetch fails or molecule not found, analysis runs normally without experimental overlay section.
- **D-12:** Both `run_analysis.py` (anharmonic) and `run_analysis_harmonic.py` (harmonic) support experimental overlay when data is available.
- **D-13:** Batch report (`batch_report.py`) includes per-molecule experimental overlay when available.

### Claude's Discretion
- JCAMP-DX parsing details (handling malformed files, unit detection)
- Legend placement and labeling for the experimental trace
- How to note "no experimental data available" in the report (if at all)
- Frequency range alignment between experimental and computed grids
- Cache file format (raw JCAMP-DX, parsed numpy arrays, or both)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Analysis code
- `mace_gaussian/analysis/analyze_spectra.py` — `SpectralAnalyzer` with `broaden_spectrum()`, spectrum plotting, FWHM parameter. Experimental trace will be added to existing plot methods.
- `mace_gaussian/analysis/html_report_generator.py` — HTML report generation. Must add experimental overlay section/toggle.
- `mace_gaussian/analysis/analysis_workflow.py` — Orchestrates analysis. NIST fetch/cache logic hooks in here.
- `mace_gaussian/analysis/batch_report.py` — Batch multi-molecule report. Must support per-molecule experimental overlay.

### Analysis entry points
- `run_analysis.py` — Anharmonic analysis entry point
- `run_analysis_harmonic.py` — Harmonic analysis entry point

### Result structure
- `mace_gaussian/results_manager.py` — Manages `comparison_results/{molecule}/` directory structure. New `experimental/` subdirectory goes here.

### Requirements
- `.planning/REQUIREMENTS.md` §Experimental Comparison — NIST-01 (fetch+cache), NIST-02 (overlay). NIST-03 deferred.

### Prior phase context
- `.planning/phases/18-spectral-quality-foundations/18-CONTEXT.md` — Lorentzian broadening decisions (D-01, D-02, D-06, D-07). Experimental overlay uses the same broadened spectrum infrastructure.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `SpectralAnalyzer.broaden_spectrum()`: Lorentzian broadening already implemented (Phase 18). Experimental overlay goes on the same frequency grid.
- `SpectralAnalyzer` plotting methods: Existing DFT/ML comparison plots — add third trace here.
- `results_manager.py`: Directory management for `comparison_results/{molecule}/` — extend with `experimental/` subdirectory pattern.

### Established Patterns
- Analysis workflow pattern: `analysis_workflow.py` orchestrates data loading, analysis, and report generation. NIST fetch is a new data loading step.
- Best-effort pattern: STATE.md already establishes "never blocks analysis" as the error handling philosophy for NIST.
- Seaborn colorblind palette: All plots use this. Experimental trace (black dashed) is deliberately outside the palette.

### Integration Points
- `analysis_workflow.py` — Where NIST data loading hooks in (before plotting)
- `analyze_spectra.py` — Where the third trace gets added to spectrum plots
- `html_report_generator.py` — Where experimental overlay section appears in the report
- `batch_report.py` — Where per-molecule experimental data appears in batch reports
- `comparison_results/{molecule}/experimental/` — New directory for cached JCAMP-DX data

</code_context>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches for JCAMP-DX parsing and NIST WebBook interaction.

</specifics>

<deferred>
## Deferred Ideas

- **NIST-03 quantitative peak comparison**: Peak position matching with MAE/RMSE metrics. User wants visual comparison first; quantitative matching is a future phase.
- **SDBS as secondary source**: The todo mentioned SDBS (Spectral Database for Organic Compounds) as a backup. Not in scope for this phase.
- **Condensed-phase fallback**: Could optionally show condensed-phase data with a shift warning. Deferred — gas-phase only for now.

### Reviewed Todos (not folded)
- "Automatic functional group peak labeling" — separate analysis enhancement, not related to NIST overlay
- "Interactive HTML spectrum viewer with Plotly" — visualization upgrade, different phase
- "JCAMP-DX spectral data export" — export of computed spectra, opposite direction from this phase

</deferred>

---

*Phase: 21-nist-experimental-overlay*
*Context gathered: 2026-04-01*
