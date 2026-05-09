# Phase 23 — Next Session TODO

## What's Done This Session
- Bug 1: Regression plot axes equal range (scaleanchor + constrain="domain")
- Bug 2: Overtones section pulls from SpectrumData labels, DFT vs ML side-by-side with Δfreq
- Bug 3: Heatmap was working; now base64-embedded for SSH download
- Feature 4: `mode="harmonic"` wired in analysis_workflow.py
- Regression plots color-coded by type (fundamental/overtone/combination) — both freq and intensity
- Per-Category Accuracy Ranking section (was "Category Awards")
- Experimental line: alpha 0.15, solid
- Combined spectrum: 10 maximally distinct colors
- Regression y-axis: simplified to "ML frequencies"/"ML intensities"
- Regression plots enlarged (550x550), legend inside plot (top-left)
- Comparisons sorted by MAE (best first in nav + sections)
- Metrics stat boxes above spectrum, timing integrated into boxes
- Wider layout (95% max-width)
- Heatmap: numbers hidden when >12 modes, x-labels rotated for large matrices
- Overtone tables: summary stats + top-15 worst errors for large molecules
- Degenerate notes collapsed into single compact line
- Timing fallback to runtime_s when gaussian_timing unavailable
- Mode count overview (fundamentals/overtones/combinations/total) in executive summary
- DegenerateGroup .get() → getattr() fix

## Remaining / Open Items

### 1. Octane fundamentals all IR-inactive
All fundamental and overtone DFT intensities are 0.0 for octane (symmetric molecule).
Only combination bands have nonzero intensity. This is correct physics, not a bug.
The intensity regression plot correctly shows only combos after the >= 0.1 filter.

### 2. Feature 5: Eigenvector overlap confidence flags per mode
In regression plots and/or per-method metrics table, show eigenvector overlap score
for each matched mode pair. Flag modes with overlap < 0.7 as low-confidence.
Data available from `comp["mode_mapping"]` and overlap matrix.

### 3. Investigate 5600-5800 cm⁻¹ intensity difference (mace_anicc_mace_ml on octane)
User noticed a large intensity difference in this overtone/combo frequency region.
Check which specific mode_id via hover tooltip. May be a real ML prediction error.

### 4. Timing data missing for pre-phase-20 calculations
Octane calcs were run before timing was implemented. Need to either:
- Recalculate octane with timing enabled
- Or accept N/A for old calcs (currently shows 0.0s as fallback)

### 5. Harmonic report still generates with new Plotly style
User noted harmonic report was replaced with new style (not what they wanted).
Decide: keep single anharmonic report or maintain separate harmonic report.
For thesis, one anharmonic report per molecule is likely sufficient.

## Files Changed
- `mace_gaussian/analysis/plotly_builders.py` — colors, regression layout, type color-coding, legend
- `mace_gaussian/analysis/_shared_css.py` — wider layout, awards CSS, nav flex-wrap
- `mace_gaussian/analysis/html_report_generator.py` — MAE sorting, metrics boxes, timing, overtones, mode overview, degenerate fix
- `mace_gaussian/analysis/analyze_spectra.py` — matched_mode_ids in match_stats
- `mace_gaussian/analysis/mode_matching.py` — heatmap text/label scaling
- `mace_gaussian/analysis/analysis_workflow.py` — mode="harmonic" wiring
