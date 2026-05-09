# Phase 23 — Next Session TODO

## What's Done This Session

### Verification gap fixes (commit `d9879cb`)
- `batch_report.py` now imports `_shared_css.build_css` and `html` stdlib;
  `_esc()` helper wraps user-controlled molecule/combo/timing strings (T-23-01)
- `analysis_workflow.py` passes `mode="harmonic" | "anharmonic"` to
  `HTMLReportGenerator` based on `use_harmonic`
- 23-VERIFICATION should now flip 5/7 → 7/7 on re-run

### Report polish (commit `98af64c`)
- Regression plots: equal-range axes, type color-coding, R²/count legend
- Overtones section, heatmap, executive summary, layout improvements
- See commit body for the full list

### Feature 5 — eigenvector overlap confidence flags (commit `e2f2ffe`)
- `analyze_spectra.match_by_mode` accepts `mode_overlaps`, returns
  `matched_mode_overlaps` aligned with matched arrays
- `plotly_builders` exposes `LOW_OVERLAP_THRESHOLD = 0.7` and a shared
  `_add_regression_traces` helper. Both regression builders accept
  `mode_overlaps`; low-overlap fundamentals render as hollow markers
  (`circle-open`) at 0.55 opacity in their own legend trace, with the
  overlap value shown on hover
- `html_report_generator` threads overlaps through and adds a
  "Low-overlap (< 0.7)" stat box (n/N fundamentals)
- `analysis_workflow` stores `mode_overlaps` in `comp` so the report has
  the data without re-running the match

### Adjacent work — `mace_polar1` dipole calculator (commit `0a014a5`)
- New `MACEPolar1DipoleCalculator` registered alongside `espaloma` and
  `mace_ml`, wired through factory / CLI / workflow defaults
- Spike scripts in `scripts/spike_polar_*` validate static dipoles vs
  B3LYP and autograd-vs-finite-difference derivative agreement
- End-to-end IR intensity comparison on water + HF: POLAR-1 and MACE4IR
  trade wins per molecule — kept as additional dipole option, not
  default replacement

## Resolved Decisions

- **Harmonic report** — keep unified Plotly style with `mode=` flag.
  Anharmonic stays the default; harmonic report is fast sanity check.
  No separate template.
- **Octane all-zero fundamentals** — correct physics for the symmetric
  alkane, not a bug. Intensity regression filter (≥ 0.1 km/mol) keeps
  combos visible.

## Remaining / Open Items (deferred)

### Octane intensity discrepancy at 5600–5800 cm⁻¹
Investigation deferred until octane is recalculated with timing
instrumentation during the v1.3 benchmark campaign. Recipe to follow
when picked up:
1. Hover the regression-plot points in 5600–5800 cm⁻¹ to read off
   `mode_id`; look up overtone vs combination in
   `comparison_results/octane/mace_anicc_mace_ml/results.json`
2. Compare the parent fundamentals (DFT vs ML) — discrepancy is either
   in the underlying ML fundamental (dipole-derivative error amplified
   into the overtone) or in the χ matrix (anharmonic coupling, energy
   side)
3. Check Gaussian freq log for Fermi resonance comments around `2νₖ`

### Pre-Phase-20 octane timing data missing
Resolved-by-deferring: octane will be re-run with timing enabled in the
benchmark campaign. Until then the report falls back to `runtime_s`.

## Files Changed (this session)
- `mace_gaussian/analysis/plotly_builders.py` — overlap helper, hollow
  markers, hover tooltips
- `mace_gaussian/analysis/analyze_spectra.py` — overlap plumbing in
  `match_by_mode`
- `mace_gaussian/analysis/html_report_generator.py` — overlap stat box,
  threading
- `mace_gaussian/analysis/analysis_workflow.py` — store mode_overlaps in
  comp; mode flag wiring
- `mace_gaussian/analysis/batch_report.py` — XSS mitigation, shared CSS
- `mace_gaussian/calculators/mace_polar1.py` — new dipole calculator
- `mace_gaussian/calculators/{__init__,factory}.py`,
  `mace_gaussian/{cli,workflow}.py` — wiring
- `tests/test_plotly_builders.py`, `tests/test_calculators.py`,
  `tests/test_batch_report.py` — coverage for above
