---
phase: 18-spectral-quality-foundations
verified: 2026-03-30T09:33:00Z
status: passed
score: 5/5 must-haves verified
gaps:
  - truth: "SPEC-02 and 18-02-PLAN completion reflected in tracking documents"
    status: resolved
    reason: "Tracking documents updated inline — REQUIREMENTS.md SPEC-02 marked [x], ROADMAP.md 18-02-PLAN.md marked [x]"
    artifacts:
      - path: ".planning/REQUIREMENTS.md"
        issue: "Line 73: SPEC-02 still marked [ ] (pending); traceability table line 175 also says Pending"
      - path: ".planning/ROADMAP.md"
        issue: "Line 76: 18-02-PLAN.md still marked [ ] (incomplete); Phase 18 top-level phase also still [ ]"
    missing:
      - "Change `- [ ] **SPEC-02**` to `- [x] **SPEC-02**` in REQUIREMENTS.md (line 73)"
      - "Change `| SPEC-02 | Phase 18 | Pending |` to `| SPEC-02 | Phase 18 | Complete |` in REQUIREMENTS.md (line 175)"
      - "Change `- [ ] 18-02-PLAN.md` to `- [x] 18-02-PLAN.md` in ROADMAP.md (line 76)"
      - "Change `- [ ] **Phase 18: Spectral Quality Foundations**` to `- [x]` in ROADMAP.md (top-level phase list)"
human_verification:
  - test: "Open any previously-run molecule's HTML report and inspect the methodology section"
    expected: "Text reads 'Broadening was applied using Lorentzian line shapes with 10 cm-1 FWHM' with intensity filtering note — NOT 'Gaussian convolution with 8 cm-1 FWHM'"
    why_human: "Requires browser rendering to confirm HTML display is visually correct"
---

# Phase 18: Spectral Quality Foundations Verification Report

**Phase Goal:** Simulated IR spectra are physically correct (Lorentzian line shapes) and intensity statistics are meaningful (zero-intensity modes excluded from regression)
**Verified:** 2026-03-30T09:33:00Z
**Status:** gaps_found (documentation tracking only — all code is complete)
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `broaden_spectrum()` produces a Lorentzian line shape (peak at center equals intensity, half-max at +/- gamma) | VERIFIED | `test_lorentzian_peak_height` and `test_lorentzian_fwhm` pass; `gamma_sq / ((freq_grid - freq)**2 + gamma_sq)` confirmed in analyze_spectra.py line 338 |
| 2 | Default FWHM is 10 cm-1 throughout the pipeline (SpectrumAnalyzer, ComparisonWorkflow, CLI) | VERIFIED | `test_default_fwhm` passes; `bandwidth_fwhm: float = 10.0` in SpectrumAnalyzer.__init__, ComparisonWorkflow.__init__, analyze_molecule(), analyze_molecule_harmonic(); no `bandwidth_fwhm=8.0` anywhere in mace_gaussian/ |
| 3 | User can pass --fwhm to both run_analysis.py and run_analysis_harmonic.py to override FWHM | VERIFIED | `python run_analysis.py --help` and `python run_analysis_harmonic.py --help` both show `--fwhm FWHM` with Lorentzian description; argparse wired in run_analysis_main() and run_analysis_harmonic_main() |
| 4 | Intensity regression R-squared and RMSE exclude modes where DFT intensity < 0.1 km/mol | VERIFIED | `int_mask = dft_int >= INTENSITY_THRESHOLD` at analyze_spectra.py line 521; `num_intensity_filtered` field on ComparisonMetrics; all 3 filtering tests pass |
| 5 | Frequency regression metrics still include all modes regardless of intensity | VERIFIED | Frequency metrics block (lines 509-517) has no `int_mask` or filter; `test_frequency_metrics_unfiltered` asserts `num_matched == 5` and `mae_freq == 5.0` for 5-mode dataset |
| 6 | HTML report methodology text says Lorentzian with the actual FWHM value used | VERIFIED | `html_report_generator.py` line 488-489: `"Lorentzian line shapes with {self.bandwidth_fwhm} cm<sup>-1</sup> FWHM"`; no "Gaussian convolution" text remains; `bandwidth_fwhm` threaded from ComparisonWorkflow to HTMLReportGenerator |
| 7 | SPEC-03 (stick spectrum) is deferred per user decision D-01 and documented | VERIFIED | D-01 defined in 18-CONTEXT.md, 18-RESEARCH.md, ROADMAP.md note (line 78); 18-02-SUMMARY.md documents deferral to Phase 23 |
| 8 | PIPE-02 (mace_polar dipole) is deferred per user request and documented | VERIFIED | Standalone todo at `.planning/todos/pending/2026-03-26-reevaluate-mace-polar-as-dipole-calculator.md` exists; deferral documented in 18-02-SUMMARY.md |
| 9 | SPEC-02 and Plan 02 completion reflected in tracking documents | FAILED | REQUIREMENTS.md line 73 still `[ ] SPEC-02`; REQUIREMENTS.md line 175 still `Pending`; ROADMAP.md line 76 `18-02-PLAN.md` still `[ ]` |

**Score:** 8/9 truths verified (code goal fully achieved; tracking document updates pending)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `tests/test_spectral_broadening.py` | Unit tests for Lorentzian broadening and intensity filtering | VERIFIED | 104 lines (min 100); 6 tests: test_lorentzian_peak_height, test_lorentzian_fwhm, test_default_fwhm, TestIntensityFiltering (3 methods) |
| `mace_gaussian/analysis/analyze_spectra.py` | Lorentzian kernel in broaden_spectrum(), int_mask in calculate_metrics(), num_intensity_filtered field | VERIFIED | `gamma = self.bandwidth_fwhm / 2.0` at line 334; `int_mask = dft_int >= INTENSITY_THRESHOLD` at line 521; `num_intensity_filtered: int` in ComparisonMetrics dataclass line 69 |
| `mace_gaussian/analysis/analysis_workflow.py` | FWHM CLI wiring and 10.0 default, argparse | VERIFIED | `import argparse` (local import in both main() fns); `bandwidth_fwhm: float = 10.0` in ComparisonWorkflow.__init__, analyze_molecule(), analyze_molecule_harmonic(); `parser.add_argument("--fwhm"` in both entry points |
| `mace_gaussian/analysis/html_report_generator.py` | Dynamic FWHM in methodology text, "Lorentzian" present | VERIFIED | `bandwidth_fwhm: float = 10.0` parameter in __init__; `self.bandwidth_fwhm` stored; f-string with `{self.bandwidth_fwhm}` in methodology text; "0.1 km/mol" filtering note present |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `mace_gaussian/analysis/analysis_workflow.py` | `mace_gaussian/analysis/analyze_spectra.py` | `bandwidth_fwhm` parameter threading | WIRED | ComparisonWorkflow passes `bandwidth_fwhm=bandwidth_fwhm` to SpectrumAnalyzer at line 76; runtime inspection confirmed |
| `run_analysis.py` | `mace_gaussian/analysis/analysis_workflow.py` | calls `run_analysis_main()` | WIRED | `--help` output confirms argparse active; `analyze_molecule(args.molecule_name, bandwidth_fwhm=args.fwhm)` at line 898 |
| `mace_gaussian/analysis/analyze_spectra.py` | `ComparisonMetrics` dataclass | `num_intensity_filtered` field | WIRED | Field declared at line 69; populated in calculate_metrics() return at line 558 and early-return at line 506 |
| `mace_gaussian/analysis/analysis_workflow.py` | `mace_gaussian/analysis/html_report_generator.py` | `bandwidth_fwhm=self.bandwidth_fwhm` in generate_html_report | WIRED | `HTMLReportGenerator(..., bandwidth_fwhm=self.bandwidth_fwhm)` at workflow.py line 789; runtime inspection confirmed |

### Data-Flow Trace (Level 4)

Not applicable — phase modifies analysis computation and CLI entry points, not data ingestion or rendering pipelines with independent data sources. The broadening function is a pure transform (no database/API/store); the intensity filter is applied inline during metric computation. Both are covered by unit tests that verify correct numeric output.

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| All 6 spectral tests pass | `micromamba run -n mace4ir_v2 pytest tests/test_spectral_broadening.py -v` | 6 passed, 0 failed | PASS |
| --fwhm in anharmonic help | `python run_analysis.py --help \| grep fwhm` | `--fwhm FWHM  Full width at half maximum for Lorentzian broadening in cm-1` | PASS |
| --fwhm in harmonic help | `python run_analysis_harmonic.py --help \| grep fwhm` | `--fwhm FWHM  Full width at half maximum for Lorentzian broadening in cm-1` | PASS |
| No old Gaussian convolution text | `grep "Gaussian convolution" mace_gaussian/analysis/html_report_generator.py` | No matches | PASS |
| No old broad_param | `grep -r "broad_param" mace_gaussian/` | No matches | PASS |
| No hardcoded bandwidth_fwhm=8.0 | `grep -r "bandwidth_fwhm=8.0" mace_gaussian/` | No matches | PASS |
| bandwidth_fwhm in ComparisonWorkflow and HTMLReportGenerator | Python inspect of both __init__ signatures | Both confirmed | PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| SPEC-01 | 18-01-PLAN.md | Lorentzian line shapes, configurable FWHM default 10 cm-1 | SATISFIED | Lorentzian kernel verified in code; 3 unit tests pass; REQUIREMENTS.md already marks as `[x]` |
| SPEC-02 | 18-02-PLAN.md | Modes < 0.1 km/mol filtered from intensity regression, retained in frequency metrics | SATISFIED (code); NOT UPDATED (docs) | `int_mask = dft_int >= 0.1` confirmed; 3 unit tests pass; but REQUIREMENTS.md line 73 still `[ ]` |
| SPEC-03 | 18-02-PLAN.md (deferred) | Stick spectrum overlay | INTENTIONALLY DEFERRED | D-01 decision documented in 18-CONTEXT.md, 18-RESEARCH.md, ROADMAP.md note; deferral to Phase 23 confirmed |
| PIPE-02 | 18-02-PLAN.md (deferred) | mace_polar dipole investigation | INTENTIONALLY DEFERRED | Standalone todo at `.planning/todos/pending/2026-03-26-reevaluate-mace-polar-as-dipole-calculator.md` confirmed |

**Orphaned requirements check:** REQUIREMENTS.md maps SPEC-02, SPEC-03, PIPE-02 to Phase 18 (traceability table lines 175-177). All three are accounted for in the plans — SPEC-02 implemented, SPEC-03 and PIPE-02 deferred with documented decisions. No orphaned requirements.

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `.planning/REQUIREMENTS.md` line 73 | `[ ] SPEC-02` — not updated to `[x]` after implementation | Info | Tracking inaccuracy; no code impact |
| `.planning/REQUIREMENTS.md` line 175 | `SPEC-02 \| Phase 18 \| Pending` — not updated to `Complete` | Info | Tracking inaccuracy; no code impact |
| `.planning/ROADMAP.md` line 76 | `[ ] 18-02-PLAN.md` — not updated to `[x]` after execution | Info | Tracking inaccuracy; no code impact |

No code-level anti-patterns found. No TODOs, FIXMEs, stubs, or hardcoded empty values in the modified files.

### Human Verification Required

#### 1. HTML Report Methodology Section

**Test:** Run analysis on a molecule that has existing results (e.g., `python run_analysis.py water`), then open the generated `analysis_results/water/report.html` in a browser.
**Expected:** Methodology section reads "Broadening was applied using Lorentzian line shapes with 10 cm-1 FWHM" and includes the sentence "Modes with DFT IR intensity below 0.1 km/mol were excluded from intensity regression metrics." The old text "Gaussian convolution with 8 cm-1 FWHM" must not appear.
**Why human:** HTML rendering and browser display cannot be verified programmatically in this environment.

### Gaps Summary

The phase goal is **fully achieved in code**: Lorentzian broadening is live and mathematically correct (verified by 3 unit tests), FWHM is threaded through the full CLI-to-SpectrumAnalyzer pipeline, intensity regression correctly excludes near-zero DFT modes (verified by 3 unit tests), and HTML reports dynamically display the Lorentzian methodology text. All 6 tests pass. No stubs, no placeholder implementations, no old Gaussian code remains.

The single gap is **documentation-only**: REQUIREMENTS.md and ROADMAP.md tracking checkboxes were not updated after Plan 02 executed. Specifically:
- REQUIREMENTS.md line 73: SPEC-02 still `[ ]` (should be `[x]`)
- REQUIREMENTS.md line 175: SPEC-02 traceability row still `Pending` (should be `Complete`)
- ROADMAP.md line 76: 18-02-PLAN.md still `[ ]` (should be `[x]`)
- ROADMAP.md Phase 18 top-level entry still `[ ]` (should be `[x]` now that both plans are done)

These are three-line edits. They do not block any downstream phase (Phase 19 depends on code behavior, not checkbox state).

---

_Verified: 2026-03-30T09:33:00Z_
_Verifier: Claude (gsd-verifier)_
