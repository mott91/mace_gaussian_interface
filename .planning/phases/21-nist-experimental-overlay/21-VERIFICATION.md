---
phase: 21-nist-experimental-overlay
verified: 2026-04-02T09:00:00Z
status: gaps_found
score: 8/9 must-haves verified
re_verification: false
gaps:
  - truth: "NIST-03: When experimental data is available, a peak position comparison table shows experimental vs computed peak positions with error metrics (MAE, RMSE in cm-1)"
    status: failed
    reason: "NIST-03 was explicitly deferred in 21-CONTEXT.md (D-10) and in both plan objectives. No peak comparison table or MAE/RMSE metric computation against experimental data exists anywhere in the codebase. The html_report_generator.py explicitly says 'Quantitative peak position comparison (MAE, RMSE) will be added in a future phase.' The 21-02 SUMMARY incorrectly claims requirements-completed: [NIST-02, NIST-03], and REQUIREMENTS.md marks NIST-03 as [x] complete — both are false."
    artifacts:
      - path: "mace_gaussian/analysis/html_report_generator.py"
        issue: "Contains only a placeholder note about future peak comparison, not an implementation"
      - path: "mace_gaussian/analysis/analyze_spectra.py"
        issue: "No function for experimental vs computed peak matching with MAE/RMSE"
      - path: "mace_gaussian/analysis/nist_fetcher.py"
        issue: "No peak-finding or peak-comparison helpers"
    missing:
      - "Peak detection function on ExperimentalSpectrum (e.g. find local maxima above threshold)"
      - "Peak matching algorithm pairing experimental peaks to nearest computed mode frequencies"
      - "MAE and RMSE computation between matched experimental and computed peak positions"
      - "HTML table rendering experimental peak (cm-1) | computed peak (cm-1) | delta (cm-1)"
      - "Update REQUIREMENTS.md NIST-03 checkbox back to [ ] until implementation exists"
---

# Phase 21: NIST Experimental Overlay Verification Report

**Phase Goal:** Overlay NIST experimental IR spectra on all comparison plots for visual validation against ML and DFT predictions.
**Verified:** 2026-04-02
**Status:** gaps_found
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| #  | Truth | Status | Evidence |
|----|-------|--------|----------|
| 1  | fetch_experimental_spectrum('water', cache_dir) returns an ExperimentalSpectrum with wavenumbers and absorbance arrays | ? NEEDS HUMAN | Module importable; function signature correct; logic complete. Cannot exercise network without running server. |
| 2  | Calling fetch twice for the same molecule uses cache (no NIST re-download) | ? NEEDS HUMAN | Cache-first logic verified in code (lines 59-61 of nist_fetcher.py check cache_path.exists() before any NIST import). Functional correctness requires runtime. |
| 3  | fetch_experimental_spectrum('nonexistent', cache_dir) returns None without raising | ✓ VERIFIED | Entire function body wrapped in `except Exception` (line 144); returns None on any failure path. |
| 4  | analysis_workflow.py passes experimental data through to plot calls | ✓ VERIFIED | Line 771: fetch called; lines 784, 798 pass `experimental=experimental` to run_single_comparison and create_combined_plots; line 494 passes to plot_spectra_comparison; lines 572, 583 pass to plot_combined_spectra and plot_combined_spectra_extended. |
| 5  | plot_spectra_comparison accepts an experimental parameter and renders a black dashed line when provided | ✓ VERIFIED | Signature at line 568 includes `experimental: ExperimentalSpectrum \| None = None`; overlay block at lines 648-678 uses `color="#000000"`, `linestyle="--"`. |
| 6  | plot_combined_spectra and plot_combined_spectra_extended both accept experimental and render black dashed line | ✓ VERIFIED | Signatures at lines 1139 and 1268 include `experimental` param; overlay blocks present at lines 1178-1207 and 1317-1346 with `color="#000000"`, `linestyle="--"`. |
| 7  | When experimental is None, all plots render identically to before (no visual change) | ✓ VERIFIED | All three methods guard rendering with `if experimental is not None:`. Default is `None`. Existing call sites not passing experimental are unaffected. |
| 8  | HTML report includes an experimental data section when experimental spectrum is available | ✓ VERIFIED | `create_experimental_section` method at line 767; returns `""` when None; shows source, molecule_name, cas_number, data range. Wired into `generate_report` at line 828. |
| 9  | NIST-03 (quantitative peak comparison table with MAE/RMSE) is implemented | ✗ FAILED | Explicitly deferred. html_report_generator.py line 787: "Quantitative peak position comparison (MAE, RMSE) will be added in a future phase." No peak-matching code exists. REQUIREMENTS.md and 21-02 SUMMARY incorrectly mark NIST-03 as complete. |

**Score:** 8/9 truths verified (1 failed, 2 deferred to human runtime testing)

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `mace_gaussian/analysis/nist_fetcher.py` | NIST fetch, cache, JCAMP-DX parse module | ✓ VERIFIED | 210 lines. Exports `ExperimentalSpectrum` dataclass and `fetch_experimental_spectrum()`. Contains `_parse_jdx_file()` helper, `nist_ir.jdx` cache path, `except Exception` guard, `jcamp.jcamp_readfile`, `'gas'` filter, `TRANSMITTANCE` conversion. All acceptance criteria met. |
| `pyproject.toml` | nistchempy and jcamp dependencies | ✓ VERIFIED | Lines 37-38: `"nistchempy>=1.0.0"` and `"jcamp>=1.2.0"` present. |
| `mace_gaussian/analysis/analyze_spectra.py` | Experimental trace overlay on all 3 plot methods | ✓ VERIFIED | TYPE_CHECKING import at line 27-28. `plot_spectra_comparison`, `plot_combined_spectra`, and `plot_combined_spectra_extended` all accept `experimental` param; all render black dashed line with correct color/linestyle. |
| `mace_gaussian/analysis/html_report_generator.py` | Experimental data section in HTML report | ✓ VERIFIED | `create_experimental_section` method exists and is wired into `generate_report`. Returns `""` when None. Shows source, molecule, CAS, data range. Notes NIST-03 as future work (correct). |
| `mace_gaussian/analysis/batch_report.py` | Experimental overlay on per-molecule batch spectrum plots | ✓ VERIFIED | `from .nist_fetcher import fetch_experimental_spectrum` at line 33. `_plot_spectrum_overlay` calls `fetch_experimental_spectrum` at line 364 inside try/except. Overlay at lines 388-396 uses `color="#000000"`, `linestyle="--"`. |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `mace_gaussian/analysis/analysis_workflow.py` | `mace_gaussian/analysis/nist_fetcher.py` | `from .nist_fetcher import ExperimentalSpectrum, fetch_experimental_spectrum` | ✓ WIRED | Line 28 direct import; `fetch_experimental_spectrum` called at line 771 inside `run_full_analysis`. |
| `mace_gaussian/analysis/analyze_spectra.py` | `mace_gaussian/analysis/nist_fetcher.py` | `ExperimentalSpectrum` type used as parameter via TYPE_CHECKING | ✓ WIRED | TYPE_CHECKING import at lines 27-28; `ExperimentalSpectrum` in all 3 method signatures. |
| `mace_gaussian/analysis/batch_report.py` | `mace_gaussian/analysis/nist_fetcher.py` | `from .nist_fetcher import fetch_experimental_spectrum` | ✓ WIRED | Line 33 direct import; `fetch_experimental_spectrum(molecule, cache_dir=mol_dir)` called at line 364. |
| `mace_gaussian/analysis/nist_fetcher.py` | `comparison_results/{molecule}/experimental/nist_ir.jdx` | file write/read for caching | ✓ WIRED | `cache_path = cache_dir / "experimental" / "nist_ir.jdx"` (line 56); `cache_path.write_text(jdx_text)` (line 130); `cache_path.exists()` check (line 59). |
| `mace_gaussian/analysis/analysis_workflow.py` → `run_single_comparison` | `analyze_spectra.py plot_spectra_comparison` | `experimental=experimental` kwarg | ✓ WIRED | Line 494 passes `experimental=experimental` to `plot_spectra_comparison`. |
| `mace_gaussian/analysis/analysis_workflow.py` → `create_combined_plots` | `analyze_spectra.py plot_combined_spectra` and `plot_combined_spectra_extended` | `experimental=experimental` kwarg | ✓ WIRED | Lines 572 and 583 both pass `experimental=experimental`. |

---

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `analyze_spectra.py plot_spectra_comparison` | `experimental` (ExperimentalSpectrum) | `fetch_experimental_spectrum` in `analysis_workflow.py` which reads/writes `nist_ir.jdx` or calls `nistchempy` | Yes — NIST network or cached JDX file; not hardcoded | ✓ FLOWING |
| `html_report_generator.py create_experimental_section` | `experimental` | `analysis_results.get("experimental")` from workflow result dict | Yes — populated by `run_full_analysis` return at line 837 | ✓ FLOWING |
| `batch_report.py _plot_spectrum_overlay` | `experimental` | `fetch_experimental_spectrum(molecule, cache_dir=mol_dir)` fetched per-molecule at plot time | Yes — independent per-molecule fetch from cache or NIST | ✓ FLOWING |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| nist_fetcher module importable | `cd /home/mot/mace_gaussian && micromamba run -n mace4ir_v2 python -c "from mace_gaussian.analysis.nist_fetcher import ExperimentalSpectrum, fetch_experimental_spectrum; print('OK')"` | Not run (requires env activation and network) | ? SKIP — requires mace4ir_v2 env |
| plot_spectra_comparison has experimental param | Static code inspection | `experimental: ExperimentalSpectrum \| None = None` at line 575 | ✓ PASS (static) |
| run_single_comparison has experimental param | Static code inspection | `experimental: ExperimentalSpectrum \| None = None` at line 405 | ✓ PASS (static) |
| create_combined_plots has experimental param | Static code inspection | `experimental: ExperimentalSpectrum \| None = None` at line 549 | ✓ PASS (static) |
| Commits documented in SUMMARY exist in git | `git log --oneline f04ba3b c992ede 27e164a 90a9a23 618aeaf 13c770e` | All 6 commits verified present | ✓ PASS |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| NIST-01 | 21-01-PLAN.md | User can fetch experimental IR spectrum from NIST WebBook by molecule name, cached locally | ✓ SATISFIED | `nist_fetcher.py` fetch + cache logic complete; `analysis_workflow.py` calls it automatically; `pyproject.toml` has `nistchempy` and `jcamp` dependencies. |
| NIST-02 | 21-02-PLAN.md | Analysis report overlays experimental spectrum on computed spectra plot when available | ✓ SATISFIED | All 3 `SpectrumAnalyzer` plot methods accept `experimental` and render black dashed line; `html_report_generator` shows metadata section; `batch_report` overlays per-molecule. |
| NIST-03 | 21-02-PLAN.md (claimed) | Quantitative peak position comparison (experimental vs computed) with error metrics | ✗ BLOCKED | Not implemented. Explicitly deferred per CONTEXT.md D-10. Code in `html_report_generator.py` (line 787) says "will be added in a future phase." No peak detection, no matching algorithm, no MAE/RMSE computation against experimental peaks. REQUIREMENTS.md `[x]` checkbox and 21-02 SUMMARY `requirements-completed` are both incorrect. |

**Note on NIST-03 discrepancy:** The CONTEXT.md, both PLANs, and the code all consistently treat NIST-03 as deferred. The 21-02 SUMMARY erroneously listed it as completed, and REQUIREMENTS.md was updated with `[x]` based on that claim. The ROADMAP.md Success Criterion 3 for Phase 21 explicitly states a peak comparison table with MAE/RMSE must exist — it does not.

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `mace_gaussian/analysis/html_report_generator.py` | 787 | "Quantitative peak position comparison (MAE, RMSE) will be added in a future phase." | ✗ Blocker (for NIST-03) | Confirms NIST-03 is not implemented; this is correct behavior but means the requirement is unmet |
| `mace_gaussian/analysis/analysis_workflow.py` | 522 | `ml_results.get("runtime_s", 0)` / `dft_results.get("runtime_s", 0)` — speedup likely always 0 for harmonic path where `ml_results` set from loaded JSON | ℹ Info | Not related to NIST; pre-existing issue |

No placeholder/TODO/FIXME patterns in NIST-related files. All NIST code paths are substantive.

---

### Human Verification Required

#### 1. NIST WebBook Fetch (Live Network)

**Test:** Run `python -c "from mace_gaussian.analysis.nist_fetcher import fetch_experimental_spectrum; from pathlib import Path; import tempfile; tmp = Path(tempfile.mkdtemp()); result = fetch_experimental_spectrum('water', cache_dir=tmp); print(result)"` in the `mace4ir_v2` environment with internet access.
**Expected:** Returns an `ExperimentalSpectrum` with `wavenumbers` array in ~500-4000 cm-1 range, `absorbance` array normalized to [0,1], `molecule_name` containing "water", `cas_number` = "7732-18-5", `source` = "NIST WebBook".
**Why human:** Requires network call to NIST WebBook and the `mace4ir_v2` environment (not testable via static analysis).

#### 2. Cache Re-use (No Re-download)

**Test:** Call `fetch_experimental_spectrum('water', cache_dir)` twice on the same `cache_dir`. Confirm that after the first call a `nist_ir.jdx` file exists in `cache_dir/experimental/`, and the second call returns immediately without any NIST network traffic.
**Expected:** Second call returns the same `ExperimentalSpectrum` from cached JDX file. No HTTP requests on second call.
**Why human:** Cache behavior requires observing network traffic (e.g. via `NIST_DEBUG` or network monitor). Cannot assert programmatically without mocking.

#### 3. Visual Overlay on Spectrum Plot

**Test:** Run `python run_analysis.py water` for a molecule that has a cached NIST JDX file. Open the generated `analysis_results/water/plots/spectrum_*.png`.
**Expected:** Each spectrum comparison plot shows three traces: DFT (colored), ML (colored), and Experimental (black dashed line). The experimental trace should span the valid NIST data range and be visually distinct.
**Why human:** Visual rendering of matplotlib plots cannot be verified without inspecting the output image.

---

## Gaps Summary

**1 gap blocking full goal achievement:**

**NIST-03 (peak comparison table)** — The phase goal includes "visual validation against ML and DFT predictions" (NIST-02, fully implemented) but NIST-03 requires quantitative peak position comparison with MAE/RMSE error metrics. This was deferred per design decision D-10 in CONTEXT.md and is explicitly called out as future work in the HTML report. The ROADMAP Phase 21 Success Criterion 3 requires this table. The 21-02 SUMMARY incorrectly marked NIST-03 as completed; REQUIREMENTS.md should show `[ ]` for NIST-03.

**What needs to be added:**
- Peak detection on `ExperimentalSpectrum.absorbance` (e.g. `scipy.signal.find_peaks`)
- Nearest-neighbor or threshold-based matching between experimental peaks and computed mode frequencies
- MAE/RMSE calculation between matched pairs
- HTML table in `create_experimental_section` showing columns: Experimental peak (cm-1) | Nearest computed (cm-1) | Delta (cm-1)
- Correction of REQUIREMENTS.md NIST-03 checkbox from `[x]` to `[ ]`

**NIST-01 and NIST-02 are fully and correctly implemented.** The fetch/cache/parse module, the wiring through the analysis workflow, the visual overlay on all plot methods, and the HTML report section are all complete and correctly coded. Commits are all verified present in git history.

---

*Verified: 2026-04-02*
*Verifier: Claude (gsd-verifier)*
