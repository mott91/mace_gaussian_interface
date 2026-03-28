# Project Research Summary

**Project:** MACE-Gaussian v1.2 Analysis Quality Overhaul
**Domain:** Computational chemistry IR spectroscopy — ML potential benchmarking pipeline
**Researched:** 2026-03-28
**Confidence:** HIGH (stack and architecture based on direct codebase analysis; features and pitfalls verified against spectroscopy literature)

## Executive Summary

This milestone adds analysis quality features to an existing, working MACE-Gaussian pipeline. The existing foundation (Python/ASE/ZMQ/Gaussian) is sound and unchanged. What v1.2 adds is a set of scientifically expected outputs: physically correct Lorentzian line shapes, degenerate-mode-aware matching, zero-intensity filtering for intensity statistics, wall-clock timing for the cost-benefit thesis argument, and experimental spectra overlay via NIST. These features convert the current "technically correct but presentation-rough" pipeline into thesis-quality output. None of the core four features require new dependencies beyond `nistchempy` and `jcamp` for the NIST experimental overlay.

The recommended build order is driven by two constraints: (1) Lorentzian broadening and zero-intensity filtering are upstream of everything visual, so they should land first, and (2) the PyVPT2 alternative VPT2 engine is the one genuinely research-grade item — it requires a QCEngine adapter that does not yet exist, carries Fermi resonance risk, and should be deferred or prototyped in isolation without blocking the other features. All four core features (broadening, filtering, degenerate mode handling, timing) have well-defined scopes and can be implemented in 2-3 focused phases.

The most critical risks are: (a) lineshape normalization confusion when Lorentzian and Gaussian broadening coexist — always store broadening type in metadata and require matching parameters for spectral comparisons; (b) degenerate modes producing misleadingly low overlap scores unless subspace overlap (trace(M^T M)/k) replaces individual dot products; and (c) NIST retrieval fragility — treat as a best-effort cache, never block analysis on it. Fix the acetic acid parser bug and the xTB dipole unit bug before building new features on top of the parser or before any VPT2 work.

## Key Findings

### Recommended Stack

Only two new packages are needed for the entire v1.2 feature set: `nistchempy==1.0.5` (unofficial NIST WebBook scraper, MIT license, updated March 2026) and `jcamp==1.3.0` (JCAMP-DX parser, MIT). Everything else — Lorentzian broadening, timing, degenerate mode handling, zero-intensity filtering — is implemented with existing numpy/scipy/stdlib. PyVPT2 is blocked by a Psi4 hard dependency and cannot consume pre-computed force constants; the McCoy Group's Psience library may be able to ingest Gaussian .fchk force constants but needs a 2-hour investigation spike to confirm.

**Core technologies:**
- `numpy` (existing): Lorentzian broadening — trivial broadcast formula, no library needed
- `nistchempy==1.0.5`: NIST WebBook scraping — only Python package for this; unofficial but active
- `jcamp==1.3.0`: JCAMP-DX parsing — mature, MIT, handles NIST's format directly
- `time.perf_counter` / `contextlib` (stdlib): Pipeline timing — already partially used in `workflow.py`; extend, don't replace
- `scipy.spatial.distance` (existing): Degenerate mode detection — threshold-based frequency clustering
- PyVPT2 / Psience (optional spike only): VPT2 alternatives — do not install until spike confirms integration path

### Expected Features

**Must have (table stakes):**
- Lorentzian broadening — every published computed IR spectrum uses it; current Gaussian broadening is physically wrong for IR line shapes
- Zero-intensity mode filtering — methane R²=0.34 intensity issue is caused by correlating IR-inactive zeros with numerical noise; any reviewer will flag this
- Degenerate mode handling — without subspace overlap scoring, mode matching is systematically pessimistic for methane, benzene, ammonia (many thesis molecules)
- Wall-clock timing (ML vs DFT) — "X% accuracy at Y% of the computational cost" is the headline ML benchmarking result; absent timing data, the thesis has no cost-benefit argument

**Should have (differentiators):**
- NIST experimental spectra overlay — comparing ML vs DFT only tests self-consistency; comparing both against experiment is a substantially stronger thesis result
- Early SLURM DFT submission — reordering the batch pipeline to submit DFT jobs immediately after geometry optimization saves hours of queue wait during the benchmark campaign
- Anharmonic pipeline overhaul — integrates all of the above into a polished HTML report with experimental overlay panel, timing breakdown, and degenerate mode confidence indicators

**Defer (v1.3+):**
- PyVPT2 integration — research-grade, requires a QCEngine/ASE harness that does not exist; novel (no prior MACE+PyVPT2 work); high risk for a thesis deadline
- Voigt profile broadening — physically more accurate for gas-phase but overkill for computed spectra comparison; pure Lorentzian is the field standard
- Interactive Plotly spectrum viewer — significant effort with no improvement to thesis results
- Delta-ML correction model — insufficient data; premature optimization
- Automated functional group peak labels — error-prone rabbit hole

### Architecture Approach

All new features integrate into the existing three-layer architecture (pipeline layer: `workflow.py`/`batch.py`/`slurm.py`; analysis layer: `analysis_workflow.py`/`analyze_spectra.py`/`mode_matching.py`/`html_report_generator.py`; data layer: `results.json`/`.fchk`/`batch_manifest.json`) without restructuring it. The key principle is additive, optional enhancement: every new parameter defaults to preserving current behavior. New modules (`utils/timing.py`, `analysis/nist_spectra.py`) isolate cross-cutting concerns; modifications to existing modules (`analyze_spectra.py`, `mode_matching.py`, `batch.py`) are backward-compatible via new function variants and optional parameters.

**Major components (new or modified):**
1. `utils/timing.py` (NEW) — `PipelineTimer` context manager; wrap existing `time.perf_counter()` calls in `workflow.py` and `batch.py`
2. `analysis/analyze_spectra.py` (MODIFY) — add `broaden_spectrum_lorentzian()` parallel to existing Gaussian broadening; add `min_dft_intensity` filter to `calculate_metrics()`
3. `analysis/mode_matching.py` (MODIFY) — add `detect_degenerate_groups()` and `compute_subspace_overlap()`; new `match_modes_with_degeneracy()` preserves existing callers
4. `analysis/nist_spectra.py` (NEW) — `NISTFetcher` class: fetch, cache (`~/.cache/mace_gaussian/nist/`), parse JCAMP-DX; returns `SpectrumData | None`; never blocks analysis on absence
5. `batch.py` (MODIFY) — restructure loop to submit DFT after per-molecule geom-opt rather than after all molecules; add file locking to manifest writes
6. `vpt2/` subpackage (NEW, deferred) — self-contained PyVPT2 engine, parallel to `gaussian/`; must produce identical `results.json` schema

### Critical Pitfalls

1. **Lorentzian/Gaussian normalization mismatch** — Lorentzian tails fall as 1/x² vs Gaussian exp(-x²); the same FWHM value produces different peak heights and spectral overlap scores. Prevention: always store broadening type + FWHM in results metadata; enforce identical broadening when comparing two spectra; add unit test checking peak height normalization.

2. **Degenerate mode subspace rotation** — Individual eigenvector dot products for degenerate modes can be 0.3-0.5 even when subspaces are identical (any orthogonal rotation within the subspace is valid). Prevention: for groups within the frequency threshold, compute `trace(M^T M) / k` as the overlap metric rather than per-vector dot products. The existing Hungarian assignment still works for 1-to-1 pairing.

3. **NIST data fragility** — NistChemPy is an unofficial HTML scraper; ~30% of molecules will have no gas-phase IR; JCAMP-DX files mix transmittance (%T) and absorbance and cm-1 vs micron axes. Prevention: treat as best-effort cache; validate units on every retrieval; normalize both spectra to [0,1] before comparison; cache raw JCAMP-DX files locally.

4. **VPT2 Fermi resonance disasters** — Standard VPT2 denominators approach zero for near-degenerate states; anharmonic frequencies can silently diverge by hundreds of cm-1. Prevention: verify PyVPT2 implements GVPT2 (resonance detection + variational correction); add sanity check flagging anharmonic corrections > 300 cm-1; test on CO2 (textbook Fermi resonance) before any larger molecule.

5. **Early SLURM concurrent manifest writes** — With per-molecule DFT submission, two concurrent processes can perform read-modify-write on the batch manifest, causing silent data loss. Prevention: use separate `ml_manifest.json` and `dft_manifest.json`, or add `fcntl.flock` locking around all manifest read-modify-write cycles.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Spectral Quality Foundations
**Rationale:** Lorentzian broadening and zero-intensity filtering are both low-complexity, have no dependencies on each other or on downstream features, and immediately fix known problems (wrong line shape physics, methane R²=0.34 intensity artifact). Everything visual builds on broadening; all intensity statistics build on filtering. Doing these first unblocks every subsequent phase. The acetic acid parser bug fix belongs here because every subsequent analysis feature depends on the parser.
**Delivers:** Physically correct IR spectra plots; meaningful intensity R² and RMSE statistics; acoh parser bug fix (prerequisite for reliable analysis on all molecules)
**Addresses:** Table stakes #1 (Lorentzian broadening), table stakes #2 (zero-intensity filtering), pitfall #12 (acoh parser bug)
**Avoids:** Pitfall #1 (normalization mismatch) — implement both lineshapes in the same file with matching parameter enforcement

### Phase 2: Degenerate Mode Handling
**Rationale:** Degenerate mode handling is medium complexity but isolated to `mode_matching.py`. It is needed before the benchmark campaign to avoid systematically pessimistic mode overlap statistics for symmetric molecules. The algorithm (subspace overlap via trace(M^T M)/k) is already specified in the brainstorming notes. It does not depend on Phase 1 but should land before benchmark runs.
**Delivers:** Correct mode matching for methane (T2, 3-fold), BH3-NH3 (E modes), benzene (E modes); accurate average overlap statistics across the benchmark molecule set
**Addresses:** Table stakes #4 (degenerate mode handling)
**Avoids:** Pitfall #3 (subspace rotation ambiguity) — implement subspace overlap, keep per-vector dot products for reporting only

### Phase 3: Wall-Clock Timing Instrumentation
**Rationale:** Timing is partially implemented in `workflow.py` already. Completing it requires adding `utils/timing.py`, wrapping existing stage calls with the context manager, and extending the HTML report and batch report with timing summary and cost-accuracy Pareto plot. Low risk. Needed before benchmark campaign results are finalized for the thesis.
**Delivers:** Per-stage timing in `results.json`; timing breakdown table and cost-accuracy scatter plot in batch HTML report; speedup factor (DFT_time / ML_time) per molecule
**Addresses:** Table stakes #3 (wall-clock timing)
**Avoids:** Pitfall #6 (GPU warmup / CUDA sync) — add warmup calls and `torch.cuda.synchronize()` before timers

### Phase 4: NIST Experimental Spectra Overlay
**Rationale:** Depends on Lorentzian broadening (Phase 1) for visually meaningful comparison. Independent of degenerate mode handling and timing. Can be developed in parallel with Phases 2-3 if resources allow, but should integrate only after Phase 1 is stable. Uses new `nistchempy` + `jcamp` dependencies.
**Delivers:** `analysis/nist_spectra.py` module with caching; experimental overlay plots in HTML report; peak position comparison table; caveats panel flagging condensed-phase data
**Addresses:** Differentiator #5 (NIST experimental overlay)
**Avoids:** Pitfall #4 (NIST fragility) — graceful None return, local cache, unit validation; Pitfall #11 (JCAMP edge cases) — use `jcamp` library, test on multiple NIST files before automating

### Phase 5: Early SLURM DFT Submission
**Rationale:** Scheduling optimization that saves hours during the 25-molecule benchmark campaign. Restructures `batch.py` loop with no new logic — just reordering existing calls. The main risk is manifest concurrency, which has a clear prevention strategy. Should land before the benchmark campaign starts.
**Delivers:** Per-molecule DFT submission immediately after geom-opt; ML calculations run while DFT queues on cluster; net benchmark campaign time reduction proportional to DFT queue wait
**Addresses:** Differentiator #7 (early SLURM submission)
**Avoids:** Pitfall #7 (concurrent manifest writes) — add `fcntl.flock`; Pitfall #13 (SSH connection limits) — batch `sacct` queries

### Phase 6: Anharmonic Pipeline Integration and Report Overhaul
**Rationale:** Integrates all Phase 1-5 features into a polished HTML report. This is the capstone phase: timing panels, experimental overlay sections, degenerate mode confidence indicators, improved regression plot annotations. Depends on all prior phases being stable.
**Delivers:** Thesis-quality HTML report with all analysis features surfaced; batch report with cost-accuracy Pareto frontier; complete end-to-end pipeline test against all benchmark molecules
**Addresses:** Differentiator #8 (anharmonic pipeline overhaul)
**Avoids:** Integration pitfall (analysis pipeline ordering problem) — design full data flow before implementing; test on water + methane + acoh before declaring done

### Phase 7 (Optional): PyVPT2 Research Integration
**Rationale:** Research-grade work. The QCEngine/ASE adapter does not exist. Fermi resonance handling and ML numerical noise in force constants are open questions. Should be prototyped as a standalone spike on water only before any pipeline integration. May slip to v1.3 without affecting thesis timeline.
**Delivers:** Proof-of-concept PyVPT2+MACE on water; comparison against Gaussian+MACE anharmonic frequencies; feasibility assessment for full integration
**Addresses:** Differentiator #6 (VPT2 alternative)
**Avoids:** Pitfall #9 (QCEngine adapter gap) — 2-hour spike to assess feasibility first; Pitfall #2 (Fermi resonance) — verify GVPT2 support; Pitfall #10 (ML numerical noise) — test step size convergence

### Phase Ordering Rationale

- Phase 1 comes first because Lorentzian broadening is a prerequisite for visually meaningful NIST comparison (Phase 4) and the acoh parser fix unblocks reliable analysis across all molecules.
- Phases 2 and 3 are independent of each other and of Phase 4; they can be developed in parallel after Phase 1 but must land before the benchmark campaign.
- Phase 4 is technically independent of Phases 2-3 but benefits from stable broadening (Phase 1).
- Phase 5 (SLURM reordering) should land before benchmark campaign launch, alongside or after Phase 3.
- Phase 6 is the capstone integration pass; doing it last ensures stable inputs from all prior phases.
- Phase 7 is isolated from all other phases by design and carries the highest risk; decouple it from the milestone timeline entirely.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 4 (NIST Overlay):** NistChemPy is unofficial scraping; real-world coverage and JCAMP-DX format quirks need validation against actual thesis molecules before committing to this in the report. Run a manual test against all 25 benchmark molecules before writing the automated pipeline.
- **Phase 7 (PyVPT2):** Needs dedicated spike before planning. Key questions unanswered: Can Psience ingest Gaussian .fchk force constants? Does PyVPT2's GVPT2 implementation handle CH stretch Fermi resonances correctly? Budget a full week, not 2-3 days.

Phases with standard patterns (skip research-phase):
- **Phase 1 (Broadening + Filtering):** Textbook spectroscopy formulas. Direct codebase changes are fully specified in ARCHITECTURE.md. No unknowns.
- **Phase 2 (Degenerate Modes):** Algorithm is fully specified (trace(M^T M)/k). Isolated to `mode_matching.py`. No dependencies.
- **Phase 3 (Timing):** stdlib only. Context manager pattern is standard. Existing `workflow.py` instrumentation shows the pattern.
- **Phase 5 (SLURM):** Reordering of existing calls. The `slurm.py` API already supports per-molecule submission.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Two new packages (nistchempy, jcamp) well-understood; all other features use existing deps. PyVPT2 exclusion confirmed by Psi4 hard dependency. |
| Features | HIGH | Lorentzian, filtering, degenerate modes, timing are textbook spectroscopy practice with clear implementation paths. NIST overlay is MEDIUM due to scraping fragility. |
| Architecture | HIGH | Based on direct codebase analysis, not external sources. Integration points are specific (file, class, method). Backward-compatibility patterns are clear. |
| Pitfalls | HIGH | Most verified against codebase (existing `time.time()` pattern, acoh xfail test, manifest atomic write pattern) and literature (Fermi resonance, GVPT2 requirement). |

**Overall confidence:** HIGH for phases 1-6. MEDIUM for Phase 7 (PyVPT2) due to QCEngine adapter unknown.

### Gaps to Address

- **Psience/McCoy Group VPT2 feasibility:** Can it ingest Gaussian .fchk cubic/quartic force constants? Needs a 2-hour hands-on spike before Phase 7 planning. If yes, more promising than PyVPT2 (no Psi4 dependency, pure Python). If no, same fundamental blocker.
- **NIST coverage for benchmark molecules:** Not all 25 benchmark molecules will have gas-phase IR in NIST. Before committing Phase 4 to the HTML report, manually check coverage for the full molecule set.
- **xTB dipole unit bug:** Listed as pending in PROJECT.md. Must be fixed before any intensity analysis involving xTB, and before Phase 7 (VPT2 with xTB dipoles). If xTB is not a thesis-critical calculator, block it with an explicit RuntimeError rather than fixing it in v1.2.
- **Acoh parser bug scope:** Commit a4384c4 documents the failure. Needs a focused investigation at the start of Phase 1 before any new parsing code is layered on top.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis: `mace_gaussian/analysis/analyze_spectra.py`, `analysis_workflow.py`, `mode_matching.py`, `batch.py`, `slurm.py`, `workflow.py` — architecture integration points
- [PyVPT2 paper (J. Chem. Phys. 162, 032501, 2025)](https://pubs.aip.org/aip/jcp/article/162/3/032501/3331711) — PyVPT2 capabilities and QCEngine interface; confirmed Psi4 hard dependency
- [NIST Standard Reference Database 35](https://www.nist.gov/srd/nist-standard-reference-database-35) — JCAMP format and coverage (~5,228 gas-phase IR spectra)
- [KTH Computational Chemistry broadening docs](https://kthpanor.github.io/echem/docs/visualize/broadening.html) — Lorentzian vs Gaussian broadening

### Secondary (MEDIUM confidence)
- [NistChemPy GitHub](https://github.com/IvanChernyshov/NistChemPy) — unofficial NIST WebBook scraper; MIT; updated March 2026
- [jcamp PyPI](https://pypi.org/project/jcamp/) — JCAMP-DX parser v1.3.0
- [MACE-OFF23 Composite IR (JCTC 2024)](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01157) — ML+IR benchmark methodology and timing reporting conventions
- [ORCA VPT2/GVPT2 Manual](https://orca-manual.mpi-muelheim.mpg.de/contents/spectroscopyproperties/vpt2.html) — Fermi resonance handling in VPT2
- [AMS Vibrational Spectroscopy docs](https://www.scm.com/doc/AMS/Vibrational_Spectroscopy.html) — degenerate mode handling reference

### Tertiary (LOW confidence)
- [Psience/VPT2 (McCoy Group)](https://github.com/McCoyGroup/Psience) — pure Python VPT2; may accept external force constants; needs hands-on investigation before relying on it
- [Theoretical IR Spectra similarity (JCTC 2020)](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00126) — spectral similarity scoring (SID/SIS); useful if quantitative ML-vs-experiment metrics are needed

---
*Research completed: 2026-03-28*
*Ready for roadmap: yes*
