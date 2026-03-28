# Feature Landscape: v1.2 Analysis Quality Overhaul

**Domain:** IR spectroscopy analysis pipeline (computational chemistry)
**Researched:** 2026-03-28
**Context:** Subsequent milestone adding analysis quality features to existing MACE-Gaussian tool

## Table Stakes

Features that the thesis and any reviewer will expect. Missing = results look amateur or untrustworthy.

### 1. Lorentzian Broadening for Simulated IR Spectra

| Aspect | Detail |
|--------|--------|
| Why Expected | Every published computed IR spectrum uses Lorentzian (or Voigt) broadening. The current Gaussian broadening is physically wrong -- IR line shapes are lifetime-broadened (Lorentzian), not Doppler-broadened (Gaussian). A reviewer will notice. |
| Complexity | Low |
| Depends On | Nothing -- drop-in replacement in `analyze_spectra.py:broaden_spectrum()` |
| Files | `mace_gaussian/analysis/analyze_spectra.py` |

**Standard practice:**
- **Lorentzian** is the correct physical model for IR absorption (homogeneous broadening from finite excited-state lifetimes). Formula: `L(v) = I * (gamma/2)^2 / ((v - v0)^2 + (gamma/2)^2)` where gamma = FWHM.
- **FWHM values:** 4 cm-1 for high-resolution gas-phase comparison, 8-10 cm-1 for general computed spectra, 20-30 cm-1 for matching condensed-phase/NIST experimental spectra. The MACE-OFF23 composite IR paper uses scaled frequencies with ~10 cm-1 broadening. NIST average bandwidth is ~24 cm-1.
- **Voigt profile (Lorentzian + Gaussian convolution)** is more physically accurate for gas-phase spectra where both pressure broadening (Lorentzian) and Doppler broadening (Gaussian) contribute. However, for computed spectra comparison, pure Lorentzian is standard and sufficient. Voigt is overkill and adds complexity without thesis value.
- **Implementation:** Replace the Gaussian kernel in `broaden_spectrum()` with Lorentzian. Keep the same FWHM parameter. Add `--broadening-type lorentzian|gaussian` CLI flag (default: lorentzian). Keep Gaussian as option for backward compatibility.
- **Current state:** Uses `exp(-B * (v - v0)^2)` Gaussian broadening with FWHM=8 cm-1. Needs to become `gamma^2 / ((v-v0)^2 + gamma^2)` Lorentzian.

**Recommendation:** Default FWHM=10 cm-1 Lorentzian for computed spectra plots; 4 cm-1 for high-res; 24 cm-1 when overlaying with NIST experimental data. Make it a parameter, not hardcoded.

**Confidence:** HIGH -- this is textbook spectroscopy, well-documented everywhere.

### 2. Zero-Intensity Mode Filtering

| Aspect | Detail |
|--------|--------|
| Why Expected | Including IR-inactive modes (zero intensity) in intensity regression is statistically meaningless. Methane's R^2=0.34 for intensity is caused by correlating zeros with numerical noise. Any reviewer doing a sanity check on the intensity statistics will flag this. |
| Complexity | Low |
| Depends On | Nothing -- post-processing filter on existing matched mode pairs |
| Files | `mace_gaussian/analysis/analyze_spectra.py`, `mace_gaussian/analysis/analysis_workflow.py` |

**Standard practice:**
- **Threshold:** 0.1 km/mol is the standard cutoff. Modes below this are effectively IR-inactive (symmetry-forbidden or near-zero transition dipole moment). Some papers use 1.0 km/mol for stricter filtering.
- **When to filter:** Match ALL modes first using eigenvector overlap (intensity-agnostic). Then for frequency regression, include all matched pairs. For intensity regression only, exclude pairs where the DFT reference intensity < threshold.
- **Reporting:** Always report how many modes were filtered and what fraction of total. A table showing "N modes matched, M filtered (intensity < 0.1 km/mol), K used for intensity metrics" is expected.
- **Effect on statistics:** R^2 for intensity will jump dramatically (the zeros are dominating the correlation). MAE for intensity will increase (you removed the easy-to-match zeros). Both changes are correct -- the filtered stats are meaningful, the unfiltered ones are not.
- **Implementation:** Add `min_dft_intensity` parameter to `ComparisonMetrics` calculation. Report both filtered and unfiltered metrics in the HTML report so the reader can see the effect.

**Confidence:** HIGH -- well-established practice, directly addresses known methane R^2 bug.

### 3. Wall-Clock Timing Comparison (ML vs DFT)

| Aspect | Detail |
|--------|--------|
| Why Expected | "X% accuracy at Y% of the computational cost" is THE key result for any ML-vs-DFT benchmarking paper. Without timing data, the thesis has accuracy comparisons but no cost-benefit analysis. |
| Complexity | Low-Medium |
| Depends On | Existing `time.perf_counter()` instrumentation in `workflow.py` (already partially there) |
| Files | `mace_gaussian/workflow.py`, `mace_gaussian/batch.py`, `mace_gaussian/analysis/batch_report.py` |

**What to measure:**
- **Per-stage wall-clock time:** geometry optimization, frequency calculation (Gaussian+ZMQ), analysis. Already partially instrumented with `time.perf_counter()`.
- **Per-molecule total time:** ML pipeline end-to-end vs DFT baseline end-to-end. The DFT time includes SLURM queue wait (report separately: "DFT compute time" vs "DFT wall-clock with queue").
- **Per-calculator-combination:** time for mace_mp/espaloma vs mace_omol/mace_ml etc.
- **Scaling with molecule size:** plot time vs number of atoms. ML should scale much better than DFT for larger molecules.
- **Key metric:** Speedup factor = DFT_time / ML_time. Report as "Nx speedup" in the batch report.

**Standard reporting format:**
- Table: molecule, N_atoms, DFT_time, ML_time, speedup
- Plot: accuracy (R^2 or RMSE) vs computational cost (Pareto frontier)
- Both are expected in ML+spectroscopy papers (ANI+VPT2, MACE-OFF23 composite IR papers include these)

**Implementation:** Capture timing in the results JSON (already partially done). Aggregate in batch_report. Add Pareto plot (accuracy vs cost) to HTML report.

**Confidence:** HIGH -- straightforward instrumentation, standard reporting.

### 4. Degenerate Mode Handling

| Aspect | Detail |
|--------|--------|
| Why Expected | Without degenerate mode awareness, the mode matching produces misleading low overlaps for high-symmetry molecules (methane, benzene, ammonia). The brainstorming session already identified this as a known problem. |
| Complexity | Medium |
| Depends On | Existing Hungarian matching in `mode_matching.py` |
| Files | `mace_gaussian/analysis/mode_matching.py`, `mace_gaussian/analysis/analyze_spectra.py` |

**Standard practice for degenerate modes:**
- **Detection:** Group modes whose DFT reference frequencies are within a tolerance (typically 5 cm-1 for exact degeneracies, up to 20 cm-1 for near-degeneracies caused by symmetry breaking in ML potentials).
- **Subspace overlap:** For a k-fold degenerate group, the individual eigenvector dot products are meaningless because any orthogonal rotation within the degenerate subspace is equally valid. The correct metric is the subspace overlap: form the overlap matrix M (k x k) between DFT and ML eigenvectors within the group, then compute `trace(M^T @ M) / k`. This gives 1.0 if the subspaces are identical regardless of eigenvector orientation.
- **Hungarian matching still works:** Keep the Hungarian assignment for 1-to-1 mode pairing (needed for frequency/intensity correlation). But report confidence as the subspace overlap for degenerate groups, not individual dot products.
- **Practical degeneracies in benchmark set:**
  - Methane (Td): T2 modes are 3-fold degenerate (bending ~1300 cm-1, stretching ~3000 cm-1)
  - BH3-NH3 (C3v): E modes are 2-fold degenerate
  - Benzene (D6h): E modes are 2-fold degenerate
  - Water, formaldehyde, acetic acid: no degeneracies (C2v or Cs)
- **Implementation:** Add `detect_degenerate_groups()` and `compute_subspace_overlap()` functions. Modify overlap heatmap to shade degenerate blocks. Add per-group confidence to HTML report.

**Confidence:** HIGH -- the brainstorming notes already have the correct algorithm (trace of M^T M). This is standard in vibrational spectroscopy literature.


## Differentiators

Features that go beyond table stakes and make the thesis stand out. Not expected, but valued.

### 5. NIST/SDBS Experimental Spectra Overlay

| Aspect | Detail |
|--------|--------|
| Value Proposition | Comparing ML vs DFT only shows how well ML reproduces DFT, not reality. Adding experimental overlay answers "how close are we to actual measured spectra?" -- a much stronger thesis result. |
| Complexity | Medium |
| Depends On | Lorentzian broadening (#1) for proper visual comparison |
| Files | New module: `mace_gaussian/analysis/experimental_spectra.py` |

**How other tools do it:**
- **GaussView/Avogadro:** Load experimental spectrum as a file (JCAMP-DX, CSV, or image), overlay on computed spectrum. No automatic fetching -- manual file import. Peak matching is visual only.
- **AMS (Amsterdam Modeling Suite):** Has built-in NIST comparison with automatic peak assignment. Gold standard but proprietary.
- **Most papers:** Manual overlay in plotting software. Very few tools automate this.

**Data sources:**
- **NIST WebBook (primary):** Gas-phase IR spectra in JCAMP-DX format. 5,228 compounds. Gas-phase is ideal because computed spectra are in vacuum (no solvent effects). Use `nistchempy` Python package for programmatic access (search by name/CAS/InChI, download JCAMP-DX).
- **SDBS (secondary):** Japanese spectral database, broader coverage but mostly condensed-phase. Use only when NIST lacks gas-phase data.
- **JCAMP-DX parsing:** Use `jcamp` Python package (`pip install jcamp`). Returns x (wavenumber), y (transmittance or absorbance) arrays.

**Implementation plan:**
1. `fetch_experimental_spectrum(molecule_name)` -- search NIST via nistchempy, download JCAMP-DX, cache locally in `experimental_spectra/` directory.
2. `parse_jcampdx(filepath)` -- convert transmittance to absorbance if needed, interpolate onto same frequency grid.
3. `overlay_experimental(computed_spectrum, experimental_spectrum)` -- plot both on same axes with proper normalization.
4. **Quantitative comparison:** Peak position matching (find peaks in both, compute position errors in cm-1). Spectral similarity score (e.g., Pearson correlation of broadened spectra, or the SID/SIS scores from Theoretical IR Spectra paper by Laurens et al. JCTC 2020).
5. **Caveats panel in report:** Flag when experimental data is condensed-phase (expected systematic shifts for O-H: -300 to -550 cm-1, C=O: -10 to -25 cm-1, C-H: negligible).

**Key gotcha:** NIST transmittance spectra need conversion to absorbance (A = -log10(T)) or plotted inverted. Most computed spectra are in absorbance/intensity units.

**Confidence:** MEDIUM -- nistchempy exists and works, but automated peak matching between computed and experimental spectra is nontrivial (different resolutions, baselines, units). The fetching and overlaying is straightforward; quantitative metrics need care.

### 6. VPT2 Alternatives (PyVPT2 Integration)

| Aspect | Detail |
|--------|--------|
| Value Proposition | Eliminates Gaussian 16 dependency for anharmonic calculations. Enables fully open-source pipeline. Nobody has published PyVPT2+MACE foundation models -- this is novel. |
| Complexity | High |
| Depends On | MACE calculators already working (they are) |
| Files | New module: `mace_gaussian/vpt2/pyvpt2_backend.py`, modifications to `workflow.py` |

**What PyVPT2 needs:**
- **Input:** Optimized molecular geometry + a "program" that can compute energies, gradients, and Hessians at displaced geometries. PyVPT2 generates the displaced geometries internally.
- **QCEngine harness:** PyVPT2 talks to quantum chemistry codes via QCEngine, which uses a standardized QCSchema (JSON-like) input/output format. To use MACE, you need to write a QCEngine "harness" -- a thin adapter that translates QCSchema requests into ASE calculator calls and returns QCSchema results.
- **What the harness must provide:** `energy(geometry)`, `gradient(geometry)`, and `hessian(geometry)` calls. MACE ASE calculators already provide all three (via `atoms.get_potential_energy()`, `atoms.get_forces()`, and finite-difference Hessian via ASE's `Vibrations` module).
- **Force constants:** PyVPT2 computes cubic and quartic force constants by numerical differentiation of Hessians at displaced geometries. This means many single-point calculations (for N atoms: O(N) displacements for cubic, O(N^2) for quartic). With ML potentials this is fast (seconds per point vs minutes for DFT).
- **Dipole derivatives:** For IR intensities, PyVPT2 needs dipole moment at each displaced geometry. MACE dipole calculator already provides this.

**Integration path:**
1. Write QCEngine harness for MACE (translates QCSchema -> ASE calculator calls)
2. Test on water: compare PyVPT2+MACE vs Gaussian+MACE (should give identical VPT2 physics, different implementation)
3. Add `--vpt2-engine gaussian|pyvpt2` flag to workflow
4. Benchmark: PyVPT2 will be slower per-displacement (Python vs Fortran) but eliminates Gaussian license requirement

**Risk:** QCEngine harness development may be nontrivial. The QCEngine project is MolSSI-maintained but adding a new harness requires understanding their schema deeply. Also, PyVPT2 is young (published Jan 2025) -- may have bugs or missing features for larger molecules.

**Prior art:** PhysNet+VPT2 (Meuwly 2021), ANI+VPT2 (2025), MACE-OFF23+VPT2 (2024). All used Gaussian or custom code for VPT2, not PyVPT2.

**Confidence:** MEDIUM -- PyVPT2 exists and is published, but nobody has done PyVPT2+MACE. The QCEngine harness is the unknown. This is research-grade work, not engineering.

### 7. Early SLURM DFT Submission (Parallel with ML Calcs)

| Aspect | Detail |
|--------|--------|
| Value Proposition | Currently DFT runs after geometry optimization. Since DFT does its own optimization, submitting DFT to cluster simultaneously with ML calcs saves wall-clock time. For a 25-molecule benchmark, this could save hours of queue wait. |
| Complexity | Medium |
| Depends On | Existing SLURM infrastructure in `slurm.py` |
| Files | `mace_gaussian/workflow.py`, `mace_gaussian/batch.py`, `mace_gaussian/slurm.py` |

**Implementation:**
- After geometry optimization (stage 1), immediately submit DFT job to SLURM before starting ML frequency calculations (stage 3).
- ML calculations proceed locally while DFT job queues and runs on cluster.
- At end of ML calculations, poll for DFT completion (or wait if not done yet).
- This is a scheduling optimization, not a new capability.

**Confidence:** HIGH -- straightforward reordering of existing pipeline stages.

### 8. Anharmonic Analysis Pipeline Overhaul

| Aspect | Detail |
|--------|--------|
| Value Proposition | The current anharmonic analysis (overtones, combination bands) has known parsing issues (acetic acid bug) and the report format needs improvement for thesis-quality output. |
| Complexity | Medium |
| Depends On | Lorentzian broadening (#1), degenerate mode handling (#4), zero-intensity filtering (#2) |
| Files | `mace_gaussian/analysis/analyze_spectra.py`, `mace_gaussian/analysis/analysis_workflow.py`, `mace_gaussian/analysis/html_report_generator.py` |

**What needs to change:**
- Fix acetic acid overtone/combination band parsing (known bug, commit a4384c4)
- Integrate all new features (Lorentzian, degenerate mode awareness, filtering) into the HTML report
- Add experimental overlay panel when NIST data available
- Add timing summary panel
- Improve regression plot annotations (show filtered vs unfiltered metrics)
- Add degenerate mode confidence indicators

**Confidence:** HIGH -- this is integration work building on the other features.


## Anti-Features

Features to explicitly NOT build for v1.2.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Voigt profile broadening | Adds complexity for negligible benefit in computed spectra. Lorentzian is sufficient and standard. Voigt matters for experimental spectral fitting, not for plotting computed stick spectra. | Use pure Lorentzian. Mention Voigt as future work if reviewer asks. |
| Automated peak assignment labels | Functional group labeling (e.g., "C-H stretch") requires a database of characteristic ranges and is error-prone for non-standard molecules. This is a rabbit hole. | Let the user interpret peaks. Include frequency range reference in docs. |
| Interactive HTML spectrum viewer (Plotly) | Nice-to-have but significant effort for a feature that doesn't improve thesis results. Static plots are fine for the paper. | Keep matplotlib PNG plots. Revisit for v2.0 if there's demand. |
| Conformer-aware spectra (Boltzmann weighting) | Requires conformer search + Boltzmann weighting. Separate research problem. Most benchmark molecules are rigid or have one dominant conformer. | Use single lowest-energy conformer. Note limitation in thesis. |
| JCAMP-DX export of computed spectra | Nobody will import your computed spectra into their instrument software. The value is in comparing, not exporting. | Export as CSV/JSON (already done). |
| Delta-ML correction model | Training a correction model (ML_corrected = ML + delta) requires more data than we have. Premature optimization. | Report raw ML vs DFT errors. Leave delta-ML for PhD work. |
| Uncertainty quantification for ML spectra | Requires committee models or MC dropout. Separate research direction (noted in brainstorming as PhD topic). | Report accuracy statistics (MAE, RMSE) as the uncertainty measure. |


## Feature Dependencies

```
Lorentzian Broadening (#1)
    |
    +---> NIST Overlay (#5)  [needs proper broadening for visual comparison]
    |
    +---> Anharmonic Pipeline Overhaul (#8)  [integrates all improvements]

Zero-Intensity Filtering (#2) ---> Anharmonic Pipeline Overhaul (#8)

Degenerate Mode Handling (#4) ---> Anharmonic Pipeline Overhaul (#8)

Wall-Clock Timing (#3) ---> Anharmonic Pipeline Overhaul (#8)  [timing panel in report]

Early SLURM (#7) --- independent, no downstream deps

VPT2/PyVPT2 (#6) --- independent, parallel research track
```

## MVP Recommendation

**Phase order for v1.2:**

1. **Lorentzian broadening** + **Zero-intensity filtering** -- both Low complexity, immediate quality improvement, unblock everything downstream.
2. **Degenerate mode handling** -- Medium complexity, fixes known methane/BH3NH3 issues, needed before benchmark campaign.
3. **Wall-clock timing** -- Low-Medium, partially implemented, needed for thesis cost-benefit argument.
4. **Early SLURM submission** -- Medium, saves real time during benchmark campaign.
5. **NIST experimental overlay** -- Medium, strongest differentiator for thesis.
6. **Anharmonic pipeline overhaul** -- integrates 1-5 into polished report.

**Defer to v1.3 or later:**
- **PyVPT2 integration (#6):** High complexity, research-grade uncertainty. Prototype as a side experiment but don't block v1.2 on it. The Gaussian-based pipeline works and is proven. PyVPT2 is the "open-source future" story, not the "thesis needs this now" story.

## Sources

- [Spectrum broadening (KTH computational chemistry)](https://kthpanor.github.io/echem/docs/visualize/broadening.html) -- Lorentzian vs Gaussian broadening explanation
- [MACE-OFF23 Composite IR (JCTC 2024)](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01157) -- ML+IR benchmark methodology and timing
- [PyVPT2 paper (J. Chem. Phys. 2025)](https://pubs.aip.org/aip/jcp/article/162/3/032501/3331711/pyVPT2-Interoperable-software-for-anharmonic) -- PyVPT2 capabilities and QCEngine interface
- [PyVPT2 GitHub](https://github.com/philipmnel/pyvpt2) -- Source code and examples
- [NistChemPy (PyPI)](https://pypi.org/project/nistchempy/) -- Python API for NIST WebBook
- [jcamp (PyPI)](https://pypi.org/project/jcamp/) -- JCAMP-DX parser
- [NIST/EPA Gas-Phase IR Database](https://www.nist.gov/publications/nist-35-nistepa-gas-phase-infrared-database-jcamp-format) -- 5,228 gas-phase IR spectra
- [Theoretical IR Spectra similarity measures (JCTC 2020)](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00126) -- Spectral similarity scoring methods
- [AMS Vibrational Spectroscopy docs](https://www.scm.com/doc/AMS/Vibrational_Spectroscopy.html) -- Reference for degenerate mode handling
