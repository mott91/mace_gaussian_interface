# Domain Pitfalls: Analysis Quality Overhaul (v1.2)

**Domain:** Computational chemistry IR spectroscopy -- adding Lorentzian broadening, VPT2 alternatives, experimental spectra comparison, degenerate mode handling, timing benchmarks, SLURM parallelism
**Researched:** 2026-03-28
**Confidence:** HIGH (most pitfalls verified against codebase + literature)

---

## Critical Pitfalls

Mistakes that cause wrong results, silent data corruption, or major rework.

### Pitfall 1: Gaussian vs Lorentzian Broadening Confusion

**What goes wrong:** The existing `analyze_spectra.py` uses Gaussian KDE broadening (FWHM 8 cm-1). Adding Lorentzian broadening alongside it creates two problems: (a) forgetting that the lineshape function normalization differs between Gaussian and Lorentzian profiles, and (b) mixing up which broadening was used when comparing ML vs DFT vs experiment.

**Why it happens:** Gaussian broadening normalizes to unit area as `1/(sigma*sqrt(2*pi)) * exp(...)`. Lorentzian normalizes as `(gamma/pi) / ((x-x0)^2 + gamma^2)` where gamma = FWHM/2. The tail behavior is drastically different -- Lorentzian tails fall off as 1/x^2 (heavy tails), Gaussian as exp(-x^2) (light tails). Using Gaussian broadening parameters with a Lorentzian formula produces spectra with wrong peak heights and exaggerated tails.

**Consequences:**
- Spectral overlap metrics (R^2, cosine similarity) become meaningless if ML spectrum uses one broadening and DFT uses another
- FWHM parameter that looks reasonable for Gaussian (8 cm-1) is too narrow for Lorentzian comparison with experiment (typical experimental FWHM is 15-30 cm-1)
- Match scores can be artificially inflated by choosing too-large FWHM, hiding real frequency errors

**Prevention:**
1. Always store the broadening type and FWHM in result metadata (JSON and plot annotations)
2. When comparing two spectra, require identical broadening parameters -- enforce this in the comparison function signature
3. Use Lorentzian for experimental comparison (matches natural linewidth physics) and keep Gaussian KDE for the existing ML-vs-DFT analysis
4. Add a unit test that verifies peak height and area normalization for both lineshape functions
5. Current code at `analyze_spectra.py:98-100` converts FWHM to Gaussian sigma -- add a parallel Lorentzian path, do not retrofit the existing Gaussian path

**Detection:** If computed spectrum peak heights change when switching broadening type but FWHM stays the same, the normalization is wrong.

**Phase:** Lorentzian broadening implementation.

---

### Pitfall 2: VPT2 Fermi Resonance Disasters

**What goes wrong:** VPT2 calculations produce anharmonic frequencies that are hundreds or thousands of cm-1 off. A CH stretching fundamental at 2900 cm-1 suddenly appears at 4700 cm-1 in the output. The calculation "succeeds" with no error -- the numbers just silently go haywire.

**Why it happens:** VPT2 is perturbation theory. When two states are nearly degenerate (e.g., a fundamental and an overtone/combination band have similar energies), the perturbative denominator approaches zero. This is a Fermi resonance. Standard VPT2 does not handle these; it requires GVPT2 (generalized VPT2) or VPT2+K to identify resonant terms, remove them from the perturbative treatment, and handle them variationally.

**Consequences:**
- If PyVPT2 is integrated without Fermi resonance detection, results for CH stretches (the most common Fermi resonance case) will be garbage
- Automated batch processing will propagate bad results into aggregate statistics without anyone noticing
- Comparing VPT2 results against Gaussian's anharmonic (which uses GVPT2) will produce misleading discrepancies that are bugs, not physics

**Prevention:**
1. When integrating PyVPT2, verify it implements GVPT2 (resonance detection + variational correction). The published paper (J. Chem. Phys. 162, 032501, 2025) describes this -- confirm in the actual code
2. Add sanity checks: any anharmonic frequency that deviates more than 300 cm-1 from its harmonic value should be flagged as suspect
3. Test with a known Fermi resonance case (CO2 is the textbook example: v1 fundamental near 2*v2 overtone)
4. Log the Martin test values (resonance detection thresholds) so the user can audit which terms were treated variationally

**Detection:** Anharmonic corrections (anharm - harm) larger than ~200 cm-1 for any single mode. Typical corrections are 20-150 cm-1.

**Phase:** VPT2 integration research and implementation.

---

### Pitfall 3: Degenerate Mode Matching Produces Artificially Low Overlaps

**What goes wrong:** The mode matching system (using eigenvector dot products via Hungarian algorithm in `mode_matching.py`) reports poor overlap (e.g., 0.3-0.5) for modes that are physically correct. The user concludes ML modes are wrong when they are actually fine.

**Why it happens:** Degenerate modes (e.g., methane's T2 bend is 3-fold degenerate) form a subspace where any orthogonal rotation of the eigenvectors within that subspace is equally valid. DFT and ML will typically choose different rotations. The dot product of individual vectors within the subspace can be low even though the subspaces are identical. The brainstorming session notes (2026-03-26) already identified this: "compute subspace overlap (trace of M^T M) instead of individual vector dot products."

**Consequences:**
- Benchmark statistics (average overlap, match rate) are systematically pessimistic for symmetric molecules
- Methane, benzene, ammonia, BH3 -- many thesis molecules have degenerate modes
- Mode matching failures cascade into frequency pairing errors, corrupting MAE/RMSE/R^2 statistics

**Prevention:**
1. Detect degenerate groups: modes within a configurable threshold (e.g., 5 cm-1 for DFT, 10 cm-1 for ML) of each other are candidates
2. For degenerate groups, compute subspace overlap: form the matrix M = V_DFT^T @ V_ML for the group, then overlap = trace(M^T @ M) / n_modes_in_group. A value near 1.0 means the subspaces match regardless of rotation
3. The existing Hungarian assignment (scipy.optimize.linear_sum_assignment, already imported in mode_matching.py) should operate on the subspace-aware overlap matrix, not raw dot products
4. Report both individual overlaps and subspace overlaps so the user can see the difference

**Detection:** Multiple modes with similar frequencies (within threshold) that all show moderate-but-not-great overlap. If you rotate the ML eigenvectors within the degenerate subspace and the overlap jumps to ~1.0, the original matching was being fooled by the rotation ambiguity.

**Phase:** Degenerate mode handling.

---

### Pitfall 4: NIST Data Retrieval Is Fragile and Inconsistent

**What goes wrong:** Automated NIST spectra comparison breaks silently because (a) NIST has no official API, (b) the unofficial NistChemPy library scrapes HTML that changes without notice, (c) spectral resolution varies wildly between compounds, and (d) some JCAMP-DX files have format quirks.

**Why it happens:** NIST Chemistry WebBook is designed for manual browser use, not programmatic access. The data was compiled from many sources over decades. Some spectra are digitized from old paper records. NistChemPy (the best Python option) parses HTML pages -- any website redesign breaks it.

**Consequences:**
- Pipeline fails on ~30% of molecules because NIST does not have gas-phase IR for them
- Retrieved spectra have different x-axis ranges (some are 400-4000, some 600-3800, some use wavelength in microns instead of cm-1)
- Resolution mismatch: NIST gas-phase spectra have ~2-4 cm-1 resolution; computed spectra are stick spectra broadened to whatever FWHM you choose. Direct pixel-wise comparison is meaningless without careful interpolation
- Transmittance vs absorbance confusion: some NIST spectra are in %T, some in absorbance. The conversion is A = -log10(T/100) but forgetting this produces inverted spectra

**Prevention:**
1. Treat NIST retrieval as a best-effort cache, not a pipeline requirement. Never let a missing NIST spectrum block the analysis
2. Store raw JCAMP-DX files locally after first download. Do not re-fetch every run
3. Validate retrieved data: check x-axis units (cm-1 vs micron), y-axis type (absorbance vs transmittance), and spectral range before using
4. Normalize both computed and experimental spectra to [0, 1] range before comparison. Use area-normalized Lorentzian broadening for computed spectra with FWHM matching the experimental resolution (~24 cm-1 per NIST average)
5. Build a manual override table: for key thesis molecules, hand-verify the NIST data once and store the verified version

**Detection:** Correlation between ML and experiment is worse than correlation between DFT and experiment for a molecule where ML-vs-DFT is good. This means the experimental comparison pipeline has a bug, not the ML model.

**Phase:** Experimental spectra comparison.

---

## Moderate Pitfalls

Mistakes that waste time or produce misleading (but not catastrophically wrong) results.

### Pitfall 5: Zero-Intensity Filtering Removes Physically Real Modes

**What goes wrong:** Filtering out "zero-intensity" modes to clean up spectra removes modes that have small but nonzero intensity. The threshold is set too aggressively, and symmetric stretch modes (which can have near-zero IR intensity due to symmetry) disappear from the analysis.

**Why it happens:** In a perfect symmetric molecule, some modes are IR-inactive (intensity exactly 0.0 by symmetry). But ML potentials break exact symmetry slightly, so these modes get small but nonzero intensities (0.01-1.0 km/mol). A threshold of "< 1 km/mol = zero" would remove them. Meanwhile, legitimate weak modes (some bending overtones) also have intensities in the 0.1-5 km/mol range.

**Prevention:**
1. Use a very conservative threshold: only filter modes with intensity exactly 0.0 or below a tiny epsilon (< 0.001 km/mol)
2. Never filter by intensity for mode matching -- match all modes, then annotate intensities
3. Report filtered modes in a separate "suppressed modes" table so the user can audit
4. For experimental comparison, modes below the experimental noise floor (~0.5 km/mol) can be suppressed from the broadened spectrum but should remain in the stick spectrum

**Phase:** Degenerate mode handling and zero-intensity filtering.

---

### Pitfall 6: GPU Warmup and CUDA Sync Make Timing Comparisons Unfair

**What goes wrong:** ML potential appears 2-5x slower than it actually is because the first inference call includes CUDA kernel JIT compilation, cuDNN autotuning, and model weight transfer to GPU. Or the opposite: timing uses CPU wall clock without `torch.cuda.synchronize()`, so GPU operations appear instant because they are still running asynchronously.

**Why it happens:** PyTorch CUDA operations are asynchronous by default. `torch.cuda.synchronize()` must be called before reading the wall clock. Additionally, the first forward pass through a model triggers CUDA context creation (~0.5-2s), kernel compilation, and cuBLAS handle creation. Subsequent calls are much faster.

**Consequences:**
- Timing comparison between ML and DFT is the headline result for the thesis. Getting it wrong by 2-5x undermines the main selling point
- If DFT timing includes SLURM queue wait time, the comparison is even more unfair (in the opposite direction)

**Prevention:**
1. Run 3-5 warmup inference calls before starting the timer. Do not count these
2. Always call `torch.cuda.synchronize()` before `time.perf_counter()` for GPU timing
3. Separate timing into: (a) model loading, (b) geometry optimization, (c) single-point + derivative calculation, (d) Gaussian external interface overhead, (e) total wall clock
4. For DFT, time only the Gaussian computation (parse the CPU/wall time from the log file), not SLURM queue wait
5. Report both "first call" and "amortized per-molecule" timings
6. Use `time.perf_counter()` not `time.time()` (the latter has poor resolution on some systems)

**Detection:** If ML timing varies by more than 50% between the first molecule and subsequent molecules in a batch, warmup was not properly handled.

**Phase:** Wall-clock timing comparison.

---

### Pitfall 7: Early SLURM Submission Creates Race Conditions with Manifest

**What goes wrong:** DFT jobs are submitted to SLURM early (parallel with local ML calculations), but the batch manifest is updated by the local process while SLURM results are being retrieved asynchronously. The manifest gets corrupted or DFT results overwrite ML results (or vice versa).

**Why it happens:** The current `batch.py` uses atomic writes (tempfile + os.replace) for the manifest, which prevents corruption from interrupts. But the design assumes sequential execution: run molecule, update manifest, next molecule. With early SLURM submission, you have two concurrent writers: (a) the local ML pipeline updating manifest as each molecule completes, and (b) the SLURM polling/retrieval updating manifest as DFT results come back. `os.replace` is atomic for a single write but does not prevent lost updates from concurrent read-modify-write cycles.

**Consequences:**
- DFT results that arrived while an ML calculation was in progress get silently dropped from the manifest on the next ML manifest write
- Manifest shows "pending" for a DFT job that actually completed, causing unnecessary re-submissions on restart
- In the worst case, the batch resumes from a stale manifest state and re-runs already-complete ML calculations

**Prevention:**
1. Use separate manifest sections for ML and DFT status, or separate manifest files entirely (ml_manifest.json + dft_manifest.json)
2. If sharing one manifest, use file-level locking (fcntl.flock on Linux) around all read-modify-write cycles
3. The current atomic write pattern in `batch.py:49-69` is good for crash safety but insufficient for concurrency. Add a lock
4. Design the early submission as a clear two-phase protocol: (a) submit all DFT jobs and record job IDs, then (b) run ML calculations, then (c) poll/retrieve DFT results. Phases (a) and (b) can overlap, but manifest writes should not

**Detection:** After a batch run, check that every DFT job ID in the SLURM manifest has a corresponding result directory with .log and .fchk files.

**Phase:** Early SLURM DFT submission.

---

### Pitfall 8: xTB Dipole Unit Bug Propagates into New Analysis

**What goes wrong:** The known xTB dipole unit bug (e*Bohr vs e*Angstrom, ~1.89x error) is not accounted for when adding new analysis features. Lorentzian-broadened spectra from xTB look reasonable because intensity scaling is relative, masking the absolute error. But VPT2 with xTB dipoles produces wrong anharmonic intensities because the derivatives are off by 1.89x.

**Why it happens:** The bug is known but not yet fixed (listed in PROJECT.md as pending). New features may be developed and tested only against MACE calculators, then silently produce wrong results when someone runs xTB.

**Prevention:**
1. Fix the xTB unit bug before integrating xTB into any new analysis pipeline. If the fix is deferred, add an explicit RuntimeError when xTB is selected with VPT2 or intensity-dependent analysis
2. Add a unit test that checks dipole magnitude against a known reference (water dipole ~1.85 Debye)
3. Document the conversion: 1 Debye = 0.3934 e*Angstrom = 0.2082 e*Bohr (a.u.)

**Phase:** Should be resolved before VPT2 integration; at minimum before the benchmark campaign (v1.3).

---

### Pitfall 9: PyVPT2 Requires QCEngine, Not ASE

**What goes wrong:** You assume PyVPT2 can take an ASE calculator directly because MACE uses ASE. But PyVPT2 interfaces through QCEngine (MolSSI's quantum chemistry engine wrapper). Bridging ASE calculators to QCEngine is nontrivial -- there is no existing adapter.

**Why it happens:** The brainstorming notes say "MACE ASE calc -> QCEngine -> PyVPT2" as if this is straightforward. It is not. QCEngine expects programs that speak the QCSchema format (energy, gradient, Hessian as JSON). MACE via ASE returns numpy arrays through the ASE Atoms interface. Someone needs to write the adapter.

**Consequences:**
- The VPT2 integration takes 2-3x longer than expected because the adapter must handle: (a) geometry format conversion (QCSchema vs ASE Atoms), (b) gradient/Hessian format and unit conventions, (c) QCEngine's process model (it wants to launch subprocesses, not call Python functions)
- Alternatively, bypass QCEngine entirely and call PyVPT2's internal functions directly -- but this couples to PyVPT2 internals that may change

**Prevention:**
1. Do a 2-hour spike: install PyVPT2, read its source for the QCEngine interface, determine exactly what adapter code is needed
2. Consider writing a thin QCEngine "harness" for MACE (a Python class that QCEngine can call). This is the intended extension point -- QCEngine has a harness system for exactly this
3. Alternative: compute the required cubic and quartic force constants using ASE's finite difference capabilities, then feed those directly to PyVPT2's VPT2 solver, bypassing the QCEngine dispatch entirely
4. Budget 1-2 weeks for VPT2 integration, not 2-3 days

**Phase:** VPT2 alternative research.

---

### Pitfall 10: VPT2 Numerical Noise from ML Potentials

**What goes wrong:** VPT2 requires cubic and quartic force constant derivatives, computed via finite differences of Hessians. ML potentials have small numerical noise in forces/energies that gets amplified through successive differencing. The resulting anharmonic frequencies fluctuate by 5-20 cm-1 between runs or with different step sizes.

**Why it happens:** VPT2 acts as a "magnifying glass" for numerical noise (per ORCA manual). DFT codes have consistent numerical precision. ML potentials have model-dependent noise that varies with atomic environment. The finite difference step size that works for DFT may be too small for ML, amplifying noise into the force constants.

**Prevention:**
1. Test multiple finite difference step sizes (0.005, 0.01, 0.02, 0.05 Angstrom) and check convergence of anharmonic frequencies
2. Use central differences (not forward differences) -- this cancels first-order numerical error
3. Compare ML VPT2 anharmonic frequencies against Gaussian's ML anharmonic frequencies (same model, different VPT2 implementation) to isolate numerical vs methodological differences
4. Consider tighter geometry optimization convergence criteria before VPT2 -- a geometry that is not at a true minimum will produce garbage force constants

**Detection:** Run VPT2 twice with slightly different step sizes. If results differ by more than 2-3 cm-1, noise is a problem.

**Phase:** VPT2 integration and validation.

---

## Minor Pitfalls

Annoying but not blocking.

### Pitfall 11: JCAMP-DX Parsing Edge Cases

**What goes wrong:** The JCAMP-DX format (used by NIST) has multiple compression schemes (AFFN, PAC, SQZ, DIF, DIFDUP) and vendor-specific extensions. A simple parser handles the common case but fails on ~10% of files.

**Prevention:** Use an existing JCAMP-DX parser library (e.g., `jcamp` on PyPI, or the parser in `scipy.io` if available) rather than writing one from scratch. Test with at least 20 different NIST files before trusting automated retrieval.

**Phase:** Experimental spectra comparison.

---

### Pitfall 12: Acetic Acid Parser Bug Interacts with New Features

**What goes wrong:** The known acetic acid (acoh) frequency parsing bug (commit a4384c4, xfail test) may not just affect regression plots. If the parser fails to correctly separate harmonic and anharmonic sections for acoh, then Lorentzian broadening, mode matching, and VPT2 comparison will all inherit the same bug.

**Prevention:** Fix the acoh parser bug before adding new analysis features that depend on parsed frequencies. At minimum, ensure the xfail test covers the specific parsing failure mode so new code does not regress differently.

**Phase:** Should be fixed early in v1.2, before new analysis features are built on top of the parser.

---

### Pitfall 13: SSH Connection Limits Under Parallel SLURM Polling

**What goes wrong:** When polling many SLURM jobs simultaneously, each poll opens an SSH connection. HPC login nodes typically limit concurrent SSH connections per user (often 10-20). Exceeding this causes "Connection refused" errors that look like network failures.

**Prevention:** The existing `_ssh_with_backoff` in `slurm.py:127` handles retries, which partially mitigates this. Additionally: (a) batch all job status queries into a single SSH command (`sacct --jobs=id1,id2,id3`), (b) rate-limit SSH connections to max 5 concurrent, (c) use SSH connection multiplexing (ControlMaster) to reuse a single TCP connection.

**Phase:** Early SLURM DFT submission.

---

### Pitfall 14: Mace_polar Calculator Returns False Silently

**What goes wrong:** mace_polar already fails for some molecules (noted in PROJECT.md context). When reevaluating it for v1.2, the failure mode is not an exception but a return of `False` from `run_frequency_calculation`. This means batch processing will silently skip mace_polar for failing molecules without logging why.

**Prevention:** Before reevaluating mace_polar, add explicit logging at the failure point. Capture the actual error (likely a model inference error for unsupported elements or geometries) and log it. Then reevaluate with better diagnostics.

**Phase:** Mace_polar dipole calculator reevaluation.

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation | Severity |
|---|---|---|---|
| Lorentzian broadening | Normalization mismatch with existing Gaussian KDE | Keep both, require matching params for comparison | Critical |
| Lorentzian broadening | FWHM too narrow for experiment comparison | Use 15-30 cm-1 for experimental overlay, document choice | Moderate |
| VPT2 (PyVPT2) | QCEngine adapter gap -- no ASE bridge exists | Budget spike time, consider direct force constant input | Critical |
| VPT2 (PyVPT2) | Fermi resonance unhandled | Verify GVPT2 support, add sanity checks on output | Critical |
| VPT2 (PyVPT2) | ML numerical noise amplified in force constants | Test step size convergence, use central differences | Moderate |
| NIST/SDBS comparison | No official API, scraping is fragile | Cache locally, treat as best-effort, use NistChemPy | Moderate |
| NIST/SDBS comparison | Transmittance vs absorbance, cm-1 vs micron | Validate units on every retrieval | Critical |
| Degenerate modes | Subspace rotation produces low dot products | Compute subspace overlap, not individual overlaps | Critical |
| Zero-intensity filter | Overly aggressive threshold removes real modes | Use near-zero threshold (< 0.001), report suppressed | Moderate |
| Timing comparison | GPU warmup inflates ML time | Warmup calls + torch.cuda.synchronize | Moderate |
| Timing comparison | Including SLURM queue time in DFT timing | Parse Gaussian log for actual compute time only | Moderate |
| Early SLURM | Concurrent manifest writes cause lost updates | Separate manifests or file locking | Critical |
| Early SLURM | SSH connection limits on login node | Batch sacct queries, use ControlMaster | Minor |
| Mace_polar reevaluation | Silent False return, no diagnostics | Add logging before reevaluation | Minor |
| Acoh parser bug | New features inherit existing parse failure | Fix parser bug early in v1.2 | Moderate |
| xTB unit bug | Wrong dipole derivatives propagate into VPT2 | Fix or block xTB before VPT2 integration | Moderate |

## Integration Pitfalls (Cross-Feature)

### The Analysis Pipeline Ordering Problem

Multiple v1.2 features modify the analysis pipeline. If added incrementally without a clear data flow design, they create spaghetti:

```
Current:  parse results -> mode match -> statistics -> KDE broadening -> plots -> HTML report
v1.2 adds: parse -> degenerate grouping -> mode match -> zero-intensity filter
                -> Lorentzian broadening -> NIST overlay -> timing annotation
                -> VPT2 comparison -> statistics -> plots -> HTML report
```

**Prevention:** Design the full v1.2 analysis data flow before implementing any single feature. Each feature should be a composable transform on a `SpectrumData` or `AnalysisResult` object, not a monolithic function that does everything. The existing `SpectrumData` dataclass and `ComparisonMetrics` dataclass in `analyze_spectra.py` are the right pattern -- extend them, do not bypass them.

### The "Works for Water, Fails for Everything" Problem

Water (3 atoms, C2v symmetry, no degeneracy, no Fermi resonance, available in NIST) is the easiest possible test case. Every new feature will work perfectly on water and then fail on the first real molecule. Methane has degeneracy. Acetic acid has parsing bugs. Benzene has many modes. Ethanol has conformational flexibility.

**Prevention:** Test every new feature on at least 3 molecules: water (easy), methane (degenerate), and one larger molecule (e.g., acetic acid or ethanol). Do not consider a feature "done" until it handles all three.

## Sources

- [Spectrum broadening -- Computational Chemistry from Laptop to HPC](https://kthpanor.github.io/echem/docs/visualize/broadening.html)
- [Lineshape Functions -- Chemistry LibreTexts](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Spectroscopy/Fundamentals_of_Spectroscopy/Lineshape_Functions)
- [iGVPT2: an interface to computational chemistry (arXiv)](https://arxiv.org/pdf/1704.02144)
- [ORCA VPT2/GVPT2 Manual](https://orca-manual.mpi-muelheim.mpg.de/contents/spectroscopyproperties/vpt2.html)
- [pyVPT2: Interoperable software for anharmonic vibrational frequency calculations (J. Chem. Phys. 2025)](https://pubs.aip.org/aip/jcp/article-abstract/162/3/032501/3331711/pyVPT2-Interoperable-software-for-anharmonic)
- [PyVPT2 GitHub repository](https://github.com/philipmnel/pyvpt2)
- [How to VPT2: Accurate and Intuitive Simulations of CH Stretching (JPCA)](https://pubs.acs.org/doi/abs/10.1021/acs.jpca.0c09526)
- [An Effective and Automated Processing of Resonances in VPT (JPCA 2022)](https://pubs.acs.org/doi/10.1021/acs.jpca.2c06460)
- [NistChemPy: Python API for NIST Chemistry WebBook (PyPI)](https://pypi.org/project/nistchempy/)
- [NistChemPy GitHub](https://github.com/IvanChernyshov/NistChemPy)
- [NIST Standard Reference Database 35 -- JCAMP Format](https://www.nist.gov/srd/nist-standard-reference-database-35)
- [NIST Quantitative IR Database](https://webbook.nist.gov/chemistry/quant-ir/)
- [Vibrational Spectroscopy -- AMS 2025.1 documentation](https://www.scm.com/doc/AMS/Vibrational_Spectroscopy.html)
- [COMSOL: Tracking Eigenmodes over Parametric Sweeps](https://www.comsol.com/blogs/tracking-eigenmodes-over-parameteric-sweeps)
- [PyTorch Benchmark -- Lei Mao](https://leimao.github.io/blog/PyTorch-Benchmark/)
- [SLURM race condition with restart (PyTorch Lightning issue)](https://github.com/Lightning-AI/pytorch-lightning/issues/12369)
