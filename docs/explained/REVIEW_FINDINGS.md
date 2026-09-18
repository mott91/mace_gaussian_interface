# Code Review Findings (2026-09-18)

Read-only review of `mace_gaussian/` (41 files, 14.3k lines). Nothing was changed.
Each finding has: where, what, evidence, why it matters, and a suggested fix.
We go through these together and decide what to do. Nothing here is fixed yet.

Severity scale:
- **H** changes numbers that could end up in the thesis, or a thesis conclusion
- **M** wrong data stored, or a failure that is hidden instead of reported
- **L** robustness, dead code, docs out of sync

Line numbers refer to the state on branch `spike/24-vpt2-psience` at commit f929741.

---

## H1. Frequencies for four of five energy models are computed at another model's minimum

**Where:** `mace_gaussian/workflow.py:747` (default energy list), `workflow.py:485-486`
(the same `optimized_atoms` is reused for every energy model), `mace_gaussian/gaussian/io.py:154`
(route is `# freq (anharm)`, no `opt`).

**What:** Stage 1 optimizes the geometry with `mace_omol` only. Stage 3 then runs
`freq(anharm)` with `mace_mp`, `mace_off`, `mace_anicc`, `mace_polar` at that same geometry.
For those four models the geometry is not a minimum of their own energy surface.

**Evidence (water, Gaussian log line 1371, "Maximum Force" vs threshold 0.000450):**

| energy model | max force (Ha/Bohr) | converged? |
|---|---|---|
| mace_omol | 0.000000 | yes |
| mace_off | 0.003245 | no |
| mace_anicc | 0.004198 | no |
| mace_mp | 0.012873 | no (28x over threshold) |

Methane with `mace_mp` (`comparison_results/methane/mace_mp_espaloma/gaussian_freq.log`):
three imaginary modes at -274 cm-1 (line 191), plus `WARNING: Unreliable CUBIC force constant`
(lines 1037-1041) and a rotor/framework inconsistency warning (line 1001).

**Why it matters:** Gaussian's harmonic frequencies at a non-stationary point are contaminated
by the gradient (Gaussian does not project it out unless `freq=projected`). VPT2 assumes you
sit at a minimum; cubic constants along the gradient direction are meaningless, which is what
the "Unreliable CUBIC" warning says. The "mace_mp fails on methane bends (-193 / -274 cm-1)"
story in `thesis/MOLECULES.md` may be partly a geometry artifact, not a PES failure. This
needs to be settled before it goes into Chapter 5.

`docs/methods.md` §2 presents the shared geometry as a deliberate choice ("differences arise
from the PES, not the geometry"). That is a defensible design, but only if the thesis says so
explicitly and shows the residual forces. The DFT twin re-optimizes at B3LYP anyway
(`dft_baseline.py:183`, route `# opt freq(anharm)`), so DFT and ML are already not at the
same geometry.

**Suggested fix:** Re-optimize with each energy model in ASE before its Gaussian run
(seconds of GPU time), and record the residual max force in `results.json`. Alternatively keep
the shared geometry but add `freq=projected` and state the choice in Chapter 4.

**How to verify:** Re-run methane with `--optimization-calculator mace_mp
--energy-calculators mace_mp`. If the -274 cm-1 modes disappear, H1 is confirmed.

**Status 2026-09-18: FIXED on branch `fix/review-2026-09`** (`workflow.py`: LBFGS on the
pair's energy model before the gjf is written; `OPT_FMAX`/`OPT_MAX_STEPS` module constants;
`calculation_parameters.reoptimization` in results.json records model, steps, convergence,
max force). Verified with MACE-MP on methane and water at MP's own minimum
(max force 5e-8 and 2e-7 eV/Å, 4 and 8 LBFGS steps):

| | DFT | MP before (OMOL geom) | MP after (own min) |
|---|---|---|---|
| methane rotations ("Low frequencies" line 1) | ~0 | -274 ×3 | -0.19 ×3 |
| methane bend (T2) harm | 1356 | 1163 | 1181 |
| methane bend (E) harm | 1579 | 1410 | 1415 |
| methane C-H sym stretch harm | 3047 | 3076 | 2988 |
| methane C-H asym stretch (T2) harm | 3162 | 3163 | 3077 |
| water bend harm / VPT2 | 1665 / 1615 | 1462 / 1445 | 1497 / 1480 |
| water sym stretch harm / VPT2 | 3799 / 3624 | 3852 / 3685 | 3696 / 3564 |
| water asym stretch harm / VPT2 | 3912 / 3722 | 3975 / 3780 | 3815 / 3655 |

Reading: the imaginary triple was the three rotations carrying the residual gradient; it is
gone. The soft bends are **real** (still -170 to -175 cm-1 vs DFT at MP's own minimum), so
the "crystal-trained model fails on molecular bends" story survives. What changes is the
stretches: at the OMOL geometry MP sat on its repulsive wall and looked *stiffer* than DFT
(+50 to +60 cm-1); at its own minimum it is *softer* than DFT by 90 to 100 cm-1. So the
"before" data mis-stated the sign of MP's stretch error. Every non-OMOL ML run in
`comparison_results/` is affected to some degree (residual forces 0.003 to 0.013 Ha/Bohr)
and must be rerun before any number goes into the thesis. `thesis/MOLECULES.md` still quotes
the old -193 figure.

The after-run `results.json` files are kept in `docs/explained/h1_check/`.

**Observation O1 (physics, not a bug): ML surfaces trip Gaussian's cubic-constant
consistency check even at a proper minimum.** `Unreliable CUBIC force constant` warnings in
the methane logs: DFT 0, mace_omol 2, mace_mp at OMOL geometry 3, mace_mp at its own minimum
8. Gaussian gets each cubic constant φᵢⱼₖ from several displacement pairs and warns when they
disagree. With DFT they agree; with ML Hessians they do not always, which means the ML
Hessian changes less smoothly between displaced geometries than DFT's, i.e. Hessian noise
relative to Gaussian's step size. This is exactly the point raised in `03_hard_questions.md`
A1 and it needs (a) a parser that counts these warnings into results.json and (b) one
step-size sensitivity test on methane before the campaign. Which constants: for MP the
warnings involve modes 1, 4, 5, 6 (the A₁ and T₂ stretches), so it is stretch-stretch
coupling, not the bends.

---

## H2. Anharmonic mode matching applies the eigenvector map in the wrong index space

**Where:** `mace_gaussian/analysis/analyze_spectra.py:411-439` (`match_by_mode` remap),
`analysis_workflow.py:417-421` (mapping keyed by checkpoint index),
`gaussian/parser.py:162` (`mode` = Gaussian's anharmonic `Mode(n)` number).

**What:** The eigenvector mapping is built from `.fchk` modes, which are stored in ascending
frequency order (index 0 = lowest). The anharmonic fundamentals come from the log's
"Fundamental Bands" table, where Gaussian numbers modes by symmetry block, not by frequency.
`match_by_mode` treats the anharmonic mode number as if it were the checkpoint index.

**Evidence (water, B3LYP):**

| source | order |
|---|---|
| `.fchk` Vib-E2 (ascending) | 1665.3, 3799.2, 3912.4 |
| log "Fundamental Bands" Mode(1), (2), (3) | 3799.2, 1665.3, 3912.4 |

Methane: Mode(1) = 3046.5 (checkpoint index 5), Modes 7-9 = 1356.2 (indices 0-2).

**Why it matters:** When the Hungarian mapping is the identity (all small molecules run so far),
nothing breaks, because F1 maps to F1. When the mapping is a real permutation, which is the
case mode matching exists for, the remapped IDs point at the wrong DFT modes in every
anharmonic report. Overtone (`O{m}_{l}`) and combination (`C{m1}_{m2}`) IDs are never
remapped at all (`analyze_spectra.py:433-435`), so they pair by raw Gaussian numbering.

The harmonic analysis (`run_analysis_harmonic.py`) is **not** affected: both sides come from
`.fchk` in the same index space.

**Suggested fix:** Each anharmonic entry carries `freq_harmonic` (parser.py:170, 195).
Translate Gaussian mode number to checkpoint index by matching `freq_harmonic` to the
`.fchk` harmonic list (they agree to 1e-3 cm-1), then apply the mapping. Do the same
translation for overtone/combination indices.

**How to verify:** Pick a molecule where the harmonic report shows a non-identity mapping
(check `Mode_Overlap` column order in `analysis_results_harmonic/<mol>/data/comparison_*.csv`)
and compare the anharmonic `comparison_*.csv` pairs by hand.

**Status 2026-09-18: FIXED on `fix/review-2026-09`.** New
`analyze_spectra.gaussian_mode_to_checkpoint_index()` ranks the anharmonic rows by
`freq_harmonic` and `extract_spectrum_data` labels fundamentals, overtones and
combinations in checkpoint index space. Old result files without `freq_harmonic` fall
back to the raw number. Tests in `tests/test_review_fixes.py`. Before/after comparison of
the analysis output: see the table at the end of this file.

Left open as **H2b**: `html_report_generator._build_anharmonicity_section` (lines 685-712)
pairs ML and DFT anharmonic rows by raw Gaussian mode number with no eigenvector mapping at
all. Same bug class; it feeds only the anharmonicity-ratio plot, not the metrics.

---

## H3. Mode overlaps use Cartesian displacements, not mass-weighted eigenvectors

**Where:** `mace_gaussian/analysis/mode_matching.py:128-162` (`compute_mode_overlap`),
`:224-250` (`create_alignment_matrix`). `masses` is returned by `gaussian/fchk.py:268` and
never used for weighting. `compute_reduced_masses` (`mode_matching.py:79`) exists but is unused.

**What:** Gaussian's `Vib-Modes` block stores Cartesian displacement vectors normalized to 1.
These are not orthogonal to each other. The true normal-mode eigenvectors are orthonormal in
mass-weighted coordinates (multiply each atom's displacement by sqrt(m)).

**Evidence (self-overlap of a calculation with itself, max off-diagonal element):**

| molecule | Cartesian | mass-weighted |
|---|---|---|
| water B3LYP | 0.052 | 9e-12 |
| methane B3LYP | 0.126 | 7e-10 |
| methane mace_omol | 0.126 | 9e-10 |

**Why it matters:** Reported "overlap" values are not projections onto an orthonormal basis, so
cross terms are inflated by up to ~0.13 for H-rich molecules. The 0.5 confidence threshold, the
degenerate subspace overlap trace(M^T M)/k (`mode_matching.py:665-693`), and the heatmaps are
all in the wrong metric. Assignments are probably still right in most cases (the diagonal
dominates), but `docs/methods.md` §6 and thesis §2.4 claim "mass-weighted normal-mode
eigenvectors", which is currently false.

**Suggested fix:** In `extract_mode_data_from_checkpoint` or at the top of `match_modes`:
`modes = modes * np.sqrt(masses)[None, :, None]`, then renormalize each mode. Two lines.

**Status 2026-09-18: FIXED on `fix/review-2026-09`** in
`mode_matching.extract_mode_data_from_checkpoint`. Test: self-overlap of the water and
methane fixture checkpoints is the identity to 1e-6.

---

## H4. `parse_final_energy` returns a thermal correction, not the energy

**Where:** `mace_gaussian/gaussian/parser.py:352` (`Energy=\s+([-\d\.]+)`).

**What:** Gaussian prints `Thermal correction to Energy=` and `Thermal correction to Gibbs Free
Energy=` with spaces after `=`. The regex matches those and takes the last one. The external
`Energy= -76.4 NIter=` lines and the DFT `SCF Done:` line are not what gets returned.

**Evidence:** Parser output for the three logs checked:

| log | returned (Ha) | actual |
|---|---|---|
| water B3LYP | 0.003705 | -76.4196 (SCF Done, line 1593) |
| water mace_omol+mace_ml | 0.003678 | about -76.4 |
| methane B3LYP | 0.027695 | about -40.5 |

Every `energy_eV` in `comparison_results/water/*/results.json` is 0.096-0.101 eV, which is
the Gibbs correction times 27.2.

**Why it matters:** Frequencies and intensities are unaffected. But all stored energies are
wrong, `mace-gaussian list` prints them, and any future "ML vs DFT energy" statement would be
garbage. Cheap to fix, embarrassing to leave.

**Suggested fix:** For external logs use `^\s*Energy=\s+(-?[\d.]+)\s+NIter`; for DFT logs use
`SCF Done:\s+E\(\w+\)\s+=\s+(-?[\d.]+)`; try both, take the last match of whichever hits.

**Status 2026-09-18: FIXED on `fix/review-2026-09`** in `parser.parse_final_energy`.
Note: `tests/test_gaussian_parser.py` had pinned the wrong value (`WATER_ENERGY_EXPECTED =
0.003705`, i.e. the Gibbs correction); that test now asserts `None` for the truncated DFT
fixture, and the real cases are in `tests/test_review_fixes.py`. Existing `results.json`
energies are still wrong until re-parsed or rerun.

---

## M1. Dipole failures are silently converted to zero intensities

**Where:** `mace_gaussian/workflow.py:221-228` (`except Exception` → zeros, run continues),
`mace_gaussian/calculators/base.py:71-72` (exception inside the finite-difference loop is
logged and the partially filled zero array is returned).

**What:** If the dipole model throws at any Gaussian call (CUDA OOM, RDKit bond perception
failure, shape mismatch), Gaussian receives zero dipole derivatives and finishes normally.
`run_frequency_calculation` returns `True`, the manifest says `complete`, and the spectrum
has zero intensities for that call, which after VPT2 assembly means partially or fully wrong
intensities with no flag in `results.json`.

**Suggested fix:** Re-raise (the runner already sends "error" to Gaussian on exception,
`runner.py:105-107`), or at minimum count fallbacks and write the count into `results.json`
so the analysis can refuse the intensities.

---

## M2. The Gaussian timeout can never fire while Gaussian is silent

**Where:** `mace_gaussian/gaussian/runner.py:82-99`, `gaussian/zmq_server.py:104-111`.

**What:** `is_calc_finished` loops internally (`while True`, `sleep(1)`) until either a ZMQ
message arrives or the process exits. The elapsed-time check in `run_gaussian_with_zmq`
only runs between messages. If Gaussian hangs inside its own code (or the helper dies before
sending), the 24 h timeout is never evaluated.

**Suggested fix:** Pass a deadline into `is_calc_finished` and return a third state, or
replace `sleep(1)` with `socket.poll(timeout=1000)` and check the clock in that loop.

---

## M3. Up to 1 second of idle latency per Gaussian external call

**Where:** `mace_gaussian/gaussian/zmq_server.py:106-111` (poll 10 ms, then `sleep(1)`).

**What:** Gaussian launches the helper, the helper sends its message, and the server notices
it only on its next loop iteration, on average 0.5 s later. VPT2 makes 6N-11 Hessian requests
(water 7, methane 19, decane 181) plus a few gradient calls.

**Why it matters:** That idle time is inside `runtime_s` and inside Gaussian's own
"Elapsed time", i.e. inside the cost-chapter numbers. Water ML: Gaussian elapsed 16.6 s for
1.6 s CPU and 7 calls. For decane it is roughly 1.5 minutes of pure sleeping in a run that
is supposed to demonstrate speed. `docs/methods.md` §4 even says IPC was chosen "to minimize
per-call latency".

**Suggested fix:** `if socket.poll(timeout=1000) != 0: return False` and drop the sleep.
One line. Then re-measure the alkane ladder before the cost chapter is written.

---

## M4. Geometry optimization: 1e-6 eV/Å target, 10000-step cap, convergence never checked

**Where:** `mace_gaussian/workflow.py:310-324` (`fmax=0.000001`, `steps=10000`),
`workflow.py:426` (`converged=True` hardcoded).

**What:** `docs/methods.md` §2 says 0.01 eV/Å. The code asks for 1e-6 eV/Å, which is at the
float64 noise floor of an ML potential, and records `converged=True` regardless of whether
LBFGS actually met it. Water met it in 8 steps. A floppy molecule may run 10000 steps, stop,
and still be recorded as converged.

**Suggested fix:** Use `opt.converged()` for the flag; pick and freeze one fmax (this is P0
item 3 in `thesis/TODO.md`); store the final max force in `results.json`.

---

## M5. Batch report metrics are computed differently from per-molecule metrics

**Where:** `mace_gaussian/analysis/batch_report.py:154-172`.

**What:** The leaderboard pairs sorted harmonic frequencies (no eigenvector matching,
no imaginary-frequency exclusion) and computes R² as 1 - SSres/SStot, while the per-molecule
reports use Hungarian matching, exclude imaginary pairs, and report Pearson r². For a molecule
with a mode swap or an imaginary mode the two reports disagree by construction.

**Suggested fix:** Have the batch report read `analysis_results_harmonic/<mol>/data/
comparison_*.csv` (already mode-matched) instead of recomputing from `results.json`.

---

## M6. Charge and multiplicity are hardcoded to 0 / 1

**Where:** `mace_gaussian/workflow.py:788-806`, `batch.py:198-199`,
`calculators/espaloma.py:55` (`DetermineBonds(mol, charge=0)`).

**What:** The pipeline signature accepts charge/spin, but the CLI never sets them and the
loader overwrites them with 0/1. Fine for the thesis panel (all neutral singlets), but the
docs should say "neutral closed-shell only" rather than imply support.

---

## L1. Espaloma "dipole derivatives" are a fixed-charge model (physics note, not a bug)

**Where:** `mace_gaussian/calculators/espaloma.py:34-76`.

Espaloma predicts charges from the molecular graph only. Displacing an atom does not change
the charges, so the finite-difference derivative in `base.py:41-79` is exactly
d(mu)/dr_i = q_i (the geometric term only). No charge flux. This should be stated in the thesis
when comparing espaloma and MACE4IR intensities, because it explains a systematic gap
independent of model quality.

## L2. `DFT_BASELINES` has four identical entries

`mace_gaussian/dft_baseline.py:52-77`. All four map to B3LYP/6-31G(d,p) and the same
directory, so runs 2-4 are skipped by `check_baseline_exists`. With `skip_if_exists=False`
the same DFT job would run four times. One entry is enough.

## L3. Dipole package path fallback points to a directory that does not exist

`mace_gaussian/calculators/mace_loader.py:32` builds `mace_gaussian/mace_dipole_pkg`;
the real directory is `<repo>/mace_dipole_pkg`. It works only because the package is
pip-installed in `mace4ir_v2`. Either fix the path (`parent.parent.parent`) or delete the
fallback.

## L4. Import-time side effects in the calculator factory

`mace_gaussian/calculators/factory.py:64` instantiates every calculator at import;
`espaloma.py:25-26` runs a live espaloma inference on NH3 as its availability test and
catches only `ImportError`. Any other exception (torch/DGL version issue, CUDA) breaks
`import mace_gaussian.calculators`, and therefore `workflow`. Same for `mace_ml.py:45`.
Consider catching `Exception` in the availability checks.

## L5. Fallback final energy uses the last displaced geometry

`mace_gaussian/workflow.py:581, 586, 589`: `mol` has been moved by every Gaussian request,
so `mol.get_potential_energy()` in the fallback is the energy at the last VPT2 displacement,
not the equilibrium. Moot once H4 is fixed, since the fallback then never triggers.

## L6. `extract_spectrum_from_fchk` picks an arbitrary log

`mace_gaussian/analysis/analyze_spectra.py:150-152` takes `glob("*.log")[0]`. Directories
have one log today; if a second one ever lands there (e.g. a SLURM-retrieved
`<mol>_freq_anharm.log` next to `gaussian_dft.log`), intensities may come from the wrong file.

## L7. Subprocess stdout pipe is never drained

`mace_gaussian/gaussian/runner.py:72-78` and `dft_baseline.py:243-249`: `stdout=PIPE` with
no reader. g16 writes to the `.log` file so this is fine in practice, but if it ever emits
more than the pipe buffer (64 KB) to stdout the process deadlocks silently.

## L8. Docs out of sync with the code (`docs/methods.md`)

| methods.md says | code does |
|---|---|
| Gaussian broadening, FWHM 8.0 | Lorentzian, FWHM 10.0 (`analyze_spectra.py:331-353`, default at `analysis_workflow.py:122`) |
| fmax 0.01 eV/Å | 1e-6 eV/Å (`workflow.py:310`) |
| mass-weighted eigenvector overlap | Cartesian (H3) |
| 2 energy models × 2 dipole models | 5 × 3 (`workflow.py:747-749`) |
| "same optimized geometry" as a feature | see H1 |

Chapter 3 of the thesis should be written from the code, not from `methods.md`.
(`docs/explained/02_software.md` from this session is written from the code.)

---

## Things I checked that are fine

For the defense, these are the "did you verify X" questions with a yes:

- **Unit conversions at the Gaussian boundary** (`io.py:94-97`, `workflow.py:135, 146, 210`,
  `base.py:79`, all three dipole calculators): energy eV→Ha, gradient eV/Å→Ha/Bohr,
  Hessian eV/Å²→Ha/Bohr², dipole e·Å→e·Bohr, dipole derivative e·Å/Å = e, polarizability
  Å³→Bohr³. All correct, CODATA 2018 constants.
- **Gaussian external file format** (`io.py:103-140`): energy+dipole line, 3N gradient
  lines, 2 polarizability lines, 3N dipole-derivative lines, lower-triangle Hessian three per
  line, Fortran `D` exponents. Matches the Gaussian External keyword spec.
- **Sign of the gradient**: `gradient = -forces` (`workflow.py:81`). Correct.
- **Hessian reshaping**: MACE returns `(3N, N, 3)`, reshaped to `(3N, 3N)` (`workflow.py:136`).
  Correct. Finite-difference fallback is symmetrized.
- **Dipole derivative layout**: `(3N, 3)` with row = (atom, cartesian), column = dipole
  component; autograd path transposes `(3, N, 3)` → `(N, 3, 3)` → `(3N, 3)`
  (`mace_loader.py:258`, `mace_polar1.py:128`). Correct.
- **Helper script argv**: Gaussian calls `script layer infile outfile msgfile fchk matel`;
  helper reads `argv[2]`, `argv[3]`. Correct.
- **VPT2 call count**: water log has exactly 7 external Hessian requests = 6·3 − 11. Correct.
- **Archive unwrapping for the dipole** (`parser.py:364-378`): correct, this was the
  July fix.
- **Intensity table selection** (`parser.py:113-149`): takes I(anharm) in km/mol, ignores
  DS(anharm) in 10⁻⁴⁰ esu²·cm². Verified on water B3LYP: 0.63 / 69.3 / 16.1 km/mol.
- **ZMQ LINGER=0 and IPC cleanup** (`zmq_server.py:53-87`): sound.
- **Manifest atomic writes** (`batch.py:49-70`): sound.
- **Transmittance→absorbance** (`nist_fetcher.py:451-458`): 2 − log10(%T). Correct.

## Test suite (run 2026-09-18, `mace4ir_v2`, 20 min 22 s)

`11 failed, 333 passed, 2 skipped`. The 11 failures are the known ones and are all test
problems, not code problems:

| tests | cause | fix |
|---|---|---|
| 8 in `tests/test_html_report.py` | `test_html_report.py:35` reads `report.html` with `read_text()` and no encoding; under pytest the default codec resolved to ASCII and the report contains `cm⁻¹` (byte 0xC2) | `read_text(encoding="utf-8")` in the test helper |
| `test_slurm.py::test_default_poll_interval_is_one_hour` | asserts 3600; code was changed to 600 (`slurm.py:33`) | update the assertion |
| `test_slurm.py::test_default_remote_base` | asserts `~/mace_gaussian_dft`; code now uses `/scratch_rune03a/mot/calculations/mace_gaussian` (`slurm.py:32`) | update the assertion |
| `test_slurm.py::test_submit_dft_jobs` | same remote-base change (`mace_gaussian_dft/water` expected in the mkdir command) | update the assertion |

**L10.** Two tests in `tests/test_workflow_calculator.py` (`TestElementGuardAtCallSites`,
the two `run_frequency_calculation` cases) run the real pipeline: they launched Gaussian on
water and rewrote `comparison_results/water/mace_anicc_*/results.json` and
`mace_off_*/results.json` with new timings during the 2026-09-18 test run. Tests must never
write into `comparison_results/`; they need a mocked `run_gaussian_with_zmq` and a
`tmp_path` output dir. (Restored with `git checkout`.)

One more thing the run showed (**L9**): seven tests in `tests/test_cli_validation.py` take
46 to 361 s each (17 of the 20 minutes) because click's `run` command validates its options
and then actually starts the pipeline, which loads MACE models and runs espaloma before the
mocked step is reached. They should stop at validation (mock `run_pipeline` before invoking)
or the suite is unusable as a pre-commit check.
