# How the software makes Gaussian believe MACE is a quantum chemistry program

This file explains the code the way you would explain it at a whiteboard: what runs, in
what order, what files exist at each moment, and where the units change. Line numbers refer
to branch `spike/24-vpt2-psience` at commit f929741. Section 6 traces `mace-gaussian run
water.xyz` from the command line to the report, with the real files and numbers.

`01_physics.md` explains why any of this is computed. `REVIEW_FINDINGS.md` lists the
places where the code does not do what this file says it should.

---

## 1. The one idea

Gaussian has a keyword, `External="<program>"`, that says: "do not run your own electronic
structure code; instead, every time you need energy and derivatives, write the geometry to a
file, run this program, and read the answer from another file." The program is
`mace_gaussian/gm_helper.py`. It is 40 lines and does no science. It forwards the two file
names to a long-running Python process that has the MACE models loaded on the GPU, waits
for "done", and exits. Gaussian then reads the output file and continues.

Two processes, one socket, a dozen file formats and unit conversions. That is the whole
bridge. Everything else in the package is orchestration around it (run the optimization,
run the DFT twin, run 15 combinations, parse, analyze, report).

Why not have `gm_helper.py` load MACE itself? Because Gaussian launches it fresh for every
one of the 6N−11 Hessian calls, and loading a foundation model takes 5 to 20 seconds.
Keeping one process alive with the model loaded and talking to it over a socket costs
milliseconds per call.

---

## 2. Package map

```
mace_gaussian/
  cli.py              click commands: run, batch, fetch, list, report, diagnose
  workflow.py         the pipeline: optimize -> DFT twin -> ML runs; the per-call callback
  batch.py            many molecules, manifest for restart
  slurm.py            ship DFT jobs to rune03 over ssh, poll sacct, scp back
  dft_baseline.py     write and run the pure-Gaussian B3LYP twin
  gm_helper.py        the 40-line script Gaussian launches
  pubchem.py          fetch a 3D structure by name
  gaussian/
    io.py             parse Gaussian's request file; write the answer file; write .gjf
    zmq_server.py     the socket server context manager + "is Gaussian done?" poll
    runner.py         launch g16, service requests until it exits, timeout, error capture
    parser.py         read frequencies, intensities, overtones, combos, dipole, timing from .log
    fchk.py           formchk wrapper; read normal modes and frequencies from .fchk
  calculators/
    base.py           dipole calculator interface + finite-difference derivatives
    mace_ml.py        MACE4IR dipole model wrapper
    mace_loader.py    loads the MACE4IR model with a pickle class-remap (see §4)
    mace_polar1.py    MACE-POLAR-1 dipole (autograd Jacobian)
    espaloma.py       graph charges -> dipole
    xtb.py            GFN2-xTB dipole (not installed in mace4ir_v2)
    factory.py        picks a dipole calculator by name
  analysis/
    analysis_workflow.py   ComparisonWorkflow: find results, match, plot, report
    analyze_spectra.py     SpectrumAnalyzer: broadening, mode-ID matching, metrics, matplotlib
    mode_matching.py       eigenvector overlap, Hungarian assignment, degenerate groups
    executive_summary.py   composite ranking + one-line verdict
    report_data.py         report_data.json + CSV export (the thesis-figure interface)
    html_report_generator.py, plotly_builders.py, _shared_css.py   the HTML report
    batch_report.py        multi-molecule leaderboard
    coverage_analysis.py   error by frequency region
    nist_fetcher.py        gas-phase spectra from NIST WebBook
  utils/
    units.py          CODATA 2018 constants
    results.py        ResultsManager: directory layout + results.json
    scratch.py        per-run scratch directories
    validation.py     prerequisite checks, device detection, version metadata
    exceptions.py     typed errors
```

Two vendored packages live next to it and are not part of the review:
`mace_ML_pkg/` (standard mace-torch, provides the energy models) and `mace_dipole_pkg/`
(a fork with dipole/polarizability heads, provides MACE4IR).

---

## 3. The three stages of one molecule

`workflow.run_pipeline` (`workflow.py:697`) runs three stages. Each depends only on stage 1's
output, which is the "cluster seam": you can run stage 2 on rune03 and stage 3 on the GPU box.

**Stage 1, geometry optimization** (`workflow.py:378-434`, `310-324`). Read the XYZ into an
ASE `Atoms`, attach the `mace_omol` calculator, run LBFGS to fmax = 1e-6 eV/Å, write
`comparison_results/<mol>/geometry_opt/optimized.xyz` and a `results.json` with energies,
step count, runtime, and software versions. Water: 8 steps, 4 s.

**Stage 2, the DFT twin** (`dft_baseline.py:280-469`). Write a normal Gaussian input,
`# opt freq(anharm) b3lyp/6-31G(d,p)`, starting from the optimized geometry, and run g16
with no external interface. Gaussian re-optimizes at B3LYP (it must; VPT2 needs its own
minimum), then does the full VPT2. Output lands in `comparison_results/<mol>/b3lyp_6-31Gdp/`.
On the cluster path (`batch.py:220-263`, `slurm.py`) the same `.gjf` is scp'd to rune03,
submitted with `templates/slurm_dft.sh`, polled with `sacct` every 10 minutes, and the
`.log/.chk/.fchk` are copied back.

**Stage 3, the ML runs** (`workflow.py:437-628`, `657-689`). For every (energy model,
dipole model) pair, by default 5 × 3 = 15:

1. Load the energy model (`workflow.calculator`, `workflow.py:335-370`), attach to a copy of
   the optimized `Atoms`.
2. Create a scratch directory `.scratch/run_<energy>_<dipole>_<timestamp>_<hex>/`.
3. Write `gaussian_freq.gjf` there with route `# freq (anharm)` and
   `# external="/abs/path/to/gm_helper.py"` (`io.py:143-176`). Absolute path, because
   Gaussian's working directory is the scratch dir. Note: no `opt`. The geometry is used as-is
   (finding H1).
4. Start the ZMQ server, launch `g16 gaussian_freq.gjf`, and service requests until g16
   exits (`runner.py:26-121`).
5. Move `.gjf`, `.log`, `.chk` to `comparison_results/<mol>/<energy>_<dipole>/`, run
   `formchk` to get `.fchk`, parse the log, write `results.json`.

Everything is best-effort per combination: one failing pair does not stop the others, and
`batch.py` records per-pair status in `batch_manifest.json` so a crashed batch restarts
where it stopped.

---

## 4. The request/response loop in detail

This is the part to be able to draw.

```
 g16 process                          Python process (workflow.run_frequency_calculation)
 ───────────                          ────────────────────────────────────────────────────
 needs E, dE/dx, d²E/dx²              GaussianZMQServer bound to ipc://.../zmq.ipc (REP socket)
 writes Gau-NNNN.EIn                  socket.poll() ... waiting
 spawns: gm_helper.py R EIn EOut ...  
     │ REQ socket connect                                  
     │ send "EIn|EOut"  ───────────────────────►  recv  → _on_request(msg)
     │ recv (blocks)                                    parse_gaussian_input(EIn)
     │                                                  update Atoms positions
     │                                                  energy, -forces from MACE
     │                                                  Hessian from MACE (autograd)
     │                                                  dipole, dμ/dx from dipole model
     │                                                  write_gaussian_output(EOut)
     │  ◄───────────────────────────────────────  send "ready"
     │ exit 0
 reads Gau-NNNN.EOut, continues
 ... repeat 6N−11 times (+ a few gradient-only calls) ...
 exits → runner sees proc.poll() != None → loop ends
```

**The request file** (`Gau-*.EIn`, parsed by `io.parse_gaussian_input`, `io.py:20-63`):

```
 natoms  deriv  charge  spin          e.g.  3  2  0  1
 Z  x  y  z  (in Bohr)  per atom
```

`deriv` is 0, 1 or 2 for energy, gradient, Hessian. Coordinates are converted Bohr → Å on
the way in.

**The response file** (`Gau-*.EOut`, written by `io.write_gaussian_output`, `io.py:66-140`),
all atomic units, Fortran `D` exponents, fixed 20.12 format:

```
 E  μx  μy  μz                       one line          (Hartree, e·Bohr)
 dE/dx dE/dy dE/dz                    N lines           (Hartree/Bohr)
 αxx αxy αyy / αxz αyz αzz            2 lines           (Bohr³, zeros unless POLAR)
 ∂μ/∂x_i  (3 components)              3N lines          (e)
 Hessian lower triangle               3 values per line (Hartree/Bohr²)
```

**The socket** (`zmq_server.py`, `gm_helper.py`): ZeroMQ REQ/REP over a Unix domain socket
file inside the scratch directory. REQ/REP is a strict ping-pong, which matches the protocol
exactly: one request, one reply, never two in flight. `LINGER=0` so the server does not hang
if Gaussian dies with a reply pending. The socket file path goes to the helper through the
`MACE_IPC_PATH` environment variable, set on the g16 subprocess and inherited by the helper.

**The poll** (`zmq_server.is_calc_finished`, `zmq_server.py:90-111`): poll the socket for
10 ms; if a message is there, return "not finished, go serve it"; else if g16 has exited,
return "finished"; else sleep one second and repeat. That one-second sleep is finding M3.

**The callback** (`workflow.run_next_calculation`, `workflow.py:231-307`) is the science
per call:

| step | code | what |
|---|---|---|
| 1 | `parse_gaussian_input` | read N, deriv, coordinates |
| 2 | `update_molecule_geometry` | set positions on the ASE `Atoms` |
| 3 | `calculate_energy_and_forces` | `atoms.get_potential_energy()`, `-atoms.get_forces()` |
| 4 | `calculate_hessian` | `calc.get_hessian(atoms)` → `(3N, N, 3)` eV/Å², reshaped to `(3N, 3N)`; finite-difference fallback for models without autograd Hessian (POLAR-1) |
| 5 | `calculate_dipole_properties` | `dipole_calc.calculate_dipole(atoms)` → e·Bohr; `calculate_dipole_derivatives(atoms)` → `(3N, 3)` in e |
| 6 | `write_gaussian_output` | convert and write |

Timing per call is logged at DEBUG (`energy=`, `hessian=`, `dipole=`).

### Loading the MACE4IR dipole model

The MACE4IR model file was pickled with class paths like `mace.modules.models.AtomicDielectricMACE`,
but the class with the right `forward()` lives in the fork at
`mace_dipole_core.modules.models`. The standard package has a class with the same name and a
different forward. `mace_loader.py` solves this by giving `torch.load` a custom `Unpickler`
whose `find_class` redirects `mace.modules.models.*` to the fork during deserialization only
(`mace_loader.py:46-78`). No `sys.modules` hacking, no global state. This was the fix for a
long-standing "wrong model silently loaded" bug and is a good defense anecdote about why
you should never trust a model that loads without error.

Dipole derivatives for MACE4IR come from `get_dielectric_derivatives()` in the fork
(autograd Jacobian, one pass instead of 6N forward passes). POLAR-1 has no such API, so
`mace_polar1.py:106-128` runs three `torch.autograd.grad` calls (one per dipole component)
against the positions. Espaloma has no gradient at all; the base class does central finite
differences with ±0.005 Å (`base.py:41-79`), which for a topology-only charge model reduces to
∂μ/∂xᵢ = qᵢ.

---

## 5. Units at every boundary

The single most important table in the project. Internal convention is ASE's: eV, Å.
Gaussian's is atomic units. Conversions live in `utils/units.py` (CODATA 2018) and are
applied exactly once each, at the file boundary.

| quantity | ASE / MACE side | Gaussian side | conversion | where |
|---|---|---|---|---|
| coordinates in | — | Bohr | × 0.529177 → Å | `io.py:59` |
| energy | eV | Hartree | ÷ 27.211386 | `io.py:94` |
| gradient | eV/Å (as −forces) | Hartree/Bohr | × 0.529177 ÷ 27.211386 | `io.py:97` |
| Hessian | eV/Å² | Hartree/Bohr² | × 0.529177² ÷ 27.211386 | `workflow.py:135, 146` |
| dipole | e·Å (all three models) | e·Bohr | ÷ 0.529177 | `mace_loader.py:214`, `mace_polar1.py:85`, `espaloma.py:73` |
| dipole derivative | (e·Å)/Å = e | e | none needed; FD path: (e·Bohr)/Å × 0.529177 | `base.py:79`, autograd paths |
| polarizability | Å³ | Bohr³ | ÷ 0.529177³ | `workflow.py:44, 210` |
| coordinates out (fchk) | — | Bohr | × 0.529177 → Å | `fchk.py:192` |
| dipole in the log archive | — | a.u. | × 2.541746 → Debye | `parser.py:405` |
| frequencies | — | cm⁻¹ | none | everywhere |
| intensities | — | km/mol (I(anharm) table, not DS) | none | `parser.py:113-149` |

All of these were checked line by line in the review and are correct. The July 2026 bugs
(intensities read from the dipole-strength table instead of km/mol; archive line-wrapping
truncating the dipole) are fixed on this branch.

---

## 6. `mace-gaussian run water.xyz`, traced

What happens, with the actual files from the July 2026 run.

**0 s.** `cli.run` validates the XYZ (3 atoms, O H H), checks that `g16` and `formchk` are on
PATH, that the MACE4IR model file exists, that `gm_helper.py` exists, reports the GPU
(RTX 2070 Super), deletes scratch directories older than 24 h, and calls `run_pipeline`.

**Stage 1 (4 s).** `mace_omol` (extra_large, float64, CUDA) loads. LBFGS runs 8 steps to
fmax 1e-6 eV/Å. Output:

```
comparison_results/water/geometry_opt/
  initial.xyz
  optimized.xyz      O 0 0 0.1181 / H 0 ±0.7627 −0.4676 ; E = −2079.866 eV
  results.json       num_steps 8, runtime 4.3 s, torch 2.4.0, mace 0.3.15, GPU name, CPU name
```

**Stage 2 (16 s on the workstation, or a SLURM job).** `gaussian_dft.gjf` with
`# opt freq(anharm) b3lyp/6-31G(d,p)`, `%NProcShared=4`, `%mem=4GB`. Gaussian optimizes
(3 steps, max force 0.000036), computes the harmonic Hessian, then 7 displaced Hessians, then
VPT2. Output directory `comparison_results/water/b3lyp_6-31Gdp/` with `.gjf`, `.log`, `.chk`,
`.fchk`, `results.json`. The log's VPT2 section is what the parser reads:

```
 Fundamental Bands            E(harm)   E(anharm)   I(harm)   I(anharm)
   1(1)                       3799.220  3624.001    1.637     0.628
   2(1)                       1665.300  1615.153   70.330    69.295
   3(1)                       3912.417  3722.273   20.228    16.098
 Overtones      1(2) 7161.6 ...  2(2) 3194.7 ...  3(2) 7346.7
 Combination Bands  2(1) 1(1) 5228.1 ... 3(1) 1(1) 7179.8 ... 3(1) 2(1) 5319.4
```

Note Gaussian's mode numbering here (1 = 3799, 2 = 1665) is not ascending; the checkpoint
stores modes ascending (1665, 3799, 3912). This matters for finding H2.

**Stage 3, first pair `mace_omol` + `mace_ml` (18 s).**

1. Scratch dir `.scratch/run_mace_omol_mace_ml_20260706_HHMMSS_xxxx/`. Inside it,
   `gaussian_freq.gjf`:
   ```
   %chk=gaussian_freq.chk
   %mem=2GB
   %NProcShared=2
   # freq (anharm)
   # external="/home/mot/mace_gaussian/mace_gaussian/gm_helper.py"

   Gaussian input generated from ASE

   0 1
   O   0.00000000  0.00000000  0.11812540
   H   0.00000000  0.76267791 -0.46756270
   H   0.00000000 -0.76267791 -0.46756270
   ```
2. `GaussianZMQServer` binds `ipc://.../zmq.ipc`. `g16 gaussian_freq.gjf` starts with
   `cwd` = scratch dir and `MACE_IPC_PATH` in its environment.
3. Gaussian writes `Gau-12345.EIn` (`3 2 0 1` then three lines of Z x y z in Bohr) and runs
   `gm_helper.py R Gau-12345.EIn Gau-12345.EOut ...`. The log shows:
   ```
   External calculation of energy, first and second derivatives.
   Running external command "/home/mot/mace_gaussian/mace_gaussian/gm_helper.py R"
   ```
4. The helper sends `"Gau-12345.EIn|Gau-12345.EOut"`. The callback runs: MACE energy
   (−2079.87 eV → −76.43 Hartree), forces, a 9 × 9 Hessian, MACE4IR dipole
   (0, 0, −0.390 e·Å → −0.738 e·Bohr), 9 × 3 dipole derivatives via autograd. The answer file
   is 1 + 3 + 2 + 9 + 15 = 30 lines. The helper gets "ready" and exits.
5. Gaussian reads the file, diagonalizes the mass-weighted Hessian, prints harmonic
   frequencies (1621.9, 3818.1, 3918.4 cm⁻¹ with intensities 80.5, 6.5, 63.3 km/mol),
   and starts VPT2. Six more external calls follow, one per ±displacement along each of the
   3 modes. Seven in total. `grep -c NIter gaussian_freq.log` = 7.
6. Gaussian prints the VPT2 tables (3643.8, 1570.3, 3739.4 cm⁻¹ fundamentals; overtones;
   combinations; `PT2 model: Deperturbed VPT2 (DVPT2)`; `No Fermi resonance found`) and the
   archive block with `Dipole=0.,0.,-0.7379...` (a.u.), and exits 0.
7. `is_calc_finished` sees the exit. Files move to
   `comparison_results/water/mace_omol_mace_ml/`. `formchk` produces `gaussian_freq.fchk`.
   `parse_gaussian_log` builds the dictionary; `ResultsManager.save_frequency_results`
   writes `results.json`:
   ```json
   {
     "energy_calculator": "mace_omol", "dipole_calculator": "mace_ml",
     "calculator_type": "ml",
     "frequencies": {
       "harmonic":   [{"freq_cm": 1621.9, "ir_intensity": 80.5}, ...],
       "anharmonic": [{"mode": 1, "freq_cm": 3643.8, "ir_intensity": 2.1, "freq_harmonic": 3818.1}, ...],
       "overtones":  [{"mode": 1, "overtone_level": 2, "freq_anharmonic": 7202.1, ...}, ...],
       "combination_bands": [{"mode1": 2, "mode2": 1, "freq_anharmonic": 5200.5, ...}, ...]
     },
     "energy_eV": 0.100,            <- wrong, finding H4 (thermal correction, not energy)
     "dipole": {"z": -1.876, "magnitude": 1.876},   (Debye)
     "runtime_s": 18.2,
     "gaussian_timing": {"total_elapsed_s": 16.6, "total_cpu_s": 1.6},
     "version_info": {...}
   }
   ```
8. Scratch dir deleted. Next pair.

Fourteen more pairs follow (`mace_omol_espaloma`, `mace_omol_mace_polar1`, `mace_mp_*`,
...). Each energy model loads once per pair (there is no sharing across pairs; a deferred
optimization). `mace_polar` + `mace_polar1` is the slow one at 58 s because its Hessian is
finite-difference and its dipole Jacobian is three backward passes per call.

**Analysis** (`python run_analysis.py water`, or `run_analysis_harmonic.py`):

1. `ComparisonWorkflow.find_dft_baseline` picks `b3lyp_6-31Gdp` (calculator_type `dft`,
   prefers a directory with an `.fchk`). `find_ml_results` lists the 15 ML directories.
2. Per ML pair, `extract_mode_mapping` reads both `.fchk` files with `force_harmonic=True`,
   builds the 3 × 3 overlap matrix, runs the Hungarian assignment, detects degenerate groups
   (none for water), and returns `{ml_idx: dft_idx}` plus overlaps (all 1.000 for water).
3. Anharmonic mode: spectra come from `results.json` (fundamentals + overtones + combos with
   IDs `F1`, `O2_2`, `C1_2`); harmonic mode: from `.fchk` plus intensities from the log.
4. `calculate_metrics` pairs by mode ID after applying the mapping, drops imaginary pairs,
   computes MAE/RMSE/R²/slope, intensity metrics on DFT modes ≥ 0.1 km/mol.
5. Plots (matplotlib PNG for the spectrum and regression, mode-overlap heatmap), a
   comparison CSV per pair, `metrics_summary.json`, then the HTML report with plotly figures,
   an executive-summary ranking, and `report_data.json` + `summary_metrics.csv` for the
   thesis figures.

Water, anharmonic, MAE over 3 fundamentals + 3 overtones + 3 combinations
(`analysis_results/water/data/metrics_summary.json`):

| energy model | MAE (cm⁻¹) | RMSE | R² |
|---|---|---|---|
| mace_anicc | 16.8 | 19.5 | 0.9999 |
| mace_off | 33.2 | 38.3 | 0.9999 |
| mace_omol | 39.0 | 44.1 | 0.9997 |
| mace_polar | 73.5 | 89.3 | 0.9995 |
| mace_mp | 129.9 | 148.9 | 0.9971 |

The frequency MAE is identical across the three dipole models for a given energy model,
which is exactly what the design predicts: frequencies come from the energy surface only.
The intensity MAE is the reverse: 6.5 km/mol for MACE4IR, 9.3 for POLAR-1, 50.6 for espaloma,
regardless of energy model. That separability is the strongest single argument that the
bridge is doing what it claims.

---

## 7. Design decisions worth defending

**Why Gaussian at all, rather than an open VPT2 code?** Because Gaussian's VPT2 is the
reference implementation everyone in the field compares against, and using it unchanged
makes the ML-vs-DFT comparison a controlled experiment. The Psience spike (branch
`spike/24-vpt2-psience`, `.planning/phases/24-vpt2-research-spike/SPIKE.md`) showed an
independent VPT2 implementation reproduces Gaussian's numbers to < 0.01 cm⁻¹ from the same
force field, so the harness is not hiding a Gaussian artifact.

**Why ZMQ over a Unix socket, not a pipe or a file?** A pipe would tie the helper's lifetime
to the server's. A file would need polling and locking. ZMQ REQ/REP gives exactly the
one-request-one-reply semantics the protocol has, with a blocking `recv` on the helper side
that turns into a clean exit when the server replies. The IPC transport avoids TCP setup per
call.

**Why a separate dipole model?** Energy models have no dipole head. MACE4IR was trained
specifically on dipoles for IR; espaloma is a cheap baseline; POLAR-1 is the "one model for
both" future. Having three lets the thesis separate frequency error (energy model) from
intensity error (dipole model).

**Why is the geometry shared across energy models?** Design choice to isolate PES curvature
from geometry differences. The review argues it should be reconsidered (H1) because the
non-optimizer models are then evaluated off their minimum, which VPT2 does not tolerate.

**Why per-run scratch directories?** Gaussian writes `Gau-*.EIn/EOut/rwf/int/d2e` files in
its working directory; running 15 pairs sequentially in one directory would let stale files
from a crashed run be picked up. Each run gets its own directory, deleted on success, kept on
failure if `--keep-scratch`.

**Why a manifest?** Fifteen Gaussian runs per molecule times 22 molecules is 330 runs.
Something will crash. `batch_manifest.json` records per-pair status atomically so a restart
skips what is done.

**Why absolute paths in the `.gjf`?** Gaussian resolves `External=` relative to its own
working directory, which is the scratch dir. This bit us early; it is in `CLAUDE.md` as a
gotcha.

---

## 8. Things the code does not do (scope statements for the thesis)

- Neutral closed-shell molecules only (charge 0, multiplicity 1 hardcoded, finding M6).
- No sharing of loaded models across pairs; each pair reloads (cost, not correctness).
- No re-optimization per energy model (H1).
- No uncertainty estimate on ML predictions (single model, no committee).
- Analysis compares to one DFT twin (B3LYP/6-31G(d,p)) and to NIST gas-phase spectra; no
  other functional or basis.
- `compare` and `export` CLI commands are stubs; analysis runs through the two
  `run_analysis*.py` scripts.
