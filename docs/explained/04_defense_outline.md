# Defense talk outline

Target: 25 minutes, about 14 content slides at 1.5 to 2 minutes each, plus backups.
Follows the narrative arc in `thesis/STORY.md`. Results slides are written for the data
the campaign will produce; where a placeholder is marked `[DATA]`, the water numbers are
given as the current stand-in so you can rehearse with real values now.

Rules of thumb used here: one message per slide, the message is the slide title, every
number on a slide is one you can defend from `03_hard_questions.md`, and no slide has more
than one plot.

---

## Content slides (25 min)

### 1. Title (0:30)

**Title:** Machine-learned foundation potentials for anharmonic IR spectroscopy

**On slide:** title, name, supervisor, date, the one-sentence pitch in small type:
"We swapped the quantum engine inside Gaussian's anharmonic machinery for ML potentials and
measured when the cheap engine is good enough."

**Say:** the sentence. Nothing else.

---

### 2. The spectrometer sees anharmonic spectra (1:30)

**Message:** The harmonic approximation is what we can afford; anharmonicity is what we observe.

**Visual:** one molecule (water or formic acid), experimental gas-phase spectrum in black,
harmonic B3LYP sticks in grey. Stretches visibly 150 cm⁻¹ too high; overtone/combination
bands present in experiment, absent in the sticks.

**Say:** Bonds are not springs. Harmonic frequencies are 3 to 5 % too high and the
harmonic model has zero intensity for every overtone and combination band. The standard
fix is to multiply by 0.96 and hope.

**Numbers to have:** water stretches harmonic 3799/3912 vs experiment 3657/3756.

---

### 3. Why nobody does VPT2 (1:30)

**Message:** Anharmonic spectra need 6N−11 Hessians, and a DFT Hessian is expensive.

**Visual:** build-up: one Hessian at the minimum; then ± displacement along each mode;
the count 1 + 2(3N−6) = 6N−11. Small table: water 7, methane 19, naphthalene 97, decane 181.
Then a bar: DFT Hessian ~ N³ each, so total ~ N⁴. "Days per molecule above ~15 atoms."

**Say:** VPT2 corrects each level with cubic and quartic force constants obtained by finite
differences of Hessians. The physics is a 1972 paper; the cost is the wall.

---

### 4. The trick: swap the engine, keep the machinery (2:00)

**Message:** Gaussian keeps doing everything except the electronic structure.

**Visual:** architecture diagram (Figure F1). Gaussian box on the left doing: normal modes,
displacements, finite differences, VPT2 algebra, resonances, intensities. An arrow labeled
"geometry in / energy, gradient, Hessian, dipole, dipole derivatives out" to a Python box
with the MACE energy model and the dipole model. The `External` keyword and the ZMQ socket
labeled on the arrow. One thing crossed out: "DFT SCF".

**Say:** Gaussian's External interface lets any program answer "what is the energy here".
We answer with a MACE model, over a socket, in atomic units. Because everything else is the
same code, every difference between ML and DFT spectra is attributable to the potential.

**Backup pointer:** B-1 (protocol details), B-2 (unit table).

---

### 5. The contenders (1:30)

**Message:** Five foundation potentials with different provenance, one dedicated dipole model.

**Visual:** table: model, training data, level of theory, expected weakness.
MACE-MP (crystals, PBE), MACE-OFF (organics, ωB97M-D3), MACE-ANICC (CHNO, CC-quality),
MACE-OMOL (100M molecules, ωB97M-V), MACE-POLAR-1 (OMOL + multipoles). Dipoles: MACE4IR
(direct), espaloma (fixed charges, baseline), POLAR-1.

**Say:** The hypothesis is that spectroscopic behaviour follows training provenance.
Crystal-trained should fail on molecular bends. CHNO-only should be excellent inside its
domain and refuse outside. Intensities need a separate model because energy models have no
dipole.

---

### 6. Trust ladder, rung one: water (2:00)

**Message:** The numbers are real: ML-VPT2 matches DFT-VPT2 to tens of cm⁻¹ on water.

**Visual:** the water table (three modes × experiment / DFT harm / DFT VPT2 / OMOL harm /
OMOL VPT2). Highlight that the anharmonic correction is the same size for ML and DFT.

**Numbers (real):**

| | exp | DFT harm | DFT VPT2 | OMOL harm | OMOL VPT2 |
|---|---|---|---|---|---|
| bend | 1595 | 1665 | 1615 | 1622 | 1570 |
| sym str | 3657 | 3799 | 3624 | 3818 | 3644 |
| asym str | 3756 | 3912 | 3722 | 3918 | 3739 |

**Say:** OMOL-VPT2 lands within 20 cm⁻¹ of DFT-VPT2 on both stretches. The cubic force
constants of a neural network agree with DFT's. Plus the harness validation: an independent
VPT2 code reproduces Gaussian to 0.01 cm⁻¹ from the same force field.

---

### 7. Rung two: CO₂ and the Fermi resonance (1:30)

**Message:** VPT2's weak point, resonances, survives the ML surface `[DATA]`.

**Visual:** CO₂ level diagram ν₁ and 2ν₂, the deperturbed dyad from DFT and from the best
ML model, experiment marked.

**Say:** When two levels are near-degenerate the perturbation series diverges and Gaussian
switches to deperturbation. Whether the ML surface lands on the same side of the threshold
decides whether the spectra are even comparable. `[DATA: result]`.

---

### 8. The map: model × functional group (2:00)

**Message:** Where each model works `[DATA]`.

**Visual:** Figure F2, heatmap of MAE (cm⁻¹) with models as rows and functional groups
(C=O, O-H, N-H, C≡N, C-H, bends, ring) as columns. Harmonic and anharmonic side by side or
toggled.

**Say:** walk one row (OMOL: uniformly good) and one column (bends: MP fails).
`[DATA: the two or three cells that carry the story]`.

---

### 9. Bending versus stretching (1:30)

**Message:** Crystal-trained models miss molecular bends; it is a three-body problem.

**Visual:** scatter of ML−DFT error vs frequency for MP and for OMOL, bends coloured.

**Say:** Bending stiffness depends on three-body angular terms. A model trained on dense
periodic solids has seen few isolated-molecule angles. `[DATA]`. Caveat you must state:
this is measured at each model's own minimum (H1 was fixed before the campaign).

**Stand-in:** water bend, MP 1462 vs DFT 1665 harmonic.

---

### 10. Intensities: the neglected half (1:30)

**Message:** A dedicated dipole model reproduces intensities; fixed charges do not.

**Visual:** one molecule, three broadened spectra stacked: DFT, OMOL+MACE4IR, OMOL+espaloma.
Same frequencies, different peak heights.

**Numbers (real, water, km/mol, bend/sym/asym):** DFT 69/0.6/16; MACE4IR 94/2.1/35;
espaloma 200/95/191.

**Say:** Espaloma charges do not move with the geometry, so its derivative has no charge
flux; it makes the silent symmetric stretch bright. MACE4IR gets the ordering and the silent
mode right.

---

### 11. The central plot (2:30)

**Message:** Does cheap-potential-plus-physics beat better-potential-without-physics?
`[DATA]`

**Visual:** Figure F3. Per-molecule (or cumulative) error against gas-phase experiment for
three treatments: scaled-harmonic B3LYP (×0.96), DFT-VPT2, best-ML-VPT2.

**Say:** This is research question 3. `[DATA: the answer, in one sentence, then the
caveat]`. Either outcome is a result: if ML-VPT2 wins, anharmonic spectroscopy is routine;
if it does not, the failure analysis says which surface property is missing.

---

### 12. Cost: where the crossover is (1:30)

**Message:** ML wall time grows linearly, DFT-VPT2 as N⁴; the crossover is at `[DATA]` atoms.

**Visual:** Figure F4. Log-log wall time vs atoms for the alkane ladder, harmonic and VPT2,
ML and DFT; DFT-VPT2 extrapolated and shaded "impractical".

**Say:** For water there is no speedup, the pipeline overhead dominates. `[DATA: decane
DFT-VPT2 extrapolated to X days vs Y minutes ML]`. Same Gaussian machinery on both sides;
hardware stated on the slide (RTX 2070 Super vs rune03 CPUs).

---

### 13. The payoff: a molecule where DFT-VPT2 is impractical (1:30)

**Message:** Full anharmonic spectrum of naphthalene in minutes, overlaid on experiment
`[DATA]`.

**Visual:** Figure F5. Experimental gas-phase spectrum, ML-VPT2 broadened spectrum on top,
wall-clock time in the corner.

**Say:** 18 atoms, 97 Hessians, `[DATA: minutes]`. PAHs are what JWST sees in the interstellar
medium; anharmonic spectra of them are a live need.

---

### 14. Where it breaks, and what is next (1:30)

**Message:** Know the limits; the door is open.

**Visual:** two-column slide. Left, failure taxonomy table: element coverage (Cl), bends for
crystal-trained models, resonance threshold crossings, large-amplitude motion (rotors),
non-stationary geometries. Right, outlook: Gaussian-free VPT2 (Psience spike, 0.01 cm⁻¹),
POLAR-1 for energy and dipole in one model, committee uncertainty, experimental-feedback
training.

**Say:** A cheap surface does not fix a wrong method; it makes the right method affordable.
The next step is to remove Gaussian from the loop entirely, which the spike showed is
feasible.

---

### 15. Conclusions (1:00)

**Message:** Four research questions, four answers.

**On slide, four lines:**
1. Accuracy vs DFT: `[DATA]` cm⁻¹ MAE for the best models, provenance-dependent.
2. Intensities: MACE4IR reproduces DFT intensities; fixed-charge models do not.
3. Anharmonic-ML vs scaled-harmonic-DFT vs experiment: `[DATA]`.
4. Cost and failure: linear scaling, crossover at `[DATA]` atoms; fails for `[list]`.

**Say:** the four lines, then "thank you".

---

## Backup slides

Ordered by how likely you are to need them.

**B-1. The External protocol in one slide.** The request file (`N deriv charge spin`, then
Z x y z in Bohr) and the response file (energy+dipole line, N gradient lines, 2
polarizability lines, 3N dipole-derivative lines, Hessian lower triangle). ZMQ REQ/REP
sequence diagram. For "how does it actually talk to Gaussian".

**B-2. Units at every boundary.** The table from `02_software.md` §5. For "how do you know
there is no unit error".

**B-3. Validation of the harness.** Three items: line-by-line unit audit; Psience
cross-check (< 0.01 cm⁻¹); DFT self-consistency check (B3LYP through the external interface
vs native `freq(anharm)` at the same geometry: water max |Δν| 0.008 cm⁻¹, formaldehyde
0.034 cm⁻¹, intensities within 0.3 %, all 3 + 6 fundamentals, overtones and combination
bands; `docs/explained/self_consistency/`, run 2026-09-18 with
`scripts/dft_self_consistency.py`). Plus the July intensity bug as an honesty anecdote (factor 3.57, found,
fixed, old data flagged).

**B-4. Mode matching.** Overlap matrix heatmap for methane (T₂ triple), Hungarian assignment,
subspace overlap trace(MᵀM)/k for degenerate groups, mass-weighting. For "how do you pair
modes".

**B-5. VPT2 formulas.** Fundamental, overtone, combination in terms of χᵢⱼ; the resonance
denominators; DVPT2. For a theory-minded committee member.

**B-6. MACE architecture.** One slide: neighbour environment, spherical-harmonic features,
ACE many-body tensor product, two message-passing layers, readout, sum, autograd for forces
and Hessian. For "explain the model".

**B-7. Training data of the five models.** Sizes, levels of theory, element coverage, citations.
For "why should MP fail on bends".

**B-8. Water full data.** All 15 combinations, harmonic and VPT2, frequencies and
intensities, plus overtones and combinations. For any "show me the numbers" on the trust
rung.

**B-9. Espaloma as a fixed-charge model.** ∂μ/∂xᵢ = qᵢ derivation in two lines; what charge
flux is; why bends and C-H stretches suffer most.

**B-10. The panel and why each molecule is there.** The table from `thesis/MOLECULES.md`.
For "is 22 molecules enough / why these".

**B-11. Protocol settings.** Optimizer fmax, Gaussian displacement step, resonance thresholds,
hardware, software versions, per-model re-optimization. For reproducibility questions.
(Fill after the P0 "freeze protocol" item.)

**B-12. Timing breakdown per call.** Energy / Hessian / dipole / overhead for one molecule
per model; where the fixed costs are. For "why is water not faster".

**B-13. Known limitations of the code.** Neutral closed-shell only; single models, no
uncertainty; one reference functional; VPT2's small-amplitude assumption.

**B-14. Reference-level caveat.** Training levels of theory vs B3LYP/6-31G(d,p); why
experiment is the arbiter.

---

## Timing check

| slides | minutes |
|---|---|
| 1-3 problem | 3.5 |
| 4-5 idea and contenders | 3.5 |
| 6-7 trust | 3.5 |
| 8-10 evidence | 5.0 |
| 11 central plot | 2.5 |
| 12-13 cost and payoff | 3.0 |
| 14-15 limits, outlook, conclusions | 2.5 |
| **total** | **23.5** |

Leaves 90 seconds of slack. If you run long, slide 9 (bending vs stretching) folds into
slide 8.

## Rehearsal notes

- Slides 2, 4, 6, 11 are the ones the committee will remember. Rehearse those with a timer.
- Every `[DATA]` gets filled from `report_data.json` / `summary_metrics.csv`, never typed by
  hand.
- The honest-limitation slide (14) is where credibility is won. Do not rush it.
- If asked anything from `03_hard_questions.md` section E, answer from that file verbatim
  the first time; it is written to be spoken.
