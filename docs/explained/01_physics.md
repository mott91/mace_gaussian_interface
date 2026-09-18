# What is being computed, and why it means something

This is the physics of MACE-Gaussian in the order a listener needs it: what an IR spectrum
is, where the numbers come from, where the machine learning goes in, and what can go wrong.
Every section ends with a one-sentence version you can say out loud.

The companion file `02_software.md` explains how the code does it. `03_hard_questions.md`
has the questions a committee will ask.

---

## 1. An IR spectrum in one paragraph

A molecule is a set of nuclei held together by electrons. The nuclei are never still; they
vibrate around their equilibrium positions. Each vibration has a natural frequency. If you
shine infrared light on the molecule, light whose frequency matches a vibration can be
absorbed, but only if that vibration changes the molecule's dipole moment. The spectrum is a
plot of absorption against frequency. Peak positions tell you the vibration frequencies.
Peak heights tell you how strongly each vibration moves charge around.

Frequencies are quoted in wavenumbers, cm⁻¹. The mid-IR window is roughly 400 to 4000 cm⁻¹.
Bends sit low (1000 to 1700), stretches of bonds to hydrogen sit high (2800 to 3800).

> **Say it:** "An IR spectrum is the list of frequencies a molecule vibrates at, weighted by
> how much each vibration shakes the charge distribution."

---

## 2. Where frequencies come from: the potential energy surface

Everything starts from the Born-Oppenheimer picture. Electrons are light and fast, nuclei are
heavy and slow, so for any fixed set of nuclear positions **R** the electrons settle into a
ground state with energy E(**R**). That function is the potential energy surface (PES).
Vibrations are the nuclei rolling around inside the bowl of E(**R**) near its minimum.

For a molecule with N atoms, **R** has 3N numbers. Near the minimum **R₀** we can Taylor-expand:

    E(R) = E₀ + 0 · (R−R₀)  +  ½ (R−R₀)ᵀ H (R−R₀)  +  cubic terms  +  quartic terms + ...

The linear term is zero because we are at a minimum (the gradient, i.e. the force, vanishes).
**H** is the Hessian, the 3N × 3N matrix of second derivatives ∂²E/∂xᵢ∂xⱼ.

If you keep only the quadratic term you have the **harmonic approximation**. The molecule is
then a set of coupled springs. To uncouple them:

1. Mass-weight: H̃ᵢⱼ = Hᵢⱼ / √(mᵢ mⱼ).
2. Diagonalize H̃. The eigenvalues λₖ give the frequencies, ωₖ = √λₖ / (2πc).
3. The eigenvectors are the **normal modes**: the pattern of atomic displacements for each
   vibration.

3N eigenvalues come out. Six are (near) zero: three translations and three rotations of the
whole molecule. The remaining 3N−6 (3N−5 for linear molecules) are the vibrations. Water has
3·3−6 = 3 modes. Methane has 9. Decane has 90.

Units matter here and they are the single most common source of bugs in this kind of code.
Gaussian works in atomic units: energy in Hartree, length in Bohr, so the Hessian is in
Hartree/Bohr². ASE and MACE work in eV and Ångström. The code converts at the boundary
(see `02_software.md` §5).

> **Say it:** "Harmonic frequencies are the square roots of the eigenvalues of the
> mass-weighted second-derivative matrix of the energy at the minimum."

### What the Hessian costs

For DFT, one Hessian at the B3LYP/6-31G(d,p) level on a small organic molecule is minutes to
hours, because it needs the second derivative of the electronic energy, which means solving
coupled-perturbed Kohn-Sham equations. For a MACE model, energy is a neural-network forward
pass, forces are one backward pass, and the Hessian is another differentiation through the
same graph. Milliseconds to seconds on a GPU. This gap of three to four orders of magnitude is
the entire reason the project exists.

---

## 3. Why harmonic is not enough

Real bonds are not springs. Stretch a bond far enough and it breaks; compress it and the
repulsion rises steeply. The potential is asymmetric (Morse-like), and the quadratic
approximation is too stiff on the stretching side. Two consequences:

1. **Harmonic frequencies are systematically too high**, by 3 to 5 % for most modes, more for
   X-H stretches. Water: harmonic B3LYP gives 3799 and 3912 cm⁻¹ for the two O-H stretches,
   experiment says 3657 and 3756. That is 140 to 160 cm⁻¹ off, which is huge for assignment.
   The usual fix is to multiply everything by an empirical scaling factor (about 0.96 for
   B3LYP). It works on average and hides the physics.

2. **A harmonic oscillator only absorbs at its fundamental.** Real spectra also show
   overtones (2νᵢ, roughly twice a fundamental, slightly less) and combination bands
   (νᵢ + νⱼ). These are weak but they are everywhere in a real spectrum, especially above
   3000 cm⁻¹ and in the fingerprint region. Harmonic theory says they have zero intensity.

> **Say it:** "The harmonic approximation is what we can afford; anharmonicity is what
> the spectrometer sees."

---

## 4. VPT2: fixing it without solving the full problem

Second-order vibrational perturbation theory keeps the cubic and quartic terms of the Taylor
expansion and treats them as a perturbation on the harmonic solution. Written in normal
coordinates Q:

    V = ½ Σ ωᵢ² Qᵢ²  +  ⅙ Σ φᵢⱼₖ QᵢQⱼQₖ  +  1/24 Σ φᵢⱼₖₗ QᵢQⱼQₖQₗ

φᵢⱼₖ are the cubic force constants (third derivatives), φᵢⱼₖₗ the quartic ones. Second-order
perturbation theory then gives closed formulas for the energy levels in terms of
**anharmonic constants** χᵢⱼ, which are algebraic combinations of the φ's and the ω's.
The results you actually use:

- fundamental:  νᵢ = ωᵢ + 2χᵢᵢ + ½ Σⱼ≠ᵢ χᵢⱼ
- overtone:     2νᵢ = 2ωᵢ + 6χᵢᵢ + Σⱼ≠ᵢ χᵢⱼ
- combination:  νᵢ + νⱼ = ωᵢ + ωⱼ + 2χᵢᵢ + 2χⱼⱼ + 2χᵢⱼ + ½ Σₖ (χᵢₖ + χⱼₖ)

The χ's are negative for stretches, which is why anharmonic fundamentals sit below harmonic
ones and overtones sit below twice the fundamental. Water again: VPT2 on B3LYP moves the
stretches from 3799/3912 down to 3624/3722, within 30 to 70 cm⁻¹ of experiment. The bend
goes 1665 → 1615 (experiment 1595). The overtone of the bend appears at 3195 cm⁻¹ with a
small but nonzero intensity. That is the physics the harmonic calculation cannot give you.

### Where the cubic and quartic constants come from

Nobody computes third and fourth analytic derivatives of the energy. Instead Gaussian
displaces the molecule along each normal mode by ±δ and recomputes the full Hessian at every
displaced geometry. Finite differences of Hessians give the cubic constants and the
semi-diagonal quartic constants (the ones VPT2 needs). Count the Hessians:

    1 at equilibrium + 2 × (3N−6) displaced = 6N − 11

Water: 7 Hessians. Methane: 19. Decane (32 atoms): 181. Naphthalene (18 atoms): 97.
With DFT, each Hessian is an expensive job, so VPT2 on anything beyond ten atoms is days.
With MACE, each is a second. That count is the number of times Gaussian calls the ML model
in this project, and you can see it in the log (`02_software.md` §6).

### Resonances, the thing that makes VPT2 dangerous

The χᵢⱼ formulas have denominators like (2ωⱼ − ωᵢ) and (ωᵢ − ωⱼ − ωₖ). When two levels
happen to be nearly degenerate, for instance an overtone 2ν₂ sitting on top of a fundamental
ν₁, the denominator goes to zero and the perturbation series blows up. This is a **Fermi
resonance**. CO₂ is the textbook case (ν₁ ≈ 2ν₂). Darling-Dennison resonances are the
analogous 2-2 case between overtones.

Gaussian's implementation (which the thesis relies on) detects these by threshold, removes
the offending terms from the perturbation sum, and diagonalizes a small matrix for the
resonant levels instead. The log says `PT2 model: Deperturbed VPT2 (DVPT2)`. For water there
are no resonances. For CO₂ and formic acid there are, and whether the ML surface produces a
resonance pattern like DFT's is one of the benchmark's stress tests.

> **Say it:** "VPT2 adds the cubic and quartic curvature of the surface by finite differences
> of Hessians, 6N−11 of them, and corrects each level analytically; resonances are handled by
> deperturbation."

---

## 5. Intensities: the other half of the spectrum

A vibration absorbs IR light in proportion to how much the dipole moment μ changes along the
mode:

    Iᵢ ∝ |∂μ/∂Qᵢ|²

∂μ/∂Q is a vector (three components of the dipole) for each mode. In Cartesian coordinates it
is a 3N × 3 tensor ∂μ/∂x, the **atomic polar tensor**, which Gaussian projects onto the
normal modes. Units end up as km/mol. Water's bend has about 70 km/mol, the symmetric stretch
only 1.6, the asymmetric stretch 20.

This is the double-harmonic approximation (harmonic potential, linear dipole). Anharmonic
intensities also need the second derivatives of the dipole, ∂²μ/∂Q², which Gaussian gets by
finite differences of ∂μ/∂x across the same displaced geometries. This is why the external
program has to hand back dipole derivatives at every one of the 6N−11 calls, not just at
equilibrium.

The important point for this thesis: **the energy model knows nothing about the dipole.**
A MACE energy model is trained on energies and forces. It predicts E(**R**) and its derivatives,
full stop. To get intensities you need a second model that predicts μ(**R**). Three are used:

| dipole model | how it gets μ | what its derivative captures |
|---|---|---|
| MACE4IR (`mace_ml`) | direct equivariant prediction of the dipole vector, trained for IR | full: charge motion and charge flux; derivative by autograd |
| MACE-POLAR-1 (`mace_polar1`) | dipole as a by-product of a multipole charge-density expansion | full; derivative by autograd (with finite-difference fallback) |
| Espaloma (`espaloma`) | partial charges from the molecular graph, μ = Σ qᵢ rᵢ | geometric term only: the charges do not change when atoms move, so ∂μ/∂xᵢ = qᵢ exactly |

That last row matters. Espaloma's intensities are a fixed-point-charge model. It cannot
represent charge flux (the redistribution of electrons during a vibration), which is often
the dominant contribution for bends and for C-H stretches. Any systematic espaloma
intensity error is a model-class limitation, not a training-data problem. Water shows it:
the symmetric stretch gets 95 km/mol from espaloma against 0.6 from DFT and 2.1 from MACE4IR.

> **Say it:** "Peak heights need the dipole derivative, which the energy model cannot give,
> so a second ML model predicts the dipole surface and the code differentiates it."

---

## 6. Where MACE goes in, and what "engine swap" means

Gaussian's VPT2 code is a machine that asks a question and consumes an answer:

- *Question:* "Here are 3N coordinates. Give me energy, gradient, Hessian, dipole, and
  dipole derivatives."
- *Answer:* those numbers, in atomic units, in a file.

Normally the answer comes from Gaussian's own DFT code. Gaussian's `External` keyword lets
any program answer instead. MACE-Gaussian is that program. Gaussian still does the geometry
bookkeeping, the normal-mode analysis, the displacements, the finite differences, the
resonance detection, the intensity assembly, and the printing. Only the electronic structure
is replaced.

Why this is a good experimental design: **every difference between the ML spectrum and the
DFT spectrum is caused by the ML potential (or dipole surface), because everything else is
literally the same code.** Most ML-potential benchmarks compare against a separately
implemented harmonic pipeline and cannot make that claim.

What a MACE model actually is, in the amount of detail a defense needs:

- An equivariant graph neural network. Atoms are nodes, neighbours within a cutoff are edges.
- Each atom builds a description of its environment from the relative positions of its
  neighbours, using spherical harmonics so that the description rotates correctly when the
  molecule rotates (equivariance). MACE's specific trick is to build many-body (not just
  pairwise) features in one step via the Atomic Cluster Expansion.
- Two message-passing layers refine those features. A small readout maps each atom's
  features to an atomic energy. The total energy is the sum. Forces and the Hessian are
  exact derivatives of that sum via automatic differentiation.
- The five models differ only in weights and training data:

| model | trained on | level of theory | what to expect |
|---|---|---|---|
| MACE-MP-0 | Materials Project crystals | PBE | poor for molecular bends, biased frequencies |
| MACE-OFF | organic molecules (SPICE) | ωB97M-D3(BJ)/def2-TZVPPD | good organic frequencies |
| MACE-ANICC | ANI-1ccx (CHNO) | CCSD(T)*/CBS-quality | very good for CHNO, no other elements |
| MACE-OMOL-0 | OMol25 (100M+ molecules) | ωB97M-V/def2-TZVPD | broad coverage, the default optimizer |
| MACE-POLAR-1 | OMol25 + polarization targets | ωB97M-V | energy plus multipoles |

Note the reference-level problem: none of these is trained on B3LYP/6-31G(d,p). Part of any
ML-vs-DFT disagreement is functional-vs-functional disagreement. That is why comparison to
gas-phase experiment is part of the benchmark and not decoration.

> **Say it:** "Gaussian keeps doing all the vibrational machinery; we swap only the thing
> that answers 'what is the energy here', so every error is attributable to the potential."

---

## 7. Comparing two spectra honestly

### Mode matching

The naive comparison sorts both frequency lists and pairs them by rank. That fails whenever
two modes are close and swap order between methods, which happens constantly. The right way
is to compare the normal-mode **eigenvectors**: mode i from DFT and mode j from ML describe
the same physical motion if their displacement patterns overlap, |⟨qᵢ|qⱼ⟩| ≈ 1. Build the
full overlap matrix, then solve the assignment problem (Hungarian algorithm) for the
one-to-one pairing with maximum total overlap. Pairs with overlap below 0.5 are flagged as
uncertain.

The overlap has to be taken in **mass-weighted** coordinates. Normal modes are only orthonormal
there. In plain Cartesian displacements, which is what Gaussian writes to the checkpoint file,
modes of the same molecule overlap each other by up to 0.13 for hydrogen-rich molecules.
(The code currently does not mass-weight; finding H3 in `REVIEW_FINDINGS.md`.)

### Degenerate modes

In symmetric molecules several modes share one frequency (methane's three C-H stretches at
3162 cm⁻¹, the T₂ set). Any rotation inside that three-dimensional subspace is an equally
valid set of eigenvectors, so one-to-one overlaps can be arbitrarily low even when the
physics is identical. The fix is to compare subspaces: group modes within 0.5 cm⁻¹ of each
other, take the block M of the overlap matrix between the two groups, and use
trace(MᵀM)/k as the group overlap. A perfectly reproduced subspace scores 1 regardless of
how the vectors are rotated inside it.

### Metrics

On matched pairs: mean absolute error and root-mean-square error in cm⁻¹, R² and slope of
ML against DFT, the same for intensities with modes below 0.1 km/mol excluded (their relative
error is noise). Imaginary frequencies (printed as negative numbers by Gaussian) are excluded
from statistics and counted separately, because an imaginary frequency means the geometry
is not a minimum, not that the frequency is negative.

### Broadening and experiment

A computed spectrum is a set of sticks. A measured one has line shapes. To overlay them the
sticks are convolved with a Lorentzian of 10 cm⁻¹ full width at half maximum. Only gas-phase
experimental spectra (NIST WebBook) are valid references, because the calculations are for
an isolated molecule; a liquid-phase spectrum has solvent shifts of tens of cm⁻¹ on any
hydrogen-bonding mode.

> **Say it:** "We pair modes by eigenvector overlap with a global assignment, treat degenerate
> sets as subspaces, and compute errors only on physically matched pairs."

---

## 8. The worked example: water

All numbers from `comparison_results/water/`, B3LYP/6-31G(d,p) as reference, run July 2026.

| mode | experiment | DFT harm | DFT VPT2 | OMOL harm | OMOL VPT2 | MP harm | MP VPT2 |
|---|---|---|---|---|---|---|---|
| bend | 1595 | 1665 | 1615 | 1622 | 1570 | 1462 | 1445 |
| sym stretch | 3657 | 3799 | 3624 | 3818 | 3644 | 3852 | 3685 |
| asym stretch | 3756 | 3912 | 3722 | 3918 | 3739 | 3975 | 3780 |

Read it like this:

- VPT2 pulls DFT toward experiment on all three modes (harmonic error 70 to 160 cm⁻¹,
  anharmonic 20 to 35).
- MACE-OMOL harmonic frequencies are within 20 cm⁻¹ of DFT harmonic. The anharmonic
  correction is about the same size, so OMOL-VPT2 lands within 20 cm⁻¹ of DFT-VPT2 on the
  stretches and 45 on the bend. That is the "the ML surface has the right curvature and the
  right cubic terms" result.
- MACE-MP's bend is 200 cm⁻¹ too low. This run was at the OMOL geometry, where MP still
  has a residual force 28 times Gaussian's convergence threshold (finding H1). Rerun at MP's
  own minimum on 2026-09-18: bend 1497 harmonic / 1480 VPT2 (still 170 cm⁻¹ too low, so
  the soft bend is real), but the stretches move from 3852/3975 to 3696/3815, i.e. from
  stiffer than DFT to about 100 cm⁻¹ softer. At the wrong geometry MP was sitting on its
  repulsive wall. Same for methane: the −274 cm⁻¹ "imaginary modes" were the three
  rotations carrying the gradient; at MP's own minimum they vanish and the C-H stretches
  drop from 3163 to 3077. The Morse-wall picture from §3 is exactly what you see.

Intensities (km/mol), bend / sym / asym:

| | DFT | OMOL + MACE4IR | OMOL + espaloma |
|---|---|---|---|
| harmonic | 70 / 1.6 / 20 | 81 / 6.5 / 63 | 198 / 100 / 196 |
| VPT2 | 69 / 0.6 / 16 | 94 / 2.1 / 35 | 200 / 95 / 191 |

MACE4IR gets the ordering and the near-silent symmetric stretch right; espaloma makes all
three modes equally bright, which is what a fixed-charge model does.

Cost on the workstation (RTX 2070 Super, i7-6800K): DFT opt+VPT2 16 s wall, 74 s CPU.
ML VPT2 17 s wall, 1.6 s Gaussian CPU. On water the ML run is not faster, because 7 Hessians
of a 3-atom molecule are trivial for DFT and the ML pipeline has fixed overhead (model load,
process launch, and up to 1 s idle per call, finding M3). The crossover comes with size:
DFT Hessian cost grows roughly as N³ to N⁴ per Hessian times 6N−11 Hessians; ML cost grows
roughly linearly per Hessian times the same count.

---

## 9. What can go wrong, physically

These are the failure modes the benchmark is designed to find. Know them cold.

- **Bends for crystal-trained models.** Bending stiffness is a three-body property. A model
  trained mostly on dense periodic solids has seen few isolated-molecule bending environments.
  Expect MACE-MP to be soft on bends. (Confirm at its own minimum first, see H1.)
- **Elements outside the training set.** MACE-ANICC knows only H, C, N, O; the code refuses
  anything else (`workflow.py:327`). Chlorine is thin in most sets. HCl is in the panel as a
  deliberate out-of-domain probe.
- **Non-stationary geometry.** VPT2 is only meaningful at a minimum of the surface being
  differentiated. Residual gradient contaminates low modes and cubic constants.
- **Resonances that differ.** If the ML surface puts 2ν₂ 20 cm⁻¹ closer to ν₁ than DFT
  does, Gaussian may switch a resonance treatment on for one and off for the other, and the
  two spectra become hard to compare mode by mode. This is not an ML failure but it looks
  like one.
- **Large-amplitude motion.** Methyl rotors, floppy chains. VPT2 assumes small displacements
  around one minimum. Toluene is in the panel only as a documented limitation case.
- **Numerical noise in the Hessian.** ML energies are smooth but not analytically exact;
  finite differences of Hessians amplify noise into cubic constants. Gaussian's step size
  (default 0.025 in reduced normal coordinates) was chosen for DFT noise levels. Should be
  frozen and stated (TODO P0).
- **Dipole model and energy model at different geometries.** The dipole model is evaluated
  at whatever geometry Gaussian sends, which is on the energy model's surface. If the two
  models disagree about the equilibrium geometry the dipole derivatives are taken slightly
  off the dipole model's own minimum. Small effect, worth one sentence in the thesis.

---

## 10. The thesis question in physical terms

Everyone with a B3LYP harmonic spectrum multiplies by 0.96 and moves on. The scaling factor
compensates two things at once: the functional's systematic stiffness and the missing
anharmonicity. The question the thesis asks is whether you do better by spending compute on
the physics (VPT2) with a cheap surface (MACE) than by spending it on a better surface with
no anharmonic physics. Water says the ML-VPT2 fundamentals land within 20 to 45 cm⁻¹ of
DFT-VPT2 and within 15 to 25 cm⁻¹ of experiment on the stretches, and additionally produce
the overtone and combination bands that scaled-harmonic cannot produce at all. Whether that
holds across 22 molecules is the benchmark.
