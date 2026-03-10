# Fundamental Concepts Track

Concepts every examiner expects you to know cold — independent of this specific project.
Quizzes are in `quizzes/fundamentals_quiz_01.md` (once created).

---

## Energy & Units

| Unit | Used for | Key value |
|---|---|---|
| Joule (J) | SI, single-molecule | ~10⁻²¹ J per bond energy |
| eV | Atomic physics, ML potentials | 1 eV = 1.6×10⁻¹⁹ J |
| kcal/mol | Chemistry, thermodynamics | Bond energies, ΔG |
| kJ/mol | Chemistry (SI) | 1 kcal = 4.18 kJ |
| cm⁻¹ | Spectroscopy | E = hcν̃; IR range 400–4000 cm⁻¹ |
| Hartree | Quantum chemistry | 1 Eₕ = 627.5 kcal/mol = 27.2 eV |

**Must-know conversions:**
- kT at 300K ≈ 25 meV ≈ 0.6 kcal/mol ≈ 200 cm⁻¹
- 1 eV = 23 kcal/mol = 96 kJ/mol = 8065 cm⁻¹
- H-bond: ~10–40 kJ/mol; C-C bond: ~350 kJ/mol; O-H stretch: ~3500 cm⁻¹

---

## Quantum Mechanics Essentials

**Zero-point energy (ZPE)**
- Origin: Heisenberg uncertainty principle forbids simultaneous precise position + momentum
- Result: lowest energy of harmonic oscillator = ½ℏω (not zero)
- Consequence: molecules vibrate even at 0K
- Observable: Casimir effect, kinetic isotope effect, tunneling rates

**Heisenberg uncertainty principle**
- ΔxΔp ≥ ℏ/2 — not a measurement limitation, a property of reality
- Position and momentum don't simultaneously have definite values

**Quantum tunneling**
- Wavefunction has finite amplitude in classically forbidden regions (barriers)
- Probability decays exponentially with barrier width and √(mass × barrier height)
- Light atoms (H, D) tunnel significantly; heavy atoms don't
- Kinetic isotope effect: H transfers faster than D due to tunneling + ZPE difference

**Cosmological constant problem**
- QFT predicts vacuum ZPE ~10¹²⁰× larger than observed — worst prediction in physics
- Why the universe's vacuum energy is small: open problem

---

## Vibrational Spectroscopy

**Normal modes**
- Collective motions of all atoms vibrating at a single frequency
- N atoms → 3N degrees of freedom → 3N−6 vibrational modes (nonlinear), 3N−5 (linear)
- From diagonalizing mass-weighted Hessian
- Orthogonal to each other (don't mix in harmonic approximation)

**IR selection rule**
- A mode is IR-active only if the dipole moment changes during the vibration
- ∂μ/∂Q ≠ 0 for IR activity
- Centrosymmetric molecules: IR-active modes are Raman-inactive (mutual exclusion rule)

**Harmonic vs anharmonic**
- Harmonic: potential is a perfect parabola (V ∝ x²); frequencies exact, no overtones
- Real bonds: asymmetric potential (can dissociate!); frequencies shift, overtones appear
- VPT2 correction: perturbation theory using 3rd + 4th derivatives

**Scaling factors**
- B3LYP/6-31G(d,p) harmonic frequencies ~5–12% too high
- Empirical fix: multiply by ~0.96
- Better fix: use anharmonic (VPT2) corrections

---

## Quantum Chemistry Basics

**Born-Oppenheimer approximation**
- Electrons are ~1800× lighter than nuclei → treat them as moving instantaneously
- Separate electronic and nuclear Schrödinger equations
- Foundation of all quantum chemistry; breaks down for excited states near degeneracy

**DFT**
- Density Functional Theory: replace many-body wavefunction with electron density ρ(r)
- Hohenberg-Kohn theorem: ground state energy is a functional of ρ alone
- Exchange-correlation functional (B3LYP, ωB97X-D...) is the approximation
- Kohn-Sham equations: non-interacting electrons in effective potential — solved self-consistently (SCF)

**Basis sets**
- Mathematical functions used to expand molecular orbitals
- 6-31G(d,p): split-valence double-zeta + polarization functions on heavy atoms (d) and H (p)
- def2-TZVP: triple-zeta with polarization — more accurate, more expensive
- Larger basis = more accurate but O(N³) cost scaling with basis size

**Common functionals and their uses**
- B3LYP: workhorse, good for structures/frequencies, underestimates dispersion
- ωB97X-D: includes dispersion, better for non-covalent interactions
- PBE: pure GGA, fast, used in solid-state

---

## Statistical Concepts (for analysis)

**MAE**: mean absolute error — average of |predicted − true|
**RMSE**: root mean square error — penalizes large errors more
**R²**: coefficient of determination — 1.0 = perfect correlation, 0 = no correlation
**KDE**: kernel density estimation — smooth histogram; used to broaden spectral peaks
