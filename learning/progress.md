# Learning Progress

## Topics Covered

_Nothing yet — this gets filled in as we go._

---

## Known Gaps

- User hasn't explained the Hessian back yet
- Taylor expansion / smooth functions — introduced, not tested back yet
- Dipole moment as a vector — introduced, not tested back yet

---

## Session Log

### 2026-02-22 — Session 1: Setup + Big Picture

**Covered:**
- Big picture: what the project does and why (IR spectra, DFT vs ML tradeoffs, the dipole gap)
- Zero-point energy: QM origin, philosophical implications, cosmological constant problem, Casimir effect
- Why molecules vibrate at 0K and how temperature affects spectra (populations, not frequencies)
- Hydrogen tunneling and kinetic isotope effect
- DFT accuracy for IR: systematic 5–12% overestimate harmonically, VPT2 fixes this
- ML accuracy: energy/force errors, amplification when differentiating twice
- Hessian: definition, why 9×9 for water, 3N−6 real modes, eigenvalue → ω²
- Eigenvectors = normal modes, orthogonality = modes don't interact (harmonic approx)
- Dipole moment, CPHF, Gaussian as VPT2 orchestrator
- Energy units: J, eV, meV, kcal/mol, kJ/mol, cm⁻¹, Hartree; kT at 300K ≈ 200 cm⁻¹
- Hot bands, kinetic isotope effect, error amplification when differentiating
- Full table of what uses autograd vs numerical FD in the pipeline
- Vibrational transitions: photon adds energy to existing vibration, v=0→v=1

**Gaps after session 1:**
- VPT2 mechanics not fully covered
- MACE architecture not covered
- Basis sets not explained
- Matrix diagonalization deferred

---

### 2026-02-22 — Session 2: Vibrational Transitions + Math Foundations

**Covered:**
- Open question answered: v=0→v=1 = amplitude increases, frequency unchanged, photon absorbed
- Introduced overtone problem: why v=0→v=2 is forbidden in harmonic, allowed in reality
- Smooth functions, Taylor expansion, cubic wiggle, flat bottom / V'(0)=0, dipole as vector
- Got through all the building blocks but did NOT yet reconnect to the overtone argument

**Stopped mid-session.** Full content + resume instructions saved at:
`learning/topics/taylor_expansion_and_overtones.md`

---

### 2026-02-23 — Session 3: Taylor Expansion Deep Dive + Multipole Expansion

**Covered:**
- Higher derivatives named: jerk, snap, crackle, pop — intuition for each level
- "Small x" clarified: you choose the expansion point; x-a is only small if x stays near that point. Not automatic — the physics has to guarantee it (small vibrations around equilibrium)
- Radius of convergence: can't make large x small by shifting; series breaks down when x-a gets large
- How Taylor was discovered: match height, slope, curvature one by one; smooth functions are locally polynomial
- k = V''(0) confirmed: spring constant is literally the curvature of the potential well
- V'(x) = -F(x): first derivative of potential is negative force. V'(0)=0 means zero net force at equilibrium
- V'''(0): encodes asymmetry of the bond — compression wall steeper than stretch wall. Why overtone frequencies are slightly less than 2× the fundamental (energy rungs squeeze together)
- Multipole expansion: monopole → dipole → quadrupole → octupole, each capturing finer charge detail
- Multipole as information reduction: tunable by distance and required precision
- When higher multipoles needed: (1) lower ones vanish by symmetry (CO₂, N₂ have no dipole), (2) close-range interactions (MD force fields), (3) precision binding calculations
- Why dipole approximation is valid for IR: photon wavelength (~10,000×) >> molecule size → photon can't resolve internal structure

**User homework before next session:**
- Watch 3Blue1Brown "Taylor series" (Essence of Calculus series)
- Search "Morse potential vs harmonic potential" for potential well visualization

**Resume checklist for next session:**
1. Check if the videos clicked — any new questions from them?
2. Walk through the overtone argument (we have all pieces, never closed the loop):
   - μ(x) = μ₀ + μ'·x + μ''·x²/2
   - Substitute x(t) = A·sin(ωt)
   - sin²(ωt) = ½(1 - cos(2ωt)) → 2ω component → overtone
3. Connect to VPT2

---

## What's Ready for a Quiz

_Topics will move here once we've covered enough ground._
