# Hard questions, and answers that hold up

Questions a committee, a supervisor, or a referee will ask, sorted by where they aim.
Each has a short answer to say and the longer reasoning to have ready. Where the honest
answer is "the current code does not do that yet", it says so and points to
`REVIEW_FINDINGS.md`.

---

## A. About the physics of doing VPT2 on an ML surface

**A1. Why should anyone trust third and fourth derivatives of a neural network?**

*Say:* "We never ask the network for third derivatives. Gaussian gets them the same way it
gets them from DFT: finite differences of Hessians at displaced geometries. The Hessian is an
exact second derivative of a smooth function by automatic differentiation. What we rely on is
that the network is smooth and accurate in a small neighbourhood of the minimum, which is
what the harmonic and anharmonic agreement with DFT tests directly."

*Have ready:* MACE uses smooth radial basis functions and polynomial cutoffs, so E(R) is C∞.
The risk is not non-differentiability but noise amplitude relative to Gaussian's displacement
step (0.025 in dimensionless normal coordinates). Water and methane show cubic corrections
matching DFT to within 10 to 20 cm⁻¹, which would not happen if the third derivatives were
noise. The open item is to freeze and state the step size (TODO P0) and, for one molecule,
show that halving it changes nothing.

**A2. Your ML models were trained on ωB97M-V, PBE, or CCSD(T)-quality data. Your reference
is B3LYP/6-31G(d,p). Isn't a disagreement just functional versus functional?**

*Say:* "Partly, yes, and the thesis says so. That is why experiment is in the benchmark as
the neutral arbiter. The ML-vs-DFT comparison tests whether the pipeline and the surface are
sane; the ML-vs-experiment comparison is the one that decides whether the method is useful."

*Have ready:* B3LYP/6-31G(d,p) harmonic frequencies are typically 3 to 5 % above
experiment; ωB97M-V/def2-TZVPD is closer. So an ML model that is faithful to its training
level should sit *below* B3LYP on stretches. Water OMOL harmonic stretches: 3818/3918 vs
B3LYP 3799/3912. Not obviously lower. That is itself a result.

**A3. You compute frequencies for four models at a geometry optimized with a fifth. VPT2 is
only valid at a minimum. How do you justify that?**

*Say:* "That was a design choice to isolate the curvature of the surface from geometry
differences, and the review this September showed it is not tenable for VPT2: MACE-MP has a
residual force 28 times Gaussian's threshold at the OMOL geometry and produces imaginary
modes on methane. The fix is to re-optimize with each model before its frequency run, which
costs seconds, and the campaign will be run that way."

*Have ready:* finding H1. Gaussian does not project the gradient out of the Hessian unless
`freq=projected`. The "MACE-MP fails on bends" story must be re-established at MP's own
minimum before it goes in the thesis. Note that the DFT twin already re-optimizes at B3LYP,
so "same geometry" was never true for the DFT comparison anyway.

**A4. What happens when the ML surface and DFT put a Fermi resonance on different sides of
Gaussian's threshold?**

*Say:* "Then Gaussian applies deperturbation to one and not the other, and those modes are
not mode-by-mode comparable. We detect it from the log (the resonance tables), report it,
and treat it as a failure mode in its own right, because for a user it is one: a small shift
in the surface changed the spectrum qualitatively."

*Have ready:* CO₂ and formic acid are in the panel as exactly this stress test. Gaussian
prints `1-2 Fermi resonances:` and `Darling-Dennison resonances:` blocks; the parser does
not read them yet. Worth one small parser addition before the campaign.

**A5. Why do you need a second model for intensities? Can't you differentiate the energy
model's charges?**

*Say:* "Energy models have no charges. MACE predicts a scalar energy per atom and nothing
else. The dipole is a separate observable and needs its own head or its own model. MACE4IR
predicts the dipole vector equivariantly and we differentiate it by autograd."

*Have ready:* POLAR-1 is the exception; it carries a multipole expansion and gives a dipole as
a by-product. Using it for both energy and dipole is in the outlook. Espaloma gives charges
but they are topology-only, so its derivatives are a fixed-charge model: no charge flux
(finding L1). That explains why espaloma makes water's silent symmetric stretch bright.

**A6. Is VPT2 the right level? Why not variational (VCI, VSCF) or a full-dimensional
treatment?**

*Say:* "VPT2 is the level that practitioners actually use, it is what the reference code
implements, and it is where the cost wall is. The question is whether cheap surfaces make the
standard method routine, not whether a better method exists. Once the surface is cheap,
variational methods become the natural next step; that is in the outlook."

**A7. Harmonic DFT scaled by 0.96 gets fundamentals within ~30 cm⁻¹ for free. Why is your
approach better?**

*Say:* "Scaling compensates two errors at once, the functional's stiffness and missing
anharmonicity, with one number tuned on a training set. It cannot give overtones or
combination bands at all, it cannot handle resonances, and it fails for exactly the modes
where anharmonicity is largest. Whether anharmonic-ML beats scaled-harmonic-DFT against
experiment is research question 3 and the central plot; the answer is empirical."

---

## B. About the pipeline and its correctness

**B1. How do you know your bridge is not introducing errors? A unit mistake would be
invisible.**

*Say:* "Three independent checks. First, every conversion is applied exactly once at the file
boundary and was audited line by line. Second, the water numbers are physically sensible in
absolute terms: an ML harmonic bend at 1622 versus DFT 1665 and experiment 1595 cannot happen
with a unit error, which would be a factor of 27 or 0.53. Third, the Psience cross-check: an
independent VPT2 code fed the same force constants reproduces Gaussian's anharmonic
frequencies to better than 0.01 cm⁻¹."

*Have ready:* The July 2026 intensity bug is a good story: intensities were read from the
dipole-strength table (10⁻⁴⁰ esu²·cm²) instead of km/mol, a factor 3.57. It was caught by
the mismatch between harmonic and anharmonic intensity scales, and it is exactly the kind of
error the audit is for. The loop on the harness itself is closed by the DFT self-consistency
check (2026-09-18): B3LYP/6-31G(d,p) run *through* the external interface, with Gaussian as
the "ML model", versus native `freq(anharm)` at the identical geometry. Water: max |Δν|
0.008 cm⁻¹ over all fundamentals, overtones and combination bands, intensities within
0.02 % (0.11 % with the finite-difference dipole derivatives the ML models use).
Formaldehyde: 0.034 cm⁻¹, 0.22 %. Tables in `docs/explained/self_consistency/`.

**B2. What exactly does Gaussian still do, and what does the ML model do?**

*Say:* "Gaussian: reads the geometry, decides which derivatives it needs, generates the
displaced geometries along normal modes, assembles cubic and quartic constants by finite
differences, runs the VPT2 algebra, detects and deperturbs resonances, projects dipole
derivatives onto modes for intensities, prints everything. The ML side: for a given geometry,
return energy, gradient, Hessian, dipole, and dipole derivatives in atomic units. Nothing
else."

**B3. How many times is the ML model called, and how does that scale?**

*Say:* "6N−11 Hessian calls: one at equilibrium plus two per mode. Water 7, methane 19,
naphthalene 97, decane 181. Each call is a forward pass, a backward pass for forces, a
second differentiation for the Hessian, and a dipole-model pass. On the RTX 2070 that is
under a second per call for anything in the panel. DFT is the same count of far more
expensive Hessians."

**B4. Your water ML run took as long as the DFT run. Where is the speedup?**

*Say:* "For a 3-atom molecule there is none, and the thesis will say so. The ML pipeline has
fixed costs: loading a foundation model, launching a helper process per call, and a polling
loop that adds up to a second per call. The crossover comes with size because a DFT Hessian
scales as N³ or worse and an ML Hessian roughly linearly. The alkane ladder measures where
that crossover is on this hardware."

*Have ready:* finding M3, the one-second polling sleep. It should be fixed before the cost
chapter is measured, otherwise the ML times carry up to 0.5 s × (6N−11) of pure idle.

**B5. Why not implement VPT2 yourself and drop Gaussian entirely?**

*Say:* "That is the outlook, and the spike showed it is feasible: Psience reproduces
Gaussian to 0.01 cm⁻¹. For the thesis the point is the controlled comparison, which needs
the same VPT2 code on both sides. A Gaussian-free pipeline would also lose the licence
constraint, which is a real practical benefit."

**B6. What is your test coverage and how do you know a parser change does not break the
numbers?**

*Say:* "There is a pytest suite of 27 files with fixture logs, including regression tests on
parsed frequencies. It has known failures that predate the branch and are tracked; CI is
disabled until they are fixed. The stronger guarantee is the `results.json` snapshots under
version control for water, which any parser change is diffed against."

---

## C. About the analysis and metrics

**C1. How do you match modes between two calculations?**

*Say:* "By eigenvector overlap, not by frequency order. We take the harmonic normal modes
from both checkpoint files, compute the full matrix of absolute dot products, and solve the
assignment problem with the Hungarian algorithm for the one-to-one pairing with maximum total
overlap. Pairs under 0.5 are flagged. Degenerate sets are compared as subspaces, because
individual vectors inside a degenerate subspace are arbitrary up to rotation."

*Have ready:* two caveats from the review. The overlaps are currently computed on Cartesian
displacements rather than mass-weighted eigenvectors, which are the ones that are actually
orthonormal (H3, two-line fix). And in the anharmonic report the mapping is applied using
Gaussian's symmetry-block mode numbering as if it were the checkpoint's ascending index (H2).
Both are known and will be fixed before the campaign's analysis is run. The harmonic report
is unaffected by H2.

**C2. Why use harmonic eigenvectors to match anharmonic frequencies?**

*Say:* "Because the harmonic normal modes are the basis in which VPT2 is formulated. The
anharmonic 'modes' Gaussian writes are perturbed mixtures. Matching on the harmonic basis
and then looking up the anharmonic frequency of each matched mode is the consistent choice."

**C3. Your R² values are all 0.999. Isn't that meaningless?**

*Say:* "R² on frequencies is nearly meaningless, yes, because the frequency range spans 400
to 4000 and any method gets the gross ordering right. It is reported because readers expect
it. The informative numbers are MAE and RMSE in cm⁻¹ on matched pairs, the slope (systematic
stiffness), and per-region errors, which is why the coverage analysis bins by frequency
region."

**C4. What about modes the ML model predicts that DFT does not, or vice versa?**

*Say:* "For fundamentals the count is fixed by 3N−6, so every mode is matched, possibly with
low overlap. Imaginary modes are excluded from statistics and counted; an imaginary frequency
is a geometry problem, not a negative wavenumber. For overtones and combinations, the count
is the same on both sides too, since Gaussian prints all of them, but weak ones are filtered
by the 0.1 km/mol threshold on the DFT side."

**C5. How do you compare to experiment when the calculation is for an isolated molecule?**

*Say:* "Only gas-phase NIST spectra, only for molecules that have one. Sticks are broadened
with a 10 cm⁻¹ Lorentzian and both curves normalized. Agreement is judged visually and by a
correlation on the broadened curves; that correlation is deliberately kept out of the headline
metrics because it is near noise for small molecules."

**C6. You have 22 molecules. Is that enough?**

*Say:* "For a functional-group heatmap with 5 models, yes: about 20 modes per functional
group per model. Every molecule is there to test one specific thing, and the panel is
documented with the reason for each. It is not a statistical survey of chemical space; it is
a structured probe of where the method works and where it breaks."

---

## D. About the ML models themselves

**D1. Explain MACE in two minutes without equations.**

*Say:* "Each atom looks at its neighbours within a cutoff and builds a description of that
environment. The description is built from vectors, not just distances, using spherical
harmonics, so it rotates correctly when the molecule rotates. MACE's specific contribution is
to build many-body correlations, three-, four-body and beyond, in a single tensor product
step (the Atomic Cluster Expansion) instead of stacking many message-passing layers. Two
layers of that, then a small network maps each atom's description to an energy. Sum the
atomic energies. Differentiate for forces. Because everything is built from smooth functions
of positions, the second derivative is also exact and smooth."

**D2. Why five models?**

*Say:* "They span training-data provenance: crystals (MP), organic molecules at a hybrid
functional (OFF), CHNO at coupled-cluster quality (ANICC), a very large diverse molecular set
(OMOL), and OMOL plus polarization targets (POLAR-1). The hypothesis is that spectroscopic
behaviour follows provenance: crystal-trained models should struggle with molecular bends,
CHNO-only models should be excellent inside their domain and refuse outside it."

**D3. Foundation models versus a potential trained for this molecule?**

*Say:* "A bespoke potential would be more accurate and would cost a DFT training set per
molecule, which defeats the purpose. The question here is whether off-the-shelf models are
good enough with zero training, because that is the regime where anharmonic spectra become
routine."

**D4. Uncertainty. How do you know when the model is extrapolating?**

*Say:* "We do not, in this version. Single models, no committee, no uncertainty estimate.
Out-of-domain behaviour is probed empirically with HCl and the element guard. Committee
disagreement as an uncertainty proxy is in the outlook and is straightforward to add because
the calculator interface already returns per-model outputs."

---

## E. The uncomfortable ones

**E1. What did you find wrong in your own code, and what did it change?**

*Say:* "Two things this summer. Intensities were being read from the wrong table by a factor
of 3.57, and a wrapped archive line truncated the dipole. Both are fixed and all older results
are flagged for rerun. The September review found four more: frequencies computed at another
model's minimum, a mode-index mismatch in the anharmonic matching, Cartesian rather than
mass-weighted overlaps, and a regex that stored the thermal correction as the energy. None
of them changed a frequency that has been reported so far, because the water and methane
mappings are the identity and energies are not used, but all four will be fixed before the
campaign. I would rather present a list of things I found and fixed than claim there were
none."

**E2. Why is a master's student's pipeline more trustworthy than a published tool?**

*Say:* "It is not more trustworthy; it is more transparent. Every number in the thesis has a
`results.json`, a Gaussian log, a checkpoint, and a version stamp behind it, and every figure
regenerates from `report_data.json`. The validation section shows the audit, the
self-consistency check, and the independent VPT2 cross-check. Trust comes from that, not from
the author."

**E3. What is the single most important limitation?**

*Say:* "VPT2 itself. It assumes small-amplitude motion around one minimum. Methyl rotors,
floppy chains, and strong resonances break it regardless of the surface. The thesis restricts
headline claims to rigid molecules and labels the rest as limitation cases. A cheap surface
does not fix a wrong method; it makes the right method affordable."

**E4. If you had one more month, what would you do?**

*Say:* "Fix the four findings, run the DFT self-consistency check through the external
interface, re-run the panel with per-model optimization and a fixed step size, and add
committee uncertainty. In that order."
