# Brainstorming Session: Research Landscape & PhD Directions
*2026-03-26*

## Key Concepts Covered

### Multi-task Loss Balancing
When training one model on multiple properties (energy + dipole + polarizability), the loss is typically a weighted sum: L = w₁L₁ + w₂L₂ + w₃L₃. Different scales/sensitivities make weights tricky. Methods: manual tuning, uncertainty weighting (learnable weights), GradNorm (gradient magnitude normalization), Pareto optimization, gradient surgery (PCGrad). Also possible: alternating optimization (train one task at a time).

### Multi-head MACE
One shared message-passing backbone, multiple small output heads. Solves the "different models see the molecule differently" problem. Used for: (a) multi-property prediction (energy head + dipole head + polarizability head), (b) committee uncertainty (multiple heads predicting same property on different data subsets → disagreement = uncertainty).

### Fine-tuning ML Potentials
Take pretrained foundation model (e.g., MACE-MP-0), adjust with small dataset (~100-1000 structures). Key: small learning rate (1e-4 to 1e-5), optionally freeze backbone layers, watch for catastrophic forgetting. MACE supports this via `--foundation_model` flag in `mace_run_train`.

### Foundation Model Definition
Trained on broad diverse data, transferable out-of-the-box, fine-tunable with minimal extra data, scales with data/model size. MACE-MP-0 and MACE4IR qualify. ANI-2x does not (only 7 elements — "transferable potential," not foundation model).

### Normal Mode Displaced Geometries
Training data that includes atoms displaced along vibrational normal modes, not just equilibrium structures. Critical for learning dipole derivatives (IR intensities) and force constant curvature (frequencies). Equilibrium-only data gives one energy/dipole point; displaced geometries teach the model how properties *change* during vibration.

### Degenerate Modes
Symmetry-equivalent vibrations with identical frequencies (e.g., methane T₂ modes are 3-fold degenerate). Eigenvectors within a degenerate subspace have arbitrary orientation — any rotation is equally valid. This causes artificially low dot products in mode matching even when modes are correct. Fix: compute subspace overlap (trace of M^T M) instead of individual vector dot products. Also: wrong PES → wrong eigenvectors on top of the rotation ambiguity.

### Why Bending Modes Are Harder Than Stretching
- Stretching = two-body (pairwise distance). Strong, clear signal. All ML potentials learn this well.
- Bending = three-body (angular). Subtler, flatter energy variation. Requires training data with angular diversity.
- mace_mp (trained on crystals) fails for methane bends (-193 cm⁻¹ error) because crystal training data lacks molecular bending diversity. mace_anicc (trained on molecular data) gets -10 cm⁻¹.

### Vacuum vs Experiment
- O-H/N-H: huge solvent shifts (-300 to -550 cm⁻¹), unpredictable, H-bonding dominated
- C=O: moderate shifts (-10 to -25 cm⁻¹), systematic, Stark effect
- C-H/C-C: negligible (<5 cm⁻¹)
- Best comparison: gas-phase NIST data or matrix isolation experiments
- For thesis: comparing ML vs DFT (both vacuum) is valid for benchmarking the ML model

## Alternative VPT2 Implementations (Replace Gaussian)

| Code | Type | Integration | Notes |
|------|------|-------------|-------|
| **PyVPT2** | Python VPT2 via QCEngine | MACE ASE calc → QCEngine → PyVPT2 | **Best option.** Published Jan 2025. Eliminates Gaussian entirely. |
| **iGVPT2** | C/C++ GVPT2, plugin system | Write bash script calling MACE | Good Fermi resonance handling |
| **PyVCI** | Python/C variational VCI | Pre-compute force constant tensors | More accurate than VPT2 for pathological cases |
| **Psience/PyVibPTn** | Python VPT2 (McCoy group) | Manual setup | Less maintained |
| **MULTIMODE** | Fortran VSCF/VCI (Bowman) | Hard, not open source | Gold standard science |

CFOUR and ORCA: no external interface, cannot inject ML potentials. Dead ends.

Prior art: PhysNet+VPT2 (2021), ANI+VPT2 (2025), MACE-OFF23+VPT2 (2024). Nobody has done PyVPT2+MACE foundation models yet.

## MACE4IR Paper Details

- "MACE4IRmol" by Bhatia, Krejci, Botti, Rinke, Marques (arXiv:2508.19118)
- Trained on QCML dataset (~10-16M structures from 33.5M total, PBE0/FHI-aims)
- v2 has ensemble-based UQ built in
- Already does MD-based IR spectra via dipole autocorrelation function (DACF)
- IS a foundation model (~80 elements)
- Opportunity: VPT2 vs DACF head-to-head comparison using same models

## Multi-property Foundation Models

- **MACE-POLAR-1** (Feb 2026): 100M structures, energy + polarizable charges. Closest to unified model.
- Gap: nobody has one model for energy + dipole + polarizability at MACE-MP-0 scale
- Reason: data scarcity (dipole/polarizability data is orders of magnitude smaller than energy data) + loss balancing difficulty

### Polarizability Datasets for Raman
- QMe14S (2025): 186K molecules, 14 elements, full tensors + Raman — best current option
- QMugs: ~2M molecules but lower QM level
- QM9S: 130K, H/C/N/O/F only
- All much smaller than QCML's 33.5M for energies — polarizability data is the bottleneck

## Active Learning for Spectroscopy

- **PALIRS** (2025, npj Comp. Mat.): Already does AL for IR with MLIPs. 100x fewer DFT calcs.
- MACE supports committee models natively (multi-head or multiple seeds)
- Practical: train 3 models, flag structures where dipole predictions disagree > 0.05 Debye, run DFT on those, retrain
- Key papers: PALIRS, Multi-head committees (J. Chem. Phys. 2025), UDD-AL (Nature Comp. Sci. 2023)

## PhD Directions (Ranked by Novelty)

### Genuinely new (not benchmarking):
1. **Experimental Feedback Loop** — use experimental IR spectra as training signal to improve ML potentials. Han & Yu (2025, Nature Comms) did differentiable MD with computed reference; nobody has used actual experimental data yet. Strongest PhD candidate.
2. **Differentiable Spectroscopy** — make entire pipeline (ML potential → displacements → VPT2 → spectrum) end-to-end differentiable. Backpropagate spectrum error to improve potential.
3. **Anharmonic Transfer Learning** — train on cheap harmonic data (millions of structures), fine-tune on expensive anharmonic data (hundreds). Learn harmonic→anharmonic mapping with minimal VPT2 data.
4. **Vacuum-to-Experiment Learned Correction** — learn systematic solvent shifts as function of molecular environment and functional group.

### Valuable but more benchmarking-flavored:
5. VPT2 vs DACF systematic comparison with same ML models
6. PyVPT2 + MACE foundation models at scale
7. UQ calibration (do ensemble error bars match actual errors?)

## Research Groups in Central Europe

### Most relevant (by topic match + proximity to Innsbruck):
- **Gonzalez / Marquetand** — U. Vienna. ML + spectroscopy + dynamics. SchNarc. Most natural fit.
- **Reiher** — ETH Zurich. Reaction mechanisms (SCINE) + vibrational spectroscopy. World-leading.
- **Ceriotti** — EPFL. ML potentials, representations, spectroscopic properties.
- **Behler** — U. Göttingen. Founded the field of ML potentials. ML dipole surfaces for IR.

### Other strong options:
- **Müller / Schütt** — TU Berlin. SchNet, PaiNN architectures. sGDML for IR/Raman.
- **Noé** — FU Berlin. Boltzmann generators, biomolecular ML. Published in Science.
- **Reuter** — Fritz Haber Institute, Berlin. ML for catalysis.
- **Tkatchenko** — U. Luxembourg. ML force fields, van der Waals.
- **Riniker** — ETH Zurich. ML for drug design, free energy calcs.
- **Dellago / Franchini** — U. Vienna. Enhanced sampling, materials.
- **MACE4IR authors** (Bhatia, Marques) — MLU Halle-Wittenberg. Direct collaboration potential.

### PhD application advice:
- Finish v1.1 benchmark campaign + one experimental comparison
- Prototype PyVPT2 integration
- Clean public GitHub repo
- Email PIs directly with 1-page summary + benchmark figure + repo link
- Strongest pitch: "I built open-source ML+anharmonic spectroscopy pipeline, benchmarked across 30 molecules and 5 foundation models, want to extend to [their direction]"

## Todos Created This Session
- Clean up temp files permanently in base directory
- NIST/SDBS experimental spectra overlay
- JCAMP-DX spectral data export
- Interactive HTML spectrum viewer (Plotly)
- Conformer-aware spectra (Boltzmann weighting)
- Automatic functional group peak labeling
- Peak assignment confidence scores
- Delta-ML correction model
- Uncertainty quantification for ML spectra
- Active learning loop for dipole model
- Wall-clock timing comparison ML vs DFT
- Automated experimental spectra comparison via NIST API
- Integrate PyVPT2 and alternative anharmonic methods
- Reevaluate MACE-POLAR-1 as dipole calculator
- Handle degenerate modes and zero-intensity filtering
