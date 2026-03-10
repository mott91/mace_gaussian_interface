> Note: This document describes the computational methodology as implemented.
> It is written for inclusion in a thesis methods section. For implementation details,
> see [docs/ARCHITECTURE.md](ARCHITECTURE.md).

---

# Computational Methods

## 1. Computational Framework Overview

Molecular infrared spectra were computed through a pipeline that combines machine learning
(ML) interatomic potentials with Gaussian 16's anharmonic frequency calculation framework.
The pipeline proceeded in four stages: geometry optimization using an ML potential,
frequency calculation with ML-derived energy and dipole quantities injected into Gaussian
via an external interface, parsing of the resulting frequency data, and statistical analysis
of the predicted spectra against DFT reference calculations.

Two types of frequency calculations were performed for each molecule: anharmonic calculations
based on second-order vibrational perturbation theory (VPT2), which constitute the primary
analysis and capture overtone and combination bands, and harmonic calculations restricted to
fundamental frequencies. Four combinations of ML models were evaluated systematically: two
potential energy models (MACE-MP-0 and MACE-OMOL-0) were each paired with two dipole models
(Espaloma and MACE4IR), yielding four independent predictions per molecule.

DFT reference calculations at the B3LYP/6-31G(d,p) level of theory were performed using
Gaussian 16 as the ground truth for frequency comparison. The goal was to assess how closely
ML-based predictions reproduce DFT anharmonic frequencies; absolute accuracy against
experimental spectra was not the primary objective of this study.

## 2. MACE Potential Energy Models

The MACE (Multi-ACE) architecture [Batatia et al., 2022] was used for potential energy
surface evaluation. MACE is an E(3)-equivariant graph neural network that employs
higher-order equivariant message passing to represent atomic environments. The equivariance
property — invariance of energy predictions to rotation, reflection, and translation of the
molecular system — ensures physical consistency of the potential energy surface without
imposing these symmetries as explicit constraints during training.

Two pre-trained MACE models were selected for this study based on differences in their
training data distributions. MACE-MP-0 [Batatia et al., 2023] was pre-trained on
approximately 150,000 diverse inorganic and organic structures drawn from the Materials
Project database, providing broad coverage across chemical space. MACE-OMOL-0 [OMol25
Authors, 2025] was pre-trained on the Open Molecules 2025 (OMol25) dataset, which focuses
on organic molecules and is expected to yield better descriptions of hydrogen-bonded and
drug-like systems. Both models were used without further fine-tuning on IR-specific data.
The inclusion of both models allowed assessment of how training data distribution affects the
accuracy of ML-predicted vibrational frequencies.

Molecular geometry optimization was performed with MACE-OMOL-0 as the default energy model,
using the Atomic Simulation Environment (ASE) [Larsen et al., 2017] L-BFGS optimizer with
convergence criteria of 0.01 eV/Å on forces. All four energy+dipole combinations then used
the same optimized geometry as starting point to ensure that frequency differences across
models arose from the potential energy surface representation rather than from differences in
equilibrium geometry.

## 3. Dipole Derivative Calculation

Infrared intensities are determined by the derivative of the molecular dipole moment with
respect to nuclear displacement. Two ML backends were compared for dipole moment prediction:
Espaloma [Wang et al., 2022] and MACE4IR.

Espaloma is a graph neural network trained to predict partial atomic charges. Dipole moments
were derived from these predicted partial charges weighted by atomic positions. This approach
is computationally inexpensive and has demonstrated good transferability across diverse
organic molecules.

MACE4IR is a custom fork of the MACE architecture trained specifically for direct dipole
moment prediction, covering 78 elements. Rather than predicting charges and inferring dipole
moments, MACE4IR directly outputs a molecular dipole vector, making it more directly suited
to IR intensity prediction. MACE4IR was trained on data selected for relevance to molecular
IR spectroscopy.

Dipole derivatives required for IR intensity calculation were computed by numerical finite
differences. Each of the N atoms in the molecule was displaced by ±δ along each Cartesian
axis, yielding 6N molecular configurations. The dipole moment was evaluated at each
displaced geometry, and finite-difference derivatives were assembled into the dipole
derivative tensor. The four combinations of energy model and dipole model (MACE-MP-0 or
MACE-OMOL-0 paired with Espaloma or MACE4IR) allowed independent assessment of the
contribution of each model choice to spectral accuracy.

## 4. ZMQ Injection Mechanism

The central methodological contribution of this work is a software bridge that enables VPT2
anharmonic frequency calculations in Gaussian 16 using ML-derived potential energy surfaces
and dipole moments. This capability is not natively available from ML potential packages,
which typically provide energy and forces but do not integrate with the anharmonic frequency
analysis infrastructure of quantum chemistry codes.

The bridge was implemented using Gaussian's `External=` keyword, which instructs Gaussian to
delegate energy, force, and dipole evaluation to a user-supplied external program rather than
performing a self-consistent field calculation. At each step of the VPT2 numerical
differentiation — which requires evaluation of the energy surface at hundreds of displaced
geometries to extract cubic and quartic force constants — Gaussian launched the helper
script `gm_helper.py` as a subprocess, passing the current geometry in Gaussian's standard
external interface format.

A main Python process was initialized before each Gaussian run. This process loaded the ML
models into GPU memory and bound a ZMQ inter-process communication (IPC) socket using a
Unix domain socket file. When Gaussian invoked `gm_helper.py` for a given geometry,
the helper script connected to this IPC socket, transmitted the geometry coordinates to the
waiting main process, and blocked until a response was received.

The main process, upon receiving a geometry, evaluated the ML energy model to obtain total
energy and atomic forces, and separately evaluated the dipole model to obtain the molecular
dipole moment. Where required by the Gaussian external interface protocol, dipole
derivatives were also computed by the main process via additional displaced geometry
evaluations. Results were written to the output file in Gaussian's expected external
interface format, and a "done" signal was sent to the helper script over the ZMQ socket.
Gaussian then read the output file and continued the frequency computation.

Unix domain sockets (IPC) were chosen over network sockets for inter-process communication
to minimize per-call latency. The VPT2 derivative evaluation for a molecule of N atoms
requires on the order of hundreds of Gaussian launches of the external script; accumulated
communication latency would be significant with higher-overhead transports.

This architecture enabled the full Gaussian VPT2 anharmonic framework — including its
resonance analysis, Fermi resonance treatment, and combination band identification — to
operate on ML-derived potential energy surfaces and dipole surfaces, without modification
of the Gaussian source code.

## 5. VPT2 Anharmonic Treatment

Harmonic approximation calculations predict vibrational frequencies from the second-order
term in the Taylor expansion of the potential energy surface around the equilibrium geometry.
While computationally straightforward, the harmonic approximation systematically
overestimates vibrational frequencies and cannot reproduce the overtone and combination band
features that are prominently observed in experimental mid-infrared spectra of organic
molecules.

Second-order vibrational perturbation theory (VPT2) [Mills, 1972] was applied using Gaussian
16's implementation to correct for anharmonicity. VPT2 expresses the vibrational Hamiltonian
as a perturbation series through fourth order in normal coordinates, requiring cubic
(third-order) and quartic (fourth-order) force constants in addition to the harmonic
(second-order) force constants. These higher-order force constants were obtained numerically
by finite differentiation of the energy gradient at geometries displaced from the equilibrium
along the normal mode coordinates.

The VPT2 treatment produced anharmonic fundamental frequencies, overtone bands (2ν_i), and
combination bands (ν_i + ν_j). These features were retained in the spectral analysis to
enable comparison with the full experimental spectral envelope. The B3LYP/6-31G(d,p) level
of theory was used for all DFT reference calculations, consistent with common practice for
mid-sized organic molecules where this functional/basis set combination provides a good
balance of accuracy and computational cost.

ML-based anharmonic calculations used the same VPT2 framework within Gaussian, with
ML-derived energy, forces, and dipole derivatives substituted for DFT quantities via the
`External=` interface described above. The Gaussian frequency code was thus used identically
for both DFT and ML calculations; only the source of the potential energy surface differed.

## 6. Mode Matching

Direct comparison of DFT and ML vibrational frequencies by frequency ordering assumes that
the mode with the i-th lowest frequency in the DFT calculation corresponds to the i-th
lowest frequency in the ML calculation. This assumption fails when two or more frequencies
are near-degenerate and their ordering is reversed between methods, or when the character of
a mode differs substantially between the DFT and ML potential energy surfaces.

Mode matching was performed using eigenvector dot product overlap to assign physically
corresponding modes between DFT and ML calculations [Roy et al., 2020]. For each pair of
DFT normal mode i and ML normal mode j, the overlap element O_ij = |**q**_i · **q**_j| was
computed, where **q** are mass-weighted normal-mode eigenvectors (Cartesian displacement
vectors). The overlap matrix O was computed from eigenvectors extracted from the formatted
checkpoint files (`.fchk`) produced by Gaussian's `formchk` utility.

The Hungarian algorithm [Kuhn, 1955] was applied to the overlap matrix to identify the
one-to-one assignment of DFT modes to ML modes that maximized the total eigenvector overlap.
This approach ensures that modes are matched by the physical character of the atomic motion
rather than by frequency ordering, yielding robust mode assignments even when frequency
orderings differ between methods. The resulting mode correspondence was used to pair DFT and
ML frequencies for regression analysis and metric computation.

Normal mode eigenvectors from harmonic calculations were used for all mode matching, as the
harmonic normal modes represent the true vibrational basis. Anharmonic frequency values were
then looked up for the matched modes and used in subsequent spectral comparisons.

## 7. Analysis and Validation

Predicted ML frequencies were compared against DFT reference frequencies using three
quantitative metrics: mean absolute error (MAE, cm⁻¹), root mean square error (RMSE, cm⁻¹),
and coefficient of determination (R²). These metrics were computed over the set of matched
fundamental, overtone, and combination band frequencies in the region 400–4000 cm⁻¹.

Regression plots of ML frequencies against DFT frequencies were generated for each of the
four ML model combinations to provide visual assessment of systematic offsets, frequency-
dependent biases, or outlier modes. The slope and intercept of the linear regression line
were reported alongside R² to characterize the nature of any systematic error.

Simulated IR spectra were generated by convolving the discrete frequency-intensity stick
spectrum with a Gaussian kernel of full-width at half-maximum (FWHM) 8.0 cm⁻¹ to represent
instrumental broadening. Spectral comparison between ML predictions and DFT reference was
performed over the standard mid-infrared range 400–4000 cm⁻¹.

All four ML combinations were evaluated independently to isolate the contributions of the
potential energy model (MACE-MP-0 vs. MACE-OMOL-0) and the dipole model (Espaloma vs.
MACE4IR) to the overall spectral accuracy. The DFT reference method (B3LYP/6-31G(d,p)) was
treated as the computational ground truth for frequency comparison; assessment of absolute
accuracy against experimental measurements was not the primary scope of this study, though
the methodology is directly applicable to such comparisons.

## References

> Full citations to be formatted per journal requirements.

- **MACE architecture:** Batatia, I., Kovács, D. P., Simm, G. N. C., Ortner, C., & Csányi, G. (2022). MACE: Higher order equivariant message passing neural networks for fast and accurate force fields. *Advances in Neural Information Processing Systems*, 35. [Author, 2022]

- **MACE-MP-0:** Batatia, I., et al. (2023). A foundation model for atomistic materials chemistry. *arXiv:2401.00096*. [Batatia et al., 2023]

- **MACE-OMOL-0 / OMol25:** Open Molecules 2025 dataset and model. Citation details to be confirmed from publication. [OMol25 Authors, 2025]

- **Espaloma:** Wang, Y., Fass, J., Stern, C. D., Luo, K., & Chodera, J. D. (2022). Espaloma-0.3.0: Machine-learning coarse-graining for molecular mechanics force fields. *Journal of Chemical Theory and Computation*. [Wang et al., 2022]

- **Gaussian 16:** Frisch, M. J., et al. (2016). *Gaussian 16, Revision C.01*. Gaussian, Inc., Wallingford CT. [Frisch et al., 2016]

- **VPT2 theory:** Mills, I. M. (1972). Vibration-rotation structure in asymmetric and symmetric top molecules. In *Molecular Spectroscopy: Modern Research* (Vol. 1, pp. 115–140). Academic Press. [Mills, 1972]

- **ASE:** Larsen, A. H., et al. (2017). The atomic simulation environment — a Python library for working with atoms. *Journal of Physics: Condensed Matter*, 29(27), 273002. [Larsen et al., 2017]
