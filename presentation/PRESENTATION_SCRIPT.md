# Presentation Script: ML-Accelerated IR Spectroscopy
## ~10 minutes | Manuel Ott | Hofer Lab

---

## Slide 1: Title Slide [30 seconds]

Good morning everyone. Today I'm going to talk about a project I've been working on that combines machine learning with quantum chemistry to dramatically speed up infrared spectroscopy calculations.

The core idea is simple: what if we could use fast machine learning models to replace the slowest parts of quantum chemistry calculations, while maintaining the same level of accuracy?

---

## Slide 2: Motivation [1.5 minutes]

So let me start with the problem. If you've ever run DFT frequency calculations, you know they're painfully slow. We're talking minutes to hours per molecule. For a single water molecule, you might wait 5 minutes. For anything larger, you're looking at hours or even days. This makes high-throughput screening basically impossible.

The bottleneck is calculating the Hessian—the matrix of energy second derivatives—and dipole derivatives. In DFT, both of these require expensive calculations at many displaced geometries.

Our solution is a hybrid approach that uses machine learning for the expensive parts. MACE calculates energies, forces, and the Hessian using automatic differentiation—this is fast because MACE can compute second derivatives analytically rather than numerically. We also use ML dipole calculators—either our custom MACE-dipole model or Espaloma—to predict dipole moments and derivatives. All of these ML quantities get injected into Gaussian through a ZMQ bridge.

So what does Gaussian actually do? Gaussian receives the ML-calculated Hessian and dipole derivatives and uses them to perform anharmonic frequency analysis—calculating overtones, combination bands, and anharmonic corrections. This is where Gaussian's quantum chemistry expertise really shines.

This hybrid approach gives us potentially 10 to 100 times speedup compared to pure DFT. The speedup comes from MACE's fast Hessian calculation via automatic differentiation, and from avoiding expensive DFT dipole calculations.

---

## Slide 3: System Architecture [2 minutes]

Here's how the system works end-to-end.

Starting at the top, we take a molecule XYZ file and run geometry optimization using MACE-OMOL. This is fast—usually under a minute.

Then the workflow splits into two parallel paths. On the left, we have the DFT baseline calculation. This uses B3LYP with a 6-31G(d,p) basis set and runs pure Gaussian—no machine learning at all. Gaussian does its own geometry optimization, calculates the Hessian with DFT, calculates dipole derivatives with DFT, and performs anharmonic analysis. This gives us our ground truth reference—everything is pure quantum chemistry.

On the right side is where the innovation happens—the ML frequency calculation path. Here's exactly how the ZMQ bridge works: We set up a ZMQ server in Python that waits for requests from Gaussian. When Gaussian launches with the "external" keyword, it calls our Python helper script through the ZMQ bridge.

For each request, our Python code receives the molecular geometry from Gaussian, then calculates four things with machine learning: energy with MACE, forces with MACE, the Hessian matrix with MACE using automatic differentiation, and dipole moments with our ML dipole calculator—either MACE-dipole or Espaloma. We also calculate dipole derivatives using finite differences. All of these ML-calculated quantities get packaged and sent back to Gaussian in its expected format.

Gaussian receives this data and uses it to perform anharmonic frequency analysis—calculating overtones, combination bands, and anharmonic corrections. So Gaussian isn't recalculating anything; it's using the ML Hessian and dipoles we provided to do the higher-level spectroscopic analysis.

Both paths converge at the analysis stage. Here we load all the results, automatically find the DFT baseline, match vibrational modes using eigenvector similarity—which I'll explain in a moment—broaden the spectra, calculate statistical metrics, and plot regressions. Finally, everything gets packaged into an HTML report.

---

## Slide 4: Modular Dipole System [45 seconds]

The dipole calculator system is fully modular. We have two main options: MACE-dipole, which is our custom-trained model, and Espaloma, which calculates dipoles from machine learning predicted partial charges.

The implementation uses a factory pattern with automatic availability checks. If one calculator isn't available, the system automatically falls back to the next option. And crucially, each dipole calculator gets automatically combined with all available energy calculators, so we can test multiple method combinations in a single run.

---

## Slide 5: Dipole Methods [1.5 minutes]

Let me explain how these two dipole calculation methods work, because the details matter.

MACE-Dipole uses an E(3)-equivariant neural network. E(3)-equivariant means the model respects rotational and translational symmetry—if you rotate the molecule, the predicted dipole vector rotates exactly the same way. This is built into the architecture through equivariant message passing layers, which are basically graph neural networks where the node features transform correctly under rotations. The input is atomic coordinates and numbers, it processes them through these layers, and directly outputs the dipole vector: mu-x, mu-y, mu-z. It's fast—about 0.1 seconds per geometry—and physically consistent, but requires a custom-trained model.

Espaloma takes a charge-based approach. It uses a graph neural network to predict partial charges for each atom, then we calculate the dipole as the sum of charge times position for all atoms. It's widely available and interpretable since you see actual charges.

Now here's an important detail: for both methods, we need dipole DERIVATIVES—how the dipole changes when atoms move—for IR intensities. We calculate these in Python using central finite differences. We displace each atom slightly forward and backward, calculate the dipole at each displaced geometry using our ML model, and approximate the derivative as the difference divided by twice the displacement. Central differences are more accurate than forward differences. Once calculated, these dipole derivatives get sent to Gaussian along with the Hessian.

---

## Slide 6: Workflow Management [1 minute]

One thing I'm really proud of is the automated workflow orchestration.

The system automatically creates a hierarchical directory structure organized by molecule name and calculator combination. Everything is created on-the-fly—no manual setup required.

It also has smart file detection. If you run a calculation, stop it halfway, and restart later, it automatically detects what's already been done and skips completed steps. It can even find the DFT baseline automatically when doing comparisons.

For results aggregation, the system loads all JSON results, converts checkpoint files to formatted checkpoint files, parses all the log and formatted checkpoint files—which can be massive—matches modes across different methods, and generates a unified HTML report with all the plots and statistics in one place.

---

## Slide 7: Parsing [1 minute]

Speaking of parsing, this is actually non-trivial. Gaussian log files are 2 to 10 megabytes of unstructured text with 10,000 to 50,000 lines.

We need to extract harmonic frequencies, anharmonic frequencies, overtones and combination bands, IR intensities, and eigenvectors for 3N minus 6 vibrational modes. The challenge is that Gaussian uses a multi-column format with three modes per block, data is scattered throughout the file, and you need complex regex patterns to find everything.

In contrast, our ML calculations output clean JSON with frequencies, intensities, and eigenvectors in a nicely structured format. Much easier to work with.

---

## Slide 8: Mode Matching [1.5 minutes]

Now here's a subtle but critical problem: frequency order mismatch.

When you run DFT, you might get frequencies like 500, 1200, 1500, 2900, 3100 wavenumbers. Run the same molecule with ML and you get 498, 3095, 1205, 2895, 1498—completely different order! You can't just compare mode 1 to mode 1 by index.

The solution is eigenvector dot product matching. Let me explain the math behind why this works. Each vibrational mode has an eigenvector—a 3N-dimensional vector where N is the number of atoms. This vector describes exactly how each atom moves in that vibration: x, y, z displacements for atom 1, then atom 2, and so on.

If two modes represent the same physical vibration—like a symmetric stretch—their eigenvectors should point in the same direction in this 3N-dimensional space. The dot product of two normalized vectors equals the cosine of the angle between them. If they're parallel, the angle is zero, cosine is one, dot product is one. If they're orthogonal—completely different motions—the angle is 90 degrees, cosine is zero, dot product is zero.

We use the absolute value because eigenvectors can be defined with opposite signs but still represent the same motion—it's just a phase convention. So we calculate the absolute dot product between each ML eigenvector and all DFT eigenvectors, find the maximum similarity, and if it's above our threshold—typically 0.8—we call it a match.

For example, ML mode at 498 matches DFT mode at 500 with a similarity of 0.996—almost perfectly parallel vectors, so definitely the same vibration. ML mode at 3095 matches DFT mode at 3100 with 0.989 similarity. These high values give us confidence we've matched the right modes.

---

## Slide 9: Results [30 seconds]

All of this analysis gets packaged into an interactive HTML report. Here's an example for water. You get regression plots showing ML versus DFT frequencies, statistical metrics, broadened spectra, and everything you need to evaluate the method's performance—all in one self-contained file you can open in any browser.

---

## Slide 10: Current Status & Next Steps [1 minute]

So where are we now?

Current status: The ZMQ bridge is fully implemented and working. We're in the middle of validation on small molecule benchmarks. The automated workflow orchestration is solid and reliable. And importantly, we can compare multiple methods side-by-side in one place.

Next steps: First, expand to larger molecules to see how well this scales. Second, implement anharmonic frequency analysis—right now we're focused on harmonic frequencies. Third, build out a proper command line interface to make it easier to use. And fourth, add ORCA support as an open-source alternative to Gaussian.

---

## Slide 11: Questions [30 seconds]

That's it! I'm happy to take any questions.

[Pause for questions]

Thank you!

---

## Total Time: ~10 minutes

**Tips for delivery:**
- Maintain eye contact, don't just read the slides
- Use the workflow diagram (slide 3) as a visual anchor—point to it when explaining the flow
- Emphasize the speed gains (10-100×) and automation aspects
- Be ready to elaborate on technical details if asked
- Have example runtimes ready (e.g., "water takes 5 minutes with DFT, 5 seconds with ML")
