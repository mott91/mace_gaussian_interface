# Presentation Script: ML-Accelerated IR Spectroscopy
## ~10 minutes | Manuel Ott | Hofer Lab

---

## Slide 1: Title Slide [30 seconds]

Good morning everyone. Today I'm going to talk about a project I've been working on that combines machine learning with quantum chemistry to dramatically speed up infrared spectroscopy calculations.

The core idea is simple: what if we could use fast machine learning models to replace the slowest parts of quantum chemistry calculations, while maintaining the same level of accuracy?

---

## Slide 2: Motivation [1 minute]

So let me start with the problem. If you've ever run DFT frequency calculations, you know they're painfully slow. We're talking minutes to hours per molecule. For a single water molecule, you might wait 5 minutes. For anything larger, you're looking at hours or even days. This makes high-throughput screening basically impossible.

Our solution is a hybrid approach. Instead of doing everything with DFT, we use machine learning to calculate dipole derivatives—the most expensive part—and inject those into Gaussian through a ZMQ bridge. This gives us potentially 10 to 100 times speedup compared to pure DFT.

The impact? We can now run rapid IR spectral predictions on hundreds or thousands of molecules, which was simply not feasible before, while maintaining DFT-level accuracy.

---

## Slide 3: System Architecture [1.5 minutes]

Here's how the system works end-to-end.

Starting at the top, we take a molecule XYZ file and run geometry optimization using MACE-OMOL. This is fast—usually under a minute.

Then the workflow splits into two parallel paths. On the left, we have the DFT baseline calculation. This uses B3LYP with a 6-31G(d,p) basis set, and importantly, it does its own built-in geometry optimization, then calculates frequencies and dipoles entirely with Gaussian. This gives us our ground truth reference.

On the right side is where the magic happens—the ML frequency calculation path. We load MACE calculators for energies and our custom dipole calculators, set up a ZMQ server, launch Gaussian, and enter a dipole loop where we feed machine learning predictions into the quantum chemistry calculation in real-time.

Both paths converge at the analysis stage at the bottom. Here we load all the results, find the DFT baseline automatically, match vibrational modes using eigenvector similarity, broaden the spectra, calculate statistical metrics, and plot regressions. Finally, everything gets packaged into an HTML report.

---

## Slide 4: Modular Dipole System [45 seconds]

The dipole calculator system is fully modular. We have two main options: MACE-dipole, which is our custom-trained model, and Espaloma, which calculates dipoles from machine learning predicted partial charges.

The implementation uses a factory pattern with automatic availability checks. If one calculator isn't available, the system automatically falls back to the next option. And crucially, each dipole calculator gets automatically combined with all available energy calculators, so we can test multiple method combinations in a single run.

---

## Slide 5: Dipole Methods [1 minute]

Let me quickly explain how these two dipole calculation methods work.

MACE-Dipole uses an E(3)-equivariant neural network that directly predicts the dipole moment vector. The input is atomic coordinates and atomic numbers, it processes them through equivariant message passing layers, and outputs mu-x, mu-y, mu-z. It's fast—about 0.1 seconds per geometry—but requires a custom-trained model.

Espaloma takes a different, charge-based approach. It first predicts partial charges q1, q2, through qn for each atom, then calculates the dipole moment as the sum of charge times position for each atom. It's widely available and interpretable since you can see the actual charges, but it's an indirect method that depends on how you assign those charges.

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

## Slide 8: Mode Matching [1 minute]

Now here's a subtle but critical problem: frequency order mismatch.

When you run DFT, you might get frequencies like 500, 1200, 1500, 2900, 3100 wavenumbers. Run the same molecule with ML and you get 498, 3095, 1205, 2895, 1498—completely different order! You can't just compare mode 1 to mode 1 by index.

The solution is eigenvector dot product matching. Each vibrational mode has an eigenvector that describes how atoms move. We calculate the similarity between ML and DFT eigenvectors using the absolute value of their dot product. This lets us compare atomic displacement patterns, find the best match for each mode, and it's robust to the frequency ordering.

For example, ML mode at 498 matches DFT mode at 500 with a similarity of 0.996. ML mode at 3095 matches DFT mode at 3100 with 0.989 similarity. Very high confidence matches.

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
