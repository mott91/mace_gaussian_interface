---
created: 2026-03-10T16:05:00.000Z
title: Uncertainty quantification for ML-predicted spectra
area: research
files: []
---

## Problem

ML predictions currently give point estimates with no error bars. There's no way to know which predicted peaks are reliable and which might be artifacts. This is a gap in the literature — nobody reports uncertainty on ML IR spectra yet.

## Solution

1. Investigate ensemble approaches: run multiple MACE models (different seeds/architectures) and measure variance across predictions.
2. Alternatively, explore MC dropout or evidential deep learning if supported by MACE.
3. Propagate uncertainty through the frequency calculation to get error bars on each predicted frequency and intensity.
4. Visualize as shaded confidence bands on the spectrum plot.
5. Quantify: do uncertainty estimates correlate with actual DFT-ML error? (calibration check)
6. Could be a strong thesis contribution if the calibration is reasonable.
