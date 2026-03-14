---
created: 2026-03-10T16:05:00.000Z
title: Delta-ML correction model for frequency predictions
area: research
files: []
---

## Problem

Pure ML frequencies have systematic errors compared to DFT. Running full DFT on every molecule is expensive. There may be a sweet spot where a small amount of DFT data can correct ML predictions cheaply.

## Solution

1. For molecules where both DFT and ML results exist, compute residuals (DFT - ML) per mode.
2. Train a lightweight model (linear regression, small NN, or GP) on these residuals using molecular/mode features.
3. For new molecules, predict ML frequencies and apply the learned correction.
4. Evaluate: how many DFT reference molecules are needed before the correction generalizes?
5. Compare corrected ML vs pure ML vs pure DFT accuracy on held-out molecules.
