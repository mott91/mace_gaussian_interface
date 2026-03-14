---
created: 2026-03-10T16:05:00.000Z
title: Active learning loop for dipole model improvement
area: research
files: []
---

## Problem

The dipole model has finite training data and may perform poorly on certain chemistries. Running DFT on every molecule to check is expensive. An active learning approach could identify the most informative molecules to compute DFT references for, maximizing model improvement per DFT dollar spent.

## Solution

1. Use uncertainty estimates (see UQ todo) or prediction disagreement to rank candidate molecules by expected informativeness.
2. Select the top-N most uncertain molecules from a candidate pool.
3. Run DFT on those molecules to get ground-truth dipole derivatives.
4. Retrain or fine-tune the dipole model with the new data.
5. Iterate: re-rank remaining candidates, select next batch, repeat.
6. Evaluate: learning curve (accuracy vs number of DFT calculations) compared to random selection baseline.
