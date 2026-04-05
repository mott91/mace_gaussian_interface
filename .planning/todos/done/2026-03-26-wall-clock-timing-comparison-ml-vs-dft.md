---
created: 2026-03-26T18:30:00.000Z
title: Wall-clock timing comparison ML pipeline vs full DFT
area: analysis
files:
  - mace_gaussian/workflow.py
---

## Problem

No systematic measurement of computational speedup from using ML dipoles vs full DFT for the frequency calculation. "We get X% of DFT accuracy in Y% of the time" is a key result for the thesis but we don't currently capture timing data.

## Solution

1. Add timing instrumentation to the pipeline stages (geometry opt, frequency calc, analysis).
2. Record wall-clock time per molecule per calculator combination.
3. Compare against DFT baseline wall-clock time for the same molecules.
4. Include in batch report: speedup factor, cost-accuracy tradeoff plot.
5. Break down by molecule size to show scaling behavior (does ML advantage grow with size?).
6. Report in thesis as a table and/or plot: accuracy (R², RMSE) vs computational cost.
