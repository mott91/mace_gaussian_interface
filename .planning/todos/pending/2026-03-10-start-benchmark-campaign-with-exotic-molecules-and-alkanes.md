---
created: 2026-03-10T14:14:33.444Z
title: Start benchmark campaign with exotic molecules and alkanes
area: general
files: []
---

## Problem

Current benchmarking is limited to a small set of molecules (water, BH3-NH3, acetic acid, etc.). A systematic benchmark campaign is needed to evaluate ML model performance across a broader chemical space — including exotic functional groups and a homologous alkane series (methane → longer chains) to probe transferability and size scaling.

## Solution

1. Define the molecule set: select exotic molecules (varied functional groups, heteroatoms) and alkane series (C1–C8 or similar).
2. Gather or generate .xyz structures for all targets.
3. Run the full pipeline for each (requires batch pipeline todo first, or run sequentially).
4. Collect and compare: frequency errors, intensity errors, mode matching quality vs. DFT baseline.
5. Summarize findings for thesis chapter.
