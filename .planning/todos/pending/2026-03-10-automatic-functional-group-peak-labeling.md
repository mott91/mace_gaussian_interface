---
created: 2026-03-10T16:05:00.000Z
title: Automatic functional group peak labeling
area: analysis
files:
  - run_analysis.py
---

## Problem

Predicted peaks are reported as raw frequencies without chemical interpretation. Manually assigning peaks to functional groups (C=O stretch, O-H bend, etc.) is tedious and error-prone, especially for larger molecules.

## Solution

1. Build a lookup table of characteristic IR frequency ranges for common functional groups.
2. Cross-reference predicted frequencies with the lookup table.
3. Use eigenvector analysis (which atoms move) to disambiguate overlapping ranges.
4. Annotate peaks on spectrum plots with functional group labels.
5. Include assignments in the analysis JSON output and any exported formats.
