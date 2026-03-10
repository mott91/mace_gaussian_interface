---
created: 2026-03-10T14:14:33.438Z
title: Finish and test anharmonic analysis workflow
area: general
files:
  - run_analysis.py
  - mode_matching.py
---

## Problem

The anharmonic analysis workflow (`run_analysis.py`) is not fully finished or robustly tested. Known issues include the acetic acid (acoh) bug with frequency matching for regression plots (commit a4384c4), and potential edge cases in overtone/combination band parsing. The workflow needs to be brought to a reliable, tested state across all tracked molecules.

## Solution

1. Audit `run_analysis.py` end-to-end for incomplete logic or fragile parsing.
2. Fix the acoh DFT regression/frequency matching bug.
3. Add tests or at minimum run the full workflow against all molecules in `analysis_results/` and verify outputs.
4. Document any remaining known limitations.
