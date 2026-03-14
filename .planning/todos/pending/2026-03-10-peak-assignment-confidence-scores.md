---
created: 2026-03-10T16:05:00.000Z
title: Peak assignment confidence scores from mode matching
area: analysis
files:
  - mode_matching.py
---

## Problem

Mode matching reports which ML mode corresponds to which DFT mode, but doesn't convey how reliable the assignment is. A near-degenerate pair with 0.51 vs 0.49 overlap is very different from a clean 0.98 match, but both are reported the same way.

## Solution

1. Propagate the eigenvector dot product magnitude as a confidence score for each mode assignment.
2. Flag low-confidence assignments (e.g., overlap < 0.7) in the output.
3. Report the second-best match and its overlap for ambiguous cases.
4. Include confidence scores in regression plots (e.g., color-code points by confidence) and analysis JSON.
