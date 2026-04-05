---
created: 2026-03-26T20:00:00.000Z
title: Handle degenerate modes and zero-intensity filtering in analysis
area: analysis
files:
  - mode_matching.py
  - mace_gaussian/analysis/analyze_spectra.py
---

## Problem

1. Degenerate modes (e.g., methane T₂ 3-fold) have arbitrary eigenvector orientation within the degenerate subspace. Individual dot products between DFT and ML eigenvectors give misleading low overlaps even when the subspaces are identical.

2. IR-inactive modes (zero DFT intensity) pollute intensity regression metrics. Methane shows R²=0.34 for intensity because we're correlating zeros with numerical noise.

## Solution

### Degenerate mode subspace overlap
1. Detect degenerate groups: modes within ~5 cm⁻¹ in the DFT reference.
2. For degenerate groups, compute subspace overlap (trace of M^T M where M is the inter-group overlap matrix) instead of individual eigenvector dot products.
3. Report per-group overlap as a single confidence score.
4. Hungarian matching should still work on individual modes but confidence reporting should be group-aware.

### Zero-intensity filtering
1. Match ALL modes first (Hungarian on eigenvectors — intensity-agnostic).
2. For frequency regression: include all matched pairs.
3. For intensity regression: exclude pairs where DFT intensity < threshold (default 0.1 km/mol).
4. Add --min-dft-intensity CLI flag to control threshold.
5. Report how many modes were filtered and why.
