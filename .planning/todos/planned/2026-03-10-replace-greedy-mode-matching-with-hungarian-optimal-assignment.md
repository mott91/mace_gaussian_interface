---
created: 2026-03-10T15:17:30.445Z
title: Replace greedy mode matching with Hungarian optimal assignment
area: general
files:
  - mace_gaussian/analysis/mode_matching.py:161-214
---

## Problem

`match_modes()` in `mode_matching.py` is greedy: for each ML mode it independently finds the highest-overlap DFT mode. Nothing prevents two ML modes from being matched to the same DFT mode. For a 10-mode molecule this can produce duplicate assignments and leave some DFT modes unmatched.

The correct approach is the Hungarian algorithm (linear sum assignment), which finds the globally optimal bijective matching — each ML mode maps to exactly one DFT mode and vice versa, maximising total overlap across all pairs simultaneously.

## Solution

Replace the greedy nested loop in `match_modes()` with `scipy.optimize.linear_sum_assignment`:

```python
from scipy.optimize import linear_sum_assignment
overlap_matrix = create_alignment_matrix(modes_calc, modes_ref)
row_ind, col_ind = linear_sum_assignment(-overlap_matrix)  # minimise negative overlap
```

`scipy` is already a dependency. The change is ~5 lines. Keep the threshold check afterwards to flag low-confidence assignments.

Note: do this alongside or after the "wire eigenvector mode matching" todo, since both touch the same code path.
