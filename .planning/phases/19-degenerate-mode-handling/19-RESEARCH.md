# Phase 19: Degenerate Mode Handling - Research

**Researched:** 2026-03-30
**Domain:** Vibrational mode matching, linear algebra (subspace overlap), molecular symmetry
**Confidence:** HIGH

## Summary

This phase adds degenerate mode awareness to the existing Hungarian-based mode matching pipeline. The core mathematical operation -- subspace overlap via trace(M^T M)/k -- is straightforward linear algebra using numpy operations already available in the codebase. The main complexity is structural: threading degenerate group information through multiple layers (mode matching -> statistics -> comparison tables -> regression plots -> heatmaps -> HTML reports).

The codebase already extracts all degenerate modes from .fchk files (`extract_spectrum_from_fchk`), and the Hungarian matching works on individual modes. The architecture decision (D-04) to keep Hungarian unchanged and add subspace overlap as a post-processing step is sound and minimizes risk. The Gaussian log parser's `seen_freqs` deduplication (parser.py line 47-79) is the one place where degenerate modes are currently lost -- but this only affects log-based extraction, not the .fchk path used for harmonic mode matching.

**Primary recommendation:** Build a `DegenerateGroupDetector` class that takes DFT reference frequencies, identifies groups within 5 cm-1, and produces a `DegenerateGroupResult` that wraps the existing `match_modes()` output with group-level overlap scores. All downstream consumers (statistics, plots, reports) read from this result object.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Degenerate groups are defined by DFT reference frequencies only (not ML side). ML modes get assigned to DFT-defined groups via Hungarian matching.
- **D-02:** Frequency threshold for grouping: fixed at 5 cm-1. No CLI flag.
- **D-03:** Parser must stop deduplicating degenerate modes. The .fchk extraction already gets all degenerate modes; the Gaussian log parser's `seen_freqs` deduplication must be bypassed for harmonic mode matching.
- **D-04:** Two-stage approach: run Hungarian on individual modes first (unchanged), then post-process degenerate groups with subspace overlap (trace(M^T M)/k).
- **D-05:** For degenerate groups, subspace overlap REPLACES individual dot product overlaps. Non-degenerate modes keep individual scores.
- **D-06:** One entry per degenerate group in statistics. T2 counts as 1 unit, not 3.
- **D-07:** Regression: each degenerate group contributes 1 data point (averaged freq/intensity).
- **D-08:** Report labels degenerate groups with type and multiplicity (e.g., "T2 (3-fold) at 1356 cm-1").
- **D-09:** Heatmaps collapse degenerate groups into single rows/columns with subspace overlap value.
- **D-10:** Regression plots use distinct marker (diamond/square) for degenerate group averages.

### Claude's Discretion
- Specific marker shapes and colors for degenerate vs non-degenerate in plots
- Fallback labeling when symmetry type can't be auto-determined
- Internal data structure for representing degenerate groups in match results
- Whether to add a summary line in the report like "3 degenerate groups detected (2 doubly, 1 triply)"

### Deferred Ideas (OUT OF SCOPE)
- "Peak assignment confidence scores from mode matching" -- separate enhancement layer
- "Investigate and use compute_reduced_masses in mode matching" -- orthogonal concern
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| MODE-05 | Degenerate modes (within 5 cm-1 threshold) are detected and grouped; subspace overlap (trace of M^T M / k) is used for matching quality | Subspace overlap math documented below; detection algorithm is frequency-proximity clustering on DFT reference; existing `create_alignment_matrix()` builds the sub-matrices needed |
| MODE-06 | Mode matching statistics correctly handle degenerate groups without double-counting | Group-aware statistics computed from `DegenerateGroupResult`; one entry per group in regression data; match count = non-degenerate + groups |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- **Code quality:** `ruff check --fix && ruff format` (line length 100), `ty check`
- **Plots:** PNG format, DPI=300, seaborn "colorblind" palette, Arial/Helvetica font. Always include R-squared and RMSE in regression plots. Standard range 400-4000 cm-1.
- **Heatmaps:** Always use `force_harmonic=True` for mode overlap heatmaps.
- **Environment:** Use `micromamba activate mace4ir_v2` (not the `.venv`).
- **Test runner:** `pytest tests/ --strict-markers -ra`

## Architecture Patterns

### Data Flow for Degenerate Mode Handling

```
                     existing (unchanged)                    new (this phase)
                     ====================                    ================

.fchk files -> extract_mode_data_from_checkpoint()
            -> match_modes() (Hungarian, individual)
            -> matches: dict[int, tuple[int|None, float]]
                                                        -> detect_degenerate_groups(dft_freqs, threshold=5.0)
                                                        -> groups: list[DegenerateGroup]
                                                        -> compute_group_overlaps(matches, modes_calc, modes_ref, groups)
                                                        -> DegenerateGroupResult
                                                             .individual_matches  (original)
                                                             .groups              (list of groups)
                                                             .effective_overlaps  (group overlap replaces individual)
                                                             .statistics()        (group-aware counts)
```

### Recommended Data Structures

```python
@dataclass
class DegenerateGroup:
    """A group of nearly-degenerate modes in the DFT reference."""
    ref_indices: list[int]          # DFT mode indices in this group
    calc_indices: list[int]         # Matched ML mode indices (from Hungarian)
    center_frequency: float         # Mean DFT frequency of group
    multiplicity: int               # len(ref_indices), e.g. 2 for E, 3 for T
    symmetry_label: str             # "E", "T2", "A" (best-effort), or "deg-2", "deg-3" fallback
    subspace_overlap: float         # trace(M^T M) / k


@dataclass
class DegenerateGroupResult:
    """Wraps Hungarian match results with degenerate group awareness."""
    individual_matches: dict[int, tuple[int | None, float]]  # Original from match_modes()
    groups: list[DegenerateGroup]
    non_degenerate_indices: list[int]  # DFT ref indices not in any group

    def effective_overlap(self, ref_idx: int) -> float:
        """Get overlap for a ref mode: group subspace overlap if degenerate, individual otherwise."""
        ...

    def statistics(self) -> dict:
        """Group-aware match statistics: total_units, avg_overlap, confident_count."""
        ...

    def regression_data(self, freqs_calc, freqs_ref, ints_calc, ints_ref) -> tuple:
        """One data point per non-degenerate mode + one per group (averaged)."""
        ...
```

### Recommended Module Structure

```
mace_gaussian/analysis/
├── mode_matching.py          # Existing: match_modes(), create_alignment_matrix(), heatmap
│                             # Add: detect_degenerate_groups(), compute_subspace_overlap()
│                             #      DegenerateGroup, DegenerateGroupResult dataclasses
├── analyze_spectra.py        # Modify: plot_regression() to use diamond markers for groups
│                             # Modify: match_by_mode() to use group-aware data
├── analysis_workflow.py      # Modify: extract_mode_mapping() to return DegenerateGroupResult
│                             #          generate_mode_overlap_heatmaps() to collapse groups
└── html_report_generator.py  # Modify: mode matching section to show group labels
```

Keep all new code in `mode_matching.py` -- it's the natural home for mode comparison logic. The other files consume results from it.

### Anti-Patterns to Avoid
- **Modifying Hungarian matching itself:** D-04 explicitly says Hungarian stays unchanged. Subspace overlap is post-processing only.
- **Grouping on ML side:** D-01 says groups are defined by DFT frequencies only. ML modes are assigned to DFT groups via the existing Hungarian matches.
- **Using parser.py for harmonic mode matching:** The .fchk path already extracts all degenerate modes. The parser deduplication fix (D-03) is for consistency/correctness, but mode matching already uses the .fchk path.

## Standard Stack

### Core (already in project)
| Library | Purpose | Why |
|---------|---------|-----|
| numpy | Subspace overlap linear algebra (trace, matrix multiply, SVD if needed) | Already used throughout mode_matching.py |
| scipy.optimize.linear_sum_assignment | Hungarian matching (unchanged) | Already imported |
| matplotlib | Heatmaps and regression plots | Already used |
| dataclasses | DegenerateGroup, DegenerateGroupResult | Python stdlib |

No new dependencies needed. All math is standard numpy operations.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Frequency clustering | Custom clustering algorithm | Simple sorted-frequency sweep with 5 cm-1 gap threshold | Modes are 1D (frequency), sorted. A single pass identifies groups. No need for scipy clustering. |
| Subspace overlap | Custom SVD-based overlap | `np.trace(M.T @ M) / k` where M is the sub-block of the alignment matrix | Direct formula from the requirement. 3 numpy operations. |
| Symmetry classification | Point group analysis library | Multiplicity-based heuristic: 2-fold="E", 3-fold="T", else="deg-{k}" | Full point group analysis is massive overkill. Multiplicity gives the useful information. |

## Common Pitfalls

### Pitfall 1: Alignment Matrix Sub-block Extraction
**What goes wrong:** Building the M matrix for subspace overlap from the wrong indices. The alignment matrix from `create_alignment_matrix()` is indexed [calc_idx, ref_idx]. Degenerate group ref_indices map to columns, matched calc_indices map to rows.
**Why it happens:** Easy to swap row/column indexing.
**How to avoid:** Extract sub-block as `alignment_matrix[np.ix_(calc_indices, ref_indices)]`. Verify shape is (k, k) for a k-fold degenerate group.
**Warning signs:** Subspace overlap > 1.0 or negative values.

### Pitfall 2: Incomplete Hungarian Match Coverage
**What goes wrong:** A degenerate group of 3 DFT modes might only have 2 ML modes matched to it (if n_calc < n_ref or one ML mode matched elsewhere).
**Why it happens:** Hungarian matching is global -- it might assign an ML mode to a non-degenerate DFT mode instead of keeping it in the group.
**How to avoid:** When building calc_indices for a group, only include ML modes that were actually matched to DFT modes in that group. If k_matched < k_group, the subspace overlap formula should use k_matched (the actual number of ML modes assigned to the group), not k_group. Document this edge case.
**Warning signs:** `len(calc_indices) != len(ref_indices)` for a group.

### Pitfall 3: Parser Deduplication Fix Breaking Anharmonic Parsing
**What goes wrong:** Removing `seen_freqs` deduplication from `parse_harmonic_frequencies()` could cause it to return duplicate entries from repeated frequency blocks in the Gaussian log (Gaussian prints frequencies both in the initial analysis and before anharmonic calculation).
**Why it happens:** The `seen_freqs` filter serves two purposes: (1) removing true duplicates from repeated blocks, and (2) accidentally collapsing degenerate modes.
**How to avoid:** Instead of removing `seen_freqs` entirely, change the deduplication key. Currently it uses `(round(freq, 4), round(intensity, 4))`. Degenerate modes have identical (freq, intensity) pairs. Solution: add a counter or positional index so each occurrence within a single frequency block is unique, while still deduplicating across repeated blocks.
**Warning signs:** Water (3 modes) returning 6 entries (doubled from repeated blocks).

### Pitfall 4: Regression Plot Point Count Mismatch
**What goes wrong:** After group aggregation, the number of data points changes. If `confident_mask` logic in `plot_regression()` still indexes by original mode count, array shapes won't match.
**Why it happens:** The confident mask is built per-mode, but group aggregation reduces the point count.
**How to avoid:** Build the confident mask from the group-aware result, not from individual mode overlaps. Each group contributes one entry to the mask (confident if subspace overlap >= 0.5).
**Warning signs:** IndexError or misaligned scatter points in regression plots.

### Pitfall 5: Heatmap Collapsing Dimension Mismatch
**What goes wrong:** Collapsing a 9x9 alignment matrix to ~4x4 for CH4 requires careful aggregation. The collapsed cell value should be the subspace overlap for the group, not a simple average of the individual cells.
**Why it happens:** Tempting to average the sub-block values, but that's mathematically wrong -- subspace overlap (trace(M^T M)/k) is the correct measure.
**How to avoid:** For the collapsed heatmap, non-degenerate cells keep their individual overlap. Group cells get the computed subspace overlap value. Build a new collapsed matrix, don't try to modify the original in-place.

## Code Examples

### Degenerate Group Detection
```python
# Source: custom implementation based on D-01, D-02
def detect_degenerate_groups(
    ref_frequencies: np.ndarray, threshold: float = 5.0
) -> list[DegenerateGroup]:
    """Identify groups of nearly-degenerate modes in DFT reference frequencies.

    Modes are grouped if consecutive frequencies differ by <= threshold.
    Groups of size 1 are non-degenerate and excluded from the result.
    """
    sorted_indices = np.argsort(ref_frequencies)
    sorted_freqs = ref_frequencies[sorted_indices]

    groups = []
    current_group = [sorted_indices[0]]

    for i in range(1, len(sorted_freqs)):
        if sorted_freqs[i] - sorted_freqs[i - 1] <= threshold:
            current_group.append(sorted_indices[i])
        else:
            if len(current_group) >= 2:
                _add_group(groups, current_group, ref_frequencies)
            current_group = [sorted_indices[i]]

    if len(current_group) >= 2:
        _add_group(groups, current_group, ref_frequencies)

    return groups


def _add_group(groups, indices, freqs):
    k = len(indices)
    label = {2: "E", 3: "T"}.get(k, f"deg-{k}")
    groups.append(DegenerateGroup(
        ref_indices=indices,
        calc_indices=[],  # Filled after Hungarian matching
        center_frequency=float(np.mean(freqs[indices])),
        multiplicity=k,
        symmetry_label=label,
        subspace_overlap=0.0,  # Computed later
    ))
```

### Subspace Overlap Computation
```python
# Source: trace(M^T M) / k formula from MODE-05
def compute_subspace_overlap(
    alignment_matrix: np.ndarray,
    calc_indices: list[int],
    ref_indices: list[int],
) -> float:
    """Compute subspace overlap for a degenerate group.

    Parameters
    ----------
    alignment_matrix : shape (n_calc, n_ref), full overlap matrix
    calc_indices : ML mode indices matched to this group
    ref_indices : DFT mode indices in this group

    Returns
    -------
    Subspace overlap in [0, 1]. Value of 1 means perfect subspace alignment.
    """
    # Extract sub-block M of shape (k_calc, k_ref)
    M = alignment_matrix[np.ix_(calc_indices, ref_indices)]
    k = min(len(calc_indices), len(ref_indices))
    if k == 0:
        return 0.0
    # trace(M^T M) / k -- measures how well calc subspace spans ref subspace
    return float(np.trace(M.T @ M) / k)
```

### Collapsed Heatmap Matrix Construction
```python
# Source: custom implementation based on D-09
def collapse_alignment_matrix(
    alignment_matrix: np.ndarray,
    groups: list[DegenerateGroup],
    freqs_calc: np.ndarray,
    freqs_ref: np.ndarray,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Collapse alignment matrix, merging degenerate groups into single entries.

    Returns collapsed matrix, calc labels, ref labels.
    """
    # Identify which ref/calc indices are in groups vs standalone
    grouped_ref = {idx for g in groups for idx in g.ref_indices}
    grouped_calc = {idx for g in groups for idx in g.calc_indices}

    # Build ordered list of "effective" rows and columns
    # Each group becomes one entry, non-grouped modes stay individual
    ref_entries = []  # (label, indices)
    for i in range(len(freqs_ref)):
        if i not in grouped_ref:
            ref_entries.append((f"{freqs_ref[i]:.0f}", [i]))
    for g in groups:
        ref_entries.append((f"{g.symmetry_label} @{g.center_frequency:.0f}", g.ref_indices))

    # Similar for calc side...
    # Build collapsed matrix where group cells = subspace overlap
```

### Distinct Markers in Regression Plot
```python
# Source: D-10 implementation pattern
# Diamond markers for degenerate group averages
ax.scatter(
    group_dft_freqs, group_ml_freqs,
    marker="D",  # Diamond
    c="#D08770",  # Warm orange (colorblind-safe from Nord palette)
    s=100, alpha=0.8, edgecolors="white", linewidth=1.5,
    zorder=4, label="Degenerate group (averaged)"
)
# Standard circles for non-degenerate (existing)
ax.scatter(
    nondeg_dft_freqs, nondeg_ml_freqs,
    c="#5E81AC",  # Existing muted blue
    s=80, alpha=0.7, edgecolors="white", linewidth=1.5,
    zorder=3,
)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Individual dot products for all modes | Still individual dot products | Current | Degenerate modes get misleadingly low scores |
| `seen_freqs` deduplication in parser | `seen_freqs` collapses degenerates | Original | CH4 shows 4 modes instead of 9 in log-parsed output |
| .fchk extraction for harmonic matching | Already preserves all modes | Phase 13.3 | Mode matching already gets all 9 CH4 modes from .fchk |

**Key insight:** The .fchk extraction path (used for eigenvector matching) already handles degenerate modes correctly -- it extracts all of them. The parser deduplication issue (D-03) only affects the log-parsing path, which is used for anharmonic spectra. For harmonic mode matching (the primary use case), the infrastructure is already correct. The parser fix is a correctness improvement, not a blocking prerequisite.

## Open Questions

1. **Symmetry label accuracy for complex molecules**
   - What we know: Multiplicity-based heuristic (2="E", 3="T") works for simple molecules (CH4 Td, NH3 C3v)
   - What's unclear: For molecules with accidental near-degeneracies (two unrelated modes within 5 cm-1), the label would be misleading
   - Recommendation: Use fallback label "deg-{k}" when modes are near-degenerate but likely accidental. Consider adding a note in the report that symmetry labels are heuristic.

2. **Subspace overlap normalization with unequal match counts**
   - What we know: Formula is trace(M^T M)/k. When all k modes are matched, this is well-defined.
   - What's unclear: If only 2 of 3 T2 modes get matched (one ML mode went elsewhere in Hungarian), should we normalize by 2 (matched) or 3 (group size)?
   - Recommendation: Normalize by k_matched (the actual matched count), and flag the group as "partially matched" in the report. This avoids penalizing the subspace overlap for a matching artifact.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `pytest tests/test_mode_matching.py -x` |
| Full suite command | `pytest tests/ --strict-markers -ra` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| MODE-05a | Degenerate group detection from frequencies | unit | `pytest tests/test_mode_matching.py::TestDegenerateDetection -x` | No - Wave 0 |
| MODE-05b | Subspace overlap computation (trace formula) | unit | `pytest tests/test_mode_matching.py::TestSubspaceOverlap -x` | No - Wave 0 |
| MODE-05c | CH4 degenerate modes produce higher overlap with subspace method | integration | `pytest tests/test_mode_matching.py::TestCH4DegenerateOverlap -x` | No - Wave 0 |
| MODE-06a | Statistics count groups as single units | unit | `pytest tests/test_mode_matching.py::TestGroupAwareStatistics -x` | No - Wave 0 |
| MODE-06b | Regression data has one point per group | unit | `pytest tests/test_mode_matching.py::TestGroupRegressionData -x` | No - Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest tests/test_mode_matching.py tests/test_gaussian_parser.py -x`
- **Per wave merge:** `pytest tests/ --strict-markers -ra`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/test_mode_matching.py::TestDegenerateDetection` -- covers MODE-05a (group detection from synthetic frequencies)
- [ ] `tests/test_mode_matching.py::TestSubspaceOverlap` -- covers MODE-05b (trace formula on synthetic matrices)
- [ ] `tests/test_mode_matching.py::TestCH4DegenerateOverlap` -- covers MODE-05c (uses ch4_dft_fchk + ch4_ml_fchk fixtures, already available)
- [ ] `tests/test_mode_matching.py::TestGroupAwareStatistics` -- covers MODE-06a
- [ ] `tests/test_mode_matching.py::TestGroupRegressionData` -- covers MODE-06b

## Sources

### Primary (HIGH confidence)
- `mace_gaussian/analysis/mode_matching.py` -- Current Hungarian matching, alignment matrix, heatmap code (full read)
- `mace_gaussian/analysis/analyze_spectra.py` -- `match_by_mode()`, `extract_spectrum_from_fchk()`, `plot_regression()` (full read)
- `mace_gaussian/analysis/analysis_workflow.py` -- `extract_mode_mapping()`, `generate_mode_overlap_heatmaps()`, comparison pipeline (full read)
- `mace_gaussian/gaussian/parser.py` -- `seen_freqs` deduplication logic (full read)
- `tests/test_mode_matching.py` -- Existing tests including Hungarian bijection tests (full read)
- `tests/test_gaussian_parser.py` -- CH4 deduplication test documenting current behavior (read)
- `tests/conftest.py` -- CH4 fixtures (ch4_dft_fchk, ch4_ml_fchk, ch4_dft_log, ch4_ml_log) confirmed available

### Secondary (MEDIUM confidence)
- Linear algebra: trace(M^T M)/k is a standard subspace overlap metric from numerical linear algebra. The trace of the product of overlap matrices measures how well one subspace spans another. For orthonormal bases this gives the sum of squared cosines of principal angles.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - no new dependencies, all numpy/scipy already in use
- Architecture: HIGH - data flow is clear, single module owns the logic, existing tests cover baseline behavior
- Pitfalls: HIGH - all identified from direct code reading, each has concrete prevention strategy
- Math: HIGH - subspace overlap formula is standard and well-understood; trace(M^T M)/k bounded in [0,1] for normalized eigenvectors

**Research date:** 2026-03-30
**Valid until:** 2026-05-30 (stable domain, no external dependency changes expected)
