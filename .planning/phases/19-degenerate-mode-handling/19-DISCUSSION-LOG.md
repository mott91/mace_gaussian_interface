# Phase 19: Degenerate Mode Handling - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md -- this log preserves the alternatives considered.

**Date:** 2026-03-30
**Phase:** 19-degenerate-mode-handling
**Areas discussed:** Degeneracy detection, Subspace overlap math, Statistics reporting, Heatmap & visuals

---

## Degeneracy Detection

### Q1: Which side should define degenerate groups?

| Option | Description | Selected |
|--------|-------------|----------|
| DFT reference only | Group modes by DFT frequencies within threshold. DFT is ground truth. | |
| Both DFT and ML independently | Detect groups on each side separately, match group-to-group. | |
| Frequency-agnostic (eigenvector clustering) | Use eigenvector similarity instead of frequency proximity. | |

**User's choice:** Free-text response -- user asked whether grouping matters outside DFT comparison context. After discussion confirming grouping is a comparison concern only, agreed on DFT reference only.
**Notes:** User raised good point about standalone ML use. Confirmed that degenerate grouping only matters when comparing eigenvectors between two calculations. Standalone spectra have correct frequencies/intensities regardless of eigenvector orientation.

### Q2: Frequency threshold for grouping?

| Option | Description | Selected |
|--------|-------------|----------|
| 5 cm-1 | Matches ROADMAP spec and MODE-05 requirement. | ✓ |
| Configurable with 5 cm-1 default | Add --degeneracy-threshold parameter. | |
| You decide | Claude picks. | |

**User's choice:** 5 cm-1 (Recommended)

### Q3: Should the parser stop deduplicating degenerate modes?

| Option | Description | Selected |
|--------|-------------|----------|
| Keep all modes | Bypass seen_freqs deduplication for harmonic mode matching. | ✓ |
| Keep current deduplication | Parser continues collapsing degenerates. | |

**User's choice:** Keep all modes (Recommended)

---

## Subspace Overlap Math

### Q4: How should subspace overlap integrate with Hungarian matching?

| Option | Description | Selected |
|--------|-------------|----------|
| Two-stage: Hungarian first, then group eval | Run Hungarian on individual modes, then post-process degenerate groups. | ✓ |
| Group-aware Hungarian | Modify cost matrix for block matching. | |
| You decide | Claude picks. | |

**User's choice:** Two-stage: Hungarian first, then group eval (Recommended)

### Q5: Should subspace overlap replace or supplement individual scores?

| Option | Description | Selected |
|--------|-------------|----------|
| Replace for degenerate groups | Group subspace overlap replaces individual dot products. Non-degenerate modes keep individual overlap. | ✓ |
| Supplement (both scores) | Keep individual overlaps AND add group-level field. | |
| You decide | Claude picks. | |

**User's choice:** Replace for degenerate groups (Recommended)

---

## Statistics Reporting

### Q6: How should degenerate groups count toward statistics?

| Option | Description | Selected |
|--------|-------------|----------|
| One entry per group | T2 counts as 1 matched unit, not 3. | ✓ |
| Count individual modes, group-weighted | Still count 3 modes for T2, each with same group overlap. | |
| You decide | Claude picks. | |

**User's choice:** One entry per group (Recommended)

### Q7: How should degenerate groups contribute regression data points?

| Option | Description | Selected |
|--------|-------------|----------|
| Average within group | Each group = 1 data point (avg DFT freq vs avg ML freq). | ✓ |
| All individual modes | Keep all 3N-6 data points. | |
| You decide | Claude picks. | |

**User's choice:** Average within group (Recommended)

### Q8: Should the report label degenerate groups?

| Option | Description | Selected |
|--------|-------------|----------|
| Yes, label groups | Label with type and multiplicity (e.g., "T2 (3-fold) at 1356 cm-1"). | ✓ |
| Minimal -- just note count | Just mention "N degenerate groups detected". | |
| You decide | Claude picks. | |

**User's choice:** Yes, label groups (Recommended)

---

## Heatmap & Visuals

### Q9: How should degenerate groups appear in heatmaps?

| Option | Description | Selected |
|--------|-------------|----------|
| Collapsed group rows/cols | Each group = one row/column, labeled with group type. | ✓ |
| Individual modes with group borders | Keep all modes, draw borders around groups. | |
| You decide | Claude picks. | |

**User's choice:** Collapsed group rows/cols (Recommended)

### Q10: Should regression plots distinguish degenerate group averages?

| Option | Description | Selected |
|--------|-------------|----------|
| Different markers | Diamond/square for group averages vs circles for non-degenerate. | ✓ |
| Same markers, no distinction | All points look the same. | |
| You decide | Claude picks. | |

**User's choice:** Different markers (Recommended)

---

## Claude's Discretion

- Specific marker shapes and colors for degenerate vs non-degenerate
- Fallback labeling when symmetry type can't be auto-determined
- Internal data structure for degenerate groups
- Summary line format for report

## Deferred Ideas

- Peak assignment confidence scores -- separate enhancement, backlog
- Reduced masses in matching -- orthogonal concern, backlog
