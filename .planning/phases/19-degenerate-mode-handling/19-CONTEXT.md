# Phase 19: Degenerate Mode Handling - Context

**Gathered:** 2026-03-30
**Status:** Ready for planning

<domain>
## Phase Boundary

Make mode matching correctly handle degenerate modes (e.g., methane T2, ammonia E) so overlap scores reflect true subspace alignment instead of arbitrary eigenvector orientation within degenerate subspaces. Correct statistics to avoid double-counting. This is purely a mode matching and reporting improvement -- no changes to frequency calculations, broadening, or the pipeline itself.

</domain>

<decisions>
## Implementation Decisions

### Degeneracy detection
- **D-01:** Degenerate groups are defined by DFT reference frequencies only (not ML side). ML modes get assigned to DFT-defined groups via Hungarian matching. This is correct because degeneracy grouping is a comparison concern -- standalone spectra don't need it.
- **D-02:** Frequency threshold for grouping: fixed at 5 cm-1 (per MODE-05 requirement). No CLI flag for configurability.
- **D-03:** Parser must stop deduplicating degenerate modes. The .fchk extraction already gets all degenerate modes; the Gaussian log parser's `seen_freqs` deduplication must be bypassed for harmonic mode matching. We need all 9 CH4 modes (not collapsed to 4) so subspace overlap can operate on full degenerate sets.

### Subspace overlap math
- **D-04:** Two-stage approach: run Hungarian on individual modes first (unchanged), then post-process degenerate groups with subspace overlap (trace(M^T M)/k). Hungarian stays simple; degenerate awareness is a reporting layer on top.
- **D-05:** For modes within a degenerate group, the group subspace overlap REPLACES individual dot product overlaps. Non-degenerate modes keep their individual overlap score unchanged.

### Statistics reporting
- **D-06:** One entry per degenerate group in statistics. A triply-degenerate T2 group counts as 1 matched unit with 1 subspace overlap score, not 3. Total mode count = (non-degenerate modes) + (degenerate groups).
- **D-07:** For frequency/intensity regression: each degenerate group contributes 1 data point (average DFT freq vs average ML freq; same for intensities). Avoids inflating R2 with correlated near-identical points.
- **D-08:** Report labels degenerate groups explicitly with type and multiplicity (e.g., "T2 (3-fold) at 1356 cm-1") in mode matching output.

### Heatmap & visuals
- **D-09:** Mode overlap heatmaps collapse degenerate groups into single rows/columns, labeled with group type (e.g., "T2 @1356"). Cell shows subspace overlap. Heatmap shrinks from 9x9 to ~4x4 for CH4.
- **D-10:** Regression plots use distinct marker style (diamond/square) for degenerate group averages vs circles for non-degenerate modes, annotated with group multiplicity.

### Claude's Discretion
- Specific marker shapes and colors for degenerate vs non-degenerate in plots
- How to label degenerate groups when symmetry type can't be auto-determined (fallback labeling)
- Internal data structure for representing degenerate groups in match results
- Whether to add a summary line in the report like "3 degenerate groups detected (2 doubly, 1 triply)"

### Folded Todos
- **"Handle degenerate modes and zero-intensity filtering in analysis"** (`.planning/todos/pending/2026-03-26-handle-degenerate-modes-and-zero-intensity-filtering.md`) -- the degenerate mode portion directly matches this phase. Zero-intensity filtering was already completed in Phase 18.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Mode matching code
- `mace_gaussian/analysis/mode_matching.py` -- Current Hungarian matching: `match_modes()` (line 164), `create_alignment_matrix()` (line 225), `compute_mode_overlap()` (line 130). This is where subspace overlap post-processing will be added.
- `mace_gaussian/analysis/analyze_spectra.py` -- `SpectralAnalyzer.match_by_mode()` (line 343), `extract_spectrum_from_fchk()` (line 109) which already extracts all degenerate modes from .fchk files.

### Analysis pipeline
- `mace_gaussian/analysis/analysis_workflow.py` -- Lines 414-434 handle harmonic vs anharmonic mode extraction. Comments already note "Extract directly from .fchk files to get ALL degenerate modes" vs "results.json (collapsed degenerates)".
- `mace_gaussian/analysis/html_report_generator.py` -- Report generation including mode matching tables and plots.

### Gaussian parsing
- `mace_gaussian/gaussian/parser.py` -- Contains `seen_freqs` deduplication logic that needs to be bypassed for degenerate mode preservation.

### Tests
- `tests/test_mode_matching.py` -- Existing mode matching tests including Hungarian assignment test (line 164).
- `tests/test_gaussian_parser.py` -- `test_ch4_harmonic_degenerate_deduplication` (line 75) documents current deduplication behavior.

### Requirements
- `.planning/REQUIREMENTS.md` section v1.2 -- MODE-05 (degenerate detection, subspace overlap) and MODE-06 (statistics without double-counting).

### Related todo
- `.planning/todos/pending/2026-03-26-handle-degenerate-modes-and-zero-intensity-filtering.md` -- Original problem description and solution sketch (folded into this phase).

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `match_modes()` in `mode_matching.py`: Hungarian algorithm already works on individual modes -- subspace overlap is a post-processing step, not a replacement.
- `create_alignment_matrix()`: Builds full NxM overlap matrix. Can be reused to build sub-matrices for degenerate groups.
- `extract_spectrum_from_fchk()`: Already extracts all degenerate modes from .fchk files (bypasses parser deduplication).
- `scipy.optimize.linear_sum_assignment`: Already imported and used.

### Established Patterns
- Mode matching returns `dict[int, tuple[int | None, float]]` mapping calc_idx -> (ref_idx, overlap). Degenerate group handling needs to extend or wrap this return format.
- Heatmap generation uses `force_harmonic=True` to get true normal modes. Degenerate collapsing should happen after eigenvector extraction.
- Regression plots already differentiate confident (filled) vs low-overlap (open circles) matches (Phase 13.3).

### Integration Points
- `match_by_mode()` in `SpectralAnalyzer` calls `match_modes()` and builds comparison DataFrames -- this is where group-level aggregation flows into statistics.
- Heatmap plotting in `analyze_spectra.py` consumes the overlap matrix directly -- collapsing needs to happen before or during heatmap rendering.
- HTML report generator reads comparison data from analysis workflow -- group labels and degenerate flags need to flow through to the report.

</code_context>

<specifics>
## Specific Ideas

- User confirmed that degenerate grouping is a comparison concern only -- standalone spectra don't need it, which validates the DFT-reference-only detection approach.
- The existing todo describes the problem well: "Individual dot products between DFT and ML eigenvectors give misleading low overlaps even when the subspaces are identical."
- CH4 is the canonical test case: 9 modes with T2 (3-fold) degeneracies that currently produce systematically pessimistic overlap scores.

</specifics>

<deferred>
## Deferred Ideas

### Reviewed Todos (not folded)
- "Peak assignment confidence scores from mode matching" -- related but a separate enhancement layer on top of matching; could be Phase 23 or backlog.
- "Investigate and use compute_reduced_masses in mode matching" -- orthogonal concern about weighting modes by mass; backlog.

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 19-degenerate-mode-handling*
*Context gathered: 2026-03-30*
