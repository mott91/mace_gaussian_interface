# Phase 2: Content - Context

**Gathered:** 2026-03-07
**Status:** Ready for planning

<domain>
## Phase Boundary

Fill every content gap in the deck so that every claim is backed by a visible artifact. This means: completing the results tables with real numbers, making the HTML report links slide concrete and visual, building a dedicated scaling slide, and rewriting the ZMQ slide as a pedagogical single-column flow. No visual polish (color-coding, fonts, slide numbers) — that is Phase 3.

</domain>

<decisions>
## Implementation Decisions

### Plot embedding (CONT-03)
- No PNG images embedded in slides — HTML report links only
- CONT-03 is satisfied by showing the reports exist and are referenced
- The HTML links slide uses terminal-style file listing (`$ ls results/`)
- Add ASCII art of H₂O and aspirin alongside/above the report links — keeps terminal aesthetic, makes the slide visually interesting without images
- Reports to reference: `analysis_results_harmonic/water/report.html`, `analysis_results_harmonic/aspirin/report.html`

### Results table (CONT-01)
- Per-molecule tables: water results on the water slide, aspirin results on the aspirin slide
- Columns: Frequency R² | Frequency MAE (cm⁻¹) | Intensity R² | Intensity MAE — no RMSE (redundant for this audience)
- Row ordering: ranked by frequency R² descending (makes the winner immediately obvious)
- Molecules covered: water + aspirin only (complete harmonic data available; acoh has known parsing bugs)
- Data source: pull from `comparison_results/{molecule}/{combo}/results.json` files

### Scaling argument (CONT-02)
- Gets a dedicated slide (not tucked into the results slides)
- Form: ASCII bar chart — molecule name + bar of `=` characters scaled to speedup, e.g.:
  - `water    (3 atoms)   |=                    |  ~1×`
  - `glycine  (10 atoms)  |========             |  ~7–10×`
  - `aspirin  (21 atoms)  |=====================|  ~18–29×`
- Use estimated speedup ranges as stated — not benchmarked to the second, ranges are honest
- Command: something like `$ python benchmark.py --scaling` to fit the terminal theme

### ZMQ slide redesign (CONT-04)
- Redesign as a single-column pedagogical flow (not left/right split)
- Structure: problem statement → solution (ZMQ bridge) → how it works → why it's non-trivial
- Update server-side reference from `gm_main.py` (no longer exists post-refactor) → `zmq_server.py`
- `gm_helper.py` (Gaussian-side ZMQ client) is still correct — keep it
- Drop function-name level detail from both the ZMQ slide and the architecture slide — focus on what the components DO, not their Python identifiers
- Architecture slide: keep the high-level pipeline shape (geometry opt → DFT → ML → analysis), remove stale function name annotations

### Claude's Discretion
- Exact ASCII art representation of H₂O and aspirin on the links slide
- Specific bar widths/scale for the ASCII speedup chart
- Exact wording of the ZMQ pedagogical narrative
- Whether to regenerate data from JSON or hardcode the key numbers (either is fine)

</decisions>

<specifics>
## Specific Ideas

- ASCII molecule art next to or above the HTML report links — user explicitly asked for H₂O and aspirin ASCII structures
- The scaling bar chart should feel like a terminal benchmark output — tidy columns, molecule name left-padded, bar right-aligned
- ZMQ slide: "problem → solution → how → why non-trivial" as the pedagogical spine
- The results tables should make the winner (mace_omol + mace_ml) obvious from rank position alone

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `content_slide(prs, command, lines)` — standard single-column slide; use for scaling slide and ZMQ redesign
- `two_col_slide(prs, command, left_lines, right_lines, split=4.3)` — still available if needed
- `slide_results_overview(prs)` and `slide_results_table(prs)` — existing results slides to be updated
- Color constants: BG, TEXT, ACCENT, GREEN, YELLOW, DIM, RED — use GREEN for winner row, YELLOW for moderate

### Established Patterns
- Every slide has `$ command` header via `add_prompt()` — keep this for new slides
- Lines are `(text, RGBColor)` tuples — section headers use ACCENT, body TEXT, highlights GREEN/YELLOW/RED
- ASCII art: `Consolas Pt(13-15)`, `space_after=Pt(0)` for tight lines

### Integration Points
- `comparison_results/{molecule}/{combo}/results.json` — source for R²/MAE numbers
- `analysis_results_harmonic/{molecule}/report.html` — paths to reference on the links slide
- `slide_zmq(prs)` at line 325 — function to rewrite
- `slide_architecture(prs)` at line 257 — function to update (remove stale function names)
- All slides assembled in `main()` — scaling slide will need a new slot inserted

### Current ZMQ diagram (to be replaced)
- Currently: left = flow with `gm_main.py` (outdated), right = engineering challenges list
- Replace with: single column pedagogical flow using `zmq_server.py`
- `gm_helper.py` is still the correct Gaussian-side script name

</code_context>

<deferred>
## Deferred Ideas

- Color-coded results table (green/yellow/red cells) — Phase 3 (polish)
- Mode overlap heatmap for aspirin — nice-to-have, post-presentation (V2-01)

</deferred>

---

*Phase: 02-content*
*Context gathered: 2026-03-07*
