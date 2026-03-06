# Phase 1: Narrative & Structure - Context

**Gathered:** 2026-03-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Restructure the slide deck so the story arc is correct and every slide is in the right place. No new content (plots, data tables) — that's Phase 2. This phase changes slide order, removes one slide, reframes two slides, and rewrites one slide.

</domain>

<decisions>
## Implementation Decisions

### VPT2 slide
- Cut the standalone VPT2 slide entirely
- Add 1-2 lines at the end of the IR theory slide: "anharmonic (VPT2) is already running, results TBD — future direction"
- Do not make anharmonic sound broken; frame it as "next step" or "we've started, not fully validated yet"

### IR theory slide (slide 2)
- Keep the IR theory slide but reframe it completely
- New angle: "what the software needs to compute frequencies and intensities" — computational framing, not IR 101
- Cover: harmonic approximation → normal modes → why you need energy + dipole derivatives for each displaced geometry
- This makes it a natural setup for the architecture slide that follows
- End with 1-2 lines on anharmonic/VPT2 as future direction

### Results slides (replacing water + aspirin slides)
- Two results slides replacing the current two:
  - Slide A (`$ ls results/`): links/references to the HTML analysis reports — water, aspirin, glycine etc. The HTML report is the artifact; show it exists and looks good
  - Slide B: method comparison table — 8 combos, harmonic, frequencies AND intensities as separate columns, ranked
- The HTML report visual is the "beautiful results" anchor
- The comparison table is the ranked summary

### Slide 10 (status/git log) rewrite
- Keep the `$ git log --oneline --graph` command but rewrite the content entirely
- Story arc: tool built → we ran the benchmark → here's what we found → here's what's still open
- What we found: frequencies work well across all combos, intensities depend strongly on the dipole model (mace_ml wins over espaloma)
- What's open: anharmonic validation, 25-molecule systematic benchmark, the thesis question
- Physics/chemistry framing — not git versions, not CI status

### ASCII diagrams
- Keep all ASCII diagrams — user wants them preserved and improved
- Do NOT replace with PNG or matplotlib images
- Architecture and ZMQ diagrams can be made cleaner/more readable in Phase 2

### Motivation slide
- Tighten: fewer bullets per section, stronger signal-to-noise
- Keep Problem/Solution/Impact structure
- Remove the verbose secondary bullets

### Slide order (final)
0. Title (neofetch) — unchanged
1. Motivation — tightened
2. IR theory → reframed as "what the software computes"
3. Architecture — unchanged
4. Dipole methods — unchanged
5. ZMQ bridge — unchanged
6. Mode matching — unchanged
7. Results A: HTML report links (`$ ls results/`)
8. Results B: method comparison table (harmonic, 8 combos)
9. Status/journey (`$ git log`) — rewritten
10. Questions — unchanged

Total: 11 slides (VPT2 removed, water+aspirin collapsed into two new results slides = net -1)

### Claude's Discretion
- Exact wording of the 1-2 anharmonic lines at end of IR slide
- Whether the IR slide becomes single-column or stays two-column
- The exact command prompt for the HTML report slide

</decisions>

<specifics>
## Specific Ideas

- Old presentation (Nov 2025) was more concise — take inspiration from its brevity
- The HTML report slides should feel like `$ cat results/water.html` or `$ ls results/` — terminal idiom
- Slide 10 spine: "tool built → benchmark run → frequencies solved, intensities model-dependent → anharmonic + 25-mol campaign next"
- User loves the terminal dark aesthetic and ASCII art — lean into it, don't soften it

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `content_slide(prs, command, lines)` — single column, takes list of (text, color) tuples
- `two_col_slide(prs, command, left_lines, right_lines, split=4.3)` — two panels
- `blank_slide(prs)`, `add_prompt(slide, command)`, `add_footer(slide)` — building blocks
- Color constants: BG, TEXT, ACCENT, GREEN, YELLOW, DIM, RED — all defined at top of file

### Established Patterns
- Every slide has a `$ command` header via `add_prompt()` — this is the terminal idiom; keep it
- Lines are `(text, RGBColor)` tuples — section headers use ACCENT, body uses TEXT, good results GREEN, warnings YELLOW, errors RED, subdued DIM
- ASCII art goes in `add_textbox` with Consolas Pt(13-15), space_after=Pt(0) for tight lines

### Integration Points
- All slides assembled in `main()` function — add/remove/reorder function calls there
- Each slide is its own function (`slide_X(prs)`) — easy to add/remove
- Output path hardcoded in `main()`: needs updating from old path to `presentation/presentation_v2.pptx`
- HTML reports at: `analysis_results_harmonic/{molecule}/report.html`

</code_context>

<deferred>
## Deferred Ideas

- Color-coded results table (green/yellow/red cells) — Phase 3 (polish)
- Embedded spectrum plot images from analysis_results_harmonic/ — Phase 2 (content)
- Slide numbers — Phase 3 (polish)
- Improved ZMQ ASCII diagram — Phase 2 (content)
- Speaking notes — Phase 3

</deferred>

---

*Phase: 01-narrative-structure*
*Context gathered: 2026-03-06*
