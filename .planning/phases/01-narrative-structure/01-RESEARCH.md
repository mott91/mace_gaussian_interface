# Phase 1: Narrative & Structure - Research

**Researched:** 2026-03-06
**Domain:** python-pptx presentation restructuring (slide reorder, content rewrite, function add/remove)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**VPT2 slide:**
- Cut the standalone VPT2 slide (`slide_vpt2`) entirely
- Add 1-2 lines at end of IR theory slide: "anharmonic (VPT2) is already running, results TBD — future direction"
- Do not make anharmonic sound broken; frame it as "next step" or "we've started, not fully validated yet"

**IR theory slide (slide 2):**
- Keep but reframe completely
- New angle: "what the software needs to compute frequencies and intensities" — computational framing, not IR 101
- Cover: harmonic approximation → normal modes → why you need energy + dipole derivatives for each displaced geometry
- End with 1-2 lines on anharmonic/VPT2 as future direction

**Results slides (replacing water + aspirin slides):**
- Two results slides replacing the current two:
  - Slide A (`$ ls results/`): links/references to the HTML analysis reports — water, aspirin, glycine, etc. The HTML report is the artifact; show it exists and looks good
  - Slide B: method comparison table — 8 combos, harmonic, frequencies AND intensities as separate columns, ranked
- The HTML report visual is the "beautiful results" anchor
- The comparison table is the ranked summary

**Slide 10 (status/git log) rewrite:**
- Keep the `$ git log --oneline --graph` command but rewrite content entirely
- Story arc: tool built → we ran the benchmark → here's what we found → here's what's still open
- What we found: frequencies work well across all combos, intensities depend strongly on the dipole model (mace_ml wins over espaloma)
- What's open: anharmonic validation, 25-molecule systematic benchmark, the thesis question
- Physics/chemistry framing — not git versions, not CI status

**ASCII diagrams:**
- Keep all ASCII diagrams — preserved and improved
- Do NOT replace with PNG or matplotlib images

**Motivation slide:**
- Tighten: fewer bullets per section, stronger signal-to-noise
- Keep Problem/Solution/Impact structure
- Remove verbose secondary bullets

**Slide order (final):**
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

### Deferred Ideas (OUT OF SCOPE)
- Color-coded results table (green/yellow/red cells) — Phase 3 (polish)
- Embedded spectrum plot images from analysis_results_harmonic/ — Phase 2 (content)
- Slide numbers — Phase 3 (polish)
- Improved ZMQ ASCII diagram — Phase 2 (content)
- Speaking notes — Phase 3
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| NARR-01 | Slide structure reviewed and reordered to match harmonic-first narrative arc (why → how → what we found → what's next) | Current slide order documented; exact diff of what moves/stays/goes identified |
| NARR-02 | Slide 10 rewritten as research journey — physics and chemistry progress, not git versions | Current slide_status content audited; spine for replacement content documented |
| NARR-03 | VPT2/anharmonic slide scoped to context + "ongoing work" framing (not a deep dive) | slide_vpt2 function identified for deletion; IR slide reframe documented |
</phase_requirements>

---

## Summary

This phase is pure Python source editing on a single file: `/home/mot/mace_gaussian_march26_presentation/mace_gaussian/presentation/generate_pptx_v2.py`. No library research is required — the python-pptx patterns are already established and working in the file. The work is entirely about which functions to keep, which to delete, which to rewrite, and what to add.

The current file has 12 slides across 12 functions called sequentially in `main()`. Phase 1 reduces this to 11 slides by: deleting `slide_vpt2`, rewriting `slide_ir_theory`, rewriting `slide_motivation` (tighten), replacing `slide_results_water` and `slide_results_aspirin` with two new functions (`slide_results_overview` and `slide_results_table`), and rewriting `slide_status`. The architecture, dipole methods, ZMQ, mode matching, title, and questions slides are untouched.

The real data for the comparison table (Slide B) is fully available in `analysis_results_harmonic/water/data/metrics_summary.json` — 8 combos with `mae_freq`, `r2_freq`, `r2_intensity`, and `speedup` fields. HTML reports exist for 8 molecules in `analysis_results_harmonic/` and 9 molecules in `analysis_results/`.

**Primary recommendation:** Edit `generate_pptx_v2.py` directly. Three tasks: (1) delete slide_vpt2 + update main(), (2) rewrite the four content-heavy functions, (3) fix the hardcoded output path.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| python-pptx | already installed | Generate .pptx from Python | Already in use, working |

No new dependencies. All tooling is in place.

**Build command:**
```bash
cd /home/mot/mace_gaussian_march26_presentation/mace_gaussian
source .venv/bin/activate
python presentation/generate_pptx_v2.py
```

---

## Architecture Patterns

### Existing File Structure
```
presentation/
├── generate_pptx_v2.py      # Single source file — all slides here
├── presentation_v2.pptx     # Generated output (overwritten on each run)
├── master_presentation_191125.pdf  # Old Nov 2025 reference (brevity inspiration)
└── PROJECT_CONTEXT.md       # Context notes
```

### Pattern: Each slide = one function
```python
def slide_name(prs):
    # Uses content_slide(), two_col_slide(), or custom layout
    # Returns slide object (or None — not all functions return)
    ...

def main():
    prs = Presentation()
    prs.slide_width  = Inches(10)
    prs.slide_height = Inches(7.5)
    slide_title(prs)
    slide_motivation(prs)
    # ... one call per slide, in order
    prs.save(out)
```

### Pattern: Line content as (text, color) tuples
```python
lines = [
    ("# Section Header", ACCENT),      # blue — section titles
    ("", DIM),                          # blank spacer
    ("  body text here", TEXT),         # grey — normal content
    ("  good result", GREEN),           # green — positive metrics
    ("  warning", YELLOW),              # yellow — caveats
    ("  error or poor result", RED),    # red — failures
    ("  subdued", DIM),                 # dim grey — secondary info
]
```

### Pattern: Two-column vs single-column
- `content_slide(prs, command, lines)` — single textbox, 8.5 inch wide, good for ASCII diagrams and long single-thread narrative
- `two_col_slide(prs, command, left_lines, right_lines, split=4.3)` — split panel, good for comparison/contrast
- For the new results table slide, a custom layout may be needed (see Pitfalls section)

### Anti-Patterns to Avoid
- **Don't use relative paths in output:** The current output path `"/home/mot/mace_gaussian/presentation/presentation_v2.pptx"` is hardcoded to a different machine path. Fix to use the current project's `presentation/` directory.
- **Don't add imports:** All needed imports are already at the top. No new python-pptx primitives are required for this phase.
- **Don't use return values of slide functions:** `main()` calls functions for side effects; functions add slides to `prs` in place.

---

## Current Slide Inventory (Before Phase 1)

| Index | Function | Command Prompt | Keep/Change/Delete |
|-------|----------|---------------|-------------------|
| 0 | `slide_title` | neofetch-style | KEEP unchanged |
| 1 | `slide_motivation` | `$ cat motivation.md` | REWRITE (tighten bullets) |
| 2 | `slide_ir_theory` | `$ man ir_spectroscopy` | REWRITE (computational reframe) |
| 3 | `slide_vpt2` | `$ cat vpt2_theory.md` | DELETE |
| 4 | `slide_architecture` | `$ python3 main.py` | KEEP unchanged |
| 5 | `slide_dipole_methods` | `$ cat dipole_methods.md` | KEEP unchanged |
| 6 | `slide_zmq` | `$ cat zmq_bridge.md` | KEEP unchanged |
| 7 | `slide_mode_matching` | `$ cat mode_matching.md` | KEEP unchanged |
| 8 | `slide_results_water` | `$ cd results/water/...` | REPLACE with slide_results_overview |
| 9 | `slide_results_aspirin` | `$ cd results/aspirin/...` | REPLACE with slide_results_table |
| 10 | `slide_status` | `$ git log --oneline --graph` | REWRITE (physics/chemistry framing) |
| 11 | `slide_questions` | `$ ./questions.sh` | KEEP unchanged |

**After phase 1:** 11 slides. `main()` drops `slide_vpt2`, swaps water/aspirin for overview/table.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Slide layout | Custom textbox math | `content_slide()` / `two_col_slide()` | Already tuned for font, spacing, footer |
| Color coding | Inline RGBColor() calls | Existing constants: GREEN, YELLOW, RED, DIM, ACCENT, TEXT | Consistency |
| Prompt header | Manual textbox | `add_prompt(slide, command)` | Same position/size/style as all other slides |

**Key insight:** All slide primitives are already established. This phase is a content rewrite, not infrastructure work.

---

## Common Pitfalls

### Pitfall 1: Output path points to wrong machine
**What goes wrong:** `main()` currently saves to `"/home/mot/mace_gaussian/presentation/presentation_v2.pptx"` — this is a different path from the current repo at `/home/mot/mace_gaussian_march26_presentation/mace_gaussian/presentation/`.
**Why it happens:** Hardcoded absolute path from original development machine.
**How to avoid:** Change to `"presentation/presentation_v2.pptx"` (relative) or the correct absolute path.
**Warning signs:** `FileNotFoundError` on save, or writing to the wrong location silently.

### Pitfall 2: Table slide using two_col_slide() is too narrow for 8-row table
**What goes wrong:** `two_col_slide()` was designed for bullet lists with ~15-18 lines per panel. A data table with 8 combos + header needs wider columns and smaller font.
**Why it happens:** Consolas Pt(13) with the standard split gives ~30 chars per column — not enough for "mace_anicc_mace_ml | 19.0 | 0.9999 | 0.72".
**How to avoid:** For the comparison table slide (Slide B), use `content_slide()` with a monospaced table formatted as a string block, or reduce font to Pt(11-12) in a custom layout. This is within Claude's discretion on layout.
**Warning signs:** Text overflow, lines wrapping mid-row, misaligned columns.

### Pitfall 3: Deleting slide_vpt2 but forgetting to update main()
**What goes wrong:** If only the function definition is deleted but the `slide_vpt2(prs)` call remains in `main()`, Python raises `NameError` at runtime.
**How to avoid:** Delete both the function body AND the call in `main()`, update the slide count comment.

### Pitfall 4: IR slide reframe loses the "dipole gap" explanation
**What goes wrong:** The current `slide_vpt2` right panel contains the critical "dipole gap" explanation ("MACE energy models: good forces, bad dipoles → need a dedicated dipole model"). If this content is deleted with the slide and not preserved elsewhere, the architectural motivation for having two separate calculator types is lost.
**How to avoid:** The dipole gap concept already exists in `slide_dipole_methods` right panel ("Energy calculators" section). Verify it's still adequately explained after VPT2 deletion. The IR slide reframe does NOT need to duplicate it — just don't lose it entirely.

### Pitfall 5: Motivation slide over-tightening removes key thesis framing
**What goes wrong:** The existing motivation slide has "Thesis question: energy surface vs dipole model quality" in the Impact section. This connects to `slide_status`. Removing it while tightening could break the narrative thread.
**How to avoid:** Keep the thesis question line even when tightening — it's load-bearing for the story arc.

### Pitfall 6: Results data mismatch — aspirin has only 4 combos, water has 8
**What goes wrong:** The aspirin `metrics_summary.json` only contains 4 combinations (mace_mp, mace_omol each paired with espaloma/mace_ml). Water has all 8 (adds mace_anicc and mace_off pairs). For the comparison table slide (Slide B), using only water data is the right call since it has the full 8-combo coverage.
**How to avoid:** Source the table from `analysis_results_harmonic/water/data/metrics_summary.json`. Label it clearly as "Water (H₂O) — 3 modes" to set context.

---

## Real Data Available for Slide B (Comparison Table)

Source: `analysis_results_harmonic/water/data/metrics_summary.json`

All 8 combinations, sorted by frequency MAE (ascending = best first):

| Combo | MAE freq (cm⁻¹) | R² freq | R² intensity | Speedup |
|-------|----------------|---------|--------------|---------|
| mace_anicc_mace_ml | 19.0 | 0.999997 | 0.719 | 1.18× |
| mace_anicc_espaloma | 19.0 | 0.999997 | 0.518 | 1.18× |
| mace_off_espaloma | 20.7 | 0.999934 | 0.516 | 1.27× |
| mace_off_mace_ml | 20.7 | 0.999934 | 0.719 | 1.18× |
| mace_omol_espaloma | 22.7 | 0.999963 | 0.521 | 0.91× |
| mace_omol_mace_ml | 22.7 | 0.999963 | 0.719 | 0.86× |
| mace_mp_mace_ml | 106.3 | 0.999998 | 0.718 | 1.22× |
| mace_mp_espaloma | 106.3 | 0.999998 | 0.514 | 1.01× |

**Key story the data tells:**
- Frequencies: all combos R² > 0.999 (excellent) except mace_mp which has large MAE despite high R² (the R²/MAE paradox with only 3 modes)
- Intensities: mace_ml dipole consistently outperforms espaloma (~0.72 vs ~0.52 R²)
- Speedup: water is too small for meaningful speedup — the story is that speedup grows with molecule size (aspirin: 18-29×)

**HTML reports available** (for Slide A "ls results/"):
- `analysis_results_harmonic/`: ammonia, aspirin, bh3_nh3, C2H6_ase, CH4_ase, co, gly, water (8 molecules)
- `analysis_results/` (anharmonic): acoh, ammonia, aspirin, bh3_nh3, C3H8_ase, CH4_ase, co, formaldehyde, water (9 molecules)

---

## Code Examples

### Deleting slide_vpt2 — diff pattern
```python
# BEFORE main():
slide_vpt2(prs)            # 3 — VPT2, anharmonicity, the dipole gap

# AFTER main() — line removed entirely. Also delete the slide_vpt2() function.
```

### IR slide anharmonic footer lines (Claude's discretion — recommended wording)
```python
# Add at end of IR slide content, after harmonic discussion:
("", DIM),
("# Anharmonic (VPT2): ongoing work", ACCENT),
("  → Already running. Overtones + combination bands.", TEXT),
("  → Results pending full validation — next direction.", DIM),
```

### Results overview slide (Slide A) — terminal ls idiom
```python
def slide_results_overview(prs):
    content_slide(prs, "$ ls analysis_results_harmonic/", [
        ("# HTML reports generated per molecule", ACCENT),
        ("", DIM),
        ("  water/report.html       3 atoms  · 3 modes  · 8 combos", GREEN),
        ("  aspirin/report.html    19 atoms  · 51 modes · 4 combos", GREEN),
        ("  gly/report.html        10 atoms  · 24 modes · 4 combos", GREEN),
        ("  ammonia/report.html     4 atoms  · 6 modes  · 4 combos", GREEN),
        ("  CH4_ase/report.html     5 atoms  · 9 modes  · 4 combos", GREEN),
        ("  C2H6_ase/report.html    8 atoms  · 18 modes · 4 combos", GREEN),
        ("  co/report.html          2 atoms  · 1 mode   · 4 combos", DIM),
        ("  bh3_nh3/report.html     6 atoms  · 12 modes · 4 combos", DIM),
        ("", DIM),
        ("  Each report: regression plots · KDE spectra · mode matching", TEXT),
        ("  Anharmonic reports also in analysis_results/", DIM),
    ])
```

### Comparison table slide (Slide B) — monospaced table
```python
def slide_results_table(prs):
    content_slide(prs, "$ cat results/water/metrics_summary.json | sort-by-freq",  [
        ("# Harmonic benchmark — water (H₂O, 3 modes)", ACCENT),
        ("  Ref: B3LYP/6-31G(d,p)   |   sorted by frequency MAE", DIM),
        ("", DIM),
        ("  Energy model   Dipole     MAE(freq)   R²(freq)   R²(intens)", DIM),
        ("  ─────────────────────────────────────────────────────────────", DIM),
        ("  mace_anicc     mace_ml      19 cm⁻¹   0.999997     0.72  ✓", GREEN),
        ("  mace_anicc     espaloma     19 cm⁻¹   0.999997     0.52", TEXT),
        ("  mace_off       mace_ml      21 cm⁻¹   0.999934     0.72  ✓", GREEN),
        ("  mace_off       espaloma     21 cm⁻¹   0.999934     0.52", TEXT),
        ("  mace_omol      mace_ml      23 cm⁻¹   0.999963     0.72  ✓", GREEN),
        ("  mace_omol      espaloma     23 cm⁻¹   0.999963     0.52", TEXT),
        ("  mace_mp        mace_ml     106 cm⁻¹   0.999998     0.72", YELLOW),
        ("  mace_mp        espaloma    106 cm⁻¹   0.999998     0.51", YELLOW),
        ("", DIM),
        ("  → Frequencies: all combos excellent except mace_mp", TEXT),
        ("  → Intensities: mace_ml dipole consistently beats espaloma", GREEN),
    ])
```

### Status slide rewrite — physics/chemistry spine
```python
# Left panel: story arc
left_lines=[
    ("# What we built", ACCENT),
    ("  ZMQ bridge: ML dipoles injected into Gaussian VPT2", TEXT),
    ("  CLI: mace-gaussian run molecule.xyz", DIM),
    ("  4 energy × 2 dipole models = 8 combinations", TEXT),
    ("", DIM),
    ("# What we ran", ACCENT),
    ("  Harmonic benchmark across 8 molecules", TEXT),
    ("  Water, aspirin, glycine, methane, ethane, NH₃, CO", DIM),
    ("", DIM),
    ("# What we found", ACCENT),
    ("  Frequencies: R² > 0.999 across all combos ✓", GREEN),
    ("  Intensities: mace_ml >> espaloma (0.72 vs 0.52)", GREEN),
    ("  Speedup scales with molecule size (aspirin: ~18×)", GREEN),
]

# Right panel: open questions
right_lines=[
    ("# Still open", ACCENT),
    ("", DIM),
    ("  Anharmonic validation:", TEXT),
    ("  → VPT2 pipeline runs; accuracy TBD", DIM),
    ("  → Overtones + combination bands", DIM),
    ("", DIM),
    ("  25-molecule benchmark campaign:", TEXT),
    ("  → Systematic across functional groups", DIM),
    ("  → Batch runner in progress", DIM),
    ("", DIM),
    ("# Thesis question", ACCENT),
    ("  Energy surface quality", TEXT),
    ("  vs dipole model quality —", TEXT),
    ("  which dominates IR accuracy?", TEXT),
    ("", DIM),
    ("  $ mace-gaussian run molecule.xyz", GREEN),
]
```

---

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|-----------------|--------|
| slide_results_water shows anharmonic data | Slide B shows harmonic data (what's actually validated) | Honest presentation — don't show unvalidated VPT2 numbers as main results |
| VPT2 as standalone slide 3 | 1-2 lines at end of IR theory slide | Removes false impression that anharmonic is complete |
| Status slide lists git tasks and CI status | Status slide tells physics/chemistry research story | Audience cares about science, not CI |

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Manual visual inspection (no automated test for .pptx content) |
| Config file | none |
| Quick run command | `python presentation/generate_pptx_v2.py` |
| Full suite command | `python presentation/generate_pptx_v2.py && echo "Saved OK"` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| NARR-01 | 11 slides in correct order | smoke | `python -c "from pptx import Presentation; p=Presentation('presentation/presentation_v2.pptx'); print(len(p.slides), 'slides')"` | ❌ Wave 0 (run after generation) |
| NARR-02 | Status slide contains research journey text (no "CI", no "131 tests") | smoke | manual inspection of generated .pptx | manual |
| NARR-03 | VPT2 slide absent; "anharmonic" appears only in IR slide footer | smoke | `python -c "from pptx import Presentation; ..."` — text search across slides | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `python presentation/generate_pptx_v2.py` (must exit 0)
- **Per wave merge:** Visual review of generated .pptx slide count and slide order
- **Phase gate:** All 11 slides present in correct order, no NameError on generation

### Wave 0 Gaps
- [ ] Fix output path in `main()` before first run — current path `/home/mot/mace_gaussian/presentation/presentation_v2.pptx` will cause `FileNotFoundError`
- [ ] No pytest fixtures needed — this is a pure generation script

---

## Open Questions

1. **IR slide: single-column or two-column?**
   - What we know: current layout is two-column (`two_col_slide`). Reframe is computational — could flow as single column narrative.
   - What's unclear: whether the reframed content fits naturally in two columns or reads better as a single flow.
   - Recommendation: Use two-column. Left: "what Gaussian needs" (harmonic derivation, normal modes, displaced geometries). Right: "what ML provides" (energy + dipole per geometry → feed into Gaussian loop). This mirrors the architecture slide's ML/DFT split and sets up the ZMQ slide naturally.

2. **Command prompt for Slide A (results overview)?**
   - Options: `$ ls analysis_results_harmonic/`, `$ cat results/index.md`, `$ ls results/`
   - Recommendation: `$ ls analysis_results_harmonic/` — most literal, truest to the terminal idiom.

---

## Sources

### Primary (HIGH confidence)
- Direct inspection of `presentation/generate_pptx_v2.py` — full source read, all functions documented
- Direct inspection of `analysis_results_harmonic/water/data/metrics_summary.json` — all 8 combo data verified
- Direct inspection of `analysis_results_harmonic/aspirin/data/metrics_summary.json` — 4 combo data verified
- Direct filesystem listing of `analysis_results_harmonic/` and `analysis_results/` — HTML report availability confirmed
- `.planning/phases/01-narrative-structure/01-CONTEXT.md` — locked decisions source

### Secondary (MEDIUM confidence)
- `.planning/REQUIREMENTS.md` — requirement IDs and descriptions
- `.planning/STATE.md` — project phase context

### Tertiary (LOW confidence)
- None

---

## Metadata

**Confidence breakdown:**
- Current slide inventory: HIGH — read directly from source file
- Real data values for table: HIGH — read directly from metrics_summary.json
- HTML report availability: HIGH — verified by filesystem listing
- Architecture patterns: HIGH — all patterns observed in existing working code
- Pitfalls: HIGH — derived from direct code inspection (output path bug confirmed)

**Research date:** 2026-03-06
**Valid until:** 2026-04-06 (stable — no external dependencies to go stale)
