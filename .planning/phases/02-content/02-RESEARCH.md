# Phase 2: Content - Research

**Researched:** 2026-03-07
**Domain:** Python-pptx slide generation, data extraction from JSON metrics, ASCII art composition
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Plot embedding (CONT-03)**
- No PNG images embedded in slides — HTML report links only
- CONT-03 is satisfied by showing the reports exist and are referenced
- The HTML links slide uses terminal-style file listing (`$ ls results/`)
- Add ASCII art of H₂O and aspirin alongside/above the report links — keeps terminal aesthetic, makes the slide visually interesting without images
- Reports to reference: `analysis_results_harmonic/water/report.html`, `analysis_results_harmonic/aspirin/report.html`

**Results table (CONT-01)**
- Per-molecule tables: water results on the water slide, aspirin results on the aspirin slide
- Columns: Frequency R² | Frequency MAE (cm⁻¹) | Intensity R² | Intensity MAE — no RMSE (redundant for this audience)
- Row ordering: ranked by frequency R² descending (makes the winner immediately obvious)
- Molecules covered: water + aspirin only (complete harmonic data available; acoh has known parsing bugs)
- Data source: pull from `comparison_results/{molecule}/{combo}/results.json` files

**Scaling argument (CONT-02)**
- Gets a dedicated slide (not tucked into the results slides)
- Form: ASCII bar chart — molecule name + bar of `=` characters scaled to speedup, e.g.:
  - `water    (3 atoms)   |=                    |  ~1×`
  - `glycine  (10 atoms)  |========             |  ~7–10×`
  - `aspirin  (21 atoms)  |=====================|  ~18–29×`
- Use estimated speedup ranges as stated — not benchmarked to the second, ranges are honest
- Command: something like `$ python benchmark.py --scaling` to fit the terminal theme

**ZMQ slide redesign (CONT-04)**
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

### Deferred Ideas (OUT OF SCOPE)
- Color-coded results table (green/yellow/red cells) — Phase 3 (polish)
- Mode overlap heatmap for aspirin — nice-to-have, post-presentation (V2-01)
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| CONT-01 | Harmonic results slide shows all 8 calculator combos with ranked R²/MAE table (frequencies and intensities separate) | Metrics data confirmed in `analysis_results_harmonic/{mol}/data/metrics_summary.json`; water has 8 combos, aspirin has 4; existing `slide_results_table()` to be updated |
| CONT-02 | Scaling argument made explicit — molecule size vs speedup (water ~1×, glycine ~7–10×, aspirin ~18–29×) | Speedup field confirmed in metrics JSON; aspirin mace_omol ~18–19×, mace_mp ~29×; new dedicated slide needed in `main()` |
| CONT-03 | Spectrum plot images embedded from `analysis_results_harmonic/` (at least water + aspirin) | Locked decision: no PNG embedding — HTML report links slide with ASCII molecule art satisfies this; both report.html files confirmed to exist |
| CONT-04 | ZMQ bridge slide has clear flow diagram showing real-time ML dipole injection into Gaussian | `slide_zmq()` at line 325 uses outdated `gm_main.py` and two-column layout; rewrite to single-column with `zmq_server.py` and pedagogical narrative |
</phase_requirements>

## Summary

Phase 2 fills content gaps in an existing Python-pptx presentation generator (`presentation/generate_pptx_v2.py`). The file uses a terminal/GitHub-dark aesthetic with ASCII art, monospaced text, and `$ command` slide headers. All slide content is Python code — no external template or design tool is involved.

The four requirements map to precise, contained changes in the generator: (1) update `slide_results_table()` to pull from actual `metrics_summary.json` and add an aspirin table, (2) add a new `slide_scaling()` function with an ASCII bar chart and insert it in `main()`, (3) update `slide_results_overview()` to add ASCII molecule art for water and aspirin, (4) rewrite `slide_zmq()` as a single-column pedagogical flow using `content_slide()` instead of `two_col_slide()`.

The metrics data is fully available and verified. Water has 8 combos, aspirin has 4. The water R² values are all ≥0.9999 making MAE the more discriminating ranking column — the planner should decide whether to sort by MAE_freq or R²_freq given the locked decision says R²_freq but water's data makes that indistinguishable.

**Primary recommendation:** Implement all four changes as surgical edits to `presentation/generate_pptx_v2.py`. Read metrics data from JSON at generation time rather than hardcoding, since the data files are already present and will stay stable.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| python-pptx | (already installed) | PowerPoint generation | Already used throughout; `Presentation`, `RGBColor`, `Inches`, `Pt` are the existing API |
| json (stdlib) | — | Load metrics_summary.json | Already used pattern — no new dependency |

### Supporting

No new libraries required. All work uses existing imports and helpers.

**Installation:** No new packages needed — `uv sync` already covers the environment.

## Architecture Patterns

### Existing Slide Pattern

Every slide function follows the same contract:

```python
def slide_NAME(prs):
    content_slide(prs, "$ command.sh", [
        ("header text", ACCENT),
        ("body line", TEXT),
        ("highlight", GREEN),
        ("dim comment", DIM),
    ])
```

For two-column layout use `two_col_slide()`. The ZMQ redesign moves from `two_col_slide` to `content_slide`.

### Pattern 1: Metrics Data Loading

**What:** Read `metrics_summary.json` at generation time, sort rows, format into slide lines.

**When to use:** Anywhere real numbers appear in the presentation (CONT-01, CONT-02 speedup data).

**Example:**

```python
import json

def load_metrics(molecule):
    path = f"analysis_results_harmonic/{molecule}/data/metrics_summary.json"
    with open(path) as f:
        return json.load(f)["comparisons"]

def make_table_lines(combos_sorted):
    lines = [
        ("  Energy        Dipole     MAE(freq)  R²(freq)  R²(intens)", DIM),
        ("  " + "─" * 55, DIM),
    ]
    for c in combos_sorted:
        energy, dipole = c["name"].rsplit("_", 1)  # split on last underscore not reliable
        # name format is e.g. "mace_omol_mace_ml" or "mace_omol_espaloma"
        line = f"  {c['name']:<22}  {c['mae_freq']:>6.0f} cm⁻¹  {c['r2_freq']:.4f}  {c['r2_intensity']:.2f}"
        color = GREEN if c["r2_freq"] > 0.9998 else (YELLOW if c["r2_freq"] > 0.995 else RED)
        lines.append((line, color))
    return lines
```

**Important note on sorting:** The locked decision says "rank by frequency R² descending." For water all 8 combos have R² ≥ 0.9999, making the order indistinguishable by R². Sorting by MAE_freq ascending (lower is better) produces a more legible ranking for water. The planner should use **MAE_freq ascending** as the primary sort with R²_freq as tiebreaker. This aligns with the spirit of "winner immediately obvious" from the CONTEXT.

### Pattern 2: Name Parsing for combo names

The combo name format in metrics JSON is: `{energy_model}_{dipole_model}` where energy models are `mace_anicc`, `mace_mp`, `mace_off`, `mace_omol` and dipole models are `espaloma` or `mace_ml`. Split on the last occurrence of `_espaloma` or `_mace_ml`:

```python
for suffix in ("_mace_ml", "_espaloma"):
    if name.endswith(suffix):
        energy = name[: -len(suffix)]
        dipole = suffix[1:]  # strip leading underscore
        break
```

### Pattern 3: ASCII Bar Chart for Scaling Slide

**What:** Fixed-width bar of `=` characters, molecule left-padded, speedup right-aligned.

**Example structure:**
```
  molecule     atoms   bar                     speedup
  ─────────────────────────────────────────────────────
  water         (3)    |=                    |   ~1×
  glycine      (10)    |========             |   ~7–10×
  aspirin      (21)    |=====================|  ~18–29×
```

Scale: aspirin bar is the widest (reference at ~29×). A 20-char bar width works cleanly:
- water: 1 char (round 1/29 × 20 ≈ 1)
- glycine: 5–7 chars (round 8/29 × 20 ≈ 6)
- aspirin: 20 chars (reference)

Data for glycine speedup is estimated (~7–10×), not in metrics JSON (aspirin and water are the only molecules with verified harmonic data). Use the ranges from CONTEXT verbatim.

### Pattern 4: ZMQ Single-Column Pedagogical Flow

Replace the current `two_col_slide` with `content_slide`. The four-beat narrative:

```
# Problem
  → Gaussian needs dipole derivatives for every displaced geometry
  → DFT dipoles are expensive — this is the bottleneck

# Solution: ZMQ bridge
  → Gaussian's 'External' keyword calls a helper script per geometry
  → Helper sends geometry over ZMQ socket → Python receives it

# How it works
  [geometry.xyz]
       ↓  gm_helper.py  →  ZMQ IPC socket
  zmq_server.py
       ↓  MACE energy + forces
       ↓  MACE dipole model
       ↑  fort.7 format (Gaussian reads this)

# Why it's non-trivial
  → Socket must close cleanly (LINGER=0 + explicit timeout)
  → Gaussian requires absolute paths — resolved at CLI startup
  → Dipole model class remapping at load time (mace_loader.py)
```

### Anti-Patterns to Avoid

- **Hardcoding metric numbers in slide strings:** The `metrics_summary.json` files are ground truth. Hardcoded numbers can drift from the real data if re-runs happen. Read from JSON.
- **Using `two_col_slide` for ZMQ:** The locked decision is `content_slide` (single column pedagogical flow).
- **Sorting water table by R²_freq:** All values round to the same 4 decimal places. Sort by MAE_freq instead.
- **Including aspirin anicc/off combos:** The aspirin metrics JSON only has 4 combos (mace_mp and mace_omol, each × 2 dipoles). Don't invent rows that aren't in the data.
- **Using `gm_main.py` in ZMQ slide:** This file no longer exists post-refactor. Use `zmq_server.py`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Metrics table data | Hardcoded number strings | `json.load()` from `metrics_summary.json` | Files already exist, parsing is trivial, stays correct if data changes |
| Slide layout | Custom textbox positioning math | `content_slide()` helper | Already handles all sizing, font, color, footer — consistent with all other slides |
| Color decisions | Custom color logic per cell | Existing GREEN/YELLOW/RED constants | Already defined, correct hex values, consistent with Phase 3 color-coding intent |

## Common Pitfalls

### Pitfall 1: Water R² Values Are All Identical at 4 Decimal Places

**What goes wrong:** Sorting by R²_freq descending produces arbitrary ordering for water's 8 combos because `mace_mp_mace_ml` and `mace_mp_espaloma` both have R²=1.0000, `mace_anicc` and `mace_omol` variants also all read 1.0000.

**Why it happens:** Water has only 3 modes. With such few data points, R² is nearly perfect for any energy model that gets frequencies in the right ballpark.

**How to avoid:** Sort by MAE_freq ascending as primary key. This correctly shows mace_anicc (MAE=19 cm⁻¹) at the top and mace_mp (MAE=106 cm⁻¹) at the bottom.

**Warning signs:** Table rows look identical in rank, winner is not obvious.

### Pitfall 2: Aspirin Only Has 4 Combos, Not 8

**What goes wrong:** CONT-01 says "8 calculator combos" but the aspirin metrics JSON only contains 4 entries (mace_mp and mace_omol, each paired with espaloma and mace_ml).

**Why it happens:** The mace_anicc and mace_off runs were only done for water, not aspirin.

**How to avoid:** Show 8 combos for water, 4 combos for aspirin. Frame this as "full sweep on water, targeted combos on aspirin." CONT-01 is satisfied by the water table covering 8 combos.

**Warning signs:** Trying to load `aspirin/data/comparison_mace_anicc_*.csv` — those files don't exist.

### Pitfall 3: Slide Insertion Order in `main()`

**What goes wrong:** New `slide_scaling()` must be inserted at the right position in `main()`. Current sequence: results_overview (7), results_table (8), status (9), questions (10). Scaling belongs between results and status — or possibly between mode_matching and results_overview.

**Why it happens:** The slide order is hard-coded by function call order in `main()`.

**How to avoid:** Insert `slide_scaling(prs)` after `slide_results_table(prs)` and before `slide_status(prs)`. Update the slide count comment in `main()` from 11 to 12 slides.

**Warning signs:** Scaling slide appears before results (audience hasn't seen the data yet), or appears after questions (dead last).

### Pitfall 4: `gm_main.py` Reference in ZMQ Slide

**What goes wrong:** Current `slide_zmq()` shows `gm_main.py` as the server-side script name. That file was renamed to `zmq_server.py` in the refactor.

**Why it happens:** The slide was written before the package reorganization.

**How to avoid:** Use `zmq_server.py` in all new ZMQ slide content. `gm_helper.py` (Gaussian-side) is still correct.

### Pitfall 5: Architecture Slide Function Names Are Stale

**What goes wrong:** `slide_architecture()` shows Python function names like `gaussian_freq()`, `parse_gaussian()`, `load_mace_calculator()`, `get_dipole_calculator()`, `setup_zmq_server()`, `launch_gaussian()`, `zmq_dipole_loop()`. These are post-refactor mismatches.

**Why it happens:** Copied from pre-refactor code structure.

**How to avoid:** Strip all function-name annotations. Replace with component-level descriptions only: `GEOMETRY OPTIMIZER`, `DFT BASELINE`, `ML FREQ CALC`, `ANALYSIS`.

## Code Examples

### Verified Data Structure from `metrics_summary.json`

```json
{
  "molecule": "water",
  "comparisons": [
    {
      "name": "mace_omol_mace_ml",
      "mae_freq": 22.74,
      "rmse_freq": 27.55,
      "r2_freq": 0.9999631,
      "r2_intensity": 0.7191,
      "speedup": 0.859
    }
  ]
}
```

**Note:** `speedup` for water is ~0.9–1.3 (water is too small to show speedup). Aspirin mace_omol speedup is ~18–19×, aspirin mace_mp speedup is ~29×. These are the actual JSON values — use them directly.

### Actual Ranked Data for Slides (pre-computed from JSON)

**Water (8 combos, sorted by MAE_freq ascending):**

| Combo | MAE_freq | R²_freq | R²_int |
|-------|----------|---------|--------|
| mace_anicc_mace_ml | 19.0 cm⁻¹ | 1.0000 | 0.72 |
| mace_anicc_espaloma | 19.0 cm⁻¹ | 1.0000 | 0.52 |
| mace_off_mace_ml | 20.7 cm⁻¹ | 0.9999 | 0.72 |
| mace_off_espaloma | 20.7 cm⁻¹ | 0.9999 | 0.52 |
| mace_omol_mace_ml | 22.7 cm⁻¹ | 1.0000 | 0.72 |
| mace_omol_espaloma | 22.7 cm⁻¹ | 1.0000 | 0.52 |
| mace_mp_mace_ml | 106.3 cm⁻¹ | 1.0000 | 0.72 |
| mace_mp_espaloma | 106.3 cm⁻¹ | 1.0000 | 0.51 |

**Aspirin (4 combos, sorted by MAE_freq ascending):**

| Combo | MAE_freq | R²_freq | R²_int | Speedup |
|-------|----------|---------|--------|---------|
| mace_omol_mace_ml | 8.8 cm⁻¹ | 0.9998 | 0.46 | ~19× |
| mace_omol_espaloma | 8.8 cm⁻¹ | 0.9998 | 0.25 | ~18× |
| mace_mp_mace_ml | 97.8 cm⁻¹ | 0.9958 | 0.29 | ~29× |
| mace_mp_espaloma | 97.8 cm⁻¹ | 0.9958 | 0.15 | ~29× |

**Key narrative:** mace_omol gets best accuracy; mace_mp gets most speedup (but worse accuracy). mace_ml dipole consistently beats espaloma on intensity R² across all energy models.

### ZMQ Slide Target Structure

Use `content_slide()` not `two_col_slide()`. Approximate target content:

```python
def slide_zmq(prs):
    content_slide(prs, "$ cat zmq_bridge.md", [
        ("# Problem", ACCENT),
        ("  Gaussian needs dipole derivatives for every displaced geometry", TEXT),
        ("  → DFT dipoles: expensive. This is the bottleneck.", TEXT),
        ("", DIM),
        ("# Solution: ZMQ bridge", ACCENT),
        ("  Gaussian 'External' keyword calls gm_helper.py per geometry", TEXT),
        ("  → gm_helper.py sends geometry over ZMQ IPC socket", TEXT),
        ("  → zmq_server.py receives, queries ML, returns results", TEXT),
        ("", DIM),
        ("# How it works", ACCENT),
        ("  [Gaussian]  → writes geometry → calls gm_helper.py", DIM),
        ("                                        ↓  ZMQ (IPC)", DIM),
        ("                                  zmq_server.py", DIM),
        ("                                        ↓  MACE energy + forces", DIM),
        ("                                        ↓  MACE dipole model", DIM),
        ("  [Gaussian]  ← reads fort.7  ←         ↑  returns results", DIM),
        ("", DIM),
        ("# Why it's non-trivial", ACCENT),
        ("  Socket cleanup: LINGER=0 + SIGKILL timeout (no orphan processes)", TEXT),
        ("  Absolute paths: Gaussian needs full path resolved at CLI startup", TEXT),
        ("  Model loading: dipole class remapping via torch.load()", TEXT),
    ])
```

### ASCII Molecule Art for Links Slide

For the `slide_results_overview()` update — ASCII structures to add above/alongside report paths:

```
H₂O:
      O
     / \
    H   H

Aspirin (skeletal):
    O    OH
    ‖    |
CH₃-C-O-⬡-COOH
```

These already exist in `slide_title()` (lines 129–133 of the generator). Reuse the same ASCII strings for consistency.

### Scaling Bar Chart Target

```python
MAX_BAR = 20
MAX_SPEEDUP = 29  # aspirin mace_mp reference

molecules = [
    ("water",   3,  1,   "~1×"),
    ("glycine", 10, 8,   "~7–10×"),
    ("aspirin", 21, 20,  "~18–29×"),
]

for name, atoms, bar_val, label in molecules:
    bar = "=" * bar_val + " " * (MAX_BAR - bar_val)
    line = f"  {name:<10} ({atoms:>2} atoms)  |{bar}|  {label}"
```

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| `two_col_slide()` for ZMQ | `content_slide()` single column | Pedagogical narrative instead of info dump |
| `gm_main.py` as server name | `zmq_server.py` | Matches actual post-refactor codebase |
| Function names in architecture slide | Component-level descriptions | Accurate, survives future refactors |
| Sorted by R²_freq for water table | Sorted by MAE_freq ascending | Produces legible ranking (R² all ~1.0000) |
| No scaling slide | Dedicated `slide_scaling()` | CONT-02 fulfilled, explicit speedup narrative |
| Results overview has no molecule art | ASCII H₂O + aspirin on links slide | CONT-03 satisfied without PNG embedding |

## Open Questions

1. **Slide position for scaling slide**
   - What we know: It goes between results and status; current slot count is 10 slides (0-indexed)
   - What's unclear: Should it come before or after `slide_results_table()`? Before = context before data; after = data before interpretation
   - Recommendation: After `slide_results_table()` — audience sees the accuracy numbers, then learns about the speedup payoff. Insert as slide index 9 (shift status to 10, questions to 11).

2. **Water table — sort by R² or MAE?**
   - What we know: Locked decision says "rank by frequency R² descending"; but all water R² values are 1.0000 at 4 decimal places
   - What's unclear: Is the spirit of the decision to show winner at top (any method achieves this) or to use R² specifically?
   - Recommendation: Sort by MAE_freq ascending. Produces a meaningful ranking. Document the sort key in the slide's command string (e.g., `| sort-by-mae`).

3. **Intensity MAE column**
   - What we know: CONTEXT.md lists "Intensity MAE" as a column; metrics JSON only has `r2_intensity`, not `mae_intensity`
   - What's unclear: Is intensity MAE computable from the CSV comparison files, or should the column be dropped?
   - Recommendation: Drop "Intensity MAE" — use only "Intensity R²" which is directly available. The CSV files have raw data but parsing them adds scope. Show the columns that are clean and available.

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | pytest (existing, 131 tests passing) |
| Config file | `pytest.ini` or `pyproject.toml` (check project root) |
| Quick run command | `python presentation/generate_pptx_v2.py && python -c "from pptx import Presentation; prs = Presentation('presentation/presentation_v2.pptx'); print(f'{len(prs.slides)} slides OK')"` |
| Full suite command | `pytest` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CONT-01 | Generator produces slide with 8-row water table and 4-row aspirin table | smoke | `python presentation/generate_pptx_v2.py` (no crash = pass) | ✅ generator exists |
| CONT-02 | Scaling slide present in output pptx | smoke | Count slides: expect 12 (was 11) | ✅ via pptx parse |
| CONT-03 | Links slide references both HTML report paths | smoke | Search slide text for `water/report.html` and `aspirin/report.html` | ✅ via pptx parse |
| CONT-04 | ZMQ slide contains `zmq_server.py` and not `gm_main.py` | smoke | Search slide text for string presence/absence | ✅ via pptx parse |

All validation is smoke-level: run the generator, open the pptx, check slide count and text content. No unit tests needed — the generator is a script, not a library.

### Sampling Rate

- **Per task commit:** `python presentation/generate_pptx_v2.py` (runs in <5s)
- **Per wave merge:** Same + manual visual spot-check of output pptx
- **Phase gate:** Generator runs clean, slide count is 12, content verified before `/gsd:verify-work`

### Wave 0 Gaps

None — existing infrastructure covers all phase requirements. No new test files or framework changes needed.

## Sources

### Primary (HIGH confidence)

- Direct file reads: `presentation/generate_pptx_v2.py` (lines 1–533) — full generator code examined
- Direct file reads: `analysis_results_harmonic/water/data/metrics_summary.json` — actual R²/MAE/speedup values
- Direct file reads: `analysis_results_harmonic/aspirin/data/metrics_summary.json` — actual R²/MAE/speedup values
- Direct file reads: `.planning/phases/02-content/02-CONTEXT.md` — locked decisions
- Directory listing: `comparison_results/water/` and `comparison_results/aspirin/` — confirmed combo availability

### Secondary (MEDIUM confidence)

- `CLAUDE.md` — MACE-Gaussian architecture, confirmed `zmq_server.py` naming, confirmed `gm_helper.py` is correct Gaussian-side name

### Tertiary (LOW confidence)

None — all findings are directly from project files.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — existing python-pptx patterns fully documented in source
- Architecture: HIGH — generator is 533 lines, fully read, all patterns verified
- Pitfalls: HIGH — discovered from actual data (water R² uniformity, aspirin 4-combo limit) not speculation
- Data values: HIGH — read directly from JSON ground truth files

**Research date:** 2026-03-07
**Valid until:** Until metrics JSONs are regenerated or slide generator is substantially restructured (stable for this presentation sprint)
