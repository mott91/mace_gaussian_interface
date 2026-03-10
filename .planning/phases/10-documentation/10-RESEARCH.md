# Phase 10: Documentation - Research

**Researched:** 2026-02-26
**Domain:** Technical writing — quickstart guide, worked example, CLI help text, thesis methods section
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Thesis methods section**
- Audience: thesis committee — physical/computational chemists who know DFT and vibrational
  spectroscopy. No need to explain basics; focus on methodology choices and what's novel.
- Format: `docs/methods.md` (Markdown, easy to copy into thesis, version controlled)
- Cover in depth:
  1. MACE architecture & model selection (why MACE-MP-0 and MACE-OMOL-0, how chosen)
  2. ZMQ injection mechanism (the novel software bridge — how dipole derivatives are
     injected into Gaussian mid-calculation)
  3. VPT2 anharmonic treatment (why anharmonic over harmonic, how Gaussian handles it)
  4. General software usage (how the pipeline is invoked, what it does end-to-end)
- Espaloma dipole model: mention as one of the dipole backends, no deep coverage needed

**Worked example**
- Molecule: water (H2O) — existing results in `comparison_results/water/` and
  `analysis_results/` can be used as reference output
- Expected output to commit: frequencies table, JSON results file, HTML analysis report
- Location: `docs/examples/water/` directory (water.xyz + expected output + README instructions)

**Quickstart guide**
- Target audience: lab members / thesis committee — Gaussian 16 access assumed.
  No Gaussian-free demo path needed.
- Goal: clone → install → run water → view HTML report (full pipeline from scratch)
- Location: `docs/quickstart.md` (standalone file, linked from README)

**CLI documentation**
- Claude's discretion on format (inline `--help` text vs separate reference page)
- Success criterion: all subcommands and flags have complete, accurate help text

### Claude's Discretion

- CLI documentation format: inline `--help` text vs separate reference page

### Deferred Ideas (OUT OF SCOPE)

None — discussion stayed within phase scope.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| DOC-01 | README updated with complete installation steps (including custom MACE packages) | Phase 9 (09-03) already implemented; Phase 10 should verify and link from quickstart |
| DOC-02 | Quickstart guide: clone → install → run water → view results | `docs/quickstart.md` — reference README installation, add run + analysis steps |
| DOC-03 | Worked example with expected output committed to repo | `docs/examples/water/` — existing `comparison_results/water/` and `analysis_results/water/` are the source of truth |
| DOC-04 | CLI help text complete for all commands and options | `mace_gaussian/cli.py` — `run`, `list`, `compare`, `diagnose` commands; `compare` and `export` are placeholders that need accurate help text |
| DOC-05 | Method description suitable for citing in thesis methods section | `docs/methods.md` — covers MACE models, ZMQ bridge, VPT2 anharmonic treatment, pipeline overview |
</phase_requirements>

---

## Summary

Phase 10 is a pure documentation phase — no new software features, no refactoring. All five requirements involve writing Markdown files and improving CLI help text. The technical infrastructure is complete from Phases 1–9; the task is to document what exists clearly enough for (a) a lab member reproducing the calculation from scratch and (b) a thesis committee reader understanding the methodology.

Three files need to be created from scratch: `docs/quickstart.md` (DOC-02), `docs/examples/water/` directory with worked example (DOC-03), and `docs/methods.md` (DOC-05). One existing file needs verification and possible extension: `README.md` for DOC-01 (Phase 9 already wrote the Installation section). CLI help text (DOC-04) is inline in `mace_gaussian/cli.py` via Click docstrings — some subcommands (`compare`, `export`) are marked "COMING SOON" placeholders that need accurate descriptions.

**Primary recommendation:** Write all four deliverables as plain Markdown files. Keep CLI docs inline in Click decorators (no separate reference page needed). Commit expected output files directly — do not attempt to re-run Gaussian calculations to generate fresh output.

---

## Standard Stack

### Core

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Markdown | N/A | All documentation files | Version-controlled, copy-paste into LaTeX/Word thesis |
| Click | already installed | CLI help text generation via docstrings | Already the project CLI framework |

No external documentation generators (MkDocs, Sphinx) — deferred to v2 (DIST-03 is out of scope).

### What Already Exists (Do Not Duplicate)

| File | Status | Phase 10 Action |
|------|--------|-----------------|
| `README.md` — Installation section | Complete (Phase 9 Plan 03) | Verify contains correct steps; add link to `docs/quickstart.md` |
| `docs/ARCHITECTURE.md` | Complete | Reference from methods.md where appropriate |
| `docs/DEVELOPMENT.md` | Complete | Reference from methods.md |
| `docs/CONVENTIONS.md` | Complete | No action needed |
| `comparison_results/water/` | Complete | Copy key files into `docs/examples/water/` |
| `analysis_results/water/` | Complete | Copy `report.html` and plots into `docs/examples/water/` |

### Supporting Reference Material (For methods.md Content)

| Source | Content Type | Use In |
|--------|-------------|--------|
| `docs/ARCHITECTURE.md` | ZMQ bridge, dual MACE architecture | methods.md §ZMQ injection |
| `mace_gaussian/cli.py` | All CLI flags and defaults | quickstart.md + DOC-04 |
| `comparison_results/water/mace_mp_espaloma/results.json` | Actual frequency output | worked example |
| `analysis_results/water/data/metrics_summary.json` | Actual R², RMSE, speedup values | worked example expected output |
| `comparison_results/water/geometry_opt/results.json` | Geometry opt results | worked example |

---

## Architecture Patterns

### Document Structure Pattern

```
docs/
├── ARCHITECTURE.md      # already exists — technical internals
├── CONVENTIONS.md       # already exists — code style
├── DEVELOPMENT.md       # already exists — developer guide
├── quickstart.md        # NEW (DOC-02) — end-to-end user guide
├── methods.md           # NEW (DOC-05) — thesis methods section
└── examples/
    └── water/
        ├── README.md    # NEW — how to reproduce this example
        ├── water.xyz    # input molecule (already at project root)
        ├── expected_output/
        │   ├── geometry_opt_results.json
        │   ├── mace_omol_mace_ml_results.json
        │   └── metrics_summary.json
        └── expected_plots/
            ├── spectrum_combined.png
            └── regression_combined.png
```

### Pattern 1: Quickstart Reference-Don't-Duplicate

**What:** `docs/quickstart.md` links to the README Installation section rather than repeating it. Adds only what README lacks: how to actually run the calculation and view results.

**When to use:** The Installation section already exists and is correct (Phase 9 verified it). Duplicating would create a maintenance hazard.

**Structure:**
```markdown
# Quickstart

## 1. Install
See [Installation](../README.md#installation) for the full setup procedure.
(Brief summary: clone → conda env → pip install 3 packages → verify with diagnose)

## 2. Run a Calculation

```bash
cd mace_gaussian
python cli.py run water.xyz --skip-dft-baseline
```

## 3. Analyze Results

```bash
python run_analysis.py water
```

## 4. View the Report
Open analysis_results/water/report.html in a browser.
```

### Pattern 2: Self-Contained Worked Example

**What:** `docs/examples/water/` is self-contained — someone can `cd` there and follow `README.md` without reading anything else.

**Structure:**
```
docs/examples/water/
├── README.md            # instructions specific to this example
├── water.xyz            # copy of the input file
└── expected_output/     # pre-committed reference output
    ├── metrics_summary.json        # R² and RMSE numbers
    ├── geometry_results.json       # geometry opt output
    └── mace_omol_mace_ml_results.json  # one representative ML result
```

**Key point:** Do NOT commit `.log`, `.chk`, `.fchk` files (large, binary, Gaussian-specific). Commit only JSON results and PNG plots.

### Pattern 3: Thesis Methods Section Structure

**What:** `docs/methods.md` reads as a thesis methods section — passive voice, past tense, no step-by-step instructions.

**Audience note:** Thesis committee knows DFT and vibrational spectroscopy. Skip basics. Cover what is novel and what choices were made.

**Sections to cover (in order):**
1. **Computational Framework Overview** — pipeline from geometry optimization to anharmonic frequencies; what tools are combined and why
2. **MACE Potential Energy Models** — MACE architecture briefly, MACE-MP-0 and MACE-OMOL-0 selection rationale (breadth vs. training set)
3. **Dipole Derivative Calculation** — Espaloma (charge model) and MACE4IR (custom dipole model); why two backends compared
4. **ZMQ Injection Mechanism** — the novel software bridge; how Gaussian's External keyword is used; IPC socket; why this approach enables VPT2 with ML dipoles
5. **VPT2 Anharmonic Treatment** — why anharmonic (overtones, combination bands visible experimentally); Gaussian's VPT2 implementation; what anharmonic corrections capture that harmonic misses
6. **Mode Matching** — eigenvector dot product method for comparing DFT vs ML vibrational modes without relying on frequency ordering
7. **Analysis and Validation** — KDE spectral broadening, R² and RMSE metrics, regression plots

### Pattern 4: CLI Help Text Via Click Docstrings

**What:** All CLI help text lives directly in Click command docstrings and `help=` parameters. No separate reference page.

**Why:** Click auto-generates `--help` output from docstrings. Single source of truth. Already implemented for `run`, `list`, `diagnose`.

**Current gaps in CLI:**
- `compare` command: docstring says "COMING SOON" and "PLACEHOLDER" — needs accurate description of actual state (not implemented, deferred to v2)
- `export` command: same issue — placeholder text in docstring
- `diagnose` command: help text is accurate but missing mention of what Gaussian/CUDA checks show
- `run` command: missing mention of `--output-dir` in the docstring example section

**Recommended approach:** Keep `compare` and `export` as-is but remove misleading "planned features" lists. Replace with honest: "This command is not yet implemented. Use `python run_analysis.py <molecule>` for analysis."

### Anti-Patterns to Avoid

- **Duplicating installation steps:** quickstart.md must reference README, not copy the 5 install steps again
- **Committing large binary output:** `.chk` files are binary, `.log` files can be many MB — commit only JSON + PNG
- **Writing methods.md as a user guide:** passive voice, past tense, methodology framing — not "to run X, do Y"
- **Over-specifying expected output:** frequency values from JSON are enough; do not commit the full `.log` file
- **Creating a separate CLI reference page:** Click's `--help` is the reference; a separate page would go stale

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CLI help output | Separate `docs/cli-reference.md` | Click `--help` auto-generation | Already correct, single source of truth |
| LaTeX thesis section | `.tex` file | `docs/methods.md` (Markdown) | User copies to thesis manually; Markdown is version-controlled |
| Documentation site | MkDocs setup | Plain Markdown files | MkDocs deferred to v2 (DIST-03 out of scope) |

**Key insight:** This project is thesis-scale. Plain Markdown files in `docs/` are the right medium. No tooling overhead.

---

## Common Pitfalls

### Pitfall 1: DOC-01 Already Done by Phase 9

**What goes wrong:** Implementing DOC-01 from scratch, duplicating the Installation section that Phase 9 already wrote correctly.

**Why it happens:** DOC-01 appears in the requirements list without noting Phase 9 completed it.

**How to avoid:** DOC-01 verification task only: run the automated check from Phase 9 (`grep "pip install -e mace_ML_pkg" README.md`) and confirm it passes. Then add a link from the quickstart to the README Installation section. No rewriting.

**Warning signs:** If a plan task says "write README Installation section," it is redundant — Phase 9 is done.

### Pitfall 2: Committing Gaussian Binary Output

**What goes wrong:** Committing `.chk`, `.fchk`, or full `.log` files as "expected output" for the worked example.

**Why it happens:** Those files exist in `comparison_results/water/` and it seems natural to include them.

**How to avoid:** The worked example expected output should be: `results.json` files (small, human-readable) and PNG plots (small, illustrative). Never `.chk` (binary, hundreds of MB), never full `.log` files (large text, system-path-dependent).

**Warning signs:** `git status` shows `.chk` or `.log` files staged for commit.

### Pitfall 3: Methods.md Written as User Guide

**What goes wrong:** `docs/methods.md` reads as step-by-step instructions ("Run `python cli.py run water.xyz`").

**Why it happens:** It's easier to write how-to than methodology rationale.

**How to avoid:** Every sentence in methods.md should be passive voice and past tense: "Molecular geometries were optimized using..." not "To optimize, run..." The audience knows how to use software — they want to know *why* these choices were made.

**Warning signs:** Imperative verbs (run, open, install) appearing in methods.md.

### Pitfall 4: Worked Example README Requires Gaussian

**What goes wrong:** `docs/examples/water/README.md` says "to reproduce: run `python cli.py run water.xyz`" without noting that Gaussian 16 is required and takes ~30 minutes.

**Why it happens:** The quickstart target audience (lab members with Gaussian) differs from someone who just wants to see the output format.

**How to avoid:** The worked example README should have two sections: (a) "Expected output" — look at these committed files without re-running anything; (b) "To reproduce" — requirements (Gaussian 16, GPU, ~45 min) and exact commands.

### Pitfall 5: CLI Docs Describe Placeholder Commands as Real

**What goes wrong:** `compare` and `export` have help text listing "planned features" as if they exist.

**Why it happens:** The placeholders have `click.echo("Coming soon")` implementations that could mislead a reader.

**How to avoid:** Remove misleading "Planned features" lists from `compare` and `export` docstrings. Replace with accurate description: "Not yet implemented. For analysis, use `python run_analysis.py <molecule>`."

---

## Code Examples

### Current CLI Commands (Accurate State)

Verified from `mace_gaussian/cli.py` — these are the implemented commands:

```bash
# Show help
python cli.py --help

# Run full calculation workflow (requires Gaussian 16 and GPU)
python cli.py run water.xyz
python cli.py run water.xyz --skip-dft-baseline
python cli.py run water.xyz --energy-calculators mace_mp,mace_omol --dipole-calculators espaloma,mace_ml
python cli.py run water.xyz --optimization-calculator mace_mp --output-dir my_results

# List results
python cli.py list
python cli.py list water

# Diagnostic check
python cli.py diagnose

# Analysis (separate entry points, not CLI subcommands)
python run_analysis.py water           # Anharmonic analysis
python run_analysis_harmonic.py water  # Harmonic-only analysis
```

### Expected Output Structure for Worked Example

From `comparison_results/water/` and `analysis_results/water/` (verified on disk):

```
docs/examples/water/
├── water.xyz
│   # 3 atoms: O + 2H
├── expected_output/
│   ├── geometry_opt_results.json
│   │   # converged: true, num_steps: 8, runtime_s: ~4.2s, calculator: mace_omol
│   ├── mace_omol_mace_ml_results.json
│   │   # 3 harmonic, 3 anharmonic fundamentals, overtones, combination bands
│   └── analysis_metrics_summary.json
│       # mace_omol_espaloma: R²=0.9997, RMSE=44 cm⁻¹
│       # mace_omol_mace_ml:  R²=0.9997, RMSE=44 cm⁻¹
└── expected_plots/
    ├── spectrum_combined.png
    └── regression_combined.png
```

### Thesis Methods Section Key Content Points

Verified facts from the codebase to include in `docs/methods.md`:

**MACE architecture:** Equivariant GNN with E(3) symmetry. MACE-MP-0 trained on ~150k diverse molecules from the Materials Project (broad coverage). MACE-OMOL-0 trained on Open Molecules 2025 dataset (organic molecules, better for drug-like structures). Choice of both allows comparison across training distribution coverage.

**ZMQ injection mechanism** (from `docs/ARCHITECTURE.md` and `mace_gaussian/gaussian/`):
- Gaussian's `External=` keyword causes Gaussian to launch `gm_helper.py` as a subprocess at each VPT2 derivative evaluation step
- `gm_helper.py` connects to a ZMQ IPC socket bound by the main Python process
- Main process calculates dipole derivatives using ML model, writes to Gaussian's outfile format
- Gaussian receives dipole data and continues the anharmonic calculation
- This enables VPT2 anharmonic treatment with ML-derived dipole surfaces, not available natively from ML potentials

**VPT2:** Gaussian's second-order vibrational perturbation theory. Uses cubic and quartic force constants derived numerically. Produces anharmonic fundamental frequencies, overtones (2ν), and combination bands (ν_i + ν_j). B3LYP/6-31G(d,p) used as DFT reference.

**Mode matching:** Eigenvector dot product overlap between DFT and ML vibrational modes. Avoids incorrect mode assignments from frequency-ordering alone (critical when ML and DFT frequency orderings differ).

**Dipole backends compared:**
- Espaloma: ML partial charge model → dipole from charge distribution
- MACE4IR: custom MACE fork trained specifically for dipole prediction (78 element coverage)

---

## State of the Art

| Old Approach | Current Approach | Impact for Phase 10 |
|--------------|------------------|---------------------|
| Documentation as afterthought | Docs-as-code in version control | Markdown in `docs/` is the right pattern |
| Separate CLI man page | Click `--help` docstrings | No separate reference page needed |
| MkDocs/Sphinx for research code | Plain Markdown for thesis-scale | MkDocs deferred to v2 per project decision |

---

## Open Questions

1. **Should `compare` and `export` be removed or kept as stubs?**
   - What we know: Both print "COMING SOON" and are not functional
   - What's unclear: Whether the user wants them visible in `--help` output (signals future plans) or removed (reduces confusion)
   - Recommendation: Keep them in `--help` but remove the misleading "Planned features" bullet list. Change docstrings to: "Not yet implemented."

2. **Does `docs/examples/water/` need the harmonic analysis output too?**
   - What we know: Both `analysis_results/water/` and `analysis_results_harmonic/water/` exist with full output
   - What's unclear: CONTEXT.md says "frequencies table, JSON results file, HTML analysis report" — anharmonic or harmonic?
   - Recommendation: Include anharmonic output (primary use case) as the worked example. Mention harmonic analysis exists via `run_analysis_harmonic.py` in the quickstart.

3. **Should the quickstart run with or without `--skip-dft-baseline`?**
   - What we know: DFT baseline takes ~15-30 min for water; ML only takes ~5 min; Gaussian required for both
   - What's unclear: CONTEXT.md says "run water → view HTML report" — full pipeline or ML-only?
   - Recommendation: Quickstart should show `--skip-dft-baseline` as the "fast path" for first-time verification, with a note that the full run includes DFT comparison.

---

## Sources

### Primary (HIGH confidence)

- Direct codebase inspection — `mace_gaussian/cli.py`, `docs/ARCHITECTURE.md`, `docs/DEVELOPMENT.md`
- `comparison_results/water/` — all JSON result files verified on disk
- `analysis_results/water/data/metrics_summary.json` — actual R²/RMSE numbers
- `.planning/phases/09-ci-cd-distribution-prep/09-03-PLAN.md` — confirms DOC-01 was implemented by Phase 9
- `README.md` — current state of installation documentation

### Secondary (MEDIUM confidence)

- `docs/CONVENTIONS.md` — confirmed passive-voice and methodology-framing preferences align with thesis writing style

### Tertiary (LOW confidence)

- None — all claims in this research are based on direct codebase inspection

---

## Metadata

**Confidence breakdown:**
- Document content mapping: HIGH — all source files inspected directly
- CLI current state: HIGH — `mace_gaussian/cli.py` read in full
- Expected output format: HIGH — JSON files in `comparison_results/water/` verified on disk
- Methods.md technical content: HIGH — architecture docs and source confirmed key claims
- DOC-01 status (already done): HIGH — Phase 9 Plan 03 and README.md both verified

**Research date:** 2026-02-26
**Valid until:** Stable — documentation phase, no external dependencies to track
