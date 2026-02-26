# Phase 10: Documentation - Context

**Gathered:** 2026-02-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Provide complete documentation for installation, usage, and method description. Covers:
a quickstart guide (end-to-end from clone to results), a worked example with committed
expected output, complete CLI help text, and a thesis methods section. No new software
features are in scope — this phase documents what exists.

</domain>

<decisions>
## Implementation Decisions

### Thesis methods section
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

### Worked example
- Molecule: water (H2O) — existing results in `comparison_results/water/` and
  `analysis_results/` can be used as reference output
- Expected output to commit: frequencies table, JSON results file, HTML analysis report
- Location: `docs/examples/water/` directory (water.xyz + expected output + README instructions)

### Quickstart guide
- Target audience: lab members / thesis committee — Gaussian 16 access assumed.
  No Gaussian-free demo path needed.
- Goal: clone → install → run water → view HTML report (full pipeline from scratch)
- Location: `docs/quickstart.md` (standalone file, linked from README)

### CLI documentation
- Claude's discretion on format (inline `--help` text vs separate reference page)
- Success criterion: all subcommands and flags have complete, accurate help text

</decisions>

<specifics>
## Specific Ideas

- Installation steps were already handled in Phase 9 (README Installation section with
  explicit `pip install` steps). Quickstart should reference or extend that section,
  not duplicate it.
- The worked example should be self-contained: someone can `cd docs/examples/water/`
  and follow README instructions there.
- Thesis methods section should read as a methods section, not a user guide — passive
  voice, past tense, methodology framing.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 10-documentation*
*Context gathered: 2026-02-26*
