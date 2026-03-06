# Roadmap: MACE-Gaussian Presentation (March 2026)

## Overview

Three phases to take the existing 12-slide deck from rough structure to presentation-ready. Phase 1 locks the narrative shape — no point polishing slides before the story is right. Phase 2 fills the content gaps that make the argument credible (real results, real plots, real diagram). Phase 3 makes it clean and speakable in 15 minutes.

## Phases

- [ ] **Phase 1: Narrative & Structure** - Lock the story arc and slide order before any content work
- [ ] **Phase 2: Content** - Add the substance: results table, plots, ZMQ diagram, scaling argument
- [ ] **Phase 3: Polish & Speaking Notes** - Make it presentable and speakable

## Phase Details

### Phase 1: Narrative & Structure
**Goal**: The slide deck tells a coherent harmonic-first story — why, how, what we found, what's next — with each slide in the right place
**Depends on**: Nothing (first phase)
**Requirements**: NARR-01, NARR-02, NARR-03
**Success Criteria** (what must be TRUE):
  1. Slides flow in the narrative arc: motivation → IR basics → ZMQ bridge → results → ongoing/next
  2. Slide 10 reads as a research journey (physics and chemistry progress) not a git version history
  3. The VPT2/anharmonic slide frames anharmonicity as context and positions it as ongoing work, not a deep dive
**Plans**: 2 plans

Plans:
- [ ] 01-01-PLAN.md — Reframe IR theory, cut VPT2 slide, add results slides (ls + table), fix slide order and output path
- [ ] 01-02-PLAN.md — Rewrite status slide as research journey; tighten motivation slide

### Phase 2: Content
**Goal**: Every claim in the deck is backed by a visible artifact — the results table is complete, plots are embedded, the ZMQ diagram is clear, and the scaling argument is explicit
**Depends on**: Phase 1
**Requirements**: CONT-01, CONT-02, CONT-03, CONT-04
**Success Criteria** (what must be TRUE):
  1. Results slide shows all 8 calculator combos ranked by R²/MAE with frequencies and intensities in separate columns
  2. A molecule-size vs speedup table or visual exists (water ~1×, glycine ~7–10×, aspirin ~18–29×)
  3. Spectrum comparison plots from `analysis_results_harmonic/` are embedded for at least water and aspirin
  4. ZMQ bridge slide contains a flow diagram showing real-time ML dipole injection at the right Gaussian hook point
**Plans**: TBD

### Phase 3: Polish & Speaking Notes
**Goal**: The deck is visually consistent and the presenter has personal notes for every slide
**Depends on**: Phase 2
**Requirements**: POLS-01, POLS-02, POLS-03, SPKR-01
**Success Criteria** (what must be TRUE):
  1. Results table cells are color-coded: green for strong performance, yellow for moderate, red for poor
  2. Every slide has a visible slide number
  3. Font, spacing, and accent colors are consistent across all 12 slides (no visual outliers)
  4. Each slide has 2–4 bullet speaking notes covering what to say and what to emphasise
**Plans**: TBD

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Narrative & Structure | 0/2 | Not started | - |
| 2. Content | 0/? | Not started | - |
| 3. Polish & Speaking Notes | 0/? | Not started | - |
