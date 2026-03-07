# Requirements: MACE-Gaussian Presentation (March 2026)

**Defined:** 2026-03-06
**Core Value:** A compelling, honest presentation that tells the story of why we built this, how the ZMQ injection works, and what the harmonic benchmark results show — clear enough for life scientists, rigorous enough for the research group.

## v1 Requirements

### Narrative & Structure

- [x] **NARR-01**: Slide structure reviewed and reordered to match harmonic-first narrative arc (why → how → what we found → what's next)
- [ ] **NARR-02**: Slide 10 rewritten as research journey — physics and chemistry progress, not git versions
- [x] **NARR-03**: VPT2/anharmonic slide scoped to context + "ongoing work" framing (not a deep dive)

### Content

- [ ] **CONT-01**: Harmonic results slide shows all 8 calculator combos with ranked R²/MAE table (frequencies and intensities separate)
- [ ] **CONT-02**: Scaling argument made explicit — molecule size vs speedup (water ~1×, glycine ~7–10×, aspirin ~18–29×)
- [ ] **CONT-03**: Spectrum plot images embedded from `analysis_results_harmonic/` (at least water + aspirin)
- [ ] **CONT-04**: ZMQ bridge slide has clear flow diagram showing real-time ML dipole injection into Gaussian

### Polish

- [ ] **POLS-01**: Results table color-coded by performance (green = good, yellow = ok, red = poor)
- [ ] **POLS-02**: Slide numbers added to all slides
- [ ] **POLS-03**: Font and spacing consistency pass across all slides

### Speaking Notes

- [ ] **SPKR-01**: Personal speaking notes written for each slide (~2–4 bullet points per slide, what to say and what to emphasise)

## v2 Requirements

### Nice-to-Have (post-presentation)

- **V2-01**: Mode overlap heatmap for aspirin embedded as image
- **V2-02**: Live demo slide showing CLI command + sample log output
- **V2-03**: Improved title ASCII art (molecule rendered in ASCII)
- **V2-04**: Makefile: `make presentation` regenerates + opens pptx

## Out of Scope

| Feature | Reason |
|---------|--------|
| BH3·NH3 results | Future-work molecule, not ready |
| Anharmonic deep dive | Ongoing work, mentioned briefly only |
| Word-for-word script | Personal reference notes only |
| Matplotlib-generated architecture diagram | ASCII is consistent with terminal theme |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| NARR-01 | Phase 1 — Narrative & Structure | Complete |
| NARR-02 | Phase 1 — Narrative & Structure | Pending |
| NARR-03 | Phase 1 — Narrative & Structure | Complete |
| CONT-01 | Phase 2 — Content | Pending |
| CONT-02 | Phase 2 — Content | Pending |
| CONT-03 | Phase 2 — Content | Pending |
| CONT-04 | Phase 2 — Content | Pending |
| POLS-01 | Phase 3 — Polish & Speaking Notes | Pending |
| POLS-02 | Phase 3 — Polish & Speaking Notes | Pending |
| POLS-03 | Phase 3 — Polish & Speaking Notes | Pending |
| SPKR-01 | Phase 3 — Polish & Speaking Notes | Pending |

**Coverage:**
- v1 requirements: 11 total
- Mapped to phases: 11
- Unmapped: 0 ✓

---
*Requirements defined: 2026-03-06*
*Last updated: 2026-03-06 after roadmap creation*
