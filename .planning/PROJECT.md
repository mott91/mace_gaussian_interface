# MACE-Gaussian Presentation (March 2026)

## What This Is

A 15-minute research group presentation on the MACE-Gaussian interface — a tool that injects ML-calculated dipole derivatives into Gaussian 16 in real-time via ZMQ to produce fast, DFT-quality IR spectra. Built programmatically in Python (python-pptx) with a terminal dark theme. The source of truth is `presentation/generate_pptx_v2.py`; the `.pptx` is a build artifact.

**Audience:** Research group + life science theoretical chemists (good background, not specialist)
**Length:** ~15 minutes
**Branch:** `presentation/v2-march2026`

## Core Value

A compelling, honest presentation that tells the story of why we built this, how the ZMQ injection works, and what the harmonic benchmark results actually show — clear enough for life scientists, rigorous enough for the research group.

## Narrative Arc

1. **Why:** We want to match DFT IR spectra fast — both frequencies and intensities
2. **How:** We built a ZMQ bridge that injects ML dipoles into Gaussian in real-time (the novel engineering)
3. **What we found:** Harmonic benchmark across molecules — best combo, scaling wins, honest limits
4. **What's next:** Anharmonic treatment is ongoing

## Key Results (Harmonic)

- **Frequencies:** Energy model dominates. mace_omol/mace_anicc → MAE ~10–27 cm⁻¹ (R² > 0.9999). mace_mp → MAE ~100+ cm⁻¹. Dipole model irrelevant for frequencies.
- **Intensities:** Dipole model dominates. mace_ml → R² 0.72–0.999 (molecule-dependent). espaloma → R² 0.09–0.74 (much worse, especially for larger molecules).
- **Winner:** mace_omol + mace_ml
- **Scaling:** ~1× (water) → ~7–10× (glycine) → ~18–29× (aspirin). Grows with molecule size.
- **Anharmonic:** Explored for water/acoh — ongoing, not the focus of this presentation.

## Deep Dives (2 slides worth of depth)

1. **ZMQ bridge** — the real-time ML injection mechanism; the novel engineering contribution
2. **Results** — systematic harmonic comparison, model ranking, scaling argument

## Current Slide Structure (12 slides, v2)

| # | Command | Content |
|---|---------|---------|
| 0 | `$ ./presentation.sh` | Neofetch-style title |
| 1 | `$ cat motivation.md` | Problem / solution / impact |
| 2 | `$ man ir_spectroscopy` | IR basics, harmonic limit |
| 3 | `$ cat vpt2_theory.md` | Anharmonicity, VPT2, dipole gap |
| 4 | `$ python3 main.py` | Full pipeline ASCII diagram |
| 5 | `$ cat dipole_methods.md` | MACE4IR vs Espaloma + energy model table |
| 6 | `$ cat zmq_bridge.md` | ZMQ flow diagram + engineering |
| 7 | `$ cat mode_matching.md` | Eigenvector dot product, water example |
| 8 | `$ cd results/water/` | 8 combos, R²/MAE/RMSE table |
| 9 | `$ cd results/aspirin/` | 51 modes, scaling, 18× speedup |
| 10 | `$ git log --graph` | Research journey — physics/chemistry/progress |
| 11 | `$ ./questions.sh` | Questions |

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] Workshopped narrative reviewed and slide structure confirmed
- [ ] Slides updated/created to match narrative (harmonic focus, ZMQ deep dive, results)
- [ ] Embedded plot images added (spectrum comparison, regression — from `analysis_results_harmonic/`)
- [ ] Scaling argument visualised in results slides
- [ ] Slide 10 rewritten as research journey (physics/chemistry/progress, not git versions)
- [ ] Speaking notes written for each slide (personal reference)
- [ ] Final polish pass: font consistency, color-coded tables, slide numbers

### Out of Scope

- BH3·NH3 results — future-work molecule, not ready
- Anharmonic deep dive — mentioned briefly as ongoing, not a focus slide
- Exact word-for-word script — personal reference notes only
- Mode overlap heatmap for aspirin — nice-to-have, only if time allows
- Live demo slide — too risky for a 15-minute slot

## Constraints

- **Timeline:** This week — days away
- **Build:** `micromamba activate mace4ir_v2 && python presentation/generate_pptx_v2.py`
- **Toolchain:** python-pptx only; no external design tools
- **Design system:** Terminal Dark (Consolas, #0D1117 bg, #58A6FF accent, #3FB950 green, #D29922 yellow)
- **Plot source:** `analysis_results_harmonic/` for harmonic plots; `analysis_results/` for anharmonic

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Harmonic-first narrative | Results are solid and complete; anharmonic is ongoing | — Pending |
| ZMQ bridge as technical deep dive | Novel engineering contribution, surprising to any audience | — Pending |
| Model comparison in results, not as narrative spine | Energy→freq, dipole→intensity is physically expected; ranking is the finding | — Pending |
| BH3·NH3 excluded | Future-work molecule, not ready to present | — Pending |

---
*Last updated: 2026-03-06 after initialization*
