# Phase 23: Anharmonic Pipeline & Report Overhaul - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-04-05
**Phase:** 23-anharmonic-pipeline-report-overhaul
**Areas discussed:** Report structure & narrative flow, Plot quality & styling, Feature integration strategy, Todo folding, Harmonic/anharmonic coherence

---

## Report Structure & Narrative Flow

| Option | Description | Selected |
|--------|-------------|----------|
| Keep current order, polish it | Overview → combined plots → per-method detail → experimental → summary. Just improve card content and styling. | |
| Executive summary first | Start with a 1-page executive summary (best method, key findings, speedup), then drill into details. Reader gets the punchline immediately. | ✓ |
| Experimental-anchored | Lead with experimental spectrum, then show how each method compares to it. | |

**User's choice:** Executive summary first
**Notes:** None

### Summary Card Content

| Option | Description | Selected |
|--------|-------------|----------|
| Metrics + recommendation | R², RMSE, speedup, best method highlighted, auto-generated verdict | ✓ |
| Metrics only | Numbers only, no verdict | |
| Visual scorecard | Traffic-light colored badges + numbers | |

**User's choice:** Metrics + recommendation

### Harmonic/Anharmonic Report Coherence

| Option | Description | Selected |
|--------|-------------|----------|
| Shared template, mode-aware | Same CSS, layout, executive summary. Anharmonic adds overtone/combination sections. Both feel like the same product. | ✓ |
| Anharmonic extends harmonic | Anharmonic report includes harmonic results as a section, then adds anharmonic analysis. | |
| Separate but styled alike | Two independent reports with matching CSS. No cross-referencing. | |

**User's choice:** Shared template, mode-aware
**Notes:** User specifically wanted to "make the two reports the same style-wise" and "make the anharm report an extension of the harm one"

---

## Plot Quality & Styling

| Option | Description | Selected |
|--------|-------------|----------|
| Plotly interactive | Zoomable, hoverable spectra as Plotly.js divs. Adds ~3MB. Supports PNG export. | |
| Static matplotlib, higher quality | Keep matplotlib, bump to 300 DPI. Thesis-ready, no JS. | |
| Both: Plotly in report, matplotlib for export | Interactive Plotly in report + separate 300 DPI matplotlib PNGs | |

**User's initial concern:** Worried about thesis figures — needs printable, incredible-looking graphs.
**Resolution:** Plotly for interactive report + structured data export (JSON/CSV) for future thesis figure generation. Thesis figures are a separate future phase.

| Option | Description | Selected |
|--------|-------------|----------|
| Plotly + data export | Interactive Plotly in report. Clean data files for future thesis figures. | ✓ |
| Plotly + basic matplotlib too | Both rendering engines now | |
| Static matplotlib only | Skip Plotly | |

**User's choice:** Plotly + data export

---

## Feature Integration Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Inline in each method section | Each section: spectrum (with exp.), regression, heatmap, timing, degenerate notes. Everything in one place. | ✓ |
| Layered: overview then detail | Summary has timing table + exp. agreement. Per-method has detailed plots. | |
| Tabbed per-method | Each method in a tab. | |

**User's choice:** Inline in each method section

### Experimental Overlay Placement

| Option | Description | Selected |
|--------|-------------|----------|
| On every method's spectrum | DFT + ML + experimental on same plot per method | ✓ |
| Once in overview + per-method | Overview combined plot + per-method | |
| Only in overview | Experimental only in combined plot | |

**User's choice:** On every method's spectrum

---

## Todo Folding

| Todo | Score | Folded? |
|------|-------|---------|
| Interactive HTML spectrum viewer with Plotly | 0.9 | ✓ (absorbed by Plotly decision) |
| JCAMP-DX spectral data export | 0.9 | No — separate utility |
| Peak assignment confidence scores | 0.9 | No — adds scope |
| Normalize intensity regression + fit line | 0.6 | No — not blocking |

---

## Visual Design

User wants to run `/gsd-ui-phase 23` after context capture for a formal visual design contract before planning.

## Claude's Discretion

- Plotly layout configuration
- Executive summary verdict generation logic
- Data export format details
- Handling molecules with no experimental data

## Deferred Ideas

- Thesis figure generator (separate future phase)
- JCAMP-DX export
- Peak confidence scores
- Intensity normalization
- Conformer-aware spectra
