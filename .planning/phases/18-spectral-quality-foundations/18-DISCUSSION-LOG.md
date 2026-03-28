# Phase 18: Spectral Quality Foundations - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-03-28
**Phase:** 18-spectral-quality-foundations
**Areas discussed:** Spectrum display, Zero-intensity threshold, FWHM & broadening, mace_polar dipole

---

## Spectrum Display

| Option | Description | Selected |
|--------|-------------|----------|
| Overlay | Stick lines and broadened curve on same axes | |
| Side-by-side panels | Left stick, right broadened | |
| Stacked panels | Top stick, bottom broadened | |

**User's choice:** None of the above — user said "don't show both, just Lorentzian for now"
**Notes:** User explicitly rejected stick spectrum display. Only Lorentzian broadened curve needed.

**Follow-up: Layout**

| Option | Description | Selected |
|--------|-------------|----------|
| Same layout (DFT above, ML below) | Replace Gaussian with Lorentzian, keep offset layout | ✓ |
| Overlay on same axis | Both on same axes, no offset | |

**User's choice:** Same layout
**Notes:** Minimal visual disruption, just swap broadening kernel

---

## Zero-Intensity Threshold

| Option | Description | Selected |
|--------|-------------|----------|
| DFT intensity as authority | Filter when DFT < 0.1 km/mol | ✓ |
| Either below threshold | Exclude if either DFT or ML < 0.1 | |
| Both below threshold | Exclude only if both < 0.1 | |

**User's choice:** DFT intensity as authority (after explanation)
**Notes:** User initially unsure ("oof, I don't know how to handle this. What is this really about?"). After explanation of why inactive modes distort regression statistics and that filtering only affects metrics (not spectra), user agreed with DFT-based filtering. Key concern: "this should be presented as an alternative to DFT" — clarified that filtering is benchmarking-only, not data loss.

---

## FWHM & Broadening

| Option | Description | Selected |
|--------|-------------|----------|
| Analysis script arg | Pass --fwhm to analysis scripts | ✓ |
| In results JSON | Store per-molecule in results | |
| You decide | Claude picks | |

**User's choice:** Analysis script argument

**Follow-up: Default value**

| Option | Description | Selected |
|--------|-------------|----------|
| 10 cm⁻¹ | Standard in literature, matches requirement | ✓ |
| 8 cm⁻¹ | Current value, sharper peaks | |

**User's choice:** 10 cm⁻¹

---

## mace_polar Dipole

| Option | Description | Selected |
|--------|-------------|----------|
| Quick investigation + skip logic | ~1hr investigation, then skip if broken | |
| Document-only | Just add error message and skip logic | |
| Deep investigation | Thorough debugging | |

**User's choice:** Defer entirely — "I want to defer this to a later time, put it in todo and move on"
**Notes:** PIPE-02 removed from Phase 18 scope by user request

---

## Claude's Discretion

- Lorentzian implementation details (normalization, numerical cutoff)
- How to label/annotate filtered modes in report text
- Whether to mention filtering count in report methodology section

## Deferred Ideas

- Stick spectrum display (SPEC-03 partial) — deferred by user, add later
- PIPE-02 mace_polar dipole investigation — deferred by user to standalone todo
- Interactive Plotly viewer — from backlog, not in scope
