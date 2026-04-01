# Phase 21: NIST Experimental Overlay - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-04-01
**Phase:** 21-nist-experimental-overlay
**Areas discussed:** Data sourcing & caching, Spectral overlay visuals, Peak matching strategy, Integration into pipeline

---

## Data Sourcing & Caching

### Fetch method

| Option | Description | Selected |
|--------|-------------|----------|
| nistchempy library | Python wrapper for NIST WebBook. Handles CAS lookup and JCAMP-DX download. | ✓ |
| Direct HTTP + scraping | Fetch JCAMP-DX files directly from NIST WebBook URLs. No extra dependency. | |
| Both with fallback | Try nistchempy first, fall back to direct HTTP if it fails. | |

**User's choice:** nistchempy library
**Notes:** Already identified as a dependency in STATE.md.

### Cache location

| Option | Description | Selected |
|--------|-------------|----------|
| Per-molecule results dir | Store in comparison_results/{molecule}/experimental/ | ✓ |
| Central cache directory | Store in ~/.mace_gaussian/nist_cache/ | |
| Project-level cache | Store in comparison_results/.nist_cache/ | |

**User's choice:** Per-molecule results dir
**Notes:** None

### Phase preference (gas vs condensed)

| Option | Description | Selected |
|--------|-------------|----------|
| Gas-phase only | Skip molecules with only condensed-phase data | ✓ |
| Gas preferred, condensed fallback | Use condensed-phase with warning if no gas-phase | |
| Fetch both, let user choose | Download all available, add selection flag | |

**User's choice:** Gas-phase only
**Notes:** None

---

## Spectral Overlay Visuals

### Overlay style

| Option | Description | Selected |
|--------|-------------|----------|
| Third trace on same plot | Add experimental as third line on existing offset layout | ✓ |
| Separate panel below | Third panel for experimental data | |
| Side-by-side figure | New separate figure | |

**User's choice:** Third trace on same plot
**Notes:** None

### Unit handling

| Option | Description | Selected |
|--------|-------------|----------|
| Normalize both to 0-1 | Both spectra normalized for visual comparison | ✓ |
| Convert transmittance to absorbance | Apply -log10(T/100) to experimental | |
| Dual y-axes | Separate axes for each unit system | |

**User's choice:** Normalize both to 0-1
**Notes:** None

### Color/style

| Option | Description | Selected |
|--------|-------------|----------|
| Black dashed | Convention in spectroscopy papers | ✓ |
| Green solid | Distinct from blue/orange | |
| Claude decides | Let implementation choose | |

**User's choice:** Black dashed
**Notes:** None

---

## Peak Matching Strategy

**User's choice:** Skip all quantitative peak matching for this phase. Visual overlay only.
**Notes:** User wants to defer NIST-03 (peak position comparison with MAE/RMSE) to a future phase. This phase focuses on NIST-01 (fetch+cache) and NIST-02 (visual overlay).

---

## Integration into Pipeline

### Trigger mechanism

| Option | Description | Selected |
|--------|-------------|----------|
| Automatic when available | Auto-fetch, silently skip on failure | ✓ |
| Opt-in via flag | Requires --experimental / --nist flag | |
| Separate fetch command | Pre-fetch via CLI command | |

**User's choice:** Automatic when available
**Notes:** None

### Script support

| Option | Description | Selected |
|--------|-------------|----------|
| Both anharmonic and harmonic | run_analysis.py and run_analysis_harmonic.py | ✓ |
| Anharmonic only | Only run_analysis.py | |
| Claude decides | Let implementation determine | |

**User's choice:** Both anharmonic and harmonic
**Notes:** None

### Batch report

| Option | Description | Selected |
|--------|-------------|----------|
| Yes, per-molecule overlay | Each molecule in batch report gets overlay if available | ✓ |
| No, single-molecule only | Keep batch report simple | |
| Claude decides | Let implementation determine | |

**User's choice:** Yes, per-molecule overlay
**Notes:** None

---

## Claude's Discretion

- JCAMP-DX parsing details (handling malformed files, unit detection)
- Legend placement and labeling for the experimental trace
- How to note "no experimental data available" in the report
- Frequency range alignment between experimental and computed grids
- Cache file format (raw JCAMP-DX, parsed numpy arrays, or both)

## Deferred Ideas

- NIST-03 quantitative peak comparison (MAE/RMSE) — future phase
- SDBS as secondary data source — not in scope
- Condensed-phase IR data fallback — gas-phase only for now
