# Phase 24 — VPT2 Research Spike: Psience as an alternative VPT2 engine

**Date**: 2026-07-06 · **Branch**: `spike/24-vpt2-psience` · **Verdict: PROOF-OF-CONCEPT ACHIEVED**

## Goal

Determine whether the Psience library (McCoy group) can replace Gaussian's VPT2
engine using MACE force constants (roadmap success criteria: run on water from
existing Gaussian output, compare against Gaussian VPT2, or document blockers).

## Result

Psience VPT2 reproduces Gaussian's VPT2 fundamentals **to < 0.01 cm⁻¹** on water
when fed the normal-mode reduced force field parsed from the Gaussian log
(`scripts/spike_psience_vpt2.py`, runs in an isolated venv — NOT mace4ir_v2):

| Surface | Mode | Harmonic | Gaussian VPT2 | Psience VPT2 |
|---|---|---|---|---|
| B3LYP/6-31G(d,p) | ν₁ sym | 3799.220 | 3624.001 | 3624.000 |
| | ν₂ bend | 1665.300 | 1615.153 | 1615.153 |
| | ν₃ asym | 3912.417 | 3722.273 | 3722.273 |
| mace_omol + mace_ml | ν₁ sym | 3818.058 | 3643.791 | 3643.791 |
| | ν₂ bend | 1621.852 | 1570.251 | 1570.251 |
| | ν₃ asym | 3918.362 | 3739.402 | 3739.401 |

The exact match (including implicit Coriolis/Watson contributions, which Psience
derives itself from the fchk geometry + modes) shows the two engines implement
the same VPT2 for this system.

## Key findings

1. **Psience's advertised one-liner (`VPTRunner.run_simple(fchk, n)`) does NOT
   work on our G16 Rev C.02 fchk files.** It yields unphysical fundamentals
   (OH stretches shift *up*: 3912 → 4084 cm⁻¹; stretch overtones land above
   2× the fundamental). Verified on both DFT and ML fchks. The McCoy group's
   own G09-era MP2 water test fchk works correctly, so this is a
   G16-format incompatibility in reading `Cartesian 3rd/4th derivatives`
   (suspected mode-block ordering — Gaussian's VPT2 numbers modes
   spectroscopically, not ascending — plus a block normalization difference;
   a pure block permutation experiment moved results toward sane but did not
   fix them).
2. **Workaround (used here): bypass the fchk derivative block entirely.**
   Parse the `CUBIC/QUARTIC FORCE CONSTANTS IN NORMAL MODES` tables (reduced,
   cm⁻¹) from the log, remap Gaussian's spectroscopic mode numbering to
   ascending order, and pass them via `potential_terms` (with
   `zero_element_warning=False`, since Gaussian truncates sub-threshold
   constants). The fchk still supplies geometry/masses/modes.
3. The diagnosis route matters for trust: Psience's internally-derived V3/V4
   from our fchk *matched Gaussian's printed tables* in magnitude — the
   corruption only shows downstream — so naive spot-checks of force constants
   would not have caught the direct-fchk bug.

## Limitations / open questions

- **Intensities not yet wired**: Psience computed harmonic intensities from the
  fchk dipole derivatives; anharmonic intensities need the dipole expansion
  (`dipole_terms`) — not attempted in this timebox.
- **Resonance handling**: water has no strong Fermi resonance; larger molecules
  will exercise Gaussian's deperturbation vs Psience's degenerate-PT defaults
  (`gaussian_resonance_handling=True` exists but was not needed here). Expect
  divergence on resonant systems until tuned.
- Log tables carry 5 decimals and threshold-truncate small constants —
  agreement is limited by print precision for bigger systems.

## Implication for the thesis (the speed story)

Today the anharmonic ML pipeline exists only because Gaussian's VPT2 driver
orchestrates the finite differencing (many ZMQ round-trips). With Psience as
the VPT2 solver, the Gaussian-free path is: MACE Hessians finite-differenced
along normal modes (2·(3N−6)+1 cheap Hessian calls) → cubic/quartic field →
Psience. For water that is seconds of compute. Recommended next step if
pursued: generate the force field directly from MACE (never printing/parsing
logs) and validate against this spike's numbers, then tackle `dipole_terms`
for anharmonic intensities.
