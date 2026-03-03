# Backlog

Living list of ideas, deferred work, bugs, and research directions. Items here are **not committed** — they inform milestone scoping. Move items into REQUIREMENTS.md when they enter active scope.

*Last updated: 2026-02-28*

---

## Bugs

- **Acoh frequency parser bug** — regression plot frequency matching fails for acetic acid DFT calculations. Xfailed in test suite (commit a4384c4). Needs root cause diagnosis and fix.

---

## Tech Debt / Cleanup

- **Update ARCHITECTURE.md and DEVELOPMENT.md** — both docs still reference pre-refactor module names and structure. Should reflect the current `mace_gaussian/` package layout.
- **E2E GPU+Gaussian pipeline test** — Phase 5 confirmation requires hardware access; currently untested in CI. Needs a hardware-gated test or documented manual procedure.
- **Implement `compare` and `export` CLI commands** — currently honest stubs. Deferred from v1.0.

---

## New Calculators / Backends

- **ORCA support for VPT2** — implement ORCA as an alternative to Gaussian 16 for anharmonic (VPT2) frequency calculations. Goal: fully open-source package that doesn't require a Gaussian license. Would be a major enabler for sharing the method with the community.
- **Alternative energy calculators** — survey what else exists beyond mace_mp and mace_omol for geometry optimization and energy evaluation.
- **Alternative intensity/dipole calculators** — survey alternatives to espaloma and mace_ml for dipole moment prediction. What else is viable?
- **Autograd-based dipole derivatives** — replace numerical finite differences (18 ZMQ calls/water, 126/aspirin) with autograd. ~10× speedup. Revisit when scaling to larger molecules.

---

## HPC / Infrastructure

- **Decouple DFT baseline calculations** — currently tightly coupled to the local pipeline. Need to support running DFT baselines on an HPC cluster (submit job, wait, collect results) and have the ML pipeline consume the results independently.
- **Full HPC workflow** — explore putting the entire pipeline (geometry opt + freq calc + analysis) on a cluster with job scheduling. May require rethinking the ZMQ/IPC architecture (localhost assumption breaks across nodes).

---

## Batch Processing & Validation Campaign

- **Batch molecule workflow** — create tooling to run a systematic set of molecules through the full pipeline:
  - Size-scaling series: same atom types, increasing molecule size — shows where ML becomes viable vs. DFT (crossover point)
  - Diversity series: molecules with very different electronic structure and functional groups — shows which calculator combinations work well for what chemistry
  - Benchmark against experimental IR data where available
- **PubChem structure fetcher** — tool to query PubChem by name/CID, download 3D structures, convert to XYZ, and feed into the batch workflow. Would make it easy to build a systematic benchmark set without manually hunting for structures.
- **Systematic benchmark design** — research session to design the benchmark: which molecule classes, what size range, what experimental data sources, how to quantify "ML viable vs DFT necessary."

---

## Periodic Systems

- **Feasibility study: periodic systems** — explore whether the workflow can be adapted for periodic systems (surfaces, crystals, MOFs). Questions to answer:
  - Does VPT2 even make sense for periodic systems? (phonons vs molecular vibrations)
  - Which MACE variants support periodic boundary conditions?
  - What replaces Gaussian for the QM reference? (VASP, CP2K, Quantum ESPRESSO)
  - What is the complexity delta vs. molecular systems?
  - Is this a thesis scope item or a future paper?

---

## Research / Exploration

- **Calculator combination brainstorm** — systematic review of what other ML models could plug into the energy/dipole/frequency slots. Literature review of recent models.
- **Other cool directions** — open-ended: what else could be interesting given the ZMQ injection architecture and the ML potential ecosystem? (e.g., uncertainty quantification on spectra, active learning for dipole models, multi-fidelity approaches)

---

## Deferred (explicitly out of scope for now, revisit later)

- PyPI distribution — clone-and-reproduce is the target for the thesis
- Web interface or GUI — CLI appropriate for research use
- Support for non-Gaussian QM packages other than ORCA (already in backlog above)
- Performance optimization for very large molecules (>100 atoms) — thesis molecules are smaller
