# Surface scans along X–H stretches (Section 4.2)

Energy and exact curvature along three X–H stretches, for B3LYP and the four MACE energy models,
each on its own minimum and its own normal mode from the 2026 campaign runs. Report:
thesis repo `docs/handoff/surface-scans-report.md`.

| script | what | where it runs |
|---|---|---|
| `common.py` | targets, mode loading from the run fchk, grids, units | – |
| `scan_mace.py` | energy + autograd Hessian on the fine (0.002 Å), wide (0.01 Å) and step grids | GPU, campaign env `mace4ir_v3_accel`, from a scratch dir |
| `scan_b3lyp.py` | the same with Gaussian 16 (analytic Hessians from `freq` jobs; single points `scf=tight`) | CPU, `g16` |
| `steps_q.py` | curvature at Gaussian's VPT2 points Q = ±0.005 … ±0.040 amu^½ Å | both |
| `zoom.py` | curvature on a 0.0002 Å grid within ±0.012 Å | GPU |
| `scan_2d.py` | methanol O–H stretch × C–O–H bend, 21 × 21, energies | both |
| `vpt2_constants.py` | k(i,i,i), k(i,i,i,i) of one mode from a Gaussian VPT2 log | – |
| `analyse.py <out dir>` | constants, wobble, energy shape, tables (`data/`), figures `surface_*` | any env with matplotlib |

Order: `scan_mace.py`, `scan_b3lyp.py`, `steps_q.py mace|b3lyp`, `zoom.py <target> <model> …`,
`scan_2d.py mace|b3lyp`, `analyse.py`. Gaussian work files go to `work/` (not committed);
all results are in `data/` (JSON per scan, `results.csv`, `tables.md`, `chi_diag_estimate.csv`).

Gaussian's VPT2 displacement is 0.010 in the mass-weighted normal coordinate Q (amu^½ Å), not
0.010 Å of Cartesian displacement: only that reproduces Gaussian's printed quartic constants on
the rough MACE-OMOL-0 surface (`steps_q.py`, all 15 cases agree to the printed digits).
