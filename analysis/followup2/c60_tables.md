# C60: harmonic T1u intensities with all three dipole models (c60_dipoles.py)

B3LYP T1u levels (cm-1): 538.3, 588.5, 1214.4, 1460.4

## Sanity check (POLAR-1 at the OMOL geometry)

- max |APT(now, CPU) - APT(fchk of the Gaussian run)| = 2.2e-08 e
- T1u from the stored APT: 75.0 / 70.3 / 0.1 / 33.7 km/mol; recomputed now: 75.0 / 70.3 / 0.1 / 33.7 (RESULTS_NOTES §12: 75.0 / 70.3 / 0.1 / 33.7)

## T1u intensities (km/mol, sum of the three components) on each set of harmonic modes

| geometry and modes | T1u freq (cm-1) | mode overlap with B3LYP T1u | MDP | POLAR-1 | MACE4IR |
|---|---|---|---|---|---|
| B3LYP | 538 / 588 / 1214 / 1460 | 1.000 / 1.000 / 1.000 / 1.000 | 101.5 / 1056.2 / 59.3 / 2356.9 | 68.0 / 94.9 / 0.0 / 23.4 | 149.8 / 367111.0 / 2797.3 / 570468.9 |
| ANI-cc | 533 / 585 / 1146 / 1334 | 0.956 / 0.960 / 1.000 / 0.999 | 38.1 / 724.9 / 87.4 / 3550.3 | 54.6 / 148.1 / 0.1 / 7.3 | 18132.9 / 333602.9 / 1245.3 / 509598.9 |
| OMOL | 551 / 612 / 1218 / 1444 | 0.995 / 0.998 / 0.993 / 0.992 | 41.1 / 1205.4 / 62.6 / 453.1 | 75.0 / 70.3 / 0.1 / 33.7 | 1255.8 / 366351.1 / 18795.7 / 472488.8 |
| OFF | 576 / 650 / 1239 / 1537 | 0.994 / 0.997 / 0.999 / 0.998 | 12.5 / 1039.9 / 128.7 / 78.4 | 96.4 / 63.9 / 0.0 / 45.3 | 1834.1 / 367351.9 / 13054.2 / 367668.1 |
| B3LYP's own APT (reference) | | | 85.1 / 33.1 / 27.1 / 35.4 | | |

Element-wise APT at the B3LYP geometry, relative Frobenius error against B3LYP (largest element in e; B3LYP 0.16): POLAR-1 0.94 (0.16), MDP 5.24 (0.75), MACE4IR 72.05 (11.48)

CPU time per APT evaluation (s): POLAR-1 14.0, MDP 13.4, MACE4IR 0.4; total 123 s
