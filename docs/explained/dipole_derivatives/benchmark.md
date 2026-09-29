# Dipole derivatives: autograd vs. finite differences

Run 2026-09-29 with `scripts/benchmark_dipole_derivatives.py`.
Model: MACE4IR (`model_1_dipole.model`), GPU: NVIDIA GeForce RTX 2070 SUPER, CPU: x86_64.
Times are the median of 6 calls (finite differences: 3), one call = what the harness does per Gaussian external call.

| molecule | atoms | FD dipole evals | FD [s] | autograd [s] | speedup | rel. diff δ=0.01 Å | rel. diff δ=0.001 Å |
|---|---|---|---|---|---|---|---|
| water | 3 | 18 | 0.279 | 0.0920 | 3x | 3.1e-04 | 3.1e-06 |
| methanol | 6 | 36 | 0.559 | 0.0909 | 6x | 1.3e-03 | 1.4e-05 |
| gly | 10 | 60 | 0.938 | 0.0916 | 10x | 5.0e-04 | 5.1e-06 |
| aspirin | 21 | 126 | 2.458 | 0.0958 | 26x | 2.0e-03 | 2.0e-05 |
| decane | 32 | 192 | 4.646 | 0.1212 | 38x | 1.3e-03 | 1.3e-05 |

Relative difference = ||FD - autograd|| / ||autograd|| (Frobenius norm over the 3N x 3 matrix).
Autograd cost is roughly flat in N; finite differences cost 6N dipole evaluations,
so the speedup grows with molecule size. The difference drops 100x for a 10x smaller
step, the O(δ²) signature of central differences converging onto the autograd value.
