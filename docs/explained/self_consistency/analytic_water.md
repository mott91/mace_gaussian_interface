# DFT self-consistency check: water (analytic)

Harness: `gaussian_b3lyp` through `run_frequency_calculation` ({'force': 14, 'freq': 7} inner Gaussian jobs). Native: `b3lyp/6-31G(d,p) freq(anharm)` at the harness geometry.

**Verdict: PASS** (max |Δν| = 0.0080 cm⁻¹ against 0.1 cm⁻¹; max relative |ΔI| = 0.02% against 1%)

Final energy: harness -76.4196339000 Ha, native -76.4196339461 Ha, Δ = 4.61e-08 Ha
Dipole: harness 2.0429325991564684 D, native 2.0429325991564684 D

## Harmonic fundamentals

| band | freq harness | freq native | Δ (cm⁻¹) | I harness | I native | ΔI (km/mol) | ok |
|---|---|---|---|---|---|---|---|
| 1 | 1665.2926 | 1665.2958 | -0.0032 | 70.3196 | 70.3207 | -0.0011 | ✓ |
| 2 | 3799.6942 | 3799.6955 | -0.0013 | 1.6382 | 1.6383 | -0.0001 | ✓ |
| 3 | 3912.8261 | 3912.8257 | +0.0004 | 20.2180 | 20.2180 | +0.0000 | ✓ |

## VPT2 fundamentals

| band | freq harness | freq native | Δ (cm⁻¹) | I harness | I native | ΔI (km/mol) | ok |
|---|---|---|---|---|---|---|---|
| 1 | 3624.4690 | 3624.4710 | -0.0020 | 0.6290 | 0.6291 | -0.0001 | ✓ |
| 2 | 1615.1440 | 1615.1480 | -0.0040 | 69.2839 | 69.2849 | -0.0010 | ✓ |
| 3 | 3722.6850 | 3722.6850 | +0.0000 | 16.0898 | 16.0895 | +0.0003 | ✓ |

## Overtones

| band | freq harness | freq native | Δ (cm⁻¹) | I harness | I native | ΔI (km/mol) | ok |
|---|---|---|---|---|---|---|---|
| 1/2 | 7162.5100 | 7162.5130 | -0.0030 | 0.7756 | 0.7756 | +0.0000 | ✓ |
| 2/2 | 3194.7370 | 3194.7450 | -0.0080 | 0.8058 | 0.8058 | -0.0000 | ✓ |
| 3/2 | 7347.5010 | 7347.5010 | +0.0000 | 0.0492 | 0.0492 | -0.0000 | ✓ |

## Combination bands

| band | freq harness | freq native | Δ (cm⁻¹) | I harness | I native | ΔI (km/mol) | ok |
|---|---|---|---|---|---|---|---|
| 2/1 | 5228.4910 | 5228.4960 | -0.0050 | 0.4065 | 0.4065 | -0.0000 | ✓ |
| 3/1 | 7180.6840 | 7180.6850 | -0.0010 | 3.6819 | 3.6818 | +0.0001 | ✓ |
| 3/2 | 5319.7560 | 5319.7590 | -0.0030 | 3.3074 | 3.3073 | +0.0001 | ✓ |
