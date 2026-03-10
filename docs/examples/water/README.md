# Water (H2O) -- Worked Example

This directory contains a complete worked example for water using the MACE-Gaussian pipeline.
The expected output files allow inspection of results without re-running the calculation.

## Input

- `water.xyz` -- 3-atom water geometry (O + 2H)

## Expected Output

Pre-committed reference results (generated with MACE-OMOL-0 energy, MACE4IR dipole):

### Geometry Optimization (`expected_output/geometry_opt_results.json`)

- Calculator: mace_omol
- Converged: yes (8 steps, ~4.2 seconds on RTX 2070)
- Final energy: -2079.866 eV

### Frequency Calculation (`expected_output/mace_omol_mace_ml_results.json`)

Water has 3 vibrational modes:

| Mode | Harmonic (cm-1) | Anharmonic (cm-1) | Assignment |
|------|-----------------|-------------------|------------|
| 1 | 1621.85 | 1570.25 | HOH bending |
| 2 | 3818.06 | 3643.79 | O-H symmetric stretch |
| 3 | 3918.36 | 3739.40 | O-H asymmetric stretch |

### Analysis Metrics (`expected_output/analysis_metrics_summary.json`)

Comparison against B3LYP/6-31G(d,p) DFT reference:

| ML combination | R² (freq) | RMSE (cm-1) |
|----------------|-----------|-------------|
| mace_omol + mace_ml | 0.9997 | 44 |
| mace_omol + espaloma | 0.9997 | 44 |
| mace_mp + espaloma | 0.9971 | 149 |
| mace_mp + mace_ml | 0.9971 | 149 |

### Plots (`expected_plots/`)

- `spectrum_combined.png` -- IR spectrum comparison (all 4 ML combinations vs DFT)
- `regression_combined.png` -- Frequency regression plot (ML vs DFT frequencies)

---

## To Reproduce

**Requirements:** Gaussian 16 (`g16`) in PATH, CUDA GPU, ~45 minutes.

```bash
# From the project root (two levels up from this directory):
cd ../../..

# Run the full pipeline (includes DFT baseline for comparison)
python cli.py run water.xyz

# Or skip DFT baseline for ML-only results (~15 minutes):
python cli.py run water.xyz --skip-dft-baseline

# After calculation, run analysis:
python run_analysis.py water

# View the report:
xdg-open analysis_results/water/report.html
```

Results are written to `comparison_results/water/` and `analysis_results/water/`.

---

## File Structure

```
docs/examples/water/
├── README.md               (this file)
├── water.xyz               (input: 3-atom water geometry)
├── expected_output/
│   ├── geometry_opt_results.json       (geometry optimization)
│   ├── mace_omol_mace_ml_results.json  (frequency calculation, best combination)
│   └── analysis_metrics_summary.json   (R², RMSE for all 4 combinations)
└── expected_plots/
    ├── spectrum_combined.png           (IR spectrum comparison)
    └── regression_combined.png         (frequency regression)
```

Note: `.chk`, `.fchk`, and `.log` files from Gaussian are not committed
(large, system-dependent). Only portable JSON and PNG files are included.
