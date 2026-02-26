# Quickstart Guide

This guide walks you through running a complete MACE-Gaussian calculation on water (H2O)
and viewing the analysis report.

**Prerequisites:** Gaussian 16 (`g16`) in PATH, CUDA-capable GPU, Python 3.10.

---

## Step 1: Install

Follow the [Installation section in README.md](../README.md#installation) to set up
the conda environment, install the custom MACE packages, and verify with `python cli.py diagnose`.

---

## Step 2: Run a calculation

The project root contains `water.xyz` as a sample input. Run the ML-only pipeline
(skipping the DFT baseline to save ~30 minutes):

```bash
python cli.py run water.xyz --skip-dft-baseline
```

This runs geometry optimization with MACE-OMOL-0, then anharmonic frequency calculations
for all four ML model combinations (mace_mp/mace_omol x espaloma/mace_ml).
Expected runtime: ~10-20 minutes on an RTX 2070 or equivalent.

Results are saved to `comparison_results/water/`.

To also run the DFT baseline for comparison (B3LYP/6-31G(d,p), ~30 additional minutes):

```bash
python cli.py run water.xyz
```

---

## Step 3: Check calculation results

```bash
python cli.py list water
```

This shows each ML combination with its harmonic and anharmonic frequency counts,
runtime, and energy.

---

## Step 4: Run analysis

```bash
python run_analysis.py water
```

This reads the frequency results, matches modes against DFT (if available), computes
R² and RMSE metrics, generates regression and spectrum comparison plots, and writes
an HTML report to `analysis_results/water/`.

For harmonic-only analysis with mode overlap heatmaps:

```bash
python run_analysis_harmonic.py water
```

Results go to `analysis_results_harmonic/water/`.

---

## Step 5: View the report

Open the report in a browser:

```bash
# Linux
xdg-open analysis_results/water/report.html

# Or navigate to it manually
ls analysis_results/water/
```

The report contains spectrum comparison plots, regression analysis, and a metrics table
(MAE, RMSE, R² for each ML combination vs DFT reference).

---

## Expected results

For reference, the committed expected output in `docs/examples/water/` shows:

| ML combination | R² (freq) | RMSE (cm-1) |
|----------------|-----------|-------------|
| mace_omol + mace_ml | 0.9997 | 44 |
| mace_omol + espaloma | 0.9997 | 44 |
| mace_mp + espaloma | 0.9971 | 149 |
| mace_mp + mace_ml | 0.9971 | 149 |

See [`docs/examples/water/`](examples/water/) for committed JSON results and plots.

---

## Troubleshooting

**`g16` not found**: Gaussian 16 must be in PATH. Check with `which g16`.

**ZMQ socket error**: If a previous run crashed, delete the `.ipc_file` in the project
root before restarting: `rm -f .ipc_file`.

**CUDA not available**: The pipeline falls back to CPU (~10x slower). Check with
`python cli.py diagnose`.

**Import errors for mace_ml or espaloma**: Ensure both custom packages are installed:
`pip install -e mace_ML_pkg && pip install -e mace_dipole_pkg`.
