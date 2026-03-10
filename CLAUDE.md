# CLAUDE.md

MACE-Gaussian Interface: bridges ML potentials (MACE) with Gaussian 16 for molecular IR spectroscopy. Uses ZMQ to inject ML-calculated dipole derivatives into Gaussian's anharmonic frequency calculations in real-time.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for detailed module architecture and data structures.
See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) for adding calculators, metrics, and testing workflows.
See [permissions.md](permissions.md) for what requires confirmation vs. auto-approve.

## Commands

```bash
# Setup
uv sync && source .venv/bin/activate
mace-gaussian diagnose                    # Check environment

# Run calculations
mace-gaussian run water.xyz               # Full workflow (recommended)
mace-gaussian run water.xyz --skip-dft-baseline
mace-gaussian run water.xyz --energy-calculators mace_mp --dipole-calculators espaloma

# Analysis
python run_analysis.py water              # Anharmonic (overtones + combinations)
python run_analysis_harmonic.py water     # Harmonic-only (fundamentals, mode matching)

# Results
mace-gaussian list                        # List all results
python convert_all_chk_files.py           # Batch .chk → .fchk for mode matching

# Code quality
ruff check --fix && ruff format           # Lint + format (line length 100)
ty check                                  # Type check
```

## Core Pipeline

```
Input (.xyz) → Geometry Optimization (MACE) → Frequency Calc (Gaussian + ML dipoles via ZMQ)
  → Anharmonic Calc (Gaussian) → Parse Results (JSON) → Mode Matching (eigenvector dot product)
  → Statistical Analysis (regression, KDE) → HTML Report
```

## Critical Gotchas

**MACE dipole model loading**: The dipole model was saved with `mace.modules.models` class paths but requires `mace_dipole_core` classes at runtime. `calculators/mace_loader.py` handles this via `torch.load(pickle_module=...)` class remapping. No `sys.modules` cleanup is needed.

**Gaussian requires absolute paths**: Never use relative paths in `.gjf` files for external scripts. The helper script path is set via `MACE_HELPER_SCRIPT_PATH` env var.

**ZMQ socket cleanup**: If `.ipc_file` exists from a previous crash, delete it manually before a new run.

**CUDA device placement**: Some calculators default to CPU. Always pass `device="cuda"` explicitly when GPU is available.

**Frequency parsing**: Overtones match `"2(X)"` pattern, combination bands match `"X(1) + Y(1)"`. Must handle harmonic/anharmonic sections separately. Some DFT methods may lack anharmonic data.

**Heatmaps**: Always use `force_harmonic=True` for mode overlap heatmaps (harmonic eigenvectors are the true normal modes).

**Acetic Acid bug** (commit a4384c4): Frequency matching for regression plots can fail with acoh DFT calculations. Check overtone/combination band parsing logic.

## Conventions

**File naming**:
- Molecules: `{name}.xyz`, Gaussian I/O: `{name}_freq_anharm.{gjf,log,chk,fchk}`
- DFT dirs: `b3lyp_6-31Gdp`, ML dirs: `mace_mp_espaloma`, `mace_omol_mace_ml`
- Analysis output: `analysis_results/` (anharmonic), `analysis_results_harmonic/` (harmonic)

**Plots**: PNG format, DPI=300, seaborn "colorblind" palette, Arial/Helvetica font. Always include R² and RMSE in regression plots. Standard range 400-4000 cm⁻¹.

## External Dependencies

- **Gaussian 16**: Must be in PATH as `g16`
- **CUDA**: Required for GPU (CPU fallback works but slow)
- **Dipole model**: `~/mace_gaussian/mace4ir_models/pretrained_models/model_1_dipole.model` (MACE4IR, 78 elements). Old narrow HCNO model at `dipole_model/model_1.model`. Override with `MACE_DIPOLE_MODEL_PATH`.
- **MACE-OMOL-0**: Requires mace-torch >= 0.3.14
