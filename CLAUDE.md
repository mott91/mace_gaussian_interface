# CLAUDE.md

MACE-Gaussian: bridges ML potentials (MACE) with Gaussian 16 for molecular IR spectroscopy via ZMQ.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for module architecture.
See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) for adding calculators and metrics.

## Commands

```bash
mace-gaussian diagnose                    # Check environment
mace-gaussian run water.xyz               # Full workflow
python run_analysis.py water              # Anharmonic analysis
python run_analysis_harmonic.py water     # Harmonic-only analysis
ruff check --fix && ruff format           # Lint + format (line length 100)
ty check                                  # Type check
```

## Gotchas

**Gaussian requires absolute paths**: Never use relative paths in `.gjf` files. Helper script path set via `MACE_HELPER_SCRIPT_PATH` env var.

**ZMQ socket cleanup**: Delete stale `.ipc_file` manually before a new run after a crash.

**Heatmaps**: Always use `force_harmonic=True` for mode overlap heatmaps (harmonic eigenvectors are the true normal modes).
