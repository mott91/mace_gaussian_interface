# External Integrations

**Analysis Date:** 2026-02-16

## APIs & External Services

**Quantum Chemistry:**
- Gaussian 16 - Commercial quantum chemistry package
  - Integration: External process via subprocess, writes `.gjf` input files, parses `.log` output
  - Real-time dipole injection: ZMQ IPC socket communication via helper script
  - Used for: Harmonic/anharmonic frequency calculations, IR intensity data

## Data Storage

**Databases:**
- None - System uses filesystem-based storage only
  - Structure: `comparison_results/` and `analysis_results/` directories
  - No SQLite, PostgreSQL, or other database backends
  - Results serialized as JSON (see `results_manager.py`)

**File Storage:**
- Local filesystem only
  - Molecular structures: `.xyz` files
  - Gaussian I/O: `.gjf` (input), `.log` (output), `.chk` (checkpoint), `.fchk` (formatted checkpoint)
  - Calculation results: JSON metadata in `results.json`
  - Analysis output: PNG plots at 300 DPI, HTML reports
  - Dipole model: Binary `.model` file (trained PyTorch checkpoint)

**Caching:**
- None - No Redis or caching layer
- Intermediate files stored in local directories for reuse across analyses

## Authentication & Identity

**Auth Provider:**
- None - No authentication required
- All operations local to filesystem/subprocess
- No API keys, tokens, or credentials needed

## Monitoring & Observability

**Error Tracking:**
- None - No external error tracking (Sentry, Rollbar, etc.)
- Local logging via Python `logging` module to stderr/stdout

**Logs:**
- Python `logging.basicConfig()` in main modules
  - Default level: INFO
  - Format: "%(asctime)s - %(levelname)s - %(message)s"
  - Output: Console only (redirected to stdout/stderr in scripts)
- No persistent log files created
- Suppressed warnings:
  - `FutureWarning` from PyTorch (weights_only=False)
  - RDKit app logging disabled via `RDLogger.DisableLog('rdApp.*')`

## CI/CD & Deployment

**Hosting:**
- None - Designed for local research use
- Runnable on HPC clusters (supports SLURM module system via ConfigArgParse)
- Single-machine GPU computation assumed

**CI Pipeline:**
- None - No GitHub Actions, GitLab CI, or similar
- Manual testing workflow documented in CLAUDE.md
  - `ruff check --fix && ruff format` for code quality
  - `ty check` for type checking
  - `pytest` for unit tests (framework installed but no CI runner configured)

## Environment Configuration

**Required env vars:**
- `MACE_DIPOLE_MODEL_PATH` - Path to trained dipole model file
  - Default: `~/mace_gaussian/dipole_model/model_1.model`
  - Format: Binary PyTorch checkpoint (`.model`)
- `MACE_HELPER_SCRIPT_PATH` - Path to ZMQ helper script for Gaussian integration
  - Default: `./gm_helper.py` (relative to working directory)
- `CUDA_VISIBLE_DEVICES` - GPU device selection (e.g., "0,1")
- `PYTHONWARNINGS` - Set to "ignore::FutureWarning" in code to suppress PyTorch warnings

**Secrets location:**
- No secrets stored (.env files not used)
- Model files assumed to be local, not versioned in git
- Dipole model path configurable via environment only

## Webhooks & Callbacks

**Incoming:**
- None - Gaussian frequency calculations do not send webhooks
- File-based handoff: Gaussian writes `.log` output, Python reads synchronously

**Outgoing:**
- None - No external API calls or notifications
- IPC communication: Python helper script (`gm_helper.py`) communicates with Gaussian via ZMQ REQ-REP socket
  - Socket file: `.ipc_file` (Unix domain socket)
  - Protocol: "infile|outfile" message format for dipole derivative injection

## External Process Communication

**Gaussian Integration:**
- Subprocess execution: `subprocess.run(['g16', gjf_file])`
- Input format: `.gjf` (Gaussian Job File) with embedded external command hook
- Output parsing: `.log` files via regex patterns (frequencies, IR intensities)
- Checkpoint conversion: `formchk` utility called via subprocess to convert `.chk` → `.fchk`

**ZMQ Real-time Communication:**
- IPC socket at `.ipc_file` (created by server in `gm_main.py`)
- Client: Helper script `gm_helper.py` (executed by Gaussian as external hook)
- Message format: `"infile|outfile"` (pipe-separated filenames)
- Response: Acknowledgment message (string)
- Pattern: REQ-REP (synchronous request-reply)
- Cleanup: `.ipc_file` must be manually deleted if process crashes

## ML Model Integration

**MACE Models:**
- Energy calculators:
  - `mace_mp` - Pre-trained MACE checkpoint
  - `mace_omol` - MACE-OMOL-0 foundation model
  - `mace_off` - OpenFF MACE variant
- Dipole calculators:
  - `mace_ml` - Custom dipole model loaded from `MACE_DIPOLE_MODEL_PATH`
  - Located: `dipole_model/model_1.model` (must exist before first run)

**Charge Calculators:**
- `espaloma` - espaloma_charge package, computes atomic charges
- Integration: ASE atoms object passed to espaloma charge function

**Model Loading:**
- Direct PyTorch model loading: `torch.load()` on binary checkpoints
- MACE calculator instantiation via `MACECalculator()` with device specification ("cuda" or "cpu")
- Module monkey-patching system in `mace_calculators.py`:
  - `load_standard_mace_calculator()` - Standard MACE (energy only)
  - `load_dipole_mace_calculator()` - Dipole-enabled MACE (swaps `mace.modules.models`)
  - `cleanup_mace_modules()` - Must be called in finally blocks to restore original modules

## Data Format Integrations

**Input:**
- XYZ molecular geometry files (`.xyz`) - ASE `read()` function
- SMILES strings - RDKit `Chem.MolFromSmiles()`

**Output:**
- JSON: `results.json` with calculation metadata and metrics
- PNG: Spectrum plots, regression plots, mode overlap heatmaps (seaborn + matplotlib)
- HTML: Dynamic reports generated by `html_report_generator.py`
- Numpy arrays: `.npy` files (optional for large datasets)
- HDF5: `.h5` format for checkpoint storage (via h5py, not currently used)

---

*Integration audit: 2026-02-16*
