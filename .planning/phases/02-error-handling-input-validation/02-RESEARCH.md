# Phase 2: Error Handling & Input Validation - Research

**Researched:** 2026-02-16
**Domain:** Python error handling, input validation, subprocess management, reproducibility metadata
**Confidence:** HIGH

## Summary

This phase transforms the MACE-Gaussian Interface from "works on the happy path" to "fails explicitly and early with context." The codebase currently has several critical gaps: parsers that return empty data silently when sections are missing (gaussian_parser.py returns `[]` for missing sections), no startup validation of external dependencies (g16, formchk, dipole model), hardcoded `device="cuda"` with no fallback detection, no version metadata in results JSON, and hardcoded `num_steps=0` in optimization results.

The nine requirements (ERR-01 through ERR-06, REPR-01 through REPR-03) map cleanly to four work areas: (1) parser error hardening in `gaussian_parser.py`, (2) startup prerequisite validation in a new `validation.py` module integrated into `cli.py` and `gm_main.py`, (3) CUDA device detection with explicit fallback logging, and (4) reproducibility metadata enrichment in `results_manager.py`. These changes are almost entirely additive -- they augment existing code with checks and metadata rather than restructuring control flow.

A key architectural decision is where to place the validation logic. The current codebase has `diagnostics.py` for environment info but it only prints output -- it does not raise errors or return structured results. The recommended approach is a new `validation.py` module with pure functions that return results or raise exceptions, keeping the diagnostic/printing concern separate from the validation/gating concern.

**Primary recommendation:** Create a `validation.py` module with prerequisite check functions, add a custom exception hierarchy rooted in a `MaceGaussianError` base class, and integrate validation into both `cli.py` (pre-command hook) and `gm_main.run_workflow()` (defensive checks). Enrich `results_manager.py` metadata dicts with version info. Fix `geometry_optimisation()` to return step count from `opt.get_number_of_steps()`.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python stdlib | 3.9+ | `shutil.which()`, `subprocess`, `pathlib`, `platform`, `sys` | All validation can be done with stdlib |
| click | (installed) | CLI parameter validation, callbacks | Already used for CLI, has built-in validation |
| torch | (installed) | `torch.cuda.is_available()` for CUDA detection | Already a dependency |
| ase | >=3.22 | `opt.get_number_of_steps()` for step tracking | Already a dependency |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| logging (stdlib) | - | Structured warning/error messages | All validation output |
| json (stdlib) | - | Version metadata serialization | Results JSON enrichment |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Custom `validation.py` | Extend `diagnostics.py` | diagnostics.py prints to stdout; validation needs to raise exceptions and return structured data -- different concerns |
| Custom exception classes | Built-in exceptions only | Custom exceptions allow callers to catch specific MACE-Gaussian errors vs. general Python errors |
| `shutil.which()` | `subprocess.run(["which", "g16"])` | `shutil.which()` is cross-platform stdlib, no subprocess overhead |

**Installation:**
```bash
# No new dependencies needed -- all stdlib + existing packages
```

## Architecture Patterns

### Recommended Project Structure
```
validation.py           # NEW: prerequisite checks, input validation
exceptions.py           # NEW: custom exception hierarchy
gaussian_parser.py      # MODIFY: raise on missing sections (ERR-01)
gm_main.py              # MODIFY: CUDA detection, step tracking, integrate validation
cli.py                  # MODIFY: pre-command validation hook
results_manager.py      # MODIFY: version metadata in save methods
diagnostics.py          # MINOR: can reference validation.py for reuse
```

### Pattern 1: Custom Exception Hierarchy
**What:** A flat exception hierarchy rooted in `MaceGaussianError` with specific subclasses for each error domain
**When to use:** Any error that is specific to this tool (not a generic Python error)
**Example:**
```python
# exceptions.py
class MaceGaussianError(Exception):
    """Base exception for all MACE-Gaussian interface errors."""
    pass

class PrerequisiteError(MaceGaussianError):
    """A required external tool or file is missing."""
    pass

class GaussianParseError(MaceGaussianError):
    """Failed to parse expected data from Gaussian output."""
    pass

class InputValidationError(MaceGaussianError):
    """Input file or parameter failed validation."""
    pass

class CUDANotAvailableWarning(UserWarning):
    """CUDA requested but not available, falling back to CPU."""
    pass
```

### Pattern 2: Prerequisite Validation Module
**What:** Pure functions that check prerequisites and return results or raise `PrerequisiteError`
**When to use:** At CLI startup (before expensive imports) and at workflow entry
**Example:**
```python
# validation.py
import shutil
from pathlib import Path
from exceptions import PrerequisiteError, InputValidationError

def check_gaussian_available() -> str:
    """Check g16 is in PATH. Returns path or raises PrerequisiteError."""
    g16_path = shutil.which("g16")
    if g16_path is None:
        raise PrerequisiteError(
            "Gaussian 16 (g16) not found in PATH. "
            "Ensure Gaussian is installed and 'module load gaussian' has been run."
        )
    return g16_path

def check_formchk_available() -> str:
    """Check formchk is in PATH. Returns path or raises PrerequisiteError."""
    formchk_path = shutil.which("formchk")
    if formchk_path is None:
        raise PrerequisiteError(
            "Gaussian formchk utility not found in PATH. "
            "Mode matching requires formchk for .chk to .fchk conversion."
        )
    return formchk_path

def check_dipole_model(model_path: str) -> Path:
    """Check dipole model file exists. Returns Path or raises PrerequisiteError."""
    p = Path(model_path)
    if not p.exists():
        raise PrerequisiteError(
            f"MACE dipole model not found at: {model_path}\n"
            f"Set MACE_DIPOLE_MODEL_PATH environment variable or place model at default location."
        )
    return p

def check_helper_script(script_path: str) -> Path:
    """Check helper script exists. Returns Path or raises PrerequisiteError."""
    p = Path(script_path)
    if not p.exists():
        raise PrerequisiteError(
            f"Helper script not found at: {script_path}\n"
            f"Set MACE_HELPER_SCRIPT_PATH environment variable."
        )
    return p

def validate_xyz_file(file_path: str) -> dict:
    """Validate XYZ file for basic sanity. Returns info dict or raises InputValidationError."""
    p = Path(file_path)
    if not p.exists():
        raise InputValidationError(f"Input file not found: {file_path}")
    if not p.suffix == '.xyz':
        raise InputValidationError(f"Expected .xyz file, got: {p.suffix}")

    # Parse and validate content
    try:
        from ase.io import read
        atoms = read(str(p))
    except Exception as e:
        raise InputValidationError(f"Cannot parse XYZ file '{file_path}': {e}") from e

    n_atoms = len(atoms)
    if n_atoms < 1:
        raise InputValidationError(f"XYZ file has 0 atoms: {file_path}")
    if n_atoms > 200:
        # Not an error, but warn -- very large molecules will be slow
        import logging
        logging.getLogger(__name__).warning(
            f"Large molecule detected: {n_atoms} atoms. "
            f"Calculations may take hours."
        )

    return {"path": str(p), "n_atoms": n_atoms, "symbols": atoms.get_chemical_symbols()}

def validate_all_prerequisites(
    check_gaussian: bool = True,
    check_formchk_tool: bool = True,
    dipole_model_path: str = None,
    helper_script_path: str = None,
) -> dict:
    """Run all prerequisite checks. Returns summary dict. Raises on first failure."""
    results = {}
    if check_gaussian:
        results["g16_path"] = check_gaussian_available()
    if check_formchk_tool:
        results["formchk_path"] = check_formchk_available()
    if dipole_model_path:
        results["dipole_model"] = str(check_dipole_model(dipole_model_path))
    if helper_script_path:
        results["helper_script"] = str(check_helper_script(helper_script_path))
    return results
```

### Pattern 3: CUDA Device Detection
**What:** Detect CUDA at startup, log result, pass device string to calculators
**When to use:** Before any calculator instantiation in `gm_main.py`
**Example:**
```python
# In gm_main.py or validation.py
import torch
import logging

logger = logging.getLogger(__name__)

def detect_device() -> str:
    """Detect best available compute device. Returns 'cuda' or 'cpu'."""
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(0)
        logger.info(f"CUDA available: {device_name}")
        return "cuda"
    else:
        logger.warning(
            "CUDA not available. Falling back to CPU. "
            "Calculations will be significantly slower. "
            "Ensure CUDA drivers and PyTorch CUDA are installed."
        )
        return "cpu"
```

### Pattern 4: Parser Error Hardening
**What:** Parsers raise `GaussianParseError` when expected sections are missing, instead of returning empty lists
**When to use:** When a caller expects data that should exist (e.g., harmonic frequencies from an anharmonic log)
**Example:**
```python
# In gaussian_parser.py - modified parse_harmonic_frequencies()
def parse_harmonic_frequencies(self) -> List[Dict[str, float]]:
    frequencies = []
    # ... existing parsing logic ...

    if not frequencies:
        raise GaussianParseError(
            f"No harmonic frequencies found in {self.log_file}. "
            f"Check that the Gaussian calculation completed normally "
            f"and included a frequency calculation."
        )
    return frequencies
```

**Important nuance:** Not ALL empty returns should become errors. The distinction:
- `parse_harmonic_frequencies()` returning empty from a freq calculation log = ERROR (expected section missing)
- `parse_anharmonic_frequencies()` returning empty from a harmonic-only log = EXPECTED (section legitimately absent)
- `parse_overtones()` returning empty = MAY be legitimate (some small molecules have no overtones parsed)

The recommendation is: `parse_harmonic_frequencies()` always raises if empty (there must be harmonic modes in any freq log). For anharmonic/overtone/combination, add an optional `strict=False` parameter that raises when `True` but returns empty when `False` (default). This preserves backward compatibility while enabling strict mode for contexts where data is expected.

### Pattern 5: Version Metadata Collection
**What:** Collect tool version, Python version, MACE version, and calculation parameters into a metadata dict
**When to use:** Every time results are saved to JSON
**Example:**
```python
# In a new function, possibly in validation.py or results_manager.py
import sys
import platform

def collect_version_metadata() -> dict:
    """Collect version info for reproducibility."""
    metadata = {
        "tool_version": "0.2.0",  # from pyproject.toml
        "python_version": sys.version,
        "platform": platform.platform(),
    }

    # MACE version
    try:
        from mace import __version__ as mace_version
        metadata["mace_version"] = mace_version
    except ImportError:
        metadata["mace_version"] = "unknown"

    # PyTorch version
    try:
        import torch
        metadata["torch_version"] = torch.__version__
        metadata["cuda_available"] = torch.cuda.is_available()
        if torch.cuda.is_available():
            metadata["cuda_version"] = torch.version.cuda
            metadata["gpu_name"] = torch.cuda.get_device_name(0)
    except ImportError:
        metadata["torch_version"] = "unknown"

    # ASE version
    try:
        import ase
        metadata["ase_version"] = ase.__version__
    except (ImportError, AttributeError):
        metadata["ase_version"] = "unknown"

    return metadata
```

### Pattern 6: Optimization Step Tracking
**What:** Track actual optimization steps from ASE LBFGS optimizer
**When to use:** In `geometry_optimisation()` and `run_geometry_optimization()`
**Example:**
```python
# Modified geometry_optimisation() in gm_main.py
def geometry_optimisation(mol, fmax=0.000001):
    ei = mol.get_potential_energy()
    print("Initial Energy: ", ei, "eV")
    opt = LBFGS(mol)
    opt.run(fmax=fmax, steps=10000)

    num_steps = opt.get_number_of_steps()  # ASE optimizer API

    ef = mol.get_potential_energy()
    print("Final Energy: ", ef, "eV")
    print(f"Optimization steps: {num_steps}")

    return mol, num_steps  # Return step count alongside molecule
```

### Anti-Patterns to Avoid
- **Catching and silencing exceptions:** Current code has `except Exception as e: logger.warning(...)` patterns that swallow errors. Phase 2 should make these explicit -- either re-raise or return a sentinel value that the caller checks.
- **Validation at the wrong layer:** Don't validate inside deep calculation functions. Validate at entry points (CLI, `run_workflow()`), then trust the data is valid in inner functions.
- **Mixing validation and diagnostics:** `diagnostics.py` prints user-friendly output. `validation.py` should raise exceptions or return structured data. Don't conflate the two.
- **Over-validating:** Don't add type checks that Python already handles. Focus on domain-specific validation (file exists, tool in PATH, atom count reasonable).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Executable path lookup | `subprocess.run(["which", "g16"])` | `shutil.which("g16")` | Stdlib, cross-platform, returns None on not-found |
| CUDA availability | Custom GPU detection | `torch.cuda.is_available()` | Already accurate, handles edge cases |
| XYZ file parsing | Custom parser | `ase.io.read()` with error wrapping | ASE handles all XYZ format variants |
| Optimizer step count | Manual counter/callback | `opt.get_number_of_steps()` | ASE built-in, always correct |
| Version string for tool | Hardcoded string | Read from `pyproject.toml` or `__version__` | Single source of truth |

**Key insight:** Nearly all validation in this phase uses stdlib (`shutil.which`, `pathlib.Path.exists()`) or existing dependency APIs (`torch.cuda.is_available()`, `opt.get_number_of_steps()`). No new libraries needed.

## Common Pitfalls

### Pitfall 1: Breaking Backward Compatibility with Parser Changes
**What goes wrong:** Making `parse_harmonic_frequencies()` raise an exception will break all callers that currently handle empty lists gracefully (e.g., `gm_main.py` line 996 catches parse errors and falls back to empty lists).
**Why it happens:** Callers were written to handle the "return empty" contract.
**How to avoid:** Add `strict` parameter defaulting to `False` for backward compatibility. Only raise when `strict=True`. Update key callers to use `strict=True` where data is expected. Alternatively, make a wrapper function `parse_harmonic_frequencies_strict()` that raises.
**Warning signs:** Tests start failing after parser changes. Existing workflow breaks on edge cases.

### Pitfall 2: Import-Time Validation Blocking Module Load
**What goes wrong:** If validation runs at import time (module level), importing `gm_main` or `cli` will fail if Gaussian isn't installed -- even for commands that don't need Gaussian (like `cli.py list`).
**Why it happens:** The current `gm_main.py` already has some module-level checks (lines 74-77) that log warnings.
**How to avoid:** Only validate at function call time, not module import time. Use lazy validation: check prerequisites at the start of `run_workflow()` or in a Click callback that runs before specific commands.
**Warning signs:** `python cli.py --help` fails because Gaussian isn't loaded.

### Pitfall 3: Subprocess Timeout Killing Gaussian Mid-Calculation
**What goes wrong:** Setting timeout too short kills a legitimate long-running Gaussian anharmonic calculation.
**Why it happens:** Anharmonic frequency calculations for large molecules (aspirin, cocaine) can take hours.
**How to avoid:** Use a generous timeout (e.g., 24 hours for ML-assisted, 72 hours for DFT). Make timeout configurable via CLI option. Log timeout value at INFO level so user knows what to expect. For the ZMQ-based ML calculations in `run_frequency_calculation()`, the timeout should be on the overall process, not per-step.
**Warning signs:** Large molecule calculations suddenly fail with timeout errors.

### Pitfall 4: Version Metadata Breaking JSON Serialization
**What goes wrong:** Adding non-serializable objects (e.g., `torch.device`, `numpy.int64`) to the metadata dict breaks `json.dump()`.
**Why it happens:** PyTorch and numpy types are not JSON-serializable by default.
**How to avoid:** Always convert to Python native types: `str()`, `int()`, `float()`, `bool()`. Test serialization in unit tests.
**Warning signs:** `TypeError: Object of type int64 is not JSON serializable` at results save time.

### Pitfall 5: geometry_optimisation() Signature Change Breaking Callers
**What goes wrong:** Changing `geometry_optimisation()` to return `(mol, num_steps)` tuple breaks the single caller that expects just `mol`.
**Why it happens:** Function signature contract change.
**How to avoid:** Update `run_geometry_optimization()` (the single caller on line 833) simultaneously. This is a contained change since `geometry_optimisation()` is only called in one place.
**Warning signs:** `TypeError: cannot unpack non-sequence Atoms` at the call site.

## Code Examples

Verified patterns from official sources and codebase analysis:

### CLI Pre-Command Validation Hook
```python
# In cli.py - validate before run command executes
@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
# ... other options ...
def run(input_file, ...):
    """Run complete calculation workflow."""
    from validation import (
        validate_all_prerequisites,
        validate_xyz_file,
        detect_device,
    )
    from exceptions import PrerequisiteError, InputValidationError

    # Validate prerequisites before expensive imports
    try:
        prereqs = validate_all_prerequisites(
            check_gaussian=True,
            check_formchk_tool=True,
            dipole_model_path=DEFAULT_MACE_DIPOLE_MODEL,
            helper_script_path=DEFAULT_HELPER_SCRIPT,
        )
        click.echo(f"Prerequisites OK: g16={prereqs['g16_path']}")
    except PrerequisiteError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # Validate input file
    try:
        xyz_info = validate_xyz_file(input_file)
        click.echo(f"Input: {xyz_info['n_atoms']} atoms")
    except InputValidationError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # Detect compute device
    device = detect_device()
    click.echo(f"Device: {device}")

    # ... proceed with workflow ...
```

### Results JSON with Version Metadata (REPR-01, REPR-02)
```python
# Enhanced save_frequency_results in results_manager.py
metadata = {
    "molecule": molecule_name,
    "energy_calculator": energy_calculator,
    "dipole_calculator": dipole_calculator,
    "calculator_type": calculator_type,
    "timestamp": None,
    "frequencies": frequencies_data,
    "energy_eV": float(energy),
    "dipole": dipole,
    "runtime_s": float(runtime),
    "files": files,
    # NEW: reproducibility metadata
    "version_info": collect_version_metadata(),
    "calculation_parameters": {
        "charge": charge,
        "multiplicity": multiplicity,
        "optimization_calculator": optimization_calculator,
        "fmax": fmax,
        "device": device,
    },
}
```

### Subprocess Timeout for Gaussian (ERR-06)
```python
# In gm_main.py run_frequency_calculation() - add overall timeout
# For ML-assisted calculations via ZMQ
import signal

GAUSSIAN_TIMEOUT_SECONDS = 86400  # 24 hours default

# In run_frequency_calculation:
proc = subprocess.Popen(["g16", gjf_temp])
start_time = time.time()

with zmq_server(".ipc_file") as socket:
    while not is_calc_finished(proc, socket):
        # Check timeout
        elapsed = time.time() - start_time
        if elapsed > GAUSSIAN_TIMEOUT_SECONDS:
            proc.kill()
            proc.wait()
            raise TimeoutError(
                f"Gaussian calculation timed out after {elapsed/3600:.1f} hours. "
                f"Consider increasing timeout or simplifying the molecule."
            )

        msg = socket.recv_string()
        # ... rest of loop ...
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `return []` on missing section | Raise `GaussianParseError` with context | This phase | Callers must handle exceptions explicitly |
| `device="cuda"` hardcoded | `detect_device()` at startup | This phase | Works on CPU-only machines |
| No version in results JSON | Full version metadata dict | This phase | Results are reproducible |
| `num_steps=0` hardcoded | `opt.get_number_of_steps()` | This phase | Accurate optimization tracking |
| Warning log on missing g16 | `PrerequisiteError` at startup | This phase | Fails fast with clear message |

**Deprecated/outdated:**
- The current pattern of module-level path checking (gm_main.py lines 74-77) that only logs warnings should be replaced with explicit validation at workflow entry points.

## Open Questions

1. **Timeout values for different molecule sizes**
   - What we know: Water takes ~20s, aspirin could take hours, DFT anharmonic can take days
   - What's unclear: Exact timeout scaling with atom count
   - Recommendation: Make timeout configurable via CLI with generous default (24h for ML, 72h for DFT). Log elapsed time periodically.

2. **How to handle `parse_anharmonic_frequencies()` returning empty**
   - What we know: Some valid logs lack anharmonic data (harmonic-only DFT runs, or ML external logs like acoh)
   - What's unclear: Whether callers always know in advance if anharmonic data should be present
   - Recommendation: Use `strict=False` default. Let callers opt into strict mode when they know data should exist.

3. **Tool version single source of truth**
   - What we know: `pyproject.toml` has `version = "0.2.0"`, `cli.py` has `version="1.0.0"` (inconsistent!)
   - What's unclear: Which is authoritative
   - Recommendation: Use `importlib.metadata.version("mace-gaussian-interface")` at runtime, or define `__version__` in a root `__init__.py`. Fix the inconsistency as part of this phase.

4. **Backward compatibility of results.json schema**
   - What we know: Existing results.json files lack version metadata. Adding new top-level keys is additive.
   - What's unclear: Whether any downstream tools parse results.json strictly
   - Recommendation: Adding keys is safe. Do not remove or rename existing keys. Analysis scripts should use `.get()` with defaults.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `gm_main.py`, `gaussian_parser.py`, `cli.py`, `results_manager.py`, `diagnostics.py`, `fchk_parser.py`, `mace_calculators.py`, `dft_baseline.py`
- [ASE optimizer documentation](https://ase-lib.org/ase/optimize.html) - `get_number_of_steps()` method confirmed
- [Python subprocess documentation](https://docs.python.org/3/library/subprocess.html) - timeout patterns
- [PyTorch CUDA documentation](https://docs.pytorch.org/docs/stable/generated/torch.cuda.is_available.html) - device detection
- `.planning/codebase/CONCERNS.md` - existing tech debt analysis

### Secondary (MEDIUM confidence)
- [Click advanced patterns](https://click.palletsprojects.com/en/stable/advanced/) - CLI validation hooks
- [Python custom exceptions guide](https://jacobpadilla.com/articles/custom-python-exceptions) - hierarchy best practices

### Tertiary (LOW confidence)
- None -- all findings verified against codebase and official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all stdlib + existing dependencies, no new packages
- Architecture: HIGH - patterns verified against codebase structure and existing test patterns from Phase 1
- Pitfalls: HIGH - identified from actual code analysis (specific line numbers, specific failure modes)
- Reproducibility metadata: HIGH - version APIs confirmed from installed packages

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (stable domain, no fast-moving dependencies)
