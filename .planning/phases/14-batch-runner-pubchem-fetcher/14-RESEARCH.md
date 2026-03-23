# Phase 14: Batch Runner & PubChem Fetcher - Research

**Researched:** 2026-03-23
**Domain:** CLI commands, PubChem REST API, batch orchestration with manifest-based restart
**Confidence:** HIGH

## Summary

Phase 14 adds two independent CLI commands to the existing Click group: `fetch` (download 3D structures from PubChem) and `batch` (run the full pipeline over a list of molecules with per-calculator failure isolation and restart). Both are well-constrained features that compose existing infrastructure -- `run_pipeline()` in workflow.py, `ResultsManager`, Click CLI patterns, and the `requests` library for HTTP.

The PubChem PUG REST API provides a single-request path from molecule name to 3D SDF: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF?record_type=3d`. This was verified live -- aspirin returns a valid 21-atom SDF that ASE reads and converts to XYZ correctly. Error cases (unknown molecule: 404 "No CID found"; known molecule without 3D conformer: 404 "No records found for the given CID(s)") return distinct error messages that the fetch command should surface.

The batch command wraps `run_pipeline()` calls in a loop with a JSON manifest for restart tracking. Per-calculator granularity (D-10) means the manifest must track each energy x dipole combination independently per molecule, not just molecule-level pass/fail.

**Primary recommendation:** Implement fetch as a thin `requests` + ASE I/O wrapper (~60 lines), and batch as a manifest-driven loop around the existing `run_pipeline()` with per-combination status tracking (~200 lines). Add `requests` to pyproject.toml dependencies.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Output location is `molecules/` subdirectory (auto-created if needed), not cwd
- **D-02:** Name-only lookup -- no CID support. `mace-gaussian fetch aspirin` -> `molecules/aspirin.xyz`
- **D-03:** If molecule name is ambiguous or not found on PubChem, fail with clear error and exit non-zero. No interactive disambiguation.
- **D-04:** Skip with warning if file already exists: "aspirin.xyz already exists, skipping. Use --force to overwrite."
- **D-05:** Uses `requests` against PubChem PUG REST API directly (not pubchempy) -- decided at v1.1 start
- **D-06:** `molecules.txt` lists paths to existing `.xyz` files, one per line. Fetch is separate -- no auto-fetch in batch.
- **D-07:** File-only input -- no directory scanning. `mace-gaussian batch molecules.txt` only.
- **D-08:** Reuses same flags as `run` command (--energy-calculators, --dipole-calculators, --optimization-calculator, --skip-dft-baseline, --output-dir, --keep-scratch). Every molecule gets the same options.
- **D-09:** Keep partial results on failure -- whatever was written to `comparison_results/{molecule}/` survives. Manifest marks it with error details.
- **D-10:** Per-calculator tracking in manifest -- each energy x dipole combination tracked separately per molecule. Restart resumes from the exact combination that failed, not the whole molecule.
- **D-11:** Simple progress counter: `[3/25] Running aspirin... done (4m 12s)`. Failed molecules show error summary inline.
- **D-12:** Results go to existing `comparison_results/` structure (same as `mace-gaussian run`). No separate `batch_results/` dir.
- **D-13:** Manifest lives at `comparison_results/batch_manifest.json` -- co-located with results it tracks.
- **D-14:** Brief summary table printed at end: molecule | status | time. Counts: N complete, M failed, K skipped. Points to manifest for details.

### Claude's Discretion
- Manifest JSON schema (exact field names, status enum values)
- PubChem API endpoint details and 3D conformer source
- Internal batch loop implementation (sequential iteration, signal handling)
- How `molecules.txt` handles blank lines, comments, whitespace

### Deferred Ideas (OUT OF SCOPE)
- Splitting fetch and batch into separate phases
- Auto-fetch integration in batch
- Per-molecule config files for different calculator options per molecule
- Directory scanning as batch input
- CID-based PubChem lookup
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| BATCH-01 | User can run `mace-gaussian fetch <molecule-name>` to download a 3D XYZ structure from PubChem | PUG REST API verified: single URL `compound/name/{name}/record/SDF?record_type=3d` returns 3D SDF; ASE reads SDF and writes XYZ. |
| BATCH-02 | User can run `mace-gaussian batch molecules.txt` to process multiple molecules sequentially through the full pipeline | `run_pipeline()` already callable per-molecule; batch wraps it in a loop with shared flags from D-08. |
| BATCH-03 | Batch run produces a per-molecule status manifest (`batch_manifest.json`) that survives interruption -- restarting skips already-complete molecules | JSON manifest with per-calculator granularity (D-10); load on start, write after each combination completes. |
| BATCH-04 | User can run `mace-gaussian batch molecules.txt --skip-dft-baseline` to run ML calculations only | `run_pipeline()` already accepts `include_dft_baselines` param; batch passes `--skip-dft-baseline` through. |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| click | (existing) | CLI framework | Already used for all commands in cli.py |
| requests | 2.27+ | HTTP client for PubChem API | Decision D-05: use requests directly, not pubchempy |
| ase | (existing) | SDF-to-XYZ conversion via `ase.io.read/write` | Already a core dependency; handles SDF format natively |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| json (stdlib) | - | Manifest serialization | Read/write batch_manifest.json |
| time (stdlib) | - | Timing each molecule run | Progress counter (D-11) |
| signal (stdlib) | - | Graceful SIGINT handling in batch loop | Optional: save manifest on Ctrl+C |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| requests for PubChem | pubchempy | Explicitly rejected (D-05); requests is lighter, more control |
| ASE for SDF->XYZ | RDKit Chem.MolToXYZBlock | ASE already in deps, simpler; RDKit also available but unnecessary |
| JSON manifest | SQLite | JSON is human-readable, inspectable, fits the scale (~25 molecules) |

**Installation:**
```bash
# requests must be added to pyproject.toml dependencies
# All other libraries are already project dependencies
```

**Note:** `requests` is available in the working environment (v2.27.1) but is NOT listed in pyproject.toml `dependencies`. It must be added.

## Architecture Patterns

### New Files
```
mace_gaussian/
├── cli.py                    # Add fetch + batch commands (existing file)
├── pubchem.py                # PubChem API client (NEW, ~60 lines)
└── batch.py                  # Batch runner + manifest (NEW, ~200 lines)
```

### Pattern 1: Fetch Command
**What:** Thin wrapper: PubChem API -> SDF bytes -> ASE read -> XYZ write
**When to use:** `mace-gaussian fetch aspirin`
**Example:**
```python
# mace_gaussian/pubchem.py
import requests
from io import StringIO
from ase.io import read, write
from pathlib import Path

PUBCHEM_3D_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF"
    "?record_type=3d"
)

def fetch_3d_structure(molecule_name: str, output_dir: Path, force: bool = False) -> Path:
    """Fetch 3D structure from PubChem and save as XYZ."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{molecule_name}.xyz"

    if output_path.exists() and not force:
        raise FileExistsError(
            f"{molecule_name}.xyz already exists, skipping. Use --force to overwrite."
        )

    url = PUBCHEM_3D_URL.format(name=molecule_name)
    response = requests.get(url, timeout=30)

    if response.status_code == 404:
        # PubChem returns distinct messages for "no CID" vs "no 3D conformer"
        raise ValueError(f"Could not fetch '{molecule_name}' from PubChem: {response.text.strip()}")
    response.raise_for_status()

    atoms = read(StringIO(response.text), format="sdf")
    write(str(output_path), atoms, format="xyz")
    return output_path
```

### Pattern 2: Manifest-Driven Batch Loop
**What:** Load manifest -> for each molecule: for each calculator combination: check manifest -> run if needed -> update manifest -> save
**When to use:** `mace-gaussian batch molecules.txt`
**Example:**
```python
# mace_gaussian/batch.py
import json
from pathlib import Path

STATUS_PENDING = "pending"
STATUS_COMPLETE = "complete"
STATUS_FAILED = "failed"
STATUS_SKIPPED = "skipped"  # for restart: already done

def load_manifest(path: Path) -> dict:
    if path.exists():
        with path.open() as f:
            return json.load(f)
    return {"molecules": {}}

def save_manifest(manifest: dict, path: Path) -> None:
    with path.open("w") as f:
        json.dump(manifest, f, indent=2)
```

### Pattern 3: Per-Calculator Granularity in Manifest
**What:** Each molecule entry in the manifest tracks individual energy x dipole combinations
**Why:** D-10 requires restart to resume from the exact failed combination
**Example manifest schema:**
```json
{
  "batch_file": "molecules.txt",
  "started": "2026-03-23T10:00:00",
  "updated": "2026-03-23T10:45:00",
  "options": {
    "energy_calculators": ["mace_mp", "mace_off"],
    "dipole_calculators": ["espaloma", "mace_ml"],
    "skip_dft_baseline": false,
    "optimization_calculator": "mace_omol"
  },
  "molecules": {
    "water": {
      "xyz_path": "molecules/water.xyz",
      "geometry_opt": "complete",
      "dft_baseline": "complete",
      "combinations": {
        "mace_mp_espaloma": {"status": "complete", "runtime_s": 245.3},
        "mace_mp_mace_ml": {"status": "complete", "runtime_s": 189.7},
        "mace_off_espaloma": {"status": "failed", "error": "mace_off element guard: molecule contains [Br]", "runtime_s": 1.2},
        "mace_off_mace_ml": {"status": "pending"}
      }
    }
  }
}
```

### Pattern 4: Reuse run_pipeline Internals, Not run_pipeline Itself
**What:** For per-calculator restart, batch cannot call `run_pipeline()` as a black box because it runs ALL combinations. Instead, batch should call the three stage functions directly: `run_geometry_optimization()`, `run_dft_baselines()`, `run_ml_calculations()` -- or more precisely, call `run_frequency_calculation()` per-combination to get per-combination restart.
**Why:** `run_pipeline()` iterates all energy x dipole combinations internally. To skip already-complete combinations, batch needs to call `run_frequency_calculation()` individually per combination, checking the manifest before each call.
**Key insight:** The three-stage pipeline structure in workflow.py is: (1) geometry opt, (2) DFT baselines, (3) ML combinations. Stages 1 and 2 are molecule-level (run once per molecule). Stage 3 is the per-combination loop. Batch should call stages 1+2 via existing functions, then loop stage 3 per-combination with manifest checks.

### Anti-Patterns to Avoid
- **Calling run_pipeline() as black box:** Loses per-calculator restart granularity. Must decompose into stage calls.
- **Global try/except around entire batch:** Loses per-molecule isolation. Each molecule must have its own try/except.
- **Writing manifest only at end:** If interrupted mid-batch, all progress is lost. Write after EACH combination completes.
- **Storing absolute paths in manifest:** Makes manifest non-portable. Store paths relative to output_dir.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| HTTP client | urllib wrapper | `requests` | D-05 explicitly chose requests; handles redirects, encoding, timeouts |
| SDF parsing | Custom SDF parser | `ase.io.read(format='sdf')` | ASE handles V2000/V3000 SDF variants, atom typing |
| XYZ writing | String formatting | `ase.io.write(format='xyz')` | Consistent with project's existing XYZ I/O |
| CLI argument parsing | argparse | Click | Project standard; reuse existing validator callbacks |
| JSON serialization | Custom format | `json.dump/load` | Manifest is simple nested dict; JSON is sufficient |

**Key insight:** This phase composes existing infrastructure. The only new external interaction is the PubChem HTTP call. Everything else (pipeline stages, results directories, CLI patterns) already exists.

## Common Pitfalls

### Pitfall 1: PubChem 3D Conformer Not Available
**What goes wrong:** Some molecules have a PubChem CID but no pre-computed 3D conformer. The API returns 404 with "No records found for the given CID(s)" (distinct from "No CID found").
**Why it happens:** PubChem only stores 3D conformers for a subset of compounds (smaller organic molecules).
**How to avoid:** Parse the 404 response body to distinguish "molecule not found" from "no 3D data available". Provide distinct error messages for each case.
**Warning signs:** Tested with buckminsterfullerene -- has CID but 404 on 3D request.

### Pitfall 2: PubChem Rate Limiting
**What goes wrong:** PubChem limits requests to 5 per second. Batch fetch of many molecules could trigger a block.
**Why it happens:** NCBI rate limits are enforced globally.
**How to avoid:** For `fetch`, this is not an issue (single request per invocation). If batch-fetch is ever added later, add `time.sleep(0.25)` between requests. Current design (D-06) separates fetch from batch, so this is low risk.
**Warning signs:** HTTP 503 responses from PubChem.

### Pitfall 3: Manifest Corruption on Interrupt
**What goes wrong:** If `json.dump` is interrupted mid-write, the manifest file is left in a corrupted state (partial JSON).
**Why it happens:** Ctrl+C during file write.
**How to avoid:** Write to a temporary file, then `os.replace()` (atomic on POSIX). This is a single-line change: `with tempfile -> json.dump -> os.replace(tmp, manifest_path)`.
**Warning signs:** `json.load` raises `JSONDecodeError` on restart.

### Pitfall 4: Batch Restart with Different Calculator Options
**What goes wrong:** User runs batch with `--energy-calculators mace_mp`, then restarts with `--energy-calculators mace_mp,mace_off`. The manifest has completions for mace_mp but not mace_off.
**How to avoid:** On restart, compare requested combinations against manifest. Combinations not in the manifest are treated as pending. Log a warning if manifest options differ from current invocation.
**Warning signs:** "Skipped" count is unexpectedly high or low.

### Pitfall 5: run_pipeline Calls All Combinations Internally
**What goes wrong:** Calling `run_pipeline()` per molecule re-runs ALL energy x dipole combinations, ignoring the manifest's per-combination status.
**Why it happens:** `run_pipeline()` -> `run_ml_calculations()` iterates the full Cartesian product.
**How to avoid:** Batch must NOT call `run_pipeline()` as a black box. It must call the three stages separately: `run_geometry_optimization()`, `run_dft_baselines()`, and `run_frequency_calculation()` (per combination, with manifest check).
**Warning signs:** "Already complete" molecules re-run all their ML calculations.

### Pitfall 6: molecules.txt Path Resolution
**What goes wrong:** If `molecules.txt` contains relative paths like `molecules/water.xyz`, they must resolve relative to the current working directory, not relative to the text file's location.
**Why it happens:** Ambiguous path semantics.
**How to avoid:** Resolve all paths relative to `os.getcwd()`. Document this behavior.
**Warning signs:** FileNotFoundError when molecules.txt is in a different directory.

## Code Examples

### PubChem API: Name to 3D SDF (Verified Live)
```python
# Source: Verified against PubChem PUG REST API 2026-03-23
import requests
from io import StringIO
from ase.io import read, write

url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/record/SDF?record_type=3d"
response = requests.get(url, timeout=30)
# Returns: 200 with SDF content (21 atoms for aspirin)

# Error cases:
# Unknown molecule: 404, body contains "No CID found"
# No 3D conformer: 404, body contains "No records found for the given CID(s)"

atoms = read(StringIO(response.text), format="sdf")
write("molecules/aspirin.xyz", atoms, format="xyz")
# Produces valid XYZ with 3D coordinates
```

### Click Command: Reuse Existing Validators
```python
# Source: existing cli.py pattern
@cli.command()
@click.argument("molecule_name")
@click.option("--force", is_flag=True, help="Overwrite existing file")
def fetch(molecule_name, force):
    """Fetch 3D structure from PubChem by molecule name."""
    from mace_gaussian.pubchem import fetch_3d_structure
    # ...
```

### Batch Command: Mirror run's Options
```python
# Source: existing cli.py run command pattern
@cli.command()
@click.argument("batch_file", type=click.Path(exists=True))
@click.option("--optimization-calculator", default="mace_omol",
              type=click.Choice(["mace_omol", "mace_off", "mace_mp", "mace_anicc", "mace_polar"]))
@click.option("--energy-calculators", default="mace_mp,mace_omol,mace_anicc,mace_off,mace_polar",
              callback=_validate_energy_calculators)
@click.option("--dipole-calculators", default="espaloma,mace_ml",
              callback=_validate_dipole_calculators)
@click.option("--skip-dft-baseline", is_flag=True)
@click.option("--output-dir", default="comparison_results", type=click.Path())
@click.option("--keep-scratch", is_flag=True, default=False)
def batch(batch_file, optimization_calculator, energy_calculators, dipole_calculators,
          skip_dft_baseline, output_dir, keep_scratch):
    """Run pipeline for multiple molecules listed in BATCH_FILE."""
    # ...
```

### Atomic Manifest Write
```python
# Source: POSIX best practice for crash-safe JSON writes
import json, os, tempfile

def save_manifest(manifest: dict, path: Path) -> None:
    """Write manifest atomically to prevent corruption on interrupt."""
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(manifest, f, indent=2)
        os.replace(tmp, str(path))
    except BaseException:
        os.unlink(tmp)
        raise
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| pubchempy library | Direct requests to PUG REST | D-05 decision | Fewer dependencies, more control over error handling |
| Molecule-level restart | Per-calculator restart (D-10) | Phase 14 design | Finer granularity, less re-work on restart |
| batch_results/ separate dir | comparison_results/ (shared with run) | D-12 decision | Consistent output structure, no duplicate results |

## Open Questions

1. **Signal handling for graceful shutdown**
   - What we know: Python's `signal.signal(SIGINT, handler)` can save manifest before exit
   - What's unclear: Whether to catch SIGINT or just rely on atomic manifest writes
   - Recommendation: Use atomic manifest writes (Pattern in Code Examples). If SIGINT arrives between manifest saves, at most one combination's status is lost -- it re-runs on restart. Signal handling adds complexity without much benefit given atomic writes.

2. **molecules.txt comment syntax**
   - What we know: D-06 says "one per line"
   - What's unclear: Whether to support `#` comments and blank lines
   - Recommendation: Support `#` comments (strip from `#` onward) and skip blank lines. This is Claude's discretion per CONTEXT.md and costs ~3 lines of code.

3. **Manifest options drift detection**
   - What we know: Manifest stores the options used for the batch run
   - What's unclear: What to do if user restarts with different options
   - Recommendation: Log a warning ("Manifest was created with different options") but proceed. New combinations are treated as pending. This handles the common case of adding calculators to a running campaign.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 7.0+ |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `pytest tests/ -x -q` |
| Full suite command | `pytest tests/ -ra` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| BATCH-01 | fetch downloads 3D XYZ from PubChem | unit (mock requests) | `pytest tests/test_pubchem.py -x` | Wave 0 |
| BATCH-02 | batch processes molecules sequentially | unit (mock pipeline) | `pytest tests/test_batch.py::test_batch_sequential -x` | Wave 0 |
| BATCH-03 | manifest tracks status, restart skips complete | unit | `pytest tests/test_batch.py::test_manifest_restart -x` | Wave 0 |
| BATCH-04 | --skip-dft-baseline passes through to pipeline | unit | `pytest tests/test_batch.py::test_skip_dft_flag -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest tests/test_pubchem.py tests/test_batch.py -x -q`
- **Per wave merge:** `pytest tests/ -ra`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/test_pubchem.py` -- covers BATCH-01 (fetch command, PubChem API mocking, error cases)
- [ ] `tests/test_batch.py` -- covers BATCH-02, BATCH-03, BATCH-04 (batch loop, manifest CRUD, restart logic, flag passthrough)

## Sources

### Primary (HIGH confidence)
- PubChem PUG REST API -- verified live 2026-03-23: name-to-3D-SDF endpoint, error responses for unknown molecules and missing 3D conformers
- ASE SDF format support -- verified via `ase.io.formats.ioformats['sdf']` and live read/write test
- Existing codebase: cli.py, workflow.py, utils/results.py, utils/scratch.py -- read in full

### Secondary (MEDIUM confidence)
- [IUPAC FAIR Chemistry Cookbook - PUG REST](https://iupac.github.io/WFChemCookbook/datasources/pubchem_pugrest1.html) -- rate limit documentation (5 req/sec)
- [PubGrep by grimme-lab](https://github.com/grimme-lab/PubGrep) -- reference implementation for PubChem structure fetching

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all libraries verified available, API tested live
- Architecture: HIGH - composing well-understood existing code patterns
- Pitfalls: HIGH - error cases verified with live API calls, manifest atomicity is standard POSIX pattern

**Research date:** 2026-03-23
**Valid until:** 2026-04-23 (stable domain, PubChem API is versioned)
