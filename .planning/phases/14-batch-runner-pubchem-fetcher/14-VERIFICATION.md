---
phase: 14-batch-runner-pubchem-fetcher
verified: 2026-03-23T13:00:00Z
status: passed
score: 8/8 must-haves verified
re_verification: false
---

# Phase 14: Batch Runner + PubChem Fetcher Verification Report

**Phase Goal:** Users can fetch 3D structures by name and run the full pipeline over a list of molecules with per-molecule failure isolation and restart safety.
**Verified:** 2026-03-23T13:00:00Z
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User runs `mace-gaussian fetch aspirin` and receives `molecules/aspirin.xyz` with valid 3D coordinates | VERIFIED | `fetch_3d_structure()` in pubchem.py calls PubChem PUG REST, converts SDF via ASE, writes XYZ; `fetch` CLI command wired in cli.py; test_successful_fetch_writes_xyz passes |
| 2 | Unknown molecule names produce a clear error and exit non-zero | VERIFIED | pubchem.py raises `ValueError` with "not found on PubChem" message on 404/"No CID found"; cli.py catches `ValueError`, prints error, calls `sys.exit(1)`; test_unknown_molecule_raises_valueerror passes |
| 3 | Existing files are skipped with a warning unless --force is passed | VERIFIED | pubchem.py raises `FileExistsError` when file exists and `force=False`; cli.py catches and calls `sys.exit(0)` (not an error); `--force` flag wires to `force=True`; tests pass |
| 4 | User runs `mace-gaussian batch molecules.txt` and each molecule is processed sequentially through geometry opt, DFT baselines, and ML combinations | VERIFIED | batch.py `run_batch()` iterates molecules; calls `run_geometry_optimization`, `run_dft_baselines`, `run_frequency_calculation` from workflow.py; `batch` CLI command wired in cli.py |
| 5 | If the batch run is interrupted and restarted, already-complete calculator combinations are skipped | VERIFIED | batch.py checks `existing.get("status") == STATUS_COMPLETE` before each combo; skipped count incremented; test_restart_skips_complete passes |
| 6 | batch_manifest.json exists after any batch run and records per-calculator status for every molecule | VERIFIED | `manifest_path = Path(output_dir) / "batch_manifest.json"`; `save_manifest` called atomically after every combination; test_manifest_written_after_run passes |
| 7 | User can pass --skip-dft-baseline to skip DFT calculations for all molecules | VERIFIED | `--skip-dft-baseline` option on batch command; `skip_dft_baseline` parameter threads through to `run_batch()`; batch.py gates `run_dft_baselines` call on `not skip_dft_baseline`; tests pass |
| 8 | A summary table is printed at the end with molecule status, timing, and counts | VERIFIED | batch.py prints "BATCH SUMMARY" table with molecule name, ok/failed counts per molecule, per-combo runtime_s, and totals; triggered unconditionally at end of run_batch() |

**Score:** 8/8 truths verified

---

### Required Artifacts

| Artifact | Expected | Lines | Status | Details |
|----------|----------|-------|--------|---------|
| `mace_gaussian/pubchem.py` | PubChem PUG REST API client | 65 (min 40) | VERIFIED | Contains `fetch_3d_structure`, `PUBCHEM_3D_URL`, `requests.get`, `raise FileExistsError`, `"No CID found"` check, `read(StringIO(...), format="sdf")` |
| `tests/test_pubchem.py` | Unit tests with mocked HTTP | 121 (min 60) | VERIFIED | 6 test functions covering success, unknown molecule, no 3D conformer, file-exists, force-overwrite, network timeout |
| `mace_gaussian/batch.py` | Batch runner with manifest-driven restart | 310 (min 150) | VERIFIED | Exports `load_manifest`, `save_manifest`, `parse_batch_file`, `run_batch`; `os.replace` atomic write; `STATUS_COMPLETE = "complete"`, `STATUS_FAILED = "failed"` |
| `tests/test_batch.py` | Unit tests for batch runner and manifest | 353 (min 100) | VERIFIED | 16 test functions covering manifest CRUD, parse, restart, failure isolation, skip-dft, summary dict, manifest written |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `mace_gaussian/cli.py` | `mace_gaussian/pubchem.py` | `from mace_gaussian.pubchem import fetch_3d_structure` | WIRED | Found at cli.py:250 inside `fetch()` command body |
| `mace_gaussian/pubchem.py` | PubChem PUG REST URL | `requests.get` | WIRED | `requests.get(url, timeout=30)` at pubchem.py:47; URL pattern at lines 9-12 |
| `mace_gaussian/batch.py` | `mace_gaussian/workflow.py` | `from .workflow import run_geometry_optimization` etc. | WIRED | Lazy imports at batch.py:194, 207, 218; all three targets (`run_geometry_optimization`, `run_dft_baselines`, `run_frequency_calculation`) confirmed in workflow.py at lines 381, 630, 440 |
| `mace_gaussian/batch.py` | `comparison_results/batch_manifest.json` | `os.replace` atomic write | WIRED | `os.replace(tmp, str(path))` at batch.py:66; manifest_path set to `Path(output_dir) / "batch_manifest.json"` at batch.py:154 |
| `mace_gaussian/cli.py` | `mace_gaussian/batch.py` | `from mace_gaussian.batch import run_batch` | WIRED | Found at cli.py:540 inside `batch()` command body |

---

### Data-Flow Trace (Level 4)

Not applicable. These are CLI command handlers and API client utilities, not components that render dynamic data from a store. The data flows are verified via behavioral tests (all 22 pass) and key link checks above.

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| All pubchem and batch tests pass | `pytest tests/test_pubchem.py tests/test_batch.py -x -v` | 22 passed, 0 failed in 2.22s | PASS |
| CLI registers `fetch` and `batch` commands | `python -c "from mace_gaussian.cli import cli; assert 'fetch' in cli.commands; assert 'batch' in cli.commands"` | OK | PASS |
| `fetch_3d_structure` importable from pubchem | `python -c "from mace_gaussian.pubchem import fetch_3d_structure"` | OK | PASS |
| `load_manifest`, `save_manifest`, `run_batch`, `parse_batch_file` importable from batch | `python -c "from mace_gaussian.batch import ..."` | OK | PASS |
| `batch --help` shows `--skip-dft-baseline`, `--keep-scratch`, `--output-dir`, `--energy-calculators`, `--dipole-calculators` | Click test runner `batch --help` | All 5 options present | PASS |
| `fetch --help` shows `--force`, `--output-dir` | Click test runner `fetch --help` | Both options present | PASS |
| No ruff lint errors in modified files | `ruff check mace_gaussian/pubchem.py mace_gaussian/batch.py mace_gaussian/cli.py` | All checks passed | PASS |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| BATCH-01 | 14-01-PLAN.md | User can run `mace-gaussian fetch <molecule-name>` to download a 3D XYZ structure from PubChem | SATISFIED | `pubchem.py` + `cli.py` fetch command; 6 unit tests pass |
| BATCH-02 | 14-02-PLAN.md | User can run `mace-gaussian batch molecules.txt` to process multiple molecules sequentially through the full pipeline | SATISFIED | `batch.py` + `cli.py` batch command; test_basic_two_molecules passes |
| BATCH-03 | 14-02-PLAN.md | Batch run produces a per-molecule status manifest (batch_manifest.json) that survives interruption — restarting skips already-complete molecules | SATISFIED | Manifest written atomically after each combo; test_restart_skips_complete passes |
| BATCH-04 | 14-02-PLAN.md | User can run `mace-gaussian batch molecules.txt --skip-dft-baseline` to run ML calculations only | SATISFIED | `--skip-dft-baseline` flag wired; test_skip_dft_baseline_true and _false pass |

**Orphaned requirements check:** BATCH-05 is listed in REQUIREMENTS.md as Phase 15 (pending). No plan in Phase 14 claimed it. This is correct — not orphaned for Phase 14.

All four Phase 14 requirement IDs (BATCH-01 through BATCH-04) are claimed, implemented, and verified.

---

### Anti-Patterns Found

None. Scan of `pubchem.py`, `batch.py`, and `cli.py` found:
- No TODO/FIXME/HACK/PLACEHOLDER comments
- No empty stub returns (`return []`, `return {}`, `return null`)
- No console.log-only handlers
- No hardcoded empty props

One noted deviation from Plan 02 (non-blocking): `cleanup_stale_scratch_dirs` was omitted from the batch CLI because `mace_gaussian/utils/scratch.py` is not present on this branch (added in Phase 13.2 which is not yet merged). This is correctly documented in 14-02-SUMMARY.md and does not affect batch runner correctness.

---

### Human Verification Required

The following behaviors require a live network connection and cannot be verified programmatically:

**1. PubChem fetch end-to-end**
**Test:** Run `mace-gaussian fetch aspirin` in the project directory
**Expected:** `molecules/aspirin.xyz` is created with a valid 3D structure (21 atoms for aspirin)
**Why human:** Requires live network access to pubchem.ncbi.nlm.nih.gov; all automated tests use mocked HTTP

**2. Batch pipeline end-to-end with a real molecule**
**Test:** Create `molecules.txt` with a path to `molecules/water.xyz`, run `mace-gaussian batch molecules.txt --skip-dft-baseline --energy-calculators mace_mp --dipole-calculators espaloma`
**Expected:** `comparison_results/batch_manifest.json` is created; summary table is printed; exit code reflects success/failure
**Why human:** Requires MACE models loaded on GPU and Gaussian 16 in PATH; full integration path cannot be unit-tested

---

### Gaps Summary

No gaps. All 8 observable truths are verified. All 4 artifacts pass all three levels (exist, substantive, wired). All 5 key links are confirmed. All 4 requirement IDs are satisfied. No blocker anti-patterns found. 22/22 unit tests pass.

---

_Verified: 2026-03-23T13:00:00Z_
_Verifier: Claude (gsd-verifier)_
