# Phase 14: Batch Runner & PubChem Fetcher - Context

**Gathered:** 2026-03-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Two independent CLI commands:
1. `mace-gaussian fetch <name>` — download 3D XYZ structure from PubChem by molecule name
2. `mace-gaussian batch molecules.txt` — run the full pipeline over a list of molecules with per-molecule failure isolation, per-calculator restart tracking, and a summary report

These are decoupled: fetch does not auto-integrate into batch. User fetches structures first, then runs batch separately.

</domain>

<decisions>
## Implementation Decisions

### Fetch Command
- **D-01:** Output location is `molecules/` subdirectory (auto-created if needed), not cwd
- **D-02:** Name-only lookup — no CID support. `mace-gaussian fetch aspirin` → `molecules/aspirin.xyz`
- **D-03:** If molecule name is ambiguous or not found on PubChem, fail with clear error and exit non-zero. No interactive disambiguation.
- **D-04:** Skip with warning if file already exists: "aspirin.xyz already exists, skipping. Use --force to overwrite."
- **D-05:** Uses `requests` against PubChem PUG REST API directly (not pubchempy) — decided at v1.1 start

### Batch File Format & Invocation
- **D-06:** `molecules.txt` lists paths to existing `.xyz` files, one per line. Fetch is separate — no auto-fetch in batch.
- **D-07:** File-only input — no directory scanning. `mace-gaussian batch molecules.txt` only.
- **D-08:** Reuses same flags as `run` command (--energy-calculators, --dipole-calculators, --optimization-calculator, --skip-dft-baseline, --output-dir, --keep-scratch). Every molecule gets the same options.

### Failure Isolation & Restart
- **D-09:** Keep partial results on failure — whatever was written to `comparison_results/{molecule}/` survives. Manifest marks it with error details.
- **D-10:** Per-calculator tracking in manifest — each energy×dipole combination tracked separately per molecule. Restart resumes from the exact combination that failed, not the whole molecule.
- **D-11:** Simple progress counter: `[3/25] Running aspirin... done (4m 12s)`. Failed molecules show error summary inline.

### Output Organization
- **D-12:** Results go to existing `comparison_results/` structure (same as `mace-gaussian run`). No separate `batch_results/` dir.
- **D-13:** Manifest lives at `comparison_results/batch_manifest.json` — co-located with results it tracks.
- **D-14:** Brief summary table printed at end: molecule | status | time. Counts: N complete, M failed, K skipped. Points to manifest for details.

### Claude's Discretion
- Manifest JSON schema (exact field names, status enum values)
- PubChem API endpoint details and 3D conformer source
- Internal batch loop implementation (sequential iteration, signal handling)
- How `molecules.txt` handles blank lines, comments, whitespace

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Existing CLI & Pipeline
- `mace_gaussian/cli.py` — Click CLI entry point; `run` command defines the flag interface that `batch` must mirror
- `mace_gaussian/workflow.py` — `run_pipeline()` function that batch will call per-molecule
- `mace_gaussian/utils/results.py` — `ResultsManager` for per-molecule directory structure

### Requirements
- `.planning/REQUIREMENTS.md` §Batch Workflow — BATCH-01 through BATCH-04 acceptance criteria

### Prior Phase Patterns
- `mace_gaussian/utils/scratch.py` — Scratch dir context manager pattern (Phase 13.2)

### PubChem API
- PubChem PUG REST API (external) — 3D conformer endpoint for structure retrieval

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `run_pipeline()` in workflow.py — complete 3-stage pipeline callable per-molecule
- `validate_xyz_file()` in utils/validation.py — input validation for each molecule's .xyz
- `ResultsManager` in utils/results.py — per-molecule directory creation and results JSON
- Click validation callbacks (`_validate_energy_calculators`, `_validate_dipole_calculators`) — reuse for batch flags

### Established Patterns
- Click `@cli.command()` with `@click.option()` decorators for CLI commands
- `scratch_dir()` context manager for temporary file isolation (Phase 13.2)
- `--skip-dft-baseline` flag already threaded through `run_pipeline()`

### Integration Points
- `cli.py` — new `fetch` and `batch` commands added to existing Click group
- `comparison_results/` — batch writes to same directory structure as `run`
- `comparison_results/batch_manifest.json` — new file for batch state tracking

</code_context>

<specifics>
## Specific Ideas

- User explicitly wants fetch and batch decoupled — fetch structures first, then batch separately
- Per-calculator tracking was chosen over per-molecule tracking for finer restart granularity
- User questioned why fetch and batch are in the same phase — they're independent features grouped as "multi-molecule tooling" for the benchmark campaign

</specifics>

<deferred>
## Deferred Ideas

- Splitting fetch and batch into separate phases (user raised the question but decided to proceed as-is)
- Auto-fetch integration in batch (user explicitly rejected — prefers manual fetch-then-batch workflow)
- Per-molecule config files for different calculator options per molecule
- Directory scanning as batch input (`mace-gaussian batch molecules/`)
- CID-based PubChem lookup

</deferred>

---

*Phase: 14-batch-runner-pubchem-fetcher*
*Context gathered: 2026-03-23*
