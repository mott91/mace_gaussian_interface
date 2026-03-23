# Phase 14: Batch Runner & PubChem Fetcher - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-03-23
**Phase:** 14-batch-runner-pubchem-fetcher
**Areas discussed:** Fetch command behavior, Batch file format & invocation, Failure isolation & restart, Output organization

---

## Fetch Command Behavior

### Output location

| Option | Description | Selected |
|--------|-------------|----------|
| Current directory | Saves `aspirin.xyz` in cwd — matches how `run` expects input files | |
| molecules/ subdirectory | Saves to `molecules/aspirin.xyz` — keeps fetched structures separate | ✓ |

**User's choice:** molecules/ subdirectory

### Error handling for ambiguous/missing names

| Option | Description | Selected |
|--------|-------------|----------|
| Fail with clear error | Print error and exit non-zero. No guessing. | ✓ |
| Show top matches and ask | Interactive disambiguation from PubChem results | |

**User's choice:** Fail with clear error

### CID support

| Option | Description | Selected |
|--------|-------------|----------|
| Name only | Keep it simple — name covers 95% of use | ✓ |
| Name + CID | Support --cid flag for unambiguous lookups | |

**User's choice:** Name only

### Overwrite behavior

| Option | Description | Selected |
|--------|-------------|----------|
| Skip with warning | Warn and suggest --force to overwrite | ✓ |
| Always overwrite | Silently replace existing file | |

**User's choice:** Skip with warning

---

## Batch File Format & Invocation

### Input format

| Option | Description | Selected |
|--------|-------------|----------|
| Names only | One name per line, auto-fetch missing | |
| Names + XYZ paths mixed | Each line is name or path | |
| XYZ paths only | Each line is path to existing .xyz | ✓ (via Other) |

**User's choice:** XYZ paths to existing files — user wants fetcher separate from batch workflow entirely
**Notes:** User said "I think I'd like the fetcher separate from the total workflow. So just invoke the fetcher to get the structures and then start the workflow separately."

### CLI options threading

| Option | Description | Selected |
|--------|-------------|----------|
| Same flags as run | Reuse all run command flags — every molecule gets same options | ✓ |
| Per-molecule config file | YAML/JSON config with different options per molecule | |

**User's choice:** Same flags as run

### Directory input

| Option | Description | Selected |
|--------|-------------|----------|
| File only | `mace-gaussian batch molecules.txt` only | ✓ |
| File or directory | Also accept `mace-gaussian batch molecules/` | |

**User's choice:** File only

**Notes:** User also questioned why fetch and batch are in the same phase — they're genuinely independent features. Kept as-is per roadmap.

---

## Failure Isolation & Restart

### Partial results on failure

| Option | Description | Selected |
|--------|-------------|----------|
| Keep partial results | Leave whatever was written, manifest marks 'failed' | ✓ |
| Clean up on failure | Remove molecule's output on failure | |

**User's choice:** Keep partial results

### Completeness definition

| Option | Description | Selected |
|--------|-------------|----------|
| All requested calcs done | Molecule complete only when ALL combos done | |
| Per-calculator tracking | Track each energy×dipole combo separately in manifest | ✓ |

**User's choice:** Per-calculator tracking — finer restart granularity

### Progress reporting

| Option | Description | Selected |
|--------|-------------|----------|
| Simple counter | `[3/25] Running aspirin... done (4m 12s)` | ✓ |
| Detailed per-stage | Report each pipeline stage per molecule | |

**User's choice:** Simple counter

---

## Output Organization

### Results directory

| Option | Description | Selected |
|--------|-------------|----------|
| Existing comparison_results/ | Same structure as `mace-gaussian run` | ✓ |
| Separate batch_results/ | Isolated batch output directory | |

**User's choice:** Existing comparison_results/

### Manifest location

| Option | Description | Selected |
|--------|-------------|----------|
| Output dir root | `comparison_results/batch_manifest.json` | ✓ |
| Project root | `batch_manifest.json` in project root | |

**User's choice:** Output dir root

### Summary output

| Option | Description | Selected |
|--------|-------------|----------|
| Brief summary table | molecule | status | time table at end | ✓ |
| No summary | Just progress lines during execution | |

**User's choice:** Brief summary table

---

## Claude's Discretion

- Manifest JSON schema details
- PubChem API endpoint specifics
- Internal batch loop implementation
- molecules.txt parsing edge cases (blank lines, comments)

## Deferred Ideas

- Splitting fetch and batch into separate phases
- Auto-fetch integration in batch
- Per-molecule config files
- Directory scanning as batch input
- CID-based PubChem lookup
