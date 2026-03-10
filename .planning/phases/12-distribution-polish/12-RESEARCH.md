# Phase 12: Distribution Polish - Research

**Researched:** 2026-02-27
**Domain:** Python packaging metadata, CI version pinning, documentation entry points
**Confidence:** HIGH

## Summary

Phase 12 is a purely cosmetic/polish phase — no new code, no new tests, no behavior changes. It closes four concrete audit findings from the v1.0 milestone audit: placeholder GitHub URLs in `pyproject.toml`, a ruff version skew between CI and dev dependencies, and `python cli.py` references in docs that should read `mace-gaussian` (the installed entry point). Success criteria 4 and 5 (REQUIREMENTS.md checkboxes and ROADMAP.md progress table) were already completed during the planning metadata update and MUST NOT appear as tasks.

The real GitHub repository is confirmed from `git remote get-url origin`: `git@github.com:mott91/mace_gaussian_interface.git`. The real owner is `mott91` and the real repo slug is `mace_gaussian_interface`. All three `[project.urls]` fields in `pyproject.toml` use `yourusername/mace-gaussian-interface` and need to be replaced with the real coordinates.

The ruff version skew is straightforward: CI pins `ruff==0.15.1` (exact) while `pyproject.toml` declares `ruff>=0.9.0` (floor). The fix is to change the CI pin to `ruff>=0.9.0` (or align both to the same concrete version). Using `>=` in CI is standard practice when the goal is to catch real lint regressions, not to freeze a specific formatter behavior. The current pin was set in Phase 9-02 for "CI reproducibility" — that rationale applies to the formatter output but creates a maintenance burden. The right fix depends on intent (see Architecture Patterns below).

The docs entry-point issue appears in four places: `docs/quickstart.md` (4 occurrences of `python cli.py`), `README.md` (8+ occurrences), and `CLAUDE.md` (5 occurrences). The installed entry point `mace-gaussian` is already wired in `pyproject.toml` via `[project.scripts]`. The fix is to prepend `mace-gaussian` usage as the primary form and demote or remove `python cli.py`.

**Primary recommendation:** One plan (12-01) with four sequential file edits: pyproject.toml URLs, ci.yml ruff pin, docs/quickstart.md entry point, README.md entry point. Treat CLAUDE.md separately since it is developer-facing and `python cli.py` may be intentionally kept there for direct-file debugging convenience.

## Current State Inventory

This is a file-editing phase. The table below captures exactly what exists and what must change.

### pyproject.toml — Placeholder URLs

| Field | Current value (placeholder) | Target value (real) |
|-------|---------------------------|---------------------|
| `Homepage` | `https://github.com/yourusername/mace-gaussian-interface` | `https://github.com/mott91/mace_gaussian_interface` |
| `Repository` | `https://github.com/yourusername/mace-gaussian-interface.git` | `https://github.com/mott91/mace_gaussian_interface.git` |
| `Issues` | `https://github.com/yourusername/mace-gaussian-interface/issues` | `https://github.com/mott91/mace_gaussian_interface/issues` |

Verified from: `git remote get-url origin` → `git@github.com:mott91/mace_gaussian_interface.git`

Note: The `[project]` `name` field is `mace-gaussian-interface` — that is the PyPI package name (a separate concept from the GitHub repo slug). It can stay as-is unless the user wants to rename it; the audit only flagged the URL fields.

### .github/workflows/ci.yml — Ruff Version Skew

| Location | Current | Issue |
|----------|---------|-------|
| `ci.yml` line 18 | `pip install "ruff==0.15.1"` | Exact pin, diverges from pyproject.toml |
| `pyproject.toml` line 42 | `"ruff>=0.9.0"` | Floor only |

Two valid fixes:

**Option A (recommended):** Change CI to `pip install "ruff>=0.9.0"` — matches pyproject.toml exactly, no skew possible. Risk: if a future ruff minor version changes formatting rules, CI output could differ from local runs. Acceptable for a research project where ruff has been stable.

**Option B:** Change pyproject.toml to `"ruff==0.15.1"` and keep CI in sync — pins both. Downside: uv lock update needed on every ruff release, more maintenance. The Phase 9-02 decision log says "ruff==0.15.1 pinned to match uv.lock for CI reproducibility." If the user wants strict pinning, Option B is correct. If they want low-maintenance, Option A.

**Planner action:** Default to Option A (match pyproject.toml floor in CI) unless the user explicitly prefers Option B. The audit labeled this LOW risk; strict pinning is not required for a thesis codebase.

### docs/quickstart.md — Entry Point References

| Line | Current | Should be |
|------|---------|-----------|
| 14 | `python cli.py diagnose` | `mace-gaussian diagnose` |
| 24 | `python cli.py run water.xyz --skip-dft-baseline` | `mace-gaussian run water.xyz --skip-dft-baseline` |
| 33 | `python cli.py run water.xyz` | `mace-gaussian run water.xyz` |
| 40 | `python cli.py list water` | `mace-gaussian list water` |

All four are in the main workflow steps (Steps 1-3). The quickstart is the user-facing guide; this is the highest-priority doc to fix.

Note: Line 14 also says "verify with `python cli.py diagnose`" — this is inside the Installation step link (`README.md#installation`), not directly in quickstart. But line 14 of quickstart.md shows this text literally. Fix it.

### README.md — Entry Point References

| Section | Occurrences of `python cli.py` |
|---------|-------------------------------|
| Step 5 (Verify installation) | 1 |
| Quick Example | 2 |
| Step-by-Step Usage | 5 (run, run --skip, run --custom, list, list water) |

Total: 8 occurrences. Each should become `mace-gaussian <subcommand>`.

The `python run_analysis.py water` references are intentional — `run_analysis.py` is NOT replaced by the entry point. Only `cli.py` calls move to `mace-gaussian`.

### CLAUDE.md — Developer-Facing References

CLAUDE.md contains 5 occurrences of `python cli.py`. These are in the `## Commands` section, which is developer/contributor documentation. Two options:

**Option A (recommended):** Update to `mace-gaussian` — consistent with how the tool is meant to be used post-install. A developer who has done `pip install -e .` will have `mace-gaussian` in PATH.

**Option B:** Leave as-is — `python cli.py` still works from the project root. Developer docs may intentionally keep the raw invocation for clarity.

**Planner action:** Update CLAUDE.md to `mace-gaussian` for consistency, but this is lower priority than quickstart.md and README.md. The planner may include this in 12-01 or defer it.

## Architecture Patterns

### Pattern 1: Single-file Edits (No Logic Changes)

All four fix areas are string substitutions in configuration or documentation files. No imports, no test changes, no logic. The plan structure should be:

```
Task 1: Fix pyproject.toml URLs (3 line replacements)
Task 2: Fix CI ruff pin (1 line replacement)
Task 3: Fix docs/quickstart.md entry point (4 line replacements)
Task 4: Fix README.md entry point (8 line replacements, leave run_analysis.py untouched)
Task 5 (optional): Fix CLAUDE.md entry point (5 line replacements)
```

Each task is ~2-5 minute execution. All can be in a single plan (12-01).

### Pattern 2: Verify After Edit

For documentation edits, verification is a grep for remaining `python cli.py` occurrences in the target file. For pyproject.toml URLs, verification is reading the updated file. No tests need to run.

```bash
# Verify pyproject.toml URL fix
grep "yourusername" pyproject.toml  # Should return nothing

# Verify docs entry point fix
grep "python cli.py" docs/quickstart.md  # Should return nothing
grep "python cli.py" README.md  # Should return nothing

# Verify CI pin fix
grep "ruff==" .github/workflows/ci.yml  # Should return nothing (or the aligned version)
```

### Anti-Patterns to Avoid

- **Renaming the pyproject.toml `[project] name` field**: The `name = "mace-gaussian-interface"` is the PyPI package name, not the GitHub repo name. Do NOT change it to match the repo slug. Leave it as-is.
- **Changing `run_analysis.py` references**: `python run_analysis.py water` is correct — that file is not replaced by an entry point. Only `cli.py` → `mace-gaussian`.
- **Creating tasks for criteria 4 and 5**: REQUIREMENTS.md checkboxes and ROADMAP.md progress table are already updated. Do NOT re-do this work.
- **Inventing the GitHub URL**: The real URL is confirmed from git remote. Do NOT guess or leave as placeholder.
- **Changing pyproject.toml `[project] name`**: Different from GitHub repo name. Keep `mace-gaussian-interface`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| URL validation | Custom URL checker | Visual inspection + grep | It is a 3-line string replacement; over-engineering is wasteful |
| Ruff version discovery | Script to query PyPI | Read pyproject.toml floor | The floor is the authoritative spec |
| Entry point detection | Script to grep all docs | Direct targeted edits | Scope is known: 3 files, known line counts |

## Common Pitfalls

### Pitfall 1: Changing the PyPI Package Name
**What goes wrong:** Editor accidentally changes `name = "mace-gaussian-interface"` to match repo slug.
**Why it happens:** The GitHub repo name `mace_gaussian_interface` looks like a natural candidate.
**How to avoid:** Only touch `[project.urls]` — not `[project]` `name` field.
**Warning signs:** `name =` value changes in the diff.

### Pitfall 2: Removing `python cli.py` From CLAUDE.md Too Aggressively
**What goes wrong:** CLAUDE.md `## Commands` section is gutted, losing developer context.
**Why it happens:** Blanket find-replace across all files.
**How to avoid:** Update CLAUDE.md `## Commands` section to show `mace-gaussian`, but verify the rest of CLAUDE.md is not affected.

### Pitfall 3: Breaking run_analysis.py References
**What goes wrong:** `python run_analysis.py water` → `mace-gaussian run_analysis water` (wrong — no such command).
**Why it happens:** Broad substitution of `python *.py`.
**How to avoid:** Only replace `python cli.py` → `mace-gaussian`, leave `python run_analysis.py` unchanged.

### Pitfall 4: Orphaned `python cli.py` In Inline Code Blocks
**What goes wrong:** grep after edit finds remaining occurrences in code block comments or examples.
**Why it happens:** Replacement missed a code block inside a different section.
**How to avoid:** Run `grep -n "python cli.py" <file>` after edits to confirm zero remaining.

### Pitfall 5: CI Pin Leaves Skew the Other Direction
**What goes wrong:** CI becomes `ruff>=0.9.0` but a future ruff major bump changes formatting in a way that fails format --check locally but not in CI (or vice versa).
**Why it happens:** Floor-only constraints are imprecise.
**How to avoid:** This is acceptable risk for a research project. Document the tradeoff in STATE.md decisions. If the user later wants stricter pinning, they can revisit.

## Code Examples

### pyproject.toml URL fix pattern

```toml
# Before (placeholder)
[project.urls]
Homepage = "https://github.com/yourusername/mace-gaussian-interface"
Repository = "https://github.com/yourusername/mace-gaussian-interface.git"
Issues = "https://github.com/yourusername/mace-gaussian-interface/issues"

# After (real)
[project.urls]
Homepage = "https://github.com/mott91/mace_gaussian_interface"
Repository = "https://github.com/mott91/mace_gaussian_interface.git"
Issues = "https://github.com/mott91/mace_gaussian_interface/issues"
```

Source: `git remote get-url origin` (HIGH confidence — directly from the repo's own config)

### CI ruff pin fix pattern (Option A)

```yaml
# Before
- name: Install ruff
  run: pip install "ruff==0.15.1"

# After (Option A — align with pyproject.toml floor)
- name: Install ruff
  run: pip install "ruff>=0.9.0"
```

### quickstart.md entry point fix pattern

```markdown
# Before
python cli.py run water.xyz --skip-dft-baseline

# After
mace-gaussian run water.xyz --skip-dft-baseline
```

Applied to all 4 occurrences in quickstart.md.

### README.md entry point fix pattern

```markdown
# Before
python cli.py run molecule.xyz

# After
mace-gaussian run molecule.xyz
```

Applied to all cli.py invocations; `python run_analysis.py water` stays unchanged.

## Open Questions

1. **CLAUDE.md entry point update**
   - What we know: CLAUDE.md has 5 `python cli.py` occurrences; it is developer-facing
   - What's unclear: Whether the user wants dev docs to show `mace-gaussian` or raw invocation
   - Recommendation: Update to `mace-gaussian` for consistency; note that `python cli.py` still works from project root

2. **ruff CI pin strategy**
   - What we know: CI pins 0.15.1, pyproject.toml has >=0.9.0; both work
   - What's unclear: Whether the user prefers strict pinning (Option B) or floor-match (Option A)
   - Recommendation: Default to Option A (floor-match); note the Phase 9-02 rationale for pinning was reproducibility, which is already served by uv.lock for local dev

3. **Author email in pyproject.toml**
   - What we know: `email = "your.email@domain.com"` is also a placeholder
   - What's unclear: Whether the user wants this fixed in Phase 12 or left for PyPI publication
   - Recommendation: Include in 12-01 if the user wants full distribution readiness; otherwise defer. The audit did not flag this explicitly so do not add it unless user requests.

## Sources

### Primary (HIGH confidence)
- Direct file read: `/home/mot/mace_gaussian/pyproject.toml` — current URL values confirmed
- Direct file read: `/home/mot/mace_gaussian/.github/workflows/ci.yml` — ruff pin confirmed as `==0.15.1`
- Direct file read: `/home/mot/mace_gaussian/docs/quickstart.md` — all 4 `python cli.py` occurrences confirmed
- Direct file read: `/home/mot/mace_gaussian/README.md` — all `python cli.py` occurrences confirmed
- `git remote get-url origin` → `git@github.com:mott91/mace_gaussian_interface.git` — real repo coordinates confirmed
- Direct file read: `/home/mot/mace_gaussian/CLAUDE.md` — 5 `python cli.py` occurrences confirmed
- Direct file read: `.planning/v1.0-MILESTONE-AUDIT.md` — all four gap descriptions confirmed

### Secondary (MEDIUM confidence)
- Phase 9-02 decision log (STATE.md line 139): "ruff==0.15.1 pinned to match uv.lock for CI reproducibility" — explains origin of pin

### Tertiary (LOW confidence)
- None. All findings are sourced from file reads or git commands.

## Metadata

**Confidence breakdown:**
- Current state inventory: HIGH — all values read directly from files
- Fix targets: HIGH — confirmed from git remote and file grep
- Ruff pin strategy: MEDIUM — two valid options exist, tradeoffs documented
- CLAUDE.md update: MEDIUM — judgment call on developer-vs-user docs

**Research date:** 2026-02-27
**Valid until:** 60 days — this is configuration/metadata, changes only if files are edited
