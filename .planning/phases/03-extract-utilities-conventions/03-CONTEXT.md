# Phase 3: Extract Utilities & Conventions - Context

**Gathered:** 2026-02-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Extract pure functions from gm_main.py into a `utils/` package and establish documented code conventions for the project. This phase creates the foundation that later extraction phases (4-7) will follow. No new functionality — only reorganization of existing code and convention documentation.

</domain>

<decisions>
## Implementation Decisions

### Claude's Discretion

User trusts Claude's expertise on all implementation decisions for this phase. The following guidance captures sensible defaults for downstream agents:

**What gets extracted:**
- Unit conversion functions (e.g., hartree↔eV, bohr↔angstrom) → `utils/units.py`
- Pure helper functions from gm_main.py that don't depend on Gaussian/MACE runtime state
- ResultsManager class → `utils/results.py` (already somewhat standalone in `results_manager.py`)
- Boundary: if a function requires MACE models, Gaussian subprocess, or ZMQ — it stays in place for now (Phases 4-6 handle those)

**Validation consolidation:**
- Phase 2 created `validation.py` and `exceptions.py` at the top level
- Move these into `utils/validation.py` and `utils/exceptions.py` respectively
- Update all imports across the codebase to use the new locations
- Keep the same API — no behavioral changes, just relocation

**Import strategy:**
- `utils/__init__.py` should re-export key items for convenience (e.g., `from utils.units import hartree_to_ev`)
- Internal code uses direct imports: `from utils.units import ...`
- No deep nesting — `utils/` is flat (units.py, validation.py, exceptions.py, results.py)

**Convention scope:**
- Document conventions in a lightweight `docs/CONVENTIONS.md` or a section in existing `docs/DEVELOPMENT.md`
- Cover: naming conventions (snake_case functions, PascalCase classes), error handling pattern (use exceptions.py hierarchy), unit conventions (internal units: cm⁻¹ for frequencies, Å for distances), import ordering (stdlib → third-party → local)
- Keep it pragmatic — conventions should reflect what the codebase already does well, not aspirational rules
- Ruff configuration in pyproject.toml enforces what can be automated

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches. User explicitly deferred all decisions to Claude's judgment.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 03-extract-utilities-conventions*
*Context gathered: 2026-02-17*
