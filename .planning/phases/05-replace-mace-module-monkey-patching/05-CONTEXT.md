# Phase 5: Replace MACE Module Monkey-Patching - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace the `sys.modules` manipulation in `mace_calculators.py` with a safe loading pattern. After this phase, standard MACE and dipole MACE coexist without runtime module swapping, `cleanup_mace_modules()` calls are removed, and the `finally` block cleanup pattern is no longer needed. The CLAUDE.md warning about monkey-patching cleanup becomes obsolete.

Current state: `mace_calculators.py` swaps `sys.modules["mace.modules.models"]` on every dipole calculation, reloads modules, and cleans up in `finally` blocks. This is fragile (global state mutation, silent corruption if cleanup is missed) but functionally correct.

</domain>

<decisions>
## Implementation Decisions

### Claude's Discretion

User trusts Claude's expertise on all implementation decisions for this phase. The following guidance captures the design intent and constraints for downstream agents:

**Loading strategy:**
- The core problem: standard MACE (`mace.modules.models`) and dipole MACE (`mace_dipole_core.modules.models`) collide in `sys.modules`. Current code hot-swaps the cache entry.
- Research should evaluate: (1) `importlib.util` isolated loading (load dipole variant without registering in `sys.modules`), (2) namespace/import tricks to avoid the collision entirely, (3) subprocess isolation as a fallback option.
- Priority order: simplest approach that works with CUDA first. `importlib.util` is the likely winner but needs verification that CUDA contexts survive isolated loading.
- Subprocess isolation is a last resort — adds latency and GPU memory serialization complexity.
- The `fake_module_from_real()` function (shallow module copy) should be eliminated, not just replaced.

**Migration boundary:**
- `mace_calculators.py` should be absorbed into the `calculators/` package (likely as `calculators/mace_loader.py` or similar) since Phase 4 already extracted the calculator classes there.
- `MACEDipoleCalculator` wrapper class moves into calculators/ — it's the natural home after Phase 4's extraction.
- `load_standard_mace_calculator()` and `load_dipole_mace_calculator()` get replaced with the new safe loading mechanism.
- `cleanup_mace_modules()` is deleted entirely — the whole point is eliminating the need for cleanup.
- All `finally` blocks that call `cleanup_mace_modules()` are removed.
- The `calculator()` function in `gm_main.py` (lines 444-467) that loads standard MACE for energy calculations should also use the new safe loader.
- Update CLAUDE.md to remove the monkey-patching warning.

**Fallback behavior:**
- No fallback to monkey-patching. Clean break — if the new loading doesn't work, it's a bug to fix, not a case to handle.
- Clear error messages if either MACE variant fails to load (which model, what went wrong, is it installed?).
- CUDA device placement must be explicit — both models on the same GPU, with clear logging of device assignment.

**Testing approach:**
- Mock the actual model loading (torch, MACE internals) to test the loading mechanism in isolation.
- Verify: standard MACE and dipole MACE can be loaded in the same process without `sys.modules` conflicts.
- Verify: no `cleanup_mace_modules()` calls remain anywhere in the codebase.
- Verify: repeated dipole calculations don't accumulate state or leak memory.
- Integration tests (actually loading models) stay marked with `@pytest.mark.gpu`.
- The existing pipeline must still work with all 4 model combinations on water.

</decisions>

<specifics>
## Specific Ideas

No specific requirements -- open to standard approaches. The user understands the problem deeply (discussed the monkey-patching mechanism, `sys.modules` cache, the swap-compute-cleanup cycle) and trusts Claude to choose the right technical solution. Priority is reliability over cleverness.

</specifics>

<deferred>
## Deferred Ideas

- Batching dipole derivative calculations (computing all 6N displacements in one GPU call) -- potential performance optimization but separate from fixing the loading mechanism. Could be a future phase or post-thesis work.

</deferred>

---

*Phase: 05-replace-mace-module-monkey-patching*
*Context gathered: 2026-02-18*
