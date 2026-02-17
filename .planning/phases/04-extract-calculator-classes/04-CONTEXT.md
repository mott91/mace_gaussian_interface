# Phase 4: Extract Calculator Classes - Context

**Gathered:** 2026-02-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Move DipoleCalculator base class, its implementations (espaloma, mace_ml, xtb), and the calculator factory from gm_main.py into a dedicated `calculators/` package. The MACE module monkey-patching problem (sys.modules manipulation, cleanup_mace_modules) is Phase 5 scope — this phase extracts the class hierarchy only, keeping existing import patterns intact.

</domain>

<decisions>
## Implementation Decisions

### Claude's Discretion

User trusts Claude's expertise on all implementation decisions for this phase. The following guidance captures sensible defaults for downstream agents:

**Calculator boundary:**
- Extract: DipoleCalculatorBase (ABC), EspalomaCalculator, MaceMLCalculator, XtbCalculator, DipoleCalculatorFactory
- The `calculator()` function in gm_main.py that handles MACE standard model loading stays in gm_main.py — it's deeply tied to the monkey-patching (Phase 5)
- Any helper that directly calls `cleanup_mace_modules()` stays in gm_main.py for now
- The calculators themselves import from mace/espaloma/xtb — those imports move with them, but the sys.modules manipulation wrapper stays behind
- If a calculator's `__init__` or `calculate()` method calls cleanup, that call stays inline for now and Phase 5 will clean it up

**Factory design:**
- Keep the factory pattern simple — a dict mapping or function that returns calculator instances
- No need to over-engineer (no plugin registry, no entry points) — there are exactly 3 calculator types
- Factory should accept calculator name as string and return instantiated calculator
- Factory lives in `calculators/__init__.py` or `calculators/factory.py`

**Testing strategy:**
- Calculator tests should mock heavy dependencies (MACE models, espaloma, xtb)
- Test the factory pattern: correct class returned for each name, error on unknown name
- Test calculator interface compliance: all implementations have the required methods
- Integration tests (actually loading models) stay marked with @pytest.mark.gpu
- Use `unittest.mock.patch` for model loading, not actual model files

**Package layout:**
- `calculators/__init__.py` — re-exports factory and base class
- `calculators/base.py` — DipoleCalculatorBase ABC
- `calculators/espaloma.py` — EspalomaCalculator
- `calculators/mace_ml.py` — MaceMLCalculator
- `calculators/xtb.py` — XtbCalculator (if it exists as a class)
- `calculators/factory.py` — DipoleCalculatorFactory
- One class per file keeps things clean and matches the pattern established in Phase 3

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches. User explicitly deferred all decisions to Claude's judgment. Priority is robustness: the pipeline must keep working with all 4 model combinations after extraction.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 04-extract-calculator-classes*
*Context gathered: 2026-02-17*
