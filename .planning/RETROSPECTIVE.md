# Project Retrospective

*A living document updated after each milestone. Lessons feed forward into future planning.*

---

## Milestone: v1.0 — Refactoring & Distribution

**Shipped:** 2026-02-28
**Phases:** 12 | **Plans:** 32 | **Timeline:** 2025-09-08 → 2026-02-28 (~173 days)
**Commits:** 199 | **GSD execution:** ~108 min active

### What Was Built
- 131-test pytest suite with fixtures, GPU/Gaussian markers, and regression baselines for water + CH4
- Proper `mace_gaussian/` Python package (7,648 LOC, 29 modules) with `mace-gaussian` CLI entry point
- Safe MACE model loading via `mace_loader.py` — pickle_module remapping replaces sys.modules monkey-patching
- `gaussian/` subpackage: LINGER=0 ZMQ race-condition fix, SIGKILL timeout, context-manager ZMQ server
- GitHub Actions CI: ruff ≥0.9.0 + pytest (unit) on every push
- User documentation: quickstart, worked water example, thesis methods section (thesis-ready prose)

### What Worked
- **Testing infrastructure first (Phase 01):** Had fixtures and passing tests before touching source code — made all 11 subsequent refactoring phases safe to execute without fear of silent regressions
- **Phase extraction order:** utilities → calculators → I/O/ZMQ → workflow orchestrator → package wiring; each phase had clean dependency boundaries with no surprises
- **Lazy import pattern (established Phase 04):** `_check_availability` wrapping and function-level lazy imports for heavy deps (espaloma, DGL) let CI tests run without GPU, and prevented module-level side effects in tests
- **sys.modules pre-mocking in tests:** Reliable technique for testing calculator code in isolation; no DGL wheel required in CI
- **Quality GSD profile:** Thorough planning documents caught several wiring decisions (LINGER=0 placement, SIGKILL vs SIGTERM, relative-vs-absolute imports at entry boundary) before they became runtime bugs

### What Was Inefficient
- **Gap-closure phases 11 and 12:** Integration wiring gaps — error types not wired to public API, env vars not resolved before validation, calculation_parameters not passed through — were only discovered during Phase 10 testing. Should have had a dedicated integration-check phase scoped after Phase 08, not after Phase 10 (documentation).
- **CI setup underestimated (Phase 09 Plan 01, 90 min):** GitHub Actions debugging is iterative and hard to time-box. Future CI phases should be planned with explicit "debug budget" time rather than treated as a routine 3-minute plan.
- **ARCHITECTURE.md and DEVELOPMENT.md skipped during refactoring:** Both docs still reference pre-refactor module names. Writing docs *during* the relevant phase would have been cheaper than retrofitting.
- **Requirements tracked only in STATE.md decisions:** Having 38 requirements buried in a decisions log made the milestone audit harder than necessary. Formal REQUIREMENTS.md from the start would have made gap-closure clearer.

### Patterns Established
- Pre-mock heavy deps at `sys.modules` level before importing test subjects (espaloma, DGL, xtb)
- Lazy imports inside `_check_availability()` or function bodies for heavy deps — never at module top-level in production code
- Relative imports everywhere in the package *except* entry boundaries (`cli.py` uses absolute `mace_gaussian.*`)
- `LINGER=0` applied *after* `socket()` and *before* `bind()` — prevents `socket.close()` from blocking forever on Gaussian crash
- `SIGKILL` (not `SIGTERM`) for subprocesses that ignore signals (Gaussian)
- `warnings.warn(..., stacklevel=2)` so warning origin points to caller, not utility internals

### Key Lessons
1. **Start every refactor with a characterization + testing phase.** Understanding what exists before touching anything is the highest-ROI phase in any refactoring milestone.
2. **Phase extraction order mirrors dependency graph.** Plan phases bottom-up: shared utilities first, then domain modules, then orchestration, then package wiring, then CI/docs.
3. **Budget for gap-closure explicitly.** Two gap-closure phases (11, 12) were needed anyway — if planned from the start (e.g., "Phase 11: Integration check"), they wouldn't feel like overrun.
4. **CI is its own phase.** Not a 5-minute plan. Give it a dedicated slot with realistic scope.
5. **Update architecture docs during the phase that changes the architecture**, not in a separate documentation phase at the end.

### Cost Observations
- Profile: quality (GSD), yolo mode, `commit_docs: true`
- 199 commits over 173 calendar days
- ~108 min active GSD execution (avg 3.4 min/plan; Phase 09 P01 was a 90-min outlier)
- Notable: yolo mode with quality profile was the right trade-off — skipped friction, thorough planning

---

## Cross-Milestone Trends

### Process Evolution

| Milestone | Phases | Plans | Key Change |
|-----------|--------|-------|------------|
| v1.0 | 12 | 32 | First GSD milestone; established full workflow from characterization → package → CI → docs |

### Cumulative Quality

| Milestone | Tests | Coverage | LOC (src) | LOC (tests) |
|-----------|-------|----------|-----------|-------------|
| v1.0 | 131 passing, 1 xfail | ~57% | 7,648 | 2,147 |

### Top Lessons (Verified Across Milestones)

1. Characterization before refactoring — understand before changing
2. Testing infrastructure phase pays compound interest across all subsequent phases
