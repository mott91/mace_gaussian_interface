# Pitfalls Research

**Research Date:** 2026-02-16
**Domain:** Refactoring scientific Python code built with AI assistance

## Pitfall 1: Refactoring Without Tests (CRITICAL)

**What goes wrong:** You refactor a module, everything looks correct, but you've silently changed behavior. Gaussian output parsing returns slightly different values. Mode matching produces different pairings. You don't notice until weeks later when results don't match your thesis data.

**Why it's especially dangerous here:**
- No comprehensive test suite exists
- Parser logic uses regex state machines that are fragile to refactoring
- Numerical code can have subtle floating-point behavior changes
- AI-generated code may have non-obvious coupling between modules

**Warning signs:**
- "It looks the same" but no test proves it
- Results change after refactoring but "only slightly"
- Tests pass but cover only happy paths

**Prevention:**
1. Before ANY refactoring, add characterization tests: run water and CH4 through the pipeline, commit the exact outputs as reference data
2. Write parser tests with exact expected outputs before touching parser code
3. Use `numpy.testing.assert_array_almost_equal` for numerical comparisons
4. Run full pipeline on water after every structural change

**Phase:** Must be addressed in Phase 1 (testing infrastructure), before any refactoring begins.

## Pitfall 2: Breaking the Working Pipeline

**What goes wrong:** Aggressive refactoring breaks the ZMQ communication, Gaussian can't find the helper script, or import paths change and the CLI silently fails. The tool was working; now it's not.

**Why it's especially dangerous here:**
- ZMQ + Gaussian subprocess is a complex integration that's hard to debug
- The external helper script path is embedded in generated .gjf files
- Module monkey-patching cleanup is order-dependent
- You can't easily test the full pipeline in CI

**Warning signs:**
- Gaussian hangs indefinitely (ZMQ socket not found)
- "Module not found" errors after moving files
- Helper script path is wrong in generated .gjf
- `.ipc_file` left behind from crashed runs

**Prevention:**
1. Keep `gm_main.py` functional throughout refactoring — extract from it, don't rewrite it
2. Each extraction: move code → update imports → verify pipeline still runs on water
3. Don't change the helper script mechanism until everything else is stable
4. Maintain backward compatibility for entry points (old scripts still work)

**Phase:** Applies to all refactoring phases. Each phase should end with a verified pipeline run.

## Pitfall 3: Over-Engineering the Architecture

**What goes wrong:** Refactoring turns into a full rewrite. Simple file moves become architecture redesigns. Three levels of abstraction are added where one function call worked fine. The project takes months instead of weeks and still isn't done.

**Why it's especially dangerous here:**
- Author is not a programmer and may not push back on complexity
- AI assistants tend to suggest "proper" patterns that add unnecessary abstraction
- Scientific code has different norms than enterprise code — simpler is better
- Thesis has a deadline

**Warning signs:**
- Abstract base classes with one implementation
- Generic factories for two options
- Configuration systems for things that never change
- "Future-proofing" for scenarios that won't happen during the thesis

**Prevention:**
1. Rule: only refactor what's actually broken or fragile (see CONCERNS.md)
2. Rule: if the existing code works and isn't fragile, move it but don't redesign it
3. Rule: three concrete duplications before abstracting
4. Ask "does this make the thesis easier?" for every change

**Phase:** Ongoing discipline. Review at each phase boundary.

## Pitfall 4: AI-Generated Code Inconsistencies

**What goes wrong:** Different AI sessions produced code with different conventions, error handling patterns, naming styles, and assumptions. Refactoring reveals that modules have incompatible assumptions about data formats, units, or error conditions.

**Why it's especially dangerous here:**
- Codebase was built across many AI sessions with varying context
- No style guide was enforced consistently
- Different parts may assume different unit conventions (Bohr vs Angstrom)
- Error handling ranges from detailed try/except to bare-metal code

**Warning signs:**
- Same concept named differently in different files
- Inconsistent unit handling (some functions expect Angstrom, others Bohr)
- Mixed error handling patterns (some raise, some return None, some log and continue)
- Duplicate utility functions in different modules

**Prevention:**
1. Establish conventions BEFORE refactoring (naming, error handling, units)
2. Add unit conversion functions with clear names (bohr_to_angstrom, not generic converters)
3. Standardize error handling: always raise with context, never return None for errors
4. Use ruff to enforce import order and basic style consistency
5. Document conventions in CONVENTIONS.md (already exists in codebase map)

**Phase:** Should be addressed in early phases alongside testing.

## Pitfall 5: Losing Reproducibility

**What goes wrong:** After refactoring, old results can't be reproduced. The tool produces slightly different outputs due to import order changes, floating-point accumulation differences, or changed default parameters.

**Why it's especially dangerous here:**
- Thesis results must be reproducible
- ML model inference can be sensitive to import order (GPU initialization)
- Gaussian external interface communicates via text files with fixed precision
- Mode matching depends on numerical thresholds

**Warning signs:**
- Results from before and after refactoring don't match exactly
- "Close enough" differences that accumulate
- Different results on different machines after refactoring

**Prevention:**
1. Commit reference outputs (water, CH4) BEFORE refactoring begins
2. After refactoring, verify outputs match within acceptable tolerance
3. Document any intentional changes to output format
4. Pin all dependency versions (uv.lock already helps)
5. Record Python version, MACE model version, and Gaussian version in output metadata

**Phase:** Phase 1 (testing) should establish the reference outputs. Every subsequent phase should verify against them.

## Pitfall 6: Ignoring the Acetic Acid Bug

**What goes wrong:** The known acetic acid parsing bug (commit a4384c4) gets forgotten during refactoring. The parser is cleaned up but the underlying regex issue persists, manifesting later with different molecules.

**Why it's especially dangerous here:**
- Bug is documented but root cause isn't fully understood
- Regex-based parsing is fragile — fixing one molecule can break another
- Anharmonic output format varies across DFT methods

**Warning signs:**
- Parser changes that aren't tested against acetic acid case
- New molecules fail in ways that look like the same bug
- Overtone/combination band counts don't match expected values

**Prevention:**
1. Add acetic acid .log file to test fixtures specifically
2. Write a failing test that reproduces the bug
3. Fix the bug as part of parser testing phase
4. Test parser against multiple Gaussian output formats (different methods, basis sets)

**Phase:** Should be addressed in the parser testing phase (early).

## Pitfall 7: Custom MACE Package Dependency Hell

**What goes wrong:** Refactoring changes how MACE packages are imported or initialized. The custom local packages (mace_torch, mace_dipole_core) have implicit dependencies on import order, CUDA initialization state, or specific module paths that break when the project structure changes.

**Why it's especially dangerous here:**
- Custom packages aren't on PyPI — can't simply reinstall
- MACE import behavior depends on which modules are loaded first
- The monkey-patching exists BECAUSE of complex import dependencies
- CUDA device placement is sensitive to import order

**Warning signs:**
- "CUDA error: device not found" after restructuring
- "Module has no attribute" errors from MACE
- Different behavior when running from different directories
- Tests pass individually but fail when run together (import order)

**Prevention:**
1. Document exact MACE import sequence before changing anything
2. Test MACE loading in isolation before and after changes
3. If replacing monkey-patching, start with lazy imports (less disruptive)
4. Always test with `pytest -x` (stop on first failure) to catch import issues early
5. Keep a rollback plan — if MACE loading breaks, revert and try a different approach

**Phase:** Address when replacing monkey-patching. This is the highest-risk refactoring task.

## Summary: Refactoring Risk Map

| Pitfall | Severity | Likelihood | Phase to Address |
|---------|----------|------------|-----------------|
| No tests before refactoring | CRITICAL | HIGH | Phase 1 |
| Breaking working pipeline | HIGH | MEDIUM | All phases |
| Over-engineering | MEDIUM | HIGH | All phases (discipline) |
| AI code inconsistencies | MEDIUM | HIGH | Phase 1-2 |
| Losing reproducibility | HIGH | MEDIUM | Phase 1 (establish baseline) |
| Acetic acid bug | MEDIUM | HIGH | Parser testing phase |
| MACE dependency hell | HIGH | MEDIUM | Monkey-patching replacement phase |

---
*Pitfalls research: 2026-02-16*
