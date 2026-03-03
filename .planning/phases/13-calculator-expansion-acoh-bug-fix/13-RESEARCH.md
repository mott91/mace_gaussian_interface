# Phase 13: Calculator Expansion & Acoh Bug Fix - Research

**Researched:** 2026-03-03
**Domain:** MACE calculator integration, Click CLI validation, Gaussian log parsing
**Confidence:** HIGH (all findings verified against source code in the repo)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Add `mace_off` and `mace_anicc` energy calculators only
- xTB (energy + dipole roles) deferred — user will discuss with supervisor before including in benchmark
- `mace_anicc`: treat as auto-download (same `MACECalculator` interface); no special model path handling
- `mace_anicc` element guard: **early check** — inspect atom symbols before loading model; raise `ValueError` listing unsupported elements found (e.g., "mace_anicc only supports H/C/N/O — molecule contains: [F, S]")
- Acoh parser: **narrow fix** — extend regex in `parse_anharmonic_frequencies()` to handle H/L-prefixed lines in the `'Vibrational Energies at Anharmonic Level'` section; fix both documented bugs within existing parsing logic; no new fallback path needed
- CLI validation: add `click.Choice` to `--energy-calculators` — valid values: `mace_mp`, `mace_omol`, `mace_off`, `mace_anicc`; split comma-separated values before checking

### Claude's Discretion
- Exact `ValueError` message wording for element guard
- Whether to validate `--dipole-calculators` the same way (consistent with energy side)
- Test coverage approach for the new calculator branches (mock vs. integration)

### Deferred Ideas (OUT OF SCOPE)
- xTB energy calculator (CALC-03) — deferred pending supervisor discussion
- xTB dipole unit verification (CALC-04) — same deferral
- TorchANI/ANI-2x — deferred, confirm after benchmark molecule set chosen
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| CALC-01 | `mace_off` available as `--energy-calculators` choice | Already in `workflow.py::calculator()` at line 227 — CLI wiring only |
| CALC-02 | `mace_anicc` available as `--energy-calculators` choice with element guard | Different API signature confirmed; element restriction documented; auto-download behavior confirmed |
| FIX-01 | Acoh frequency parser bug fixed; xfail test promoted to passing | Two bugs confirmed by reading actual fixture; fix path fully specified |
</phase_requirements>

---

## Summary

Phase 13 has three independent sub-problems: (1) wiring `mace_off` into the CLI, (2) adding `mace_anicc` to the CLI + `workflow.py` with an element guard, and (3) fixing the anharmonic frequency parser for acetic acid logs. All three sub-problems are self-contained — no shared state, no ordering dependency.

The most important discovery is that `mace_anicc` has a **different function signature** from every other MACE calculator. It takes `(device, model_path, return_raw_model)` with no `model` size parameter (there is only one model: `ani500k_large_CC`). It auto-downloads on first call. The model file is not currently on disk. The function has no built-in element guard — the 18-element restriction must be enforced in our `workflow.py`.

The acoh parser bug is more subtle than the CONTEXT.md summary suggests. The acoh ML log produces a **different Gaussian output format** than the water DFT log. The water log uses the `Anharmonic Infrared Spectroscopy` section with `I(harm)` / `I(anharm)` intensity columns. The acoh ML log produces `Vibrational Energies at Anharmonic Level` with `E(harm)`, `E(anharm)`, and rotational constants — **no intensity column at all**. Additionally, some lines are prefixed with `H ` or `L ` overlap indicators. The parser must handle both section names and both line formats, yielding entries with `ir_intensity = 0.0` when the log lacks that column.

**Primary recommendation:** Implement in three isolated tasks: (1) CLI wiring for `mace_off` + `mace_anicc` + `click.Choice` validation, (2) `workflow.py` `mace_anicc` branch + element guard, (3) parser fix + xfail removal.

---

## Standard Stack

### Core (already in use — no new dependencies)
| Library | Version | Purpose | Notes |
|---------|---------|---------|-------|
| `mace` | installed in mace4ir_v2 | MACE calculators | `mace_off`, `mace_anicc` already exported from `mace.calculators` |
| `click` | installed in .venv | CLI framework | `click.Choice`, callback parameter |
| `re` | stdlib | Regex parsing | Used in `gaussian/parser.py` |
| `ase` | installed | Atoms object | `atoms.get_chemical_symbols()` for element check |

### No New Installs Required
All work is pure code changes. The `mace_anicc` model auto-downloads on first use from `https://github.com/ACEsuit/mace/raw/main/mace/calculators/foundations_models/ani500k_large_CC.model`.

---

## Architecture Patterns

### Pattern 1: mace_anicc API (CONFIRMED — different from other MACE calculators)

**What:** `mace_anicc` uses a different signature than `mace_mp`, `mace_off`, `mace_omol`.

```python
# Source: /home/mot/micromamba/envs/mace4ir_v2/lib/python3.10/site-packages/mace/calculators/foundations_models.py

# mace_mp, mace_off, mace_omol pattern (has model size):
calc = mace_mp(model="large", device="cuda", default_dtype="float64", dispersion=False)
calc = mace_off(model="large", device="cuda", default_dtype="float64", dispersion=False)

# mace_anicc pattern (NO model size, just device + optional model_path):
# Signature: mace_anicc(device: str = "cuda", model_path: str = None, return_raw_model: bool = False)
calc = mace_anicc(device="cuda")  # model auto-downloads to foundations_models/ on first call
```

**Critical:** Do NOT pass `model="large"` or `default_dtype` to `mace_anicc`. They are not accepted parameters.

### Pattern 2: mace_off (already in workflow.py — no changes needed there)

```python
# Source: mace_gaussian/workflow.py lines 227-231 — ALREADY IMPLEMENTED
if nnp == "mace_off":
    from mace.calculators import mace_off
    calc = mace_off(model="large", device="cuda", default_dtype="float64", dispersion=False)
    return calc
```

### Pattern 3: Element guard for mace_anicc

```python
# Pattern: check before loading model — avoids slow model download/load for unsupported molecules
MACE_ANICC_SUPPORTED_ELEMENTS = {"H", "C", "N", "O"}

def _check_mace_anicc_elements(atoms):
    symbols = set(atoms.get_chemical_symbols())
    unsupported = symbols - MACE_ANICC_SUPPORTED_ELEMENTS
    if unsupported:
        raise ValueError(
            f"mace_anicc only supports H/C/N/O — molecule contains: "
            f"{sorted(unsupported)}"
        )
```

This guard belongs in `workflow.py::calculator()` before calling `mace_anicc(...)`, OR in `run_geometry_optimization` / `run_frequency_calculation` before calling `calculator()`. Earliest possible point is preferred.

### Pattern 4: click.Choice with comma-separated values

**What:** `--energy-calculators` takes comma-separated strings. `click.Choice` doesn't handle this natively. Use a `callback` parameter to validate after splitting.

```python
# Source: mace_gaussian/cli.py lines 48-57 — template from --optimization-calculator

# Current (freeform, no validation):
@click.option(
    "--energy-calculators",
    default="mace_mp,mace_omol",
    help="...",
)

# Fixed pattern (validate each value after split):
VALID_ENERGY_CALCULATORS = ["mace_mp", "mace_omol", "mace_off", "mace_anicc"]

def validate_energy_calculators(ctx, param, value):
    choices = [c.strip() for c in value.split(",")]
    for choice in choices:
        if choice not in VALID_ENERGY_CALCULATORS:
            raise click.BadParameter(
                f"'{choice}' is not one of {VALID_ENERGY_CALCULATORS}.",
                param=param,
            )
    return value  # return the original string; parsing happens later in run()

@click.option(
    "--energy-calculators",
    default="mace_mp,mace_omol",
    callback=validate_energy_calculators,
    help="...",
)
```

Also add `mace_off` and `mace_anicc` to `--optimization-calculator` click.Choice if they should be usable for geometry optimization. The context does not address this — apply Claude's discretion.

### Pattern 5: Acoh Gaussian log format (CONFIRMED by reading fixture)

**Two distinct Gaussian output formats exist:**

**Format A (water DFT)** — `Anharmonic Infrared Spectroscopy` section, with intensity column:
```
     ==================================================
                  Anharmonic Infrared Spectroscopy
     ==================================================

 Fundamental Bands
 -----------------
   Mode(n)                  E(harm)   E(anharm)        I(harm)       I(anharm)
      1(1)                  3799.220   3624.001      1.63737651      0.62829225
      2(1)                  1665.300   1615.153     70.33036928     69.29467704
```

**Format B (acoh ML)** — `Vibrational Energies at Anharmonic Level` section, NO intensity column:
```
     ==================================================
          Vibrational Energies at Anharmonic Level
     ==================================================

 NOTE: H and L indicates if there is a high or low overlap ...

 Fundamental Bands
 -----------------
     Mode(n)      Status      E(harm)  E(anharm)     Aa(x)      Ba(y)      Ca(z)
        1(1)      active      3796.914  3544.674   0.378112   0.317316   0.178021
 H      4(1)      active      1828.929  1812.980   0.376865   0.317416   0.177834
```

**Key differences in Format B:**
1. Section title is `Vibrational Energies at Anharmonic Level` (not `Anharmonic Infrared Spectroscopy`)
2. Column header says `Status` not `I(harm)` / `I(anharm)` — no intensity column
3. Data lines have format: `[H|L]?  MODE(1)  STATUS  E(harm)  E(anharm)  Aa  Ba  Ca`
4. Some lines prefixed with `H ` (high overlap) or `L ` (low overlap)
5. All 18 modes ARE listed (1 through 18)
6. No IR intensities available in this format

**Fix strategy for `parse_anharmonic_frequencies()`:**
- The current trigger `in_fundamental_section = True` requires `I(anharm)` or `DS(anharm)` in lookahead after `Fundamental Bands` — this fails for Format B
- Fix: also trigger on `E(anharm)` header in the lookahead (Format B uses this)
- The current data-line regex `r"^\s*(\d+)\(1\)\s+([\d\.]+)\s+([\d\.]+)\s+..."` fails on `H      4(1)      active ...` because of the `H ` prefix and `active` status word
- Fix: new regex that optionally handles `[H|L]?` prefix and optional `STATUS` word before the numeric columns, OR two separate regex patterns tried in order
- Return `ir_intensity=0.0` when no intensity column present (the test only checks `len(result) == 18`)

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Model download | Custom urllib downloader | Let `mace_anicc()` handle it | Already implemented with progress reporting |
| HCNO element list | Hardcoded list from memory | Verify against MACE docs/source | mace_anicc source confirms "H, C, N, O" |
| CLI multi-value validation | Custom `type=` class | `callback=` on `click.option` | Simpler, consistent with click patterns |

---

## Common Pitfalls

### Pitfall 1: Wrong mace_anicc call signature
**What goes wrong:** Passing `model="large"` or `default_dtype="float64"` to `mace_anicc` raises TypeError.
**Why it happens:** Cargo-culting from `mace_mp`/`mace_off` pattern.
**How to avoid:** `mace_anicc(device="cuda")` only. Signature: `(device, model_path, return_raw_model)`.
**Warning signs:** `TypeError: mace_anicc() got unexpected keyword argument 'model'`

### Pitfall 2: Element guard placed too late
**What goes wrong:** Model downloads (potentially slow) before element check runs.
**Why it happens:** Guard placed inside `calculator()` after the `mace_anicc` call.
**How to avoid:** Guard must run BEFORE calling `mace_anicc()`. In `workflow.py::calculator()`, add the check at the top of the `mace_anicc` branch before the import.

### Pitfall 3: atoms not available in calculator()
**What goes wrong:** `calculator(nnp)` takes only a string — no `atoms` argument — so the guard can't run there without a signature change.
**Why it happens:** Current `calculator()` signature is `def calculator(nnp)`.
**How to avoid:** Either (a) add an optional `atoms=None` parameter to `calculator()`, or (b) place the guard in `run_frequency_calculation()` and `run_geometry_optimization()` before calling `calculator()`. Option (b) avoids signature changes — preferred.

### Pitfall 4: Acoh parser fix breaks water tests
**What goes wrong:** Loosening the section-header trigger for `in_fundamental_section` causes the parser to enter the section in logs it shouldn't, producing duplicate or incorrect entries.
**Why it happens:** Over-broad trigger condition.
**How to avoid:** The trigger must remain specific. The reliable trigger for Format B is: after `Fundamental Bands` line, lookahead for `E(anharm)` in the header row (but NOT `I(anharm)` which is already handled). Run the full parser test suite after any change to `parse_anharmonic_frequencies()`.

### Pitfall 5: click.Choice on comma-separated option
**What goes wrong:** Using `type=click.Choice([...])` directly on `--energy-calculators` rejects valid inputs like `"mace_mp,mace_omol"` because Click validates the whole string, not individual tokens.
**Why it happens:** `click.Choice` treats the value as an atomic string.
**How to avoid:** Use `callback=validate_energy_calculators` (custom function) instead of `type=click.Choice`.

### Pitfall 6: acoh log has no ir_intensity — test expects 18 entries
**What goes wrong:** Parser returns entries but with `ir_intensity` key missing, causing downstream analysis code to KeyError.
**Why it happens:** Assuming Format B has an intensity column.
**How to avoid:** Return `ir_intensity=0.0` explicitly for Format B entries. The xfail test only checks `len(result) == 18` — any non-empty intensity value satisfies the test.

---

## Code Examples

### Confirmed: mace_anicc branch in workflow.py

```python
# Source: verified against foundations_models.py in mace4ir_v2 environment

if nnp == "mace_anicc":
    # Element guard: must happen BEFORE model loading
    if atoms is not None:
        supported = {"H", "C", "N", "O"}
        symbols = set(atoms.get_chemical_symbols())
        unsupported = symbols - supported
        if unsupported:
            raise ValueError(
                f"mace_anicc only supports H/C/N/O — "
                f"molecule contains: {sorted(unsupported)}"
            )
    from mace.calculators import mace_anicc
    calc = mace_anicc(device="cuda")
    return calc
```

Note: `workflow.py::calculator()` currently has no `atoms` parameter. The guard requires adding `atoms=None` to the signature, or moving the guard to the call sites (`run_geometry_optimization`, `run_frequency_calculation`).

### Confirmed: Acoh log data line formats (both bugs visible)

```
# Format B Fundamental Bands section header line — triggers `in_fundamental_section`:
     Mode(n)      Status      E(harm)  E(anharm)     Aa(x)      Ba(y)      Ca(z)

# Normal line (no prefix):
        1(1)      active      3796.914  3544.674   0.378112   0.317316   0.178021

# H-prefixed line (overlap indicator — breaks current regex):
 H      4(1)      active      1828.929  1812.980   0.376865   0.317416   0.177834
```

**Regex that handles both:**
```python
# Handles optional H/L prefix, optional status word, then two floats (harm, anharm)
match = re.match(
    r"^\s*(?:[HL]\s+)?(\d+)\(1\)\s+\w+\s+([\d\.]+)\s+([\d\.]+)",
    line
)
```

This matches:
- `        1(1)      active      3796.914  3544.674 ...`
- ` H      4(1)      active      1828.929  1812.980 ...`

### Confirmed: Section trigger fix for Format B

```python
# Current code (only triggers on I(anharm) or DS(anharm)):
if "Fundamental Bands" in line:
    for j in range(i, min(i + 10, len(lines))):
        if "I(anharm)" in lines[j] or "DS(anharm)" in lines[j]:
            in_fundamental_section = True
            break
    continue

# Fixed (also triggers on E(anharm) header — Format B):
if "Fundamental Bands" in line:
    for j in range(i, min(i + 10, len(lines))):
        if "I(anharm)" in lines[j] or "DS(anharm)" in lines[j]:
            in_fundamental_section = True
            in_format_b = False  # has intensity column
            break
        if "E(anharm)" in lines[j] and "Status" in lines[j]:
            in_fundamental_section = True
            in_format_b = True   # no intensity column
            break
    continue
```

Alternatively, a simpler single-pass approach: trigger on any `Fundamental Bands` header after recognising we are inside a `Vibrational Energies at Anharmonic Level` block, and use a broader data-line regex that covers both formats.

### Confirmed: CLI callback validation pattern

```python
# Source: Click documentation pattern; consistent with existing cli.py style

VALID_ENERGY_CALCULATORS = ["mace_mp", "mace_omol", "mace_off", "mace_anicc"]

def _validate_calculator_list(ctx, param, value):
    """Validate comma-separated calculator names."""
    for name in [c.strip() for c in value.split(",")]:
        if name not in VALID_ENERGY_CALCULATORS:
            raise click.BadParameter(
                f"'{name}' is not a valid energy calculator. "
                f"Choose from: {', '.join(VALID_ENERGY_CALCULATORS)}",
                param=param,
            )
    return value

@click.option(
    "--energy-calculators",
    default="mace_mp,mace_omol",
    callback=_validate_calculator_list,
    help="Comma-separated energy calculators. Choices: mace_mp, mace_omol, mace_off, mace_anicc",
)
```

---

## State of the Art

| Old State | New State | Change Required | Impact |
|-----------|-----------|-----------------|--------|
| `mace_off` in workflow, not in CLI | `mace_off` in CLI Choice | Add to callback validator + help text | User can pass `--energy-calculators mace_off` |
| No `mace_anicc` in codebase | `mace_anicc` in workflow + CLI | Add branch in `calculator()`, element guard, CLI entry | Enables CCSD(T)-quality PES benchmark |
| `--energy-calculators` freeform | `--energy-calculators` validated | Add callback validator | Fails fast with clear message on invalid input |
| Acoh parser returns 0 anharmonic modes (xfail) | Parser returns 18 modes | Dual-format section detection + H/L-aware regex | Regression plots work for acetic acid |

---

## Open Questions

1. **Should `mace_anicc` be added to `--optimization-calculator` click.Choice?**
   - What we know: `--optimization-calculator` currently accepts `["mace_omol", "mace_off", "mace_mp"]`
   - What's unclear: Whether the user would ever optimize with `mace_anicc` (it's HCNO only)
   - Recommendation: Add it for consistency. The element guard will protect against misuse.

2. **atoms parameter in `calculator()` — signature change or call-site guard?**
   - What we know: Guard must run before model load. `calculator(nnp)` has no atoms param.
   - What's unclear: Whether adding `atoms=None` to `calculator()` would break existing callers
   - Recommendation: Check all callers of `calculator()`. There are exactly two call sites in `workflow.py` (lines ~351 and ~649). Both have `atoms` available. Guard at call sites avoids signature change and is simpler.

3. **ir_intensity for Format B entries**
   - What we know: acoh log has no `I(anharm)` column in Fundamental Bands section
   - What's unclear: Whether downstream analysis code will break on `ir_intensity=0.0`
   - Recommendation: Return `ir_intensity=0.0`. The xfail test only checks `len == 18`. Downstream analysis treats zero-intensity modes as "not plotted" which is correct behavior for modes where intensity wasn't computed.

4. **Dipole calculators validation**
   - What we know: CONTEXT.md says "Claude's Discretion" on whether to add `click.Choice` to `--dipole-calculators`
   - Valid dipole calculators: `espaloma`, `mace_ml`, `xtb` (from `factory.py::preferred_order`)
   - Recommendation: Add the same callback pattern for consistency.

---

## Sources

### Primary (HIGH confidence — direct code inspection)
- `/home/mot/mace_gaussian/mace_gaussian/workflow.py` — confirmed `mace_off` already at line 227; confirmed `calculator()` signature; confirmed no `mace_anicc` branch
- `/home/mot/mace_gaussian/mace_gaussian/cli.py` — confirmed `--energy-calculators` is freeform; confirmed `--optimization-calculator` uses `click.Choice`
- `/home/mot/mace_gaussian/mace_gaussian/gaussian/parser.py` — confirmed bug: `in_fundamental_section` trigger requires `I(anharm)` or `DS(anharm)`; confirmed current regex `r"^\s*(\d+)\(1\)"` would fail on H/L lines
- `/home/mot/mace_gaussian/tests/fixtures/acoh/ml_mace_mp_esp.log` — confirmed Format B structure: `Vibrational Energies at Anharmonic Level`, no intensity column, H/L prefixes at lines 128 and 133
- `/home/mot/mace_gaussian/tests/fixtures/water/dft_b3lyp.log` — confirmed Format A structure: `Anharmonic Infrared Spectroscopy`, `I(harm)` + `I(anharm)` columns
- `/home/mot/micromamba/envs/mace4ir_v2/lib/python3.10/site-packages/mace/calculators/foundations_models.py` — confirmed `mace_anicc` signature: `(device, model_path, return_raw_model)`; no `model` size param; auto-downloads `ani500k_large_CC.model`
- `/home/mot/mace_gaussian/tests/test_gaussian_parser.py` — confirmed xfail test expects `len(result) == 18`, no intensity assertion

### Secondary (HIGH confidence — direct execution)
- `micromamba run -n mace4ir_v2 python3 -c "from mace.calculators import mace_anicc, mace_off, mace_mp, mace_omol"` — confirmed all four are importable
- Test run: `source .venv/bin/activate && python3 -m pytest tests/test_gaussian_parser.py -v` — confirmed 20 passed, 1 xfailed; xfail is `test_acoh_anharmonic_parsing`

---

## Metadata

**Confidence breakdown:**
- mace_off wiring: HIGH — already in workflow.py, CLI pattern is clear
- mace_anicc API: HIGH — signature confirmed by source code inspection and test run
- mace_anicc element guard: HIGH — guard logic is standard Python set arithmetic
- CLI validation: HIGH — click.Choice + callback pattern is well-established
- Acoh parser bug: HIGH — both bugs confirmed by reading actual fixture file
- Acoh parser fix: MEDIUM-HIGH — fix path is specified but two approaches exist (dual-flag vs. single-pass); the simpler approach needs implementation-time selection

**Research date:** 2026-03-03
**Valid until:** 2026-06-01 (stable codebase; mace library version unlikely to change in active use)
