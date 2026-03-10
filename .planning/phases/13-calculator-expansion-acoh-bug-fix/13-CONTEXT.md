# Phase 13: Calculator Expansion & Acoh Bug Fix - Context

**Gathered:** 2026-03-03
**Status:** Ready for planning

<domain>
## Phase Boundary

Wire two new energy calculators (`mace_off`, `mace_anicc`) into the CLI and `workflow.py`, fix the acetic acid anharmonic frequency parser so the xfail test passes, and add `click.Choice` validation to `--energy-calculators`.

`mace_off` is already implemented in `workflow.py::calculator()` — CLI wiring only.
`mace_anicc` needs a `workflow.py` branch + CLI wiring + early element guard.
xTB (CALC-03, CALC-04) is explicitly out of scope for this phase — deferred pending supervisor discussion.

</domain>

<decisions>
## Implementation Decisions

### Calculator scope
- Add `mace_off` and `mace_anicc` only
- xTB (energy + dipole roles) deferred — user will discuss with supervisor before including in benchmark

### mace_anicc model loading
- Treat as auto-download like `mace_mp` / `mace_off` (same `MACECalculator` interface)
- No special model path handling — if it downloads on first use, great; if not, we discover it during testing

### mace_anicc element guard
- **Early check** — inspect atom symbols before loading the model
- Raise a clear `ValueError` with message listing the unsupported elements found (e.g., "mace_anicc only supports H/C/N/O — molecule contains: [F, S]")
- Do not let it reach the calculator and crash cryptically

### Acoh parser fix
- **Narrow fix** — extend the regex in `parse_anharmonic_frequencies()` to handle H/L-prefixed lines in the `'Vibrational Energies at Anharmonic Level'` section
- The two documented bugs: (1) section header mismatch when ML log lacks `'Anharmonic Infrared Spectroscopy'`, (2) lines prefixed with `H`/`L` overlap indicators that break `^\s*(\d+)\(1\)` match
- Fix both within the existing parsing logic (no new fallback path needed)
- Promote xfail test to passing (remove `@pytest.mark.xfail`)

### CLI validation
- Add `click.Choice` to `--energy-calculators` — valid values: `mace_mp`, `mace_omol`, `mace_off`, `mace_anicc`
- Currently freeform; `--optimization-calculator` already uses `click.Choice` — make consistent
- Note: `--energy-calculators` accepts comma-separated values, so validation requires splitting before checking — handle appropriately

### Claude's Discretion
- Exact `ValueError` message wording for element guard
- Whether to validate `--dipole-calculators` the same way (consistent with energy side)
- Test coverage approach for the new calculator branches (mock vs. integration)

</decisions>

<specifics>
## Specific Ideas

- Thesis framing: the core question is "can ML accelerate IR calculations a lot" — energy-vs-dipole quality is a supporting subplot, not the headline. mace_anicc adds value to that subplot (extends energy quality axis to CCSD(T)/CBS level for HCNO).
- mace_anicc is scientifically interesting because: if CCSD(T)-quality PES (mace_anicc) doesn't beat PBE-level (mace_mp) on IR accuracy, the dipole model is clearly the bottleneck — clean thesis result.

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `workflow.py::calculator()` — has `mace_mp`, `mace_off`, `mace_omol` branches; add `mace_anicc` branch following identical pattern: `from mace.calculators import mace_anicc; calc = mace_anicc(model="large", ...)`
- `mace_gaussian/utils/validation.py` — check for existing element validation utilities before writing new guard logic
- `calculators/base.py` — `DipoleCalculatorBase._check_availability()` pattern shows how to handle optional deps gracefully

### Established Patterns
- All MACE energy calculators use: `model="large"`, `device="cuda"`, `default_dtype="float64"`, `dispersion=False`
- CLI options use `click.Choice` for constrained choices (`--optimization-calculator` is the template)
- `--energy-calculators` is comma-separated string → split → iterate; validation must account for this

### Integration Points
- `mace_gaussian/cli.py` — add `mace_anicc` to `click.Choice` for `--energy-calculators`; currently `--optimization-calculator` has `["mace_omol", "mace_off", "mace_mp"]`
- `mace_gaussian/workflow.py::calculator()` — add `mace_anicc` branch (lines ~219–242)
- `mace_gaussian/gaussian/parser.py::parse_anharmonic_frequencies()` — the regex fix target (lines ~100–165)
- `tests/test_gaussian_parser.py` — remove `@pytest.mark.xfail` from `test_acoh_anharmonic_parsing` after fix

</code_context>

<deferred>
## Deferred Ideas

- **xTB energy calculator** (CALC-03) — deferred pending supervisor discussion about whether semiempirical baseline belongs in the benchmark and thesis narrative
- **xTB dipole unit verification** (CALC-04) — same deferral; unit bug (`xtb.py` comment says "e*Bohr" but no conversion applied; espaloma divides by BOHR_TO_ANGSTROM) still unresolved
- **TorchANI/ANI-2x** — already deferred in requirements; confirm after benchmark molecule set is chosen

</deferred>

---

*Phase: 13-calculator-expansion-acoh-bug-fix*
*Context gathered: 2026-03-03*
