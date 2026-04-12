---
phase: 23-anharmonic-pipeline-report-overhaul
reviewed: 2026-04-12T00:00:00Z
depth: standard
files_reviewed: 15
files_reviewed_list:
  - mace_gaussian/analysis/_shared_css.py
  - mace_gaussian/analysis/plotly_builders.py
  - mace_gaussian/analysis/report_data.py
  - mace_gaussian/analysis/executive_summary.py
  - mace_gaussian/analysis/html_report_generator.py
  - mace_gaussian/analysis/analysis_workflow.py
  - mace_gaussian/analysis/batch_report.py
  - pyproject.toml
  - tests/conftest.py
  - tests/test_batch_report.py
  - tests/test_executive_summary.py
  - tests/test_html_report.py
  - tests/test_plotly_builders.py
  - tests/test_report_data.py
  - tests/test_shared_css.py
findings:
  critical: 1
  warning: 4
  info: 4
  total: 9
status: issues_found
---

# Phase 23: Code Review Report

**Reviewed:** 2026-04-12
**Depth:** standard
**Files Reviewed:** 15
**Status:** issues_found

## Summary

Phase 23 delivers the report overhaul: shared CSS, Plotly figure builders, structured data export, executive summary ranking, and an updated `HTMLReportGenerator`. The new code is well-structured and the test suite is thorough with RED-test discipline applied consistently.

One critical security issue exists: `batch_report.py` interpolates user-controlled strings (molecule names, combo names, hardware strings) directly into HTML without escaping -- the same T-23-01 XSS class that `html_report_generator.py` correctly mitigated with `_esc()`. The test `test_batch_report_escapes_molecule_name` is already written and red for this exact reason.

Four warnings cover a logic bug in the speedup guard, a dead private method, a reachable `ZeroDivisionError`, and a state mutation hazard in `generate_report`. Four info items cover duplication, a TODO comment, and minor type annotation gaps.

---

## Critical Issues

### CR-01: XSS -- batch_report.py interpolates unescaped user-controlled strings into HTML

**File:** `mace_gaussian/analysis/batch_report.py:671-677, 516-518, 665-667`

**Issue:** `_build_timing_html` interpolates `r['molecule']`, `r['combo']`, `ml_hw`, and `dft_hw` directly into `<td>` elements (lines 671, 672, 676, 677). `_generate_html` does the same for the per-molecule spectrum overlay heading `mol_name` (lines 516, 518). `dft_hw` is assembled from `dft_cpu`, `dft_node`, and `dft_cpus` without escaping (lines 665-667). All these values ultimately derive from filesystem directory names or JSON-loaded strings, which an attacker who can write result files (or who controls the directory name -- as demonstrated by the red test) can use to inject `<script>` tags. The test `test_batch_report_escapes_molecule_name` already red-flags this gap.

**Fix:**
```python
# At top of batch_report.py
import html as _html

def _esc(s: object) -> str:
    return _html.escape(str(s), quote=True)

# In _build_timing_html, replace:
f"<td>{r['molecule']}</td>"
f"<td>{r['combo']}</td>"
f"<td style='font-size:0.8em'>{ml_hw}</td>"
f"<td style='font-size:0.8em'>{dft_hw}</td>"

# With:
f"<td>{_esc(r['molecule'])}</td>"
f"<td>{_esc(r['combo'])}</td>"
f"<td style='font-size:0.8em'>{_esc(ml_hw)}</td>"
f"<td style='font-size:0.8em'>{_esc(dft_hw)}</td>"

# In _generate_html / spectrum loop, replace:
f"<h3>{mol_name}</h3>"
f'alt="Spectrum overlay for {mol_name}">'

# With:
f"<h3>{_esc(mol_name)}</h3>"
f'alt="Spectrum overlay for {_esc(mol_name)}">'

# In _build_leaderboard_html, replace:
f"<td>{row['combo']}</td>"
# With:
f"<td>{_esc(row['combo'])}</td>"
```

---

## Warnings

### WR-01: ZeroDivisionError when dft_freqs contains a zero frequency

**File:** `mace_gaussian/analysis/analysis_workflow.py:297`

**Issue:** `create_comparison_table` computes `"Percent_Error": 100 * (ml_freq - dft_freq) / dft_freq` using numpy arrays. If any DFT frequency is 0.0 (which can happen for imaginary/transitional modes that survive as 0), this produces `inf` or `nan` silently for numpy arrays but will propagate corrupt values into the DataFrame and JSON export without warning. The result is silent data corruption in downstream CSV outputs.

**Fix:**
```python
with np.errstate(divide="ignore", invalid="ignore"):
    pct_err = np.where(
        dft_freq != 0,
        100 * (ml_freq - dft_freq) / dft_freq,
        np.nan,
    )
df_data["Percent_Error"] = pct_err
```

### WR-02: Dead method -- `_create_timing_hardware_section` is defined but never called

**File:** `mace_gaussian/analysis/html_report_generator.py:517`

**Issue:** `_create_timing_hardware_section` is defined as a private method (lines 517-600) but is never invoked anywhere in `HTMLReportGenerator`. The timing information is instead rendered inline as a simple `timing-block` div in `_create_comparison_section` (lines 362-370). The dead method has different, more detailed logic (two-row table with hardware columns) that diverges from what actually appears in reports. This creates a maintenance hazard: the dead method may be incorrectly assumed to be active when debugging timing output.

**Fix:** Either delete the dead method and add a comment referencing the inline timing block, or wire it into `_create_comparison_section` to replace the simpler timing div.

### WR-03: Speedup divisor guard in batch_report is too narrow -- zero ml_runtime_s passes through

**File:** `mace_gaussian/analysis/batch_report.py:113-116`

**Issue:** The speedup calculation is:
```python
row["speedup"] = (
    dft_gauss_s / row["ml_runtime_s"]
    if row.get("ml_runtime_s")
    else 0
)
```
`row.get("ml_runtime_s")` is falsy for both `None` and `0.0`, so a zero ML runtime correctly yields speedup=0. However, `dft_gauss_s` is taken from `dft_timing.get("total_elapsed_s", 0)` while `row["ml_runtime_s"]` is `ml_data.get("runtime_s", 0)` (line 171) -- these are **different timing fields**. The speedup therefore compares Gaussian wall time (DFT) against Python pipeline runtime (ML), not like-for-like. This produces misleadingly large speedup numbers and was already noted as a concern in the timing section of `_build_timing_html`. The batch report's `Speedup` column in the leaderboard is then misleading.

**Fix:** Use `ml_gauss_s` (Gaussian elapsed) for the speedup denominator, with a guard for zero:
```python
row["speedup"] = (
    dft_gauss_s / ml_gauss_s
    if ml_gauss_s and ml_gauss_s > 0
    else 0
)
```
where `ml_gauss_s` is extracted from `ml_gauss.get("total_elapsed_s", 0)` (already available as `ml_gauss_s` at line 172).

### WR-04: `generate_report` mutates the caller's `analysis_results` dict in place

**File:** `mace_gaussian/analysis/html_report_generator.py:84-94`

**Issue:** `generate_report` calls `_enrich_with_broadened_spectra` and `_compute_and_attach_experimental_agreement`, which write `_ml_broadened`, `_dft_broadened`, `freq_grid`, `_exp_on_grid`, and `experimental_agreement` directly onto the caller-supplied dict and its nested `comparisons` dicts. It then writes `executive_summary` onto the top-level dict. The caller's `fake_analysis_results` fixture in the test suite already pre-populates `experimental_agreement` (0.92 and 0.85), which gets silently overwritten. This is a hidden side-effect that can cause subtle test ordering issues and makes the function unsafe to call twice on the same dict.

**Fix:** Either document the mutation contract explicitly in the docstring, or make a shallow copy of `analysis_results` and deep-copy `comparisons` before mutating:
```python
import copy
analysis_results = dict(analysis_results)
analysis_results["comparisons"] = copy.deepcopy(analysis_results["comparisons"])
```

---

## Info

### IN-01: `batch_report.py` defines a private `_build_css()` that duplicates `_shared_css.build_css()`

**File:** `mace_gaussian/analysis/batch_report.py:752`

**Issue:** `batch_report.py` has its own `_build_css()` function (lines 752-834) that is entirely separate from `_shared_css.build_css()`. The Phase 23 D-04 design goal is that both report generators share identical CSS. Currently they do not -- the batch report uses a different visual theme (dark nav bar, `#2c3e50` headers) that diverges from the shared module. The `_shared_css.py` module even has a comment section labeled "Batch report specific (from batch_report.py)" at line 304, suggesting the intent was to consolidate. This is not a bug today but will cause maintenance drift.

**Fix:** Import `build_css` from `_shared_css` in `batch_report.py` and delete the local `_build_css`, or explicitly document that batch reports intentionally use a different visual theme.

### IN-02: TODO comment marks unfinished overtone integration as permanently deferred

**File:** `mace_gaussian/analysis/html_report_generator.py:606-609`

**Issue:** The `_create_overtones_section` method has a TODO comment saying overtone data is not yet surfaced. The stub path is intentional and test-covered, but the TODO comment leaves no tracking reference (no issue number, no phase reference) for when this will be resolved.

**Fix:** Replace the comment with a reference to the tracking item:
```python
# TODO(phase-24): wire in structured overtone records once analysis_workflow
# emits comparison["overtones"] -- tracked in 2026-03-26 pending todo.
```

### IN-03: `build_spectrum_figure` and `build_combined_spectrum_figure` return `object` instead of `plotly.graph_objects.Figure`

**File:** `mace_gaussian/analysis/plotly_builders.py:29, 97, 162`

**Issue:** All three builder functions are typed as returning `object` to avoid importing plotly at module level. This is a valid optimization (the docstring explains "Plotly is imported function-locally to keep CLI startup fast"), but it means callers have no type information. The `_fig_to_div` method in `html_report_generator.py` already accepts `Any`, so there is no runtime impact. However, it removes `ty check` coverage from all call sites.

**Fix:** Use `TYPE_CHECKING` to provide the type annotation at zero import cost:
```python
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import plotly.graph_objects as go

def build_spectrum_figure(...) -> go.Figure: ...
```

### IN-04: `test_aggregate_results_with_real_data` and `test_generate_batch_report_creates_html` depend on the `comparison_results/` directory existing at test time

**File:** `tests/test_batch_report.py:16-37`

**Issue:** Two tests call `aggregate_results("comparison_results")` and `generate_batch_report(results_dir="comparison_results", ...)` using a relative path that assumes the working directory contains a populated `comparison_results/` tree at test time. These tests will pass on the developer's machine but silently fail (or raise `ValueError`) in CI environments where that directory is absent. The test for `generate_batch_report` would raise `ValueError: No comparison results found` without a meaningful failure message if the directory is missing.

**Fix:** Add a `pytest.importorskip`-style guard or a `skipif` condition:
```python
import pytest
_RESULTS_DIR = Path("comparison_results")

@pytest.mark.skipif(
    not _RESULTS_DIR.is_dir(),
    reason="comparison_results/ not present -- integration test only",
)
def test_aggregate_results_with_real_data(): ...
```

---

_Reviewed: 2026-04-12_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
