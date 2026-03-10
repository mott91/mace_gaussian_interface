---
status: testing
phase: 01-testing-infrastructure-characterization
source: [01-01-SUMMARY.md, 01-02-SUMMARY.md, 01-03-SUMMARY.md, 01-04-SUMMARY.md]
started: 2026-02-16T15:20:00Z
updated: 2026-02-16T15:20:00Z
---

## Current Test

number: 1
name: Run full test suite
expected: |
  Run `python -m pytest tests/ -v` from the project root (with venv activated).
  Should show ~45 tests: most PASSED, 1 xfail (acoh), 2 skipped (marker placeholders).
  No errors, no collection failures, no import errors.
awaiting: user response

## Tests

### 1. Run full test suite
expected: Run `python -m pytest tests/ -v`. Should show ~45 tests: most PASSED, 1 xfail (acoh), 2 skipped (marker placeholders). No errors or collection failures.
result: [pending]

### 2. Custom markers filter correctly
expected: Run `python -m pytest tests/ -m "not gpu and not gaussian" -v`. The GPU and Gaussian placeholder tests should be deselected. All remaining tests should pass.
result: [pending]

### 3. Coverage report runs
expected: Run `python -m pytest tests/ --cov --cov-report=term-missing -q`. Should show coverage percentages for core modules (gaussian_parser, fchk_parser, mode_matching). No errors from coverage plugin.
result: [pending]

### 4. Test fixtures are not gitignored
expected: Run `git check-ignore tests/fixtures/water/dft_b3lyp.log tests/fixtures/water/dft_b3lyp.fchk tests/fixtures/CH4_ase/dft_b3lyp.log`. Should return empty (no output) — meaning none of these fixture files are gitignored.
result: [pending]

### 5. uv.lock is tracked by git
expected: Run `git ls-files uv.lock`. Should print `uv.lock` (confirming it's tracked). Previously it was gitignored.
result: [pending]

### 6. Acoh xfail documents parser bug
expected: Run `python -m pytest tests/test_gaussian_parser.py -v -k acoh`. Should show the acoh test as XFAIL (expected failure), not FAILED. The test reason should mention the missing "Anharmonic Infrared Spectroscopy" section.
result: [pending]

## Summary

total: 6
passed: 0
issues: 0
pending: 6
skipped: 0

## Gaps

[none yet]
