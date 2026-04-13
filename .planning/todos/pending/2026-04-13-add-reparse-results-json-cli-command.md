---
created: 2026-04-13T16:15:00.000Z
title: Add CLI command to re-parse results.json from existing Gaussian logs
area: cli
files:
  - mace_gaussian/cli.py
  - mace_gaussian/gaussian/parser.py
---

## Problem

When the Gaussian parser is improved (e.g. dual-format anharmonic parsing), existing
results.json files become stale. Currently the only way to update them is to re-run
the full calculation or manually script the re-parse.

Octane DFT had all anharmonic fundamental intensities = 0.0 because results.json was
generated before the parser fix — the log file had the correct data all along.

## Proposed Solution

Add a CLI command like `mace-gaussian reparse <molecule>` that:
1. Finds all results.json under `comparison_results/<molecule>/`
2. Re-parses the corresponding `.log` files with the current parser
3. Updates only the frequency/intensity data in results.json, preserving metadata
4. Optionally accepts `--dry-run` to show what would change without writing

Should handle both DFT and ML log formats.
