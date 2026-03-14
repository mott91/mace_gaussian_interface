---
created: 2026-03-10T16:00:00.000Z
title: Clean up temp files permanently in base directory
area: general
files: []
---

## Problem

Temporary files accumulate in the project base directory over time (e.g., from Gaussian runs, ZMQ IPC sockets, checkpoint files, intermediate outputs). These clutter the working directory and can cause confusion or stale-data issues.

## Solution

1. Identify all temp/intermediate files that end up in the base directory (`.ipc_file`, `*.chk`, `*.rwf`, scratch logs, etc.).
2. Either move them to a dedicated temp/scratch directory or ensure they are cleaned up after each run.
3. Add cleanup logic to the workflow (e.g., post-run cleanup step) or a dedicated `mace-gaussian clean` CLI command.
4. Update `.gitignore` if needed to prevent any remaining temp files from being tracked.
