---
created: 2026-03-10T14:14:33.442Z
title: Build batch calculation pipeline for cluster
area: tooling
files: []
---

## Problem

Currently, calculations are run one molecule at a time via `mace-gaussian run`. For a benchmark campaign (exotic molecules, alkane series) a batch pipeline is needed that can submit and manage many jobs on a HPC cluster (likely SLURM-based). This includes job templating, status tracking, and result collection.

## Solution

1. Design a batch job schema: list of molecules, per-job config overrides, output dir layout.
2. Write a SLURM job template for `mace-gaussian run` with appropriate resource requests.
3. Build a submit script that generates and submits jobs for all molecules in a list.
4. Add a status/collection script to aggregate completed results into the analysis pipeline.
5. Test on the cluster with a small molecule set first.
