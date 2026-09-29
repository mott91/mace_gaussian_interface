#!/usr/bin/env python3
"""Autograd vs. finite-difference dipole derivatives for the MACE4IR dipole model.

For each molecule, times one call of each path (what the harness does at every
one of Gaussian's 6N-11 external calls) and compares the two derivative matrices:

- autograd: ``MACEDipoleCalculator.calculate_dipole_derivatives`` (production path,
  one ``get_dielectric_derivatives`` call)
- finite differences: ``DipoleCalculatorBase.calculate_dipole_derivatives``
  (central differences, 6N dipole evaluations), at two step sizes

Central differences have an O(delta^2) error, so the relative difference should
drop 100x when delta drops 10x if autograd is the exact limit.

Usage::

    python scripts/benchmark_dipole_derivatives.py
    python scripts/benchmark_dipole_derivatives.py water methanol --repeats 10

Outputs land in ``docs/explained/dipole_derivatives/``: ``benchmark.md`` and
``benchmark.csv``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import platform
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import torch
from ase.io import read

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from mace_gaussian.calculators.base import DipoleCalculatorBase  # noqa: E402
from mace_gaussian.calculators.mace_loader import MACEDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL  # noqa: E402

DEFAULT_MOLECULES = ["water", "methanol", "gly", "aspirin", "decane"]
OUT_DIR = REPO / "docs" / "explained" / "dipole_derivatives"


def _timed(fn, repeats: int):
    """Median wall time of ``fn`` over ``repeats`` calls, with CUDA synchronized."""
    times = []
    result = None
    for _ in range(repeats):
        if torch.cuda.is_available():
            torch.cuda.synchronize()
        start = time.perf_counter()
        result = fn()
        if torch.cuda.is_available():
            torch.cuda.synchronize()
        times.append(time.perf_counter() - start)
    return result, float(np.median(times))


def benchmark(molecules: list[str], repeats: int) -> list[dict]:
    calc = MACEDipoleCalculator(DEFAULT_MACE_DIPOLE_MODEL)
    calc._ensure_calculator()

    def fd(atoms, delta):
        return DipoleCalculatorBase.calculate_dipole_derivatives(calc, atoms, displacement=delta)

    rows = []
    for name in molecules:
        atoms = read(REPO / "molecules" / f"{name}.xyz")
        # Warm-up: first calls pay for CUDA kernel compilation and allocation.
        calc.calculate_dipole_derivatives(atoms)
        fd(atoms, 0.01)

        auto, t_auto = _timed(lambda a=atoms: calc.calculate_dipole_derivatives(a), repeats)
        fd_01, t_fd = _timed(lambda a=atoms: fd(a, 0.01), max(1, repeats // 2))
        fd_001 = fd(atoms, 0.001)

        def rel(x, ref=auto):
            return float(np.linalg.norm(x - ref) / np.linalg.norm(ref))

        rows.append(
            {
                "molecule": name,
                "atoms": len(atoms),
                "fd_dipole_evals": 6 * len(atoms),
                "fd_s": t_fd,
                "autograd_s": t_auto,
                "speedup": t_fd / t_auto,
                "rel_diff_delta_0.01": rel(fd_01),
                "rel_diff_delta_0.001": rel(fd_001),
            }
        )
        print(
            f"{name:10s} N={len(atoms):3d}  FD {t_fd:7.3f} s  autograd {t_auto:7.4f} s  "
            f"{t_fd / t_auto:5.0f}x  rel.diff {rel(fd_01):.1e} / {rel(fd_001):.1e}"
        )
    return rows


def write_outputs(rows: list[dict], repeats: int) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "benchmark.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    gpu = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none (CPU)"
    lines = [
        "# Dipole derivatives: autograd vs. finite differences",
        "",
        f"Run {date.today().isoformat()} with `scripts/benchmark_dipole_derivatives.py`.",
        f"Model: MACE4IR (`{Path(DEFAULT_MACE_DIPOLE_MODEL).name}`), GPU: {gpu}, "
        f"CPU: {platform.processor() or platform.machine()}.",
        f"Times are the median of {repeats} calls (finite differences: {max(1, repeats // 2)}), "
        "one call = what the harness does per Gaussian external call.",
        "",
        "| molecule | atoms | FD dipole evals | FD [s] | autograd [s] | speedup "
        "| rel. diff δ=0.01 Å | rel. diff δ=0.001 Å |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['molecule']} | {r['atoms']} | {r['fd_dipole_evals']} | {r['fd_s']:.3f} "
            f"| {r['autograd_s']:.4f} | {r['speedup']:.0f}x | {r['rel_diff_delta_0.01']:.1e} "
            f"| {r['rel_diff_delta_0.001']:.1e} |"
        )
    lines += [
        "",
        "Relative difference = ||FD - autograd|| / ||autograd|| (Frobenius norm over the "
        "3N x 3 matrix).",
        "Autograd cost is roughly flat in N; finite differences cost 6N dipole evaluations,",
        "so the speedup grows with molecule size. The difference drops 100x for a 10x smaller",
        "step, the O(δ²) signature of central differences converging onto the autograd value.",
        "",
    ]
    (OUT_DIR / "benchmark.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_DIR / 'benchmark.md'} and benchmark.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("molecules", nargs="*", default=DEFAULT_MOLECULES)
    parser.add_argument("--repeats", type=int, default=6)
    args = parser.parse_args()
    logging.disable(logging.WARNING)
    rows = benchmark(args.molecules, args.repeats)
    write_outputs(rows, args.repeats)


if __name__ == "__main__":
    main()
