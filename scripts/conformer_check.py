#!/usr/bin/env python3
"""Conformer check over every run: ML final geometry vs the B3LYP final geometry.

Writes ``validation_results/conformer_check.csv`` and lists the runs that ended up in a
different conformer (see ``mace_gaussian/analysis/conformer_check.py``).

Usage::

    python scripts/conformer_check.py [molecule ...]
"""

from __future__ import annotations

import csv
import sys
from dataclasses import asdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from mace_gaussian.analysis.conformer_check import check_run  # noqa: E402

BASE = REPO / "comparison_results"
OUT = REPO / "validation_results" / "conformer_check.csv"


def dft_fchk(mol_dir: Path) -> Path | None:
    for d in sorted(mol_dir.glob("b3lyp*")):
        for name in ("gaussian_dft.fchk", f"{mol_dir.name}_freq_anharm.fchk"):
            if (d / name).exists():
                return d / name
    return None


def main(argv: list[str]) -> int:
    mols = argv or sorted(p.name for p in BASE.iterdir() if p.is_dir())
    rows = []
    for mol in mols:
        ref = dft_fchk(BASE / mol)
        if ref is None:
            continue
        for run in sorted((BASE / mol).glob("mace_*/gaussian_freq.fchk")):
            try:
                r = check_run(ref, run)
            except Exception as exc:
                print(f"  skip {mol}/{run.parent.name}: {exc}")
                continue
            rows.append({"molecule": mol, "run": run.parent.name, **asdict(r)})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if not r["same_conformer"]]
    print(
        f"{len(rows)} runs checked, {len(bad)} in a different conformer -> {OUT.relative_to(REPO)}"
    )
    for r in bad:
        print(
            f"  {r['molecule']:18s} {r['run']:24s} torsion off by {r['max_torsion_dev']:5.0f} deg"
            f" (atoms {r['worst_torsion']}), RMSD heavy {r['rmsd_heavy']:.2f} A"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
