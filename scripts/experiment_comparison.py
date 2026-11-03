#!/usr/bin/env python3
"""Fundamentals against gas-phase experiment, for every molecule in band_origins.json.

Experimental bands are put on B3LYP modes with ``master_table.assign_band_origins``
(explicit ``modes`` where the file gives them, otherwise nearest B3LYP-VPT2 frequency).
Every ML value comes from the eigenvector assignment in ``bands.csv``. A band that stands
for several modes (degenerate) is compared with the mean of those modes.

Usage::

    python scripts/experiment_comparison.py BANDS_CSV OUT_DIR
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mace_gaussian.analysis.master_table import (
    assign_band_origins,
    load_band_origins,
)

MODELS = ["mace_anicc", "mace_off", "mace_omol", "mace_polar"]
SCALE_B3LYP = 0.9614  # common B3LYP/6-31G(d) harmonic scale factor (Scott & Radom 1996)


def main(bands_csv: str, out_dir: str) -> None:
    with Path(bands_csv).open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["band_type"] == "fundamental"]
    by_mol = defaultdict(list)
    for r in rows:
        if r["dipole_model"] == "mace_ml":
            by_mol[r["molecule"]].append(r)

    out_rows = []
    for mol in sorted(by_mol):
        entry = load_band_origins(mol)
        if entry is None:
            continue
        dft = {int(r["i"]): r for r in by_mol[mol]}
        ref = {i: float(r["dft_anharmonic"]) for i, r in dft.items() if r["dft_anharmonic"]}
        taken = assign_band_origins(entry["bands"], ref)
        groups = defaultdict(list)  # band label -> DFT modes
        info = {}
        for k, b in taken.items():
            groups[b["label"]].append(k)
            info[b["label"]] = b
        for label, ks in groups.items():
            b = info[label]
            row = {"molecule": mol, "label": label, "description": b.get("description", ""),
                   "exp": b["freq_cm"], "b3lyp_modes": "+".join(map(str, sorted(ks))),
                   "assignment": b["assignment"], "note": b.get("note", "")}  # fmt: skip
            row["B3LYP_vpt2"] = np.mean([float(dft[k]["dft_anharmonic"]) for k in ks])
            row["B3LYP_harm_scaled"] = SCALE_B3LYP * np.mean(
                [float(dft[k]["dft_harmonic"]) for k in ks]
            )
            for m in MODELS:
                vals = [
                    float(r["ml_anharmonic"])
                    for r in rows
                    if r["molecule"] == mol
                    and r["energy_model"] == m
                    and r["dipole_model"] == "mace_ml"
                    and int(r["i"]) in ks
                    and r["ml_anharmonic"]
                ]
                row[m] = np.mean(vals) if len(vals) == len(ks) else ""
            out_rows.append(row)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "experiment_fundamentals.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        for r in out_rows:
            w.writerow({k: (f"{v:.1f}" if isinstance(v, float) else v) for k, v in r.items()})

    print(
        f"{'method':18s} {'n':>3s} {'MAE':>6s} {'median':>7s} {'bias':>6s} {'max':>6s}  molecules"
    )
    for m in ["B3LYP_vpt2", "B3LYP_harm_scaled", *MODELS]:
        err = [(float(r[m]) - float(r["exp"]), r["molecule"]) for r in out_rows if r[m] != ""]
        e = np.array([x for x, _ in err])
        mols = sorted({mm for _, mm in err})
        print(f"{m:18s} {len(e):3d} {np.abs(e).mean():6.1f} {np.median(np.abs(e)):7.1f} "
              f"{e.mean():+6.1f} {np.abs(e).max():6.1f}  {len(mols)}")  # fmt: skip


if __name__ == "__main__":
    main(*sys.argv[1:3])
