#!/usr/bin/env python3
"""Campaign fundamentals against the in-house VCI references (CCSD(T)-F12 surfaces).

Pairs every VCI fundamental with a B3LYP normal mode by the *character* of the mode, not by
frequency order. The mapping below was set on 2026-11-03 from the B3LYP eigenvectors
(bond-stretch amplitudes, in-plane/out-of-plane share, symmetry under the molecule's mirror
plane or C2 axis; see thesis docs/handoff/followup-1-report.md, item B). Each ML model's
value for that B3LYP mode comes from the eigenvector assignment in ``bands.csv``.

Usage::

    python scripts/vci_comparison.py VCI_CSV BANDS_CSV OUT_CSV
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

# VCI label -> B3LYP harmonic mode index(es), 1-based ascending frequency (bands.csv "i").
# A tuple is a degenerate pair; its values are averaged.
# fmt: off
MODE_MAP: dict[str, dict[str, tuple[int, ...]]] = {
    "methanol": {
        # mirror plane C-O-H(O): a' = nu1-nu8, a'' = nu9-nu12 (Dinu 2024 numbering)
        "nu12": (1,), "nu8": (2,), "nu7": (3,), "nu11": (4,), "nu6": (5,), "nu5": (6,),
        "nu10": (7,), "nu4": (8,), "nu3": (9,), "nu9": (10,), "nu2": (11,), "nu1": (12,),
    },
    "formic_acid_trans": {
        "nu7": (1,), "nu9": (2,), "nu8": (3,), "nu6": (4,), "nu5": (5,), "nu4": (6,),
        "nu3": (7,), "nu2": (8,), "nu1": (9,),
    },
    "carbon_dioxide": {"nu2": (1, 2), "nu1": (3,), "nu3": (4,)},
    "carbonic_acid": {
        # C2 about C=O: a = nu1-nu6 (nu6 = a2, the OH oop "anti"), b = nu7-nu12.
        # nu10 (C=O in-plane bend) is mode 3 and nu12 (OH oop "sym") is mode 4: the reverse of
        # a plain frequency-order pairing (the comparison of 2026-10-07 already had it this way).
        "nu6": (1,), "nu5": (2,), "nu10": (3,), "nu12": (4,), "nu11": (5,), "nu4": (6,),
        "nu9": (7,), "nu3": (8,), "nu8": (9,), "nu2": (10,), "nu7": (11,), "nu1": (12,),
    },
}
# fmt: on
MODELS = ["mace_anicc", "mace_off", "mace_omol", "mace_polar"]
DIPOLE = "mace_ml"  # frequencies do not depend on the dipole model


def _f(x: str) -> float | None:
    return float(x) if x not in ("", None) else None


def main(vci_csv: str, bands_csv: str, out_csv: str) -> None:
    with Path(vci_csv).open() as fh:
        vci = list(csv.DictReader(fh))
    # (molecule, model or "B3LYP", B3LYP mode i) -> (harmonic, vpt2)
    val: dict[tuple, tuple] = {}
    with Path(bands_csv).open() as fh:
        bands = list(csv.DictReader(fh))
    for r in bands:
        if r["band_type"] != "fundamental" or r["dipole_model"] != DIPOLE:
            continue
        if r["molecule"] not in MODE_MAP:
            continue
        i = int(r["i"])
        val[(r["molecule"], r["energy_model"], i)] = (_f(r["ml_harmonic"]), _f(r["ml_anharmonic"]))
        val[(r["molecule"], "B3LYP", i)] = (_f(r["dft_harmonic"]), _f(r["dft_anharmonic"]))

    rows = []
    for v in vci:
        idx = MODE_MAP[v["molecule"]][v["label"]]
        row = {k: v[k] for k in ("molecule", "label", "description", "vci", "gas_exp", "harm_cc")}
        row["b3lyp_mode"] = "+".join(str(i) for i in idx)
        row["flag"] = v["resonance"]
        for m in ["B3LYP", *MODELS]:
            got = [val.get((v["molecule"], m, i)) for i in idx]
            if any(g is None for g in got):
                row[f"{m}_h"] = row[f"{m}_vpt2"] = ""
                continue
            row[f"{m}_h"] = f"{np.mean([g[0] for g in got]):.1f}"
            row[f"{m}_vpt2"] = f"{np.mean([g[1] for g in got]):.1f}"
        rows.append(row)

    with Path(out_csv).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def stats(ref_key: str, suffix: str, keep) -> dict[str, tuple]:
        out = {}
        for m in ["B3LYP", *MODELS]:
            err = [
                float(r[f"{m}_{suffix}"]) - float(r[ref_key])
                for r in rows
                if keep(r) and r[ref_key] and r[f"{m}_{suffix}"]
            ]
            if err:
                e = np.array(err)
                out[m] = (len(e), np.abs(e).mean(), np.median(np.abs(e)), e.mean(), np.abs(e).max())
        return out

    for title, ref, suf, keep in [
        ("VPT2 vs VCI, clean fundamentals (flagged rows left out)", "vci", "vpt2",
         lambda r: not r["flag"]),
        ("harmonic vs CCSD(T)-F12 harmonic (methanol torsion left out)", "harm_cc", "h",
         lambda r: not (r["molecule"] == "methanol" and r["label"] == "nu12")),
    ]:  # fmt: skip
        print(f"\n{title}")
        print(f"{'method':12s} {'n':>3s} {'MAE':>6s} {'median':>7s} {'bias':>6s} {'max':>6s}")
        for m, (n, mae, med, bias, mx) in stats(ref, suf, keep).items():
            print(f"{m:12s} {n:3d} {mae:6.1f} {med:7.1f} {bias:+6.1f} {mx:6.1f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
