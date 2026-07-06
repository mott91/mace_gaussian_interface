#!/usr/bin/env python3
"""Phase 24 spike: run Psience VPT2 from a Gaussian freq=anharm log/fchk pair.

Proof-of-concept that Psience can replace Gaussian's VPT2 engine: the
normal-mode reduced force field (cubic/quartic tables) is read from the
Gaussian log, remapped from Gaussian's spectroscopic mode numbering to
ascending-frequency order, and fed to Psience via `potential_terms`.
The fchk supplies geometry, masses and modes (for the kinetic/Coriolis
side). Result reproduces Gaussian's own VPT2 fundamentals to <0.01 cm-1
on water (DFT and ML surfaces alike).

NOTE: Psience CANNOT read our G16 fchk "Cartesian 3rd/4th derivatives"
directly — it produces unphysical fundamentals (stretches shift up).
Validated against the McCoy group's own G09-era MP2 test fchk, which
works. Suspected G16 mode-block ordering/normalization difference; the
log-table route below sidesteps it entirely.

Runs in the throwaway Psience venv, NOT mace4ir_v2:
    <psience_venv>/bin/python scripts/spike_psience_vpt2.py \
        comparison_results/water/mace_omol_mace_ml
"""

import itertools
import re
import sys
from pathlib import Path

import numpy as np


def parse_fundamentals(log_text: str) -> list[tuple[int, float, float]]:
    """(mode, harmonic, anharmonic) rows from the Fundamental Bands table."""
    block = re.search(
        r"Anharmonic Infrared Spectroscopy.*?Fundamental Bands.*?\n(.*?)\n\s*\n",
        log_text,
        re.S,
    )
    rows = []
    for line in block.group(1).splitlines():
        m = re.match(r"\s*(\d+)\(1\)\s+(-?[\d\.]+)\s+(-?[\d\.]+)", line)
        if m:
            rows.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return rows


def parse_force_field(log_text: str):
    """Reduced (cm-1) cubic and quartic normal-mode constants from the log."""
    cubic, quartic = [], []
    cub_block = re.search(r"CUBIC FORCE CONSTANTS IN NORMAL MODES(.*?)Num\. of 3rd", log_text, re.S)
    for line in cub_block.group(1).splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?[\d\.]+)", line)
        if m:
            cubic.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4))))
    quart_block = re.search(
        r"QUARTIC FORCE CONSTANTS IN NORMAL MODES(.*?)Num\. of 4th", log_text, re.S
    )
    for line in quart_block.group(1).splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?[\d\.]+)", line)
        if m:
            i, j, k, length = (int(m.group(x)) for x in range(1, 5))
            quartic.append((i, j, k, length, float(m.group(5))))
    return cubic, quartic


def main() -> None:
    calc_dir = Path(sys.argv[1])
    log_path = calc_dir / "gaussian_freq.log"
    if not log_path.exists():
        log_path = next(calc_dir.glob("*anharm*.log"))
    fchk_path = log_path.with_suffix(".fchk")
    if not fchk_path.exists():
        # fchk basename may differ from the log's (e.g. DFT baseline dirs);
        # any fchk carrying the anharmonic derivative block will do
        fchk_path = next(
            p
            for p in sorted(calc_dir.glob("*.fchk"))
            if "Cartesian 3rd/4th derivatives" in p.read_text(errors="ignore")
        )

    log_text = log_path.read_text(errors="ignore")
    fundamentals = parse_fundamentals(log_text)
    cubic, quartic = parse_force_field(log_text)

    # Gaussian table numbering (spectroscopic) -> ascending-frequency index
    modes = sorted(fundamentals, key=lambda r: r[1])
    g2a = {gmode: a for a, (gmode, _, _) in enumerate(modes)}
    n = len(modes)
    freqs = np.array([h for _, h, _ in modes])

    from McUtils.Data import UnitsData
    from Psience.VPT2 import VPTRunner

    w2h = UnitsData.convert("Wavenumbers", "Hartrees")
    V2 = np.diag(freqs) * w2h
    V3 = np.zeros((n, n, n))
    V4 = np.zeros((n, n, n, n))
    for i, j, k, fi in cubic:
        for idx in set(itertools.permutations((g2a[i], g2a[j], g2a[k]))):
            V3[idx] = fi * w2h
    for i, j, k, length, fi in quartic:
        for idx in set(itertools.permutations((g2a[i], g2a[j], g2a[k], g2a[length]))):
            V4[idx] = fi * w2h

    VPTRunner.run_simple(
        str(fchk_path), 1, potential_terms=[V2, V3, V4], zero_element_warning=False
    )

    print("\n=== Gaussian VPT2 reference (same log) ===")
    print(f"{'mode':>4} {'harmonic':>12} {'anharmonic':>12}")
    for gmode, harm, anharm in fundamentals:
        print(f"{gmode:>4} {harm:>12.4f} {anharm:>12.4f}")


if __name__ == "__main__":
    main()
