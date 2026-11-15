"""Cubic and quartic constants k(i,i,i), k(i,i,i,i) of one mode from a Gaussian VPT2 log.

Gaussian prints, in the anharmonic section, QUADRATIC / CUBIC / QUARTIC FORCE CONSTANTS IN NORMAL
MODES with three columns: FI (reduced, cm-1), k (attoJ amu^-n/2 Angstrom^-n) and K (atomic units).
The anharmonic section numbers the modes by decreasing frequency; the mode is located by its
harmonic wavenumber in the QUADRATIC table, not by assuming the order.
"""

from __future__ import annotations

import re
from pathlib import Path


def _block(text: str, title: str) -> list[list[str]]:
    i = text.index(title)
    rows, started = [], False
    for ln in text[i:].splitlines()[1:]:
        f = ln.split()
        if f and all(re.fullmatch(r"-?\d+", x) for x in f[: len(f) - 3]) and len(f) >= 5 and not ln.strip().startswith(":"):
            started = True
            rows.append(f)
        elif started and (not f or not re.fullmatch(r"-?\d+", f[0])):
            break
    return rows


def constants(log: Path, freq: float) -> dict:
    t = log.read_text()
    quad = _block(t, "QUADRATIC FORCE CONSTANTS IN NORMAL MODES")
    n, best = None, 1e9
    for r in quad:
        if r[0] == r[1] and abs(float(r[2]) - freq) < best:
            n, best = r[0], abs(float(r[2]) - freq)
    assert best < 0.5, (log, freq, best)
    out = {"anharm_index": int(n), "freq": float(next(r for r in quad if r[0] == r[1] == n)[2])}
    for title, size, name in (("CUBIC FORCE CONSTANTS IN NORMAL MODES", 3, "cubic"),
                              ("QUARTIC FORCE CONSTANTS IN NORMAL MODES", 4, "quartic")):  # fmt: skip
        row = next((r for r in _block(t, title) if r[:size] == [n] * size), None)
        out[name + "_FI"] = float(row[size]) if row else 0.0  # reduced, cm-1 (absent = below threshold)
        out[name + "_k"] = float(row[size + 1]) if row else 0.0  # attoJ amu^-n/2 A^-n
    return out


if __name__ == "__main__":
    import sys

    print(constants(Path(sys.argv[1]), float(sys.argv[2])))
