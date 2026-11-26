#!/usr/bin/python3
"""Gaussian External helper that replays the reference point from a finished job's fchk.

Used for the link-717 core-count test (thesis docs/handoff/l717-cores-prompt.md). A
Freq=(Anharm,ReadAnharm) Geom=AllCheck job with DataSrc=Chk still asks External once for the
energy and derivatives at the reference geometry. This helper answers that call with the
values stored in the campaign fchk (energy, gradient, dipole, dipole derivatives, Hessian;
polarizability zero, as in the campaign) and refuses any other geometry, so no new energies
enter the test.

Environment:
    REPLAY_FCHK  path of the campaign fchk
    REPLAY_LOG   file to which each call is appended (optional)

Plain python3, no third-party packages.
"""

import os
import sys
import time

TOL_BOHR = 1e-5


def read_fchk(path):
    """Return {label: value or list} for the fields we need."""
    want = {
        "Current cartesian coordinates",
        "Total Energy",
        "Cartesian Gradient",
        "Dipole Moment",
        "Dipole Derivatives",
        "Cartesian Force Constants",
    }
    out = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        label = line[:43].strip()
        if label in want and label not in out:
            rest = line[43:].split()
            if "N=" in rest:
                n = int(rest[-1])
                vals = []
                i += 1
                while len(vals) < n:
                    vals.extend(float(x) for x in lines[i].split())
                    i += 1
                out[label] = vals
                continue
            out[label] = float(rest[-1])
        i += 1
    missing = want - out.keys()
    if missing:
        raise SystemExit(f"fchk lacks {sorted(missing)}")
    return out


def fmt(vals):
    return "".join(f"{v:20.12E}" for v in vals).replace("E", "D") + "\n"


def main():
    infile, outfile = sys.argv[2], sys.argv[3]
    fchk = read_fchk(os.environ["REPLAY_FCHK"])
    log = os.environ.get("REPLAY_LOG")

    with open(infile) as f:
        lines = f.readlines()
    natoms, deriv = (int(x) for x in lines[0].split()[:2])
    coords = []
    for line in lines[1 : natoms + 1]:
        coords.extend(float(x) for x in line.split()[1:4])
    ref = fchk["Current cartesian coordinates"]
    if len(coords) != len(ref):
        sys.stderr.write("replay_helper: atom count differs from reference\n")
        return 1
    dev = max(abs(a - b) for a, b in zip(coords, ref))

    if log:
        with open(log, "a") as f:
            f.write(f"{time.strftime('%T')} natoms={natoms} deriv={deriv} maxdev_bohr={dev:.2e}\n")
    if dev > TOL_BOHR:
        sys.stderr.write(f"replay_helper: geometry differs from reference by {dev:.2e} bohr\n")
        return 1

    n3 = 3 * natoms
    with open(outfile, "w") as f:
        f.write(fmt([fchk["Total Energy"], *fchk["Dipole Moment"]]))
        g = fchk["Cartesian Gradient"]
        for a in range(natoms):
            f.write(fmt(g[3 * a : 3 * a + 3]))
        f.write(fmt([0.0, 0.0, 0.0]))
        f.write(fmt([0.0, 0.0, 0.0]))
        dd = fchk["Dipole Derivatives"]
        for k in range(n3):
            f.write(fmt(dd[3 * k : 3 * k + 3]))
        if deriv >= 2:
            fc = fchk["Cartesian Force Constants"]
            assert len(fc) == n3 * (n3 + 1) // 2
            for k in range(0, len(fc), 3):
                f.write(fmt(fc[k : k + 3]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
