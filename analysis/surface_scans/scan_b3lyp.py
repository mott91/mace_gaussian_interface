"""B3LYP scans along the target X-H stretches with Gaussian 16 (campaign settings).

Same grids as the MACE scans (common.grid_points) on the B3LYP minimum and normal mode of the
campaign run. Level and settings as the campaign baselines: B3LYP/6-31G(d,p), int=ultrafine,
nosymm (input orientation, so the geometry frame of the campaign fchk is kept).
  - fine grid, the Gaussian step points and Q = 0: `freq` jobs -> analytic Hessian (fchk), energy;
    Gaussian runs these with tight SCF convergence;
  - wide grid: single points with scf=tight (energies only).
The curvature is k = d^T H d from the analytic Hessians, so it is not formed from differences of
SCF energies.

Work files go to analysis/surface_scans/work/ (not committed); results to data/b3lyp_<target>.json.
Resumable: finished jobs are parsed again, not rerun.   python scan_b3lyp.py [workers]
"""

from __future__ import annotations

import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

WORK = c.HERE / "work"
NPROC = 2
ROUTE_FREQ = "# freq b3lyp/6-31G(d,p) nosymm int=ultrafine"
ROUTE_SP = "# b3lyp/6-31G(d,p) nosymm int=ultrafine scf=tight"


def write_input(path: Path, mode: c.Mode, q: float, freq: bool) -> None:
    lines = [f"%chk={path.stem}.chk", "%mem=2GB", f"%NProcShared={NPROC}", ROUTE_FREQ if freq else ROUTE_SP, "",
             f"{mode.key} Q={q:+.10f}", "", "0 1"]  # fmt: skip
    lines += [f"{c.SYMBOL[z]} {a:16.10f} {b:16.10f} {cc:16.10f}" for z, (a, b, cc) in zip(mode.z, mode.geometry(q))]
    path.write_text("\n".join(lines) + "\n\n")


def run(path: Path, freq: bool) -> float:
    """Run one job unless it already terminated normally; returns seconds spent."""
    log = path.with_suffix(".log")
    if log.exists() and "Normal termination" in log.read_text()[-2000:] and (not freq or path.with_suffix(".fchk").exists()):
        return 0.0
    t = time.time()
    with path.open() as fin, log.open("w") as fout:
        subprocess.run(["g16"], stdin=fin, stdout=fout, cwd=path.parent, check=True)
    if freq:
        subprocess.run(["formchk", path.stem + ".chk", path.stem + ".fchk"], cwd=path.parent, check=True,
                       capture_output=True)  # fmt: skip
    path.with_suffix(".chk").unlink(missing_ok=True)
    return time.time() - t


def parse(path: Path, mode: c.Mode, freq: bool) -> tuple[float, float | None]:
    log = path.with_suffix(".log").read_text()
    assert "Normal termination" in log[-2000:], path
    e = float([ln for ln in log.splitlines() if "SCF Done:" in ln][-1].split("=")[1].split()[0]) * c.HARTREE_AJ
    if not freq:
        return e, None
    a = c.fchk_arrays(path.with_suffix(".fchk"), ["Cartesian Force Constants"])
    n3 = 3 * len(mode.z)
    h = np.zeros((n3, n3))
    h[np.tril_indices(n3)] = a["Cartesian Force Constants"]
    h = h + h.T - np.diag(np.diag(h))
    h *= c.HARTREE_AJ / c.BOHR_A**2
    d = mode.d.ravel()
    return e, float(d @ h @ d)


def main(workers: int = 2):
    for key, _mol, _i in c.TARGETS:
        t0 = time.time()
        mode = c.load_mode(key, "b3lyp")
        grids = {"zero": np.array([0.0]), **c.grid_points(mode)}
        wd = WORK / key
        wd.mkdir(parents=True, exist_ok=True)
        jobs = []
        for name, qs in grids.items():
            freq = name != "wide"
            for j, q in enumerate(qs):
                p = wd / f"{name}_{j:03d}.gjf"
                write_input(p, mode, q, freq)
                jobs.append((p, freq))
        with ThreadPoolExecutor(workers) as ex:
            spent = list(ex.map(lambda pj: run(*pj), jobs))
        out = {"key": key, "method": "b3lyp", "freq_run": mode.freq, "own_mode": mode.own_mode, "dr_dq": mode.dr_dq,
               "dnorm": mode.dnorm, "purity": mode.purity, "grids": {}, "seconds": sum(spent),
               "wall_seconds": time.time() - t0}  # fmt: skip
        for name, qs in grids.items():
            freq = name != "wide"
            res = [parse(wd / f"{name}_{j:03d}.gjf", mode, freq) for j in range(len(qs))]
            out["grids"][name] = {"q": qs.tolist(), "dr": mode.bond(qs).tolist(), "E": [r[0] for r in res],
                                  "k": [r[1] for r in res] if freq else None}  # fmt: skip
        w0 = float(c.k_to_wavenumber(out["grids"]["zero"]["k"][0]))
        out["freq_scan"] = w0
        assert abs(w0 - mode.freq) < 0.5, (key, w0, mode.freq)
        c.save(f"b3lyp_{key}.json", out)
        print(f"{key:12s} b3lyp run {mode.freq:8.2f} scan {w0:8.2f} cm-1  cpu-wall {out['seconds']:.0f} s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
