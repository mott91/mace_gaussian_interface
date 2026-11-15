"""MACE scans along the target X-H stretches: energy and exact curvature at every grid point.

For each target and energy model: the model's own minimum and normal mode (common.load_mode), the
three grids of common.grid_points, and at each point the energy and the analytic (autograd) Hessian
from the pipeline (mace_gaussian.workflow.calculator / calculate_hessian, double precision, as in
the campaign). The curvature along the mode is k = d^T H d. Check: k at Q = 0 must give the
harmonic wavenumber of the campaign run.

Run in the campaign environment (mace4ir_v3_accel), from a scratch directory:
    python analysis/surface_scans/scan_mace.py            -> data/mace_<target>_<model>.json
"""

from __future__ import annotations

import logging
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import common as c  # noqa: E402

sys.path.insert(0, str(c.ROOT))
logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")
from ase import Atoms  # noqa: E402

from mace_gaussian import workflow as wf  # noqa: E402

TO_AJ_A2 = c.HARTREE_AJ / c.BOHR_A**2  # pipeline Hessian is in Hartree / Bohr^2


def main(models=c.MODELS):
    for model in models:
        t0 = time.time()
        calc = wf.calculator(model)
        t_load = time.time() - t0
        for key, _mol, _i in c.TARGETS:
            t1 = time.time()
            mode = c.load_mode(key, model)
            grids = c.grid_points(mode)
            dflat = mode.d.ravel()
            out = {"key": key, "method": model, "freq_run": mode.freq, "own_mode": mode.own_mode,
                   "dr_dq": mode.dr_dq, "dnorm": mode.dnorm, "purity": mode.purity, "grids": {}}  # fmt: skip
            for name, qs in {"zero": np.array([0.0]), **grids}.items():
                e, k = [], []
                for q in qs:
                    a = Atoms(numbers=mode.z, positions=mode.geometry(q))
                    a.info.update(charge=0.0, spin=1.0)
                    a.calc = calc
                    e.append(a.get_potential_energy() * c.EV_AJ)
                    h = np.asarray(wf.calculate_hessian(a, calc, len(a))) * TO_AJ_A2
                    k.append(float(dflat @ h @ dflat))
                out["grids"][name] = {"q": qs.tolist(), "dr": mode.bond(qs).tolist(), "E": e, "k": k}
            w0 = float(c.k_to_wavenumber(out["grids"]["zero"]["k"][0]))
            out["freq_scan"] = w0
            out["seconds"] = time.time() - t1 + t_load / len(c.TARGETS)
            assert abs(w0 - mode.freq) < 0.5, (key, model, w0, mode.freq)
            c.save(f"mace_{key}_{model}.json", out)
            print(f"{key:12s} {model:11s} run {mode.freq:8.2f} scan {w0:8.2f} cm-1  {out['seconds']:.0f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or c.MODELS)
