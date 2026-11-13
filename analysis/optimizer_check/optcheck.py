"""Optimizer check on methane with MACE-OMOL-0: ASE L-BFGS vs Gaussian Berny.

Both paths start from the same geometry (ASE's built-in CH4, not a minimum of any model).
A: ASE L-BFGS to fmax 1e-4 eV/A (pipeline setting), then Gaussian `freq` on that geometry.
B: Gaussian `opt freq`, MACE-OMOL-0 as external engine (Berny finds the minimum itself).
C: as B with `opt=verytight`.
Harmonic frequencies are compared. Dipoles: MACE4IR (no effect on frequencies).
"""

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from ase.build import molecule
from ase.optimize import LBFGS

sys.path.insert(0, "/home/mot/mace_gaussian")
from mace_gaussian.gaussian.io import DEFAULT_HELPER_SCRIPT, ase_to_gjf  # noqa: E402
from mace_gaussian.gaussian.parser import parse_gaussian_log  # noqa: E402
from mace_gaussian.gaussian.runner import run_gaussian_with_zmq  # noqa: E402
from mace_gaussian.workflow import OPT_FMAX, OPT_MAX_STEPS, calculator, run_next_calculation  # noqa: E402

OUT = Path(__file__).resolve().parent
calc = calculator("mace_omol")


def gaussian(tag, mol, keyword):
    d = OUT / tag
    d.mkdir(exist_ok=True)
    route = f'# {keyword}\n# external="{DEFAULT_HELPER_SCRIPT}"'
    ase_to_gjf(mol, f"{tag}.gjf", route=route, output_dir=d)

    def on_request(msg):
        run_next_calculation(mol, msg, calc, dipole_method="mace_ml", calculate_derivatives=True)
        return "ready"

    with tempfile.TemporaryDirectory(prefix="mg-") as sock:
        ipc = str(Path(sock) / "zmq.ipc")
        env = os.environ.copy()
        env["MACE_IPC_PATH"] = ipc
        run_gaussian_with_zmq(f"{tag}.gjf", on_request=on_request, ipc_file=ipc, cwd=str(d), env=env)
    p = parse_gaussian_log(str(d / f"{tag}.log"))
    return sorted(x["freq_cm"] for x in p["harmonic"])


rng = np.random.default_rng(20261111)
starts = {"near": molecule("CH4")}
d = molecule("CH4")
d.positions += rng.uniform(-0.05, 0.05, d.positions.shape)
starts["distorted"] = d
res = {}
for sname, start in starts.items():
    start.info["charge"], start.info["spin"] = 0.0, 1.0
    a = start.copy()
    a.calc = calc
    opt = LBFGS(a, logfile=str(OUT / f"lbfgs_{sname}.log"))
    conv = opt.run(fmax=OPT_FMAX, steps=OPT_MAX_STEPS)
    f_a = gaussian(f"{sname}_A_ase_lbfgs", a, "freq")
    res[f"{sname}_A_ase_lbfgs"] = {"freq": f_a, "lbfgs_steps": opt.get_number_of_steps(), "converged": bool(conv),
                                    "fmax_eV_A": float(np.linalg.norm(a.get_forces(), axis=1).max())}
    for tag, kw in [("B_berny", "opt freq"), ("C_berny_verytight", "opt=verytight freq")]:
        m = start.copy()
        m.calc = calc
        f = gaussian(f"{sname}_{tag}", m, kw)
        res[f"{sname}_{tag}"] = {"freq": f, "max_abs_diff_vs_A": float(np.max(np.abs(np.array(f) - np.array(f_a))))}

(OUT / "optcheck.json").write_text(json.dumps(res, indent=2))
for k, v in res.items():
    print(k, [round(x, 4) for x in v["freq"]], {kk: vv for kk, vv in v.items() if kk != "freq"})
