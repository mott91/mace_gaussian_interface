"""Follow-up 2, B7 (local prep): everything the cluster analysis needs that does not depend on the
new jobs, in bundle.npz + bundle.json, so analyse.py runs with numpy alone on the cluster.

Per molecule (methanol, ethanol, isopropanol, 1-butanol, SO2):
- the B3LYP/6-31G(d,p) campaign minimum, masses, APT and normal modes (b1_b2_analyse.normal_modes,
  the B1/B3/B4 modes), and the harmonic frequencies and intensities from the campaign log (as B4);
- the dipole-model APTs at that geometry (b1_apt_<model>.npz; POLAR-1 field response with the B2
  sign, -stored, field 0.001);
- the NIST region integrals (b4_experiment.read_jdx / exp_region, the B4 code);
- the B4 harmonic table values (b4_tables.md), for the reproduction check;
- methanol: the four B3 level APTs (b3_analyse.level_apt) and B3 wall times, for the consistency
  check and the cost table;
- the campaign ML runs (results.json, gaussian_freq.log): runtime, number of external calls
  (one Hessian + APT each), hardware.

    /home/mot/micromamba/envs/mace4ir_v3_accel/bin/python analysis/followup2/b7_levels_exp/prep_bundle.py
"""

from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
F2 = HERE.parent
ROOT = F2.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(F2))
logging.disable(logging.CRITICAL)
import b3_analyse as b3  # noqa: E402
import b4_experiment as b4  # noqa: E402
from b1_b2_analyse import CAMP, fchk, normal_modes  # noqa: E402

from mace_gaussian.gaussian.parser import parse_gaussian_log  # noqa: E402

MOLS = ["methanol", "ethanol", "isopropanol", "1-butanol", "sulfur_dioxide"]
DIP = {"mace_polar1": "POLAR-1", "polar1_field": "POLAR-1 (field)", "mace_mdp": "MDP",
       "mace_ml": "MACE4IR"}  # fmt: skip
ENERGY = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]
B3_LEVELS = ["b3lyp_avtz", "pbe0_tzvp", "wb97mv_tzvpd", "wb97md3_tzvppd"]
B3_DIRS = {"b3lyp_avtz": "g_b3lyp_avtz_fix", "pbe0_tzvp": "g_pbe0_tzvp_fix",
           "wb97mv_tzvpd": "o_wb97mv_fix", "wb97md3_tzvppd": "o_wb97md3bj_fix"}  # fmt: skip


def b4_table_values():
    """Harmonic columns of b4_tables.md: {mol: {region label: {column: value}}}."""
    out, mol, cols = {}, None, None
    for line in (F2 / "b4_tables.md").read_text().split("\n"):
        if line.startswith("## ") and "Summary" not in line:
            mol = line[3:].split(" (")[0]
            out[mol] = {}
        elif mol and line.startswith("| region"):
            cols = [c.strip() for c in line.strip("|").split("|")]
        elif mol and line.startswith("| ") and not line.startswith("|---") and cols:
            c = [v.strip() for v in line.strip("|").split("|")]
            if c[0] in ("method",):
                continue
            out[mol][c[0]] = {k: float(v) for k, v in zip(cols[1:], c[1:]) if v}
        elif line.startswith("## Summary"):
            mol = None
    return out


def ml_runs(mol):
    rows = []
    for e in ENERGY:
        for d in ("mace_polar1", "mace_mdp", "mace_ml"):
            p = CAMP / "comparison_results" / mol / f"{e}_{d}"
            if not (p / "results.json").exists():
                continue
            r = json.loads((p / "results.json").read_text())
            log = (p / "gaussian_freq.log").read_text()
            ncall = log.count("Running external command")
            nproc = int(re.search(r"%NProcShared=(\d+)", (p / "gaussian_freq.gjf").read_text()).group(1))
            v = r["version_info"]
            rows.append(dict(energy=e, dipole=d, runtime_s=r["runtime_s"],
                             gaussian_elapsed_s=r["gaussian_timing"]["total_elapsed_s"],
                             gaussian_cpu_s=r["gaussian_timing"]["total_cpu_s"], n_calls=ncall,
                             nproc_gaussian=nproc, cpu=v["cpu_model"], gpu=v.get("gpu_name"),
                             torch=v["torch_version"], mace=v["mace_version"]))  # fmt: skip
    return rows


def main():
    arrays, meta = {}, {"mols": MOLS, "dipoles": DIP, "regions": {}, "exp": {}, "b4": b4_table_values(),
                        "ml_runs": {}, "b3_methanol_wall_s": {}, "campaign_b3lyp": {}}  # fmt: skip
    ml = {d: np.load(F2 / f"b1_apt_{'mace_polar1' if d == 'polar1_field' else d}.npz") for d in DIP}
    for mol in MOLS:
        f = fchk(mol)
        freqs, disp = normal_modes(f["x"], f["m"], f["hess"])
        g = parse_gaussian_log(str(CAMP / mol / f"{mol}_freq_anharm.log"))["harmonic"]
        arrays[f"{mol}/x"] = f["x"]
        arrays[f"{mol}/m"] = f["m"]
        arrays[f"{mol}/apt_b3lyp_631gdp"] = f["apt"]
        arrays[f"{mol}/freqs"] = freqs
        arrays[f"{mol}/disp"] = disp
        arrays[f"{mol}/log_freq"] = np.array([r["freq_cm"] for r in g])
        arrays[f"{mol}/log_int"] = np.array([r["ir_intensity"] for r in g])
        for d in DIP:
            arrays[f"{mol}/apt_{d}"] = (-ml[d][f"{mol}/apt_field_0.001"] if d == "polar1_field"
                                        else ml[d][f"{mol}/apt"])  # fmt: skip
        x, s = b4.read_jdx(b4.EXP / f"{mol}.jdx")
        regs = b4.REGIONS.get(mol, b4.REGIONS["default"])
        meta["regions"][mol] = [[lab, lo, hi] for lab, lo, hi in regs]
        meta["exp"][mol] = {lab: b4.exp_region(x, s, lo, hi) for lab, lo, hi in regs}
        meta["ml_runs"][mol] = ml_runs(mol)
        log = (CAMP / mol / f"{mol}_freq_anharm.log").read_text()
        el = re.findall(r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+) seconds", log)
        meta["campaign_b3lyp"][mol] = dict(
            elapsed_s=[86400 * int(a) + 3600 * int(b) + 60 * int(c) + float(e) for a, b, c, e in el],
            nproc=int(re.search(r"%NProcShared=(\d+)", log).group(1)),
            route=re.search(r"\n (#.*)\n", log).group(1).strip())  # fmt: skip
    # methanol: B3 level APTs (the B3 code path) and wall times
    f0 = fchk("methanol")
    for lev in B3_LEVELS:
        a, _dx = b3.level_apt(lev, "methanol", f0)
        arrays[f"methanol/b3_apt_{lev}"] = a
        d = b3.B3 / B3_DIRS[lev] / "methanol"
        t0, t1 = (int((d / n).read_text()) for n in ("start_epoch.txt", "end_epoch.txt"))
        meta["b3_methanol_wall_s"][lev] = t1 - t0
    np.savez(HERE / "bundle.npz", **arrays)
    (HERE / "bundle.json").write_text(json.dumps(meta, indent=1))
    print("bundle:", len(arrays), "arrays;", {m: len(meta["ml_runs"][m]) for m in MOLS}, "ML runs")
    print("exp:", json.dumps(meta["exp"], indent=None)[:600])
    print("b3 methanol wall:", meta["b3_methanol_wall_s"])
    print("campaign b3lyp:", {m: (v["nproc"], v["elapsed_s"]) for m, v in meta["campaign_b3lyp"].items()})


if __name__ == "__main__":
    main()
