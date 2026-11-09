"""Ch. 5 timing breakdown: where the wall time of each MACE run goes.

For every completed ML run of campaigns/2026 the wall time (results.json runtime_s) is split into
- energy model: the server's per-call time minus the dipole part (energy, forces, Hessian, and
  the small file I/O of each call);
- dipole model: dipole, dipole derivatives (and polarizability) per call;
- Gaussian's own work: Gaussian's elapsed time minus the summed call times (its own links, mainly
  the VPT2 analysis in l717, plus process start-up);
- Python outside Gaussian: runtime minus Gaussian's elapsed time (re-optimisation on the MACE
  surface, writing the input, parsing the log).

Sources (read only): the batch logs campaigns/2026/logs/ml_{acids,alcohols,big,inorganic}.log.
- stdout lines (in order): "FREQUENCY CALCULATION: <energy> + <dipole>", one "Gaussian call
  completed in X s" per server call (mace_gaussian/workflow.py:319, timer around the whole call),
  "Completed in N seconds";
- logging lines (in order, time-stamped to the ms): "Using dipole calculator: <name>" at the start
  and "Dipole calculated" at the end of each call's dipole step (workflow.py:195, :226).
The two streams are buffered separately, so they interleave out of order in the file; each stream
is in order by itself. Dipole steps are therefore assigned to runs by count (run k takes the next
n_k dipole steps), and the dipole name on every step is checked against the run.
Each log block is matched to the run's results.json by its runtime (to 0.1 s).

Output: timing_runs.csv (one row per run) and timing_summary.md next to this script.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAMP = ROOT / "campaigns" / "2026"
LOGS = ["ml_acids.log", "ml_alcohols.log", "ml_big.log", "ml_inorganic.log"]
GEO = Path.home() / "thesis" / "data" / "geometries"
LAB = {"mace_anicc": "ANI-cc", "mace_omol": "OMOL", "mace_off": "OFF", "mace_polar": "POLAR",
       "mace_polar1": "POLAR-1", "mace_mdp": "MDP", "mace_ml": "MACE4IR"}  # fmt: skip
DIPOLE_NAME = {"mace_ml": "mace_ml", "mace_polar1": "mace_polar1", "mace_mdp": "mace_mdp"}
TS = re.compile(r"(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) - INFO - (Using dipole calculator: (\S+)|Dipole calculated)")


STDOUT = re.compile(
    r"\[\d+/\d+\] Running (?P<mol>\S+)\.\.\."
    r"|FREQUENCY CALCULATION: (?P<energy>\S+) \+ (?P<dipole>\S+)"
    r"|Gaussian call completed in (?P<call>[\d.]+)s"
    r"|^  Completed in (?P<runtime>[\d.]+) seconds"
)


def parse(log: Path):
    """Two independent passes: a stdout line and a logging line can share one physical line."""
    text = log.read_text(errors="replace")
    runs, dip = [], []
    mol, cur, start = None, None, None
    for line in text.splitlines():
        for m in STDOUT.finditer(line):
            if m["mol"]:
                mol = m["mol"]
            elif m["energy"]:
                cur = dict(molecule=mol, energy_model=m["energy"], dipole_model=m["dipole"], calls=[], runtime=None)
                runs.append(cur)
            elif m["call"] and cur is not None:
                cur["calls"].append(float(m["call"]))
            elif m["runtime"] and cur is not None:
                cur["runtime"] = float(m["runtime"])
    for m in TS.finditer(text):
        t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f").timestamp()
        if m.group(3):
            assert start is None, (log.name, "dipole step started twice")
            start = (t, m.group(3))
        else:
            assert start is not None, (log.name, "dipole step ended without a start")
            dip.append((start[1], t - start[0]))
            start = None
    # assign dipole steps to runs by count, checking the dipole model name
    k = 0
    for r in runs:
        n = len(r["calls"])
        steps = dip[k : k + n]
        k += n
        assert all(name == DIPOLE_NAME[r["dipole_model"]] for name, _ in steps), (log.name, r["molecule"])
        r["dipole_times"] = [d for _, d in steps]
    assert k == len(dip), (log.name, k, len(dip))
    return runs


def atoms(mol: str) -> int:
    return int((GEO / f"{mol}.xyz").read_text().split()[0])


def main():
    rows = []
    for log in LOGS:
        for r in parse(CAMP / "logs" / log):
            if r["runtime"] is None:
                continue  # not completed (POLAR beta-glucose, refused ANI-cc runs)
            res = CAMP / "comparison_results" / r["molecule"] / f"{r['energy_model']}_{r['dipole_model']}" / "results.json"
            j = json.loads(res.read_text())
            if abs(j["runtime_s"] - r["runtime"]) > 0.11:
                raise SystemExit(f"runtime mismatch {res}: log {r['runtime']} vs {j['runtime_s']:.2f}")
            g = j["gaussian_timing"]["total_elapsed_s"]
            model, dipole = sum(r["calls"]), sum(r["dipole_times"])
            rows.append(dict(
                molecule=r["molecule"], atoms=atoms(r["molecule"]), energy_model=r["energy_model"],
                dipole_model=r["dipole_model"], calls=len(r["calls"]), wall_s=j["runtime_s"],
                gaussian_elapsed_s=g, model_s=model, energy_s=model - dipole, dipole_s=dipole,
                gaussian_own_s=g - model, python_s=j["runtime_s"] - g,
                first_call_s=r["calls"][0], median_call_s=float(np.median(r["calls"])),
                median_dipole_s=float(np.median(r["dipole_times"])),
            ))  # fmt: skip
    df = pd.DataFrame(rows)
    for c in ("energy", "dipole", "gaussian_own", "python"):
        df[f"{c}_%"] = 100 * df[f"{c}_s"] / df.wall_s
    df["model_%"] = df["energy_%"] + df["dipole_%"]
    df.to_csv(HERE / "timing_runs.csv", index=False)

    out = ["# Ch. 5 timing breakdown (generated by timing_breakdown.py)", "",
           f"{len(df)} completed ML runs. Shares of the wall time of each run, in %.", ""]  # fmt: skip
    share_cols = ["energy_%", "dipole_%", "gaussian_own_%", "python_%"]

    def md(t: pd.DataFrame, nd=0) -> str:
        cols = list(t.columns)
        lines = ["| " + " | ".join([t.index.name or ""] + [str(c) for c in cols]) + " |", "|" + "---|" * (len(cols) + 1)]
        for idx, row in t.iterrows():
            idx = " / ".join(map(str, idx)) if isinstance(idx, tuple) else str(idx)
            lines.append("| " + " | ".join([idx] + [f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in row]) + " |")
        return "\n".join(lines)

    # 1. per energy model, summed over all runs (time-weighted) and as a mean of run shares
    tot = df.groupby("energy_model")[["wall_s", "energy_s", "dipole_s", "gaussian_own_s", "python_s"]].sum()
    t = pd.DataFrame({c.replace("_s", " %"): 100 * tot[c] / tot.wall_s for c in ["energy_s", "dipole_s", "gaussian_own_s", "python_s"]})
    t.insert(0, "wall h", tot.wall_s / 3600)
    t.index = [LAB[i] for i in t.index]
    t.index.name = "energy model (time-weighted, all runs)"
    out += ["## Per energy model, all runs (time-weighted)", "", md(t, 1), ""]
    # 2. per dipole model
    tot = df.groupby("dipole_model")[["wall_s", "energy_s", "dipole_s", "gaussian_own_s", "python_s"]].sum()
    t = pd.DataFrame({c.replace("_s", " %"): 100 * tot[c] / tot.wall_s for c in ["energy_s", "dipole_s", "gaussian_own_s", "python_s"]})
    t.insert(0, "median dipole s/call", df.groupby("dipole_model").median_dipole_s.median())
    t.index = [LAB[i] for i in t.index]
    t.index.name = "dipole model (time-weighted, all runs)"
    out += ["## Per dipole model, all runs (time-weighted)", "", md(t, 2), ""]
    # 3. per molecule (by size) and energy model: mean over the three dipole runs
    g = df.groupby(["atoms", "molecule", "energy_model"])
    t = g[share_cols + ["wall_s"]].mean()
    t["calls"] = g.calls.first()
    t["wall_s"] = t.wall_s / 60
    t = t.rename(columns={"wall_s": "wall min", "energy_%": "energy", "dipole_%": "dipole",
                          "gaussian_own_%": "Gaussian", "python_%": "Python"})  # fmt: skip
    t = t.reset_index().assign(energy_model=lambda x: x.energy_model.map(LAB)).set_index(["atoms", "molecule", "energy_model"])
    t.index.name = None
    out += ["## Per molecule and energy model (mean over the three dipole runs; shares in %, wall in min)", "",
            md(t[["calls", "wall min", "energy", "dipole", "Gaussian", "Python"]], 1), ""]  # fmt: skip
    # 4. size bins
    bins = [0, 4, 6, 12, 15, 30]
    labels = ["3-4 atoms", "5-6", "8-12", "14-15", "17-24"]
    df["size"] = pd.cut(df.atoms, bins, labels=labels)
    s = df.groupby(["size", "energy_model"], observed=True)[["wall_s", "energy_s", "dipole_s", "gaussian_own_s", "python_s"]].sum()
    t = pd.DataFrame({c.replace("_s", " %"): 100 * s[c] / s.wall_s for c in ["energy_s", "dipole_s", "gaussian_own_s", "python_s"]})
    t = t.reset_index().assign(energy_model=lambda x: x.energy_model.map(LAB)).set_index(["size", "energy_model"])
    out += ["## By molecule size and energy model (time-weighted)", "", md(t, 1), ""]
    (HERE / "timing_summary.md").write_text("\n".join(out))
    print("\n".join(out[:30]))


if __name__ == "__main__":
    main()
