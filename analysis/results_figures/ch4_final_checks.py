"""Numbers for the figures-ch4-final report (thesis docs/handoff/figures-ch4-final-report.md).

  glycine   item 3: character of five bands from the B3LYP normal modes; positions and intensities
            (harmonic, VPT2) for B3LYP and the four energy models with MACE-POLAR-1 dipoles; VPT2
            intensity relative to the strongest B3LYP VPT2 fundamental (the scale of Fig. 4.12)
  cores     item 5a: Gaussian's 'Job cpu time' against its own wall time (timing_runs.csv
            gaussian_own_s) for every campaign run, and cpu/elapsed of the C60 runs
  speedup   item 5b: B3LYP time / model time per molecule (times as in fig_cost.py)
  logratio  item 5c: Table 4.11, mean |log10(calc/exp)| to four digits and 10 to that power

    python analysis/results_figures/ch4_final_checks.py
"""

import csv
import importlib.util
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def secs(line: str) -> float:
    d, h, m, s = map(float, re.search(r"(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+) seconds", line).groups())
    return ((d * 24 + h) * 60 + m) * 60 + s


def glycine() -> None:
    from mace_gaussian.analysis.mode_matching import extract_mode_data_from_checkpoint

    ts = _load("thesis_stats", ROOT / "analysis/followup3/thesis_stats.py")
    fchk = ROOT / "campaigns/2026/comparison_results/glycine/b3lyp_6-31Gdp/glycine_freq_anharm.fchk"
    modes, freqs, *_ = extract_mode_data_from_checkpoint(str(fchk), force_harmonic=True)
    names = ["O(H)", "O(=C)", "N", "Ca", "C", "H(Ca)", "H(Ca)", "H(N)", "H(N)", "H(O)"]
    modes = np.asarray(modes).reshape(len(freqs), -1, 3)
    lab = {18: "NH2 scissoring", 19: "C=O stretch", 22: "N-H stretch (sym)", 23: "N-H stretch (asym)", 24: "O-H stretch"}
    for i in lab:
        a = (modes[i - 1] ** 2).sum(1)
        a /= a.sum()
        print(f"B3LYP mode {i:2d} {freqs[i - 1]:7.1f} cm-1: "
              + ", ".join(f"{n} {v:.2f}" for v, n in sorted(zip(a, names), reverse=True)[:3]))  # fmt: skip
    _, _, fund, harm, _ = ts.load()
    f = fund[(fund.molecule == "glycine") & (fund.dipole_model == "mace_polar1")].set_index(["energy_model", "i"])
    h = harm[(harm.molecule == "glycine") & (harm.dipole_model == "mace_polar1")].set_index(["energy_model", "i"])
    peak = f.loc["mace_anicc"].dft_intensity.max()
    rows = []
    for i, name in lab.items():
        r, hr = f.loc[("mace_anicc", i)], h.loc[("mace_anicc", i)]
        rows.append(("B3LYP", name, r.dft_harmonic, hr.dft_int, r.dft_anharmonic, r.dft_intensity))
        for m in ts.E_MODELS:
            if (m, i) in f.index:
                r, hr = f.loc[(m, i)], h.loc[(m, i)]
                rows.append((m, name, r.ml_harmonic, hr.ml_int, r.ml_anharmonic, r.ml_intensity))
    d = pd.DataFrame(rows, columns=["method", "band", "harm_cm", "harm_int", "vpt2_cm", "vpt2_int"])
    d["vpt2_int/B3LYP peak"] = d.vpt2_int / peak
    print(f"strongest B3LYP VPT2 fundamental {peak:.1f} km/mol")
    print(d.round(4).to_string(index=False))


def cores() -> None:
    runs = list(csv.DictReader((ROOT / "analysis/ch5_timing/timing_runs.csv").read_text().splitlines()))
    ratio = []
    for r in runs:
        log = ROOT / f"campaigns/2026/comparison_results/{r['molecule']}/{r['energy_model']}_{r['dipole_model']}/gaussian_freq.log"
        cpu = sum(secs(x) for x in re.findall(r"Job cpu time:.*", log.read_text()))
        ratio.append((int(r["atoms"]), cpu / float(r["gaussian_own_s"])))
    big = [q for n, q in ratio if n >= 14]
    print(f"campaign (%NProcShared=2): cpu / Gaussian's own wall time, median {np.median([q for _, q in ratio]):.2f} "
          f"(all {len(ratio)} runs); >= 14 atoms {min(big):.2f}-{max(big):.2f}")  # fmt: skip
    for log in sorted((ROOT / "campaigns/c60/comparison_results/c60").glob("*/gaussian_freq.log")):
        t = log.read_text()
        nproc = re.search(r"%NProcShared=(\d+)", t).group(1)
        cpu = sum(secs(x) for x in re.findall(r"Job cpu time:.*", t)) / 3600
        el = sum(secs(x) for x in re.findall(r"Elapsed time:.*", t)) / 3600
        print(f"C60 {log.parent.name:24s} nproc {nproc}  cpu {cpu:5.2f} h  elapsed {el:5.2f} h")


def speedup() -> None:
    fc = _load("fig_cost", Path(__file__).resolve().parent / "fig_cost.py")
    t = fc.timings()
    s = 1 / t[fc.E].div(t.b3lyp, axis=0)
    s.insert(0, "b3lyp_min", t.b3lyp)
    s.insert(0, "atoms", t.atoms)
    print(s.sort_values(["atoms", "b3lyp_min"]).round(2).to_string())


def logratio() -> None:
    rows = list(csv.DictReader((ROOT / "analysis/followup2/b7_levels_exp/results/b7_cells.csv").read_text().splitlines()))
    for c in ["b3lyp_631gdp", "b3lyp_avtz", "pbe0_tzvp", "wb97mv_tzvpd", "wb97md3_tzvppd", "mace_polar1", "mace_mdp",
              "mace_ml", "polar1_field"]:  # fmt: skip
        m = np.mean([abs(math.log10(float(r[c]) / float(r["exp"]))) for r in rows])
        print(f"{c:16s} n {len(rows)}  mean|log10| {m:.4f}  10^ {10**m:.4f}")


if __name__ == "__main__":
    for f in (glycine, cores, speedup, logratio):
        print(f"\n## {f.__name__}")
        f()
