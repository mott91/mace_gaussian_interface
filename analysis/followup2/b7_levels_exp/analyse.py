"""Follow-up 2, B7 (analysis): higher levels and dipole models against NIST gas-phase intensities.

numpy only (runs as the dependent batch job on the cluster, or locally). Inputs:
- bundle.npz / bundle.json (prep_bundle.py): B3LYP/6-31G(d,p) geometry, modes, APT and log
  intensities; dipole-model APTs; NIST region integrals; B4 values; campaign ML run times;
- <root>/<level>/<mol>/ (this follow-up's jobs) and <b3root>/<level>/methanol/ (B3, reused).

Method, as B3 and B4:
- every APT is at the B3LYP/6-31G(d,p) minimum and is projected on the B3LYP/6-31G(d,p) modes;
- the region value of a method is the sum over the fundamentals whose B3LYP harmonic frequency
  falls in the region (the reference row uses Gaussian's own harmonic log values, as B4);
- regions: fingerprint 575-1550, C-H 2600-3300, O-H 3300-3975; SO2 nu1 1050-1250, nu3 1250-1450.

    python analyse.py [root] [b3root]      (defaults: this directory, ../b3_reference)
Writes <root>/results/b7_tables.md, b7_cells.csv, b7_bands.csv, b7_cost.csv (+ b7_calc_exp.png).
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
B3ROOT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parent / "b3_reference"
OUT = ROOT / "results"
KMMOL = 974.8801  # km/mol per (e^2 / amu)
NEW = ["ethanol", "isopropanol", "1-butanol", "sulfur_dioxide"]
LEVELS = {"b3lyp_avtz": "B3LYP/aug-cc-pVTZ", "pbe0_tzvp": "PBE0/def2-TZVP",
          "wb97mv_tzvpd": "ωB97M-V/def2-TZVPD", "wb97md3_tzvppd": "ωB97M-D3(BJ)/def2-TZVPPD"}  # fmt: skip
DIRS = {"b3lyp_avtz": "g_b3lyp_avtz_fix", "pbe0_tzvp": "g_pbe0_tzvp_fix",
        "wb97mv_tzvpd": "o_wb97mv_fix", "wb97md3_tzvppd": "o_wb97md3bj_fix"}  # fmt: skip
NUMERICAL = {"wb97mv_tzvpd"}  # ORCA NumFreq (no analytic Hessian with VV10)
REF = "b3lyp_631gdp"
DIP = {"mace_polar1": "MACE-POLAR-1 (charge)", "mace_mdp": "MACE-MDP", "mace_ml": "MACE4IR",
       "polar1_field": "MACE-POLAR-1 (field response)"}  # fmt: skip
TRAINED = {"mace_polar1": "wb97mv_tzvpd", "polar1_field": "wb97mv_tzvpd",
           "mace_mdp": "wb97md3_tzvppd", "mace_ml": "pbe0_tzvp"}  # fmt: skip
LAB = {REF: "B3LYP/6-31G(d,p) (reference)", **LEVELS, **DIP}
ROWS = [REF, *LEVELS, "mace_polar1", "mace_mdp", "mace_ml", "polar1_field"]


def fchk_any(path, key):
    t = path.read_text()
    i = t.index(key + " ")
    n = int(t[i : t.index("\n", i)].split("N=")[1])
    return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])


def orca_hess(path):
    lines = path.read_text().split("\n")
    i = lines.index("$dipole_derivatives")
    n3 = int(lines[i + 1])
    apt = np.array([[float(v) for v in lines[i + 2 + k].split()] for k in range(n3)])
    j = lines.index("$atoms")
    na = int(lines[j + 1])
    x = np.array([[float(v) for v in lines[j + 2 + k].split()[2:5]] for k in range(na)])
    return apt, x


def run_dir(lev, mol):
    base = B3ROOT if mol == "methanol" else ROOT
    return base / DIRS[lev] / mol


def wall(d):
    try:
        return int((d / "end_epoch.txt").read_text()) - int((d / "start_epoch.txt").read_text())
    except (OSError, ValueError):
        return None


def program_time(d, mol):
    """(elapsed s from the program's own output, normal termination)."""
    log, out = d / f"{mol}.log", d / f"{mol}.out"
    if log.exists():
        t = log.read_text()
        el = re.findall(r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+) seconds", t)
        return (sum(86400 * int(a) + 3600 * int(b) + 60 * int(c) + float(e) for a, b, c, e in el),
                "Normal termination" in t)  # fmt: skip
    if out.exists():
        t = out.read_text()
        m = re.search(r"TOTAL RUN TIME:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+(\d+) seconds\s+(\d+) msec", t)
        s = (86400 * int(m.group(1)) + 3600 * int(m.group(2)) + 60 * int(m.group(3)) + int(m.group(4))
             + int(m.group(5)) / 1000) if m else None  # fmt: skip
        return s, "ORCA TERMINATED NORMALLY" in t
    return None, False


def level_apt(lev, mol, x0):
    """APT (3N, 3) at the reference geometry and the frame mismatch (Bohr), or (None, reason)."""
    d = run_dir(lev, mol)
    if lev.startswith("wb97"):
        p = d / f"{mol}.hess"
        if not p.exists():
            return None, "no .hess"
        apt, x = orca_hess(p)
        # ORCA's analytic Freq translates the molecule; the APT of a neutral molecule doesn't change
        return apt, float(np.abs((x - x.mean(0)) - (x0 - x0.mean(0))).max())
    p = d / f"{mol}.fchk"
    if not p.exists():
        return None, "no .fchk"
    apt = fchk_any(p, "Dipole Derivatives").reshape(-1, 3)
    x = fchk_any(p, "Current cartesian coordinates").reshape(-1, 3)
    return apt, float(np.abs(x - x0).max())


def intens(apt, disp):
    return KMMOL * np.sum((apt.T @ disp) ** 2, 0)


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def fmt(v, nd=1):
    return "–" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{nd}f}"


def hms(s):
    if s is None:
        return "–"
    s = float(s)
    return f"{s:.0f} s" if s < 120 else f"{s / 60:.1f} min" if s < 7200 else f"{s / 3600:.2f} h"


def main():
    OUT.mkdir(exist_ok=True)
    bz = np.load(ROOT / "bundle.npz")
    meta = json.loads((ROOT / "bundle.json").read_text())
    mols = meta["mols"]
    out = ["# B7 tables (generated by analyse.py)\n"]
    w = out.append
    checks, missing, frames, bands, cells, apts_all = [], [], [], [], [], {}

    for mol in mols:
        x0, disp, freqs = bz[f"{mol}/x"], bz[f"{mol}/disp"], bz[f"{mol}/freqs"]
        apts = {REF: bz[f"{mol}/apt_{REF}"]}
        for lev in LEVELS:
            a, dx = level_apt(lev, mol, x0)
            if a is None:
                missing.append(f"{lev}/{mol} ({dx})")
                continue
            frames.append((lev, mol, dx))
            apts[lev] = a
        for d in DIP:
            apts[d] = bz[f"{mol}/apt_{d}"]
        apts_all[mol] = apts
        ints = {k: intens(a, disp) for k, a in apts.items()}
        # B3LYP/6-31G(d,p) harmonic as B4: Gaussian's log values by log frequency
        lf, li = bz[f"{mol}/log_freq"], bz[f"{mol}/log_int"]
        o, gi = np.argsort(freqs), li[np.argsort(lf)]
        checks.append(("6-31G(d,p) modes reproduce the log (max rel. intensity diff, max freq diff cm-1)",
                       mol, float(np.max(np.abs(gi - ints[REF][o]) / np.maximum(gi, 1))),
                       float(np.max(np.abs(np.sort(lf) - freqs[o])))))  # fmt: skip
        if mol == "methanol":
            for lev in LEVELS:
                if lev in apts:
                    b = bz[f"methanol/b3_apt_{lev}"]
                    checks.append(("methanol APT = B3 APT (max abs diff, e)", lev,
                                   float(np.abs(apts[lev] - b).max()),
                                   float(np.max(np.abs(intens(b, disp) - ints[lev])))))  # fmt: skip
        # the 6-31G(d,p) timing rerun must reproduce the campaign APT
        p = ROOT / "g_b3lyp_631gdp_harm" / mol / f"{mol}.fchk"
        if p.exists():
            a = fchk_any(p, "Dipole Derivatives").reshape(-1, 3)
            checks.append(("6-31G(d,p) harmonic rerun APT vs campaign (rel. Frobenius)", mol, rel(a, apts[REF]), None))
        for k in range(len(freqs)):
            bands.append(dict(mol=mol, mode=k + 1, freq=float(freqs[k]), **{c: float(v[k]) for c, v in ints.items()}))
        for lab, lo, hi in meta["regions"][mol]:
            if "(all)" in lab:
                continue
            e = meta["exp"][mol][lab]
            sel = (freqs >= lo) & (freqs <= hi)
            vals = {c: float(v[sel].sum()) for c, v in ints.items()}
            vals[REF] = float(sum(i for f, i in zip(lf, li) if lo <= f <= hi))
            vals["n_bands"] = int(sel.sum())
            cells.append(dict(mol=mol, region=lab, exp=e, **vals))
            b4 = meta["b4"].get(mol, {}).get(lab)
            if b4:
                pairs = [("experiment", e), ("B3LYP harm", vals[REF]),
                         ("POLAR-1 harm (B3LYP modes)", vals["mace_polar1"]), ("MDP harm (B3LYP modes)", vals["mace_mdp"]),
                         ("MACE4IR harm (B3LYP modes)", vals["mace_ml"]), ("POLAR-1 (field) harm", vals["polar1_field"])]  # fmt: skip
                checks.append(("B4 table reproduced (max abs diff, km/mol; table has 1 decimal)", f"{mol} {lab}",
                               max(abs(b4[k] - v) for k, v in pairs if k in b4), None))  # fmt: skip

    have = [lev for lev in LEVELS if all(lev in apts_all[m] for m in mols)]
    rows = [r for r in ROWS if r in (REF, *have, *DIP)]
    w(f"Molecules: {', '.join(mols)}. Levels complete for all molecules: {', '.join(LEVELS[v] for v in have) or 'none'}.\n")
    if missing:
        w("**Missing outputs:** " + ", ".join(missing) + "\n")

    # ---------------- checks
    w("## 0. Checks\n")
    w("Frame mismatch, max |x_level - x_B3LYP/6-31G(d,p)| (Bohr; ORCA centred): " + ", ".join(
        f"{LEVELS[lev]} {max(d for lv, _m, d in frames if lv == lev):.1e}" for lev in LEVELS
        if any(lv == lev for lv, _m, _d in frames)) + "\n")  # fmt: skip
    w("| check | item | value | second value |")
    w("|---|---|---|---|")
    for c, item, v, v2 in checks:
        w(f"| {c} | {item} | {v:.2e} | {'' if v2 is None else f'{v2:.2e}'} |")

    # ---------------- item 2: against experiment
    ncell = len(cells)
    w(f"\n## 1. Against experiment: summary over the {ncell} molecule × region values\n")
    w("Harmonic fundamentals, region sums by B3LYP/6-31G(d,p) harmonic frequency; experiment NIST "
      "(B4). Within ±25 % as a count.\n")  # fmt: skip
    w("| method | median calc/exp | mean abs(log10(calc/exp)) | within ±25 % | mean log10(calc/exp) |")
    w("|---|---|---|---|---|")
    summ = {}
    for r in rows:
        q = np.array([c[r] / c["exp"] for c in cells])
        summ[r] = q
        w(f"| {LAB[r]} | {np.median(q):.2f} | {np.mean(np.abs(np.log10(q))):.3f} | "
          f"{int(np.sum(np.abs(q - 1) < 0.25))}/{len(q)} | {np.mean(np.log10(q)):+.3f} |")  # fmt: skip

    groups = {"fingerprint": lambda c: "fingerprint" in c["region"], "C–H": lambda c: "(C-H)" in c["region"],
              "O–H": lambda c: "(O-H)" in c["region"], "SO₂ ν₁/ν₃": lambda c: "nu" in c["region"]}  # fmt: skip
    w("\nBy region type: median calc/exp (mean abs log10)\n")
    w("| method | " + " | ".join(f"{g} ({sum(1 for c in cells if f(c))})" for g, f in groups.items()) + " |")
    w("|---" * (len(groups) + 1) + "|")
    for r in rows:
        cc = []
        for f in groups.values():
            q = np.array([c[r] / c["exp"] for c in cells if f(c)])
            cc.append(f"{np.median(q):.2f} ({np.mean(np.abs(np.log10(q))):.3f})" if len(q) else "–")
        w(f"| {LAB[r]} | " + " | ".join(cc) + " |")
    w("\nWithout the O–H region (the reference's 1.8–2.5× deficit there): mean abs log10 / within ±25 %\n")
    w("| method | values | mean abs(log10) | within ±25 % |")
    w("|---|---|---|---|")
    for r in rows:
        q = np.array([c[r] / c["exp"] for c in cells if "(O-H)" not in c["region"]])
        w(f"| {LAB[r]} | {len(q)} | {np.mean(np.abs(np.log10(q))):.3f} | {int(np.sum(np.abs(q - 1) < 0.25))}/{len(q)} |")

    w("\n### Per region (km/mol)\n")
    w("| molecule | region | bands | experiment | " + " | ".join(LAB[r] for r in rows) + " |")
    w("|---" * (4 + len(rows)) + "|")
    for c in cells:
        w(f"| {c['mol']} | {c['region']} | {c['n_bands']} | {c['exp']:.1f} | " + " | ".join(fmt(c[r]) for r in rows) + " |")
    w("\n### Per region, calc/exp\n")
    w("| molecule | region | " + " | ".join(LAB[r] for r in rows) + " |")
    w("|---" * (2 + len(rows)) + "|")
    for c in cells:
        w(f"| {c['mol']} | {c['region']} | " + " | ".join(f"{c[r] / c['exp']:.2f}" for r in rows) + " |")

    if have:
        w("\n### Can the four levels be told apart? Spread among the levels per value\n")
        w("max/min of the four levels' region values, against the distance of their mean from experiment.\n")
        w("| molecule | region | levels max/min | levels mean / exp | reference / exp |")
        w("|---|---|---|---|---|")
        spread = []
        for c in cells:
            v = np.array([c[lev] for lev in have])
            spread.append(v.max() / v.min())
            w(f"| {c['mol']} | {c['region']} | {v.max() / v.min():.3f} | {v.mean() / c['exp']:.2f} | {c[REF] / c['exp']:.2f} |")
        w(f"\nMedian max/min over the {ncell} values: {np.median(spread):.3f} (range {min(spread):.3f}–{max(spread):.3f}).")
        # paired comparison between levels: how often is level A closer to experiment than B?
        w("\nPairwise: number of the values where the row level is closer to experiment than the column level "
          "(abs log10), and the mean difference of abs log10 (row − column).\n")  # fmt: skip
        w("| | " + " | ".join(LEVELS[b] for b in have) + " |")
        w("|---" * (len(have) + 1) + "|")
        for a in have:
            la = np.abs(np.log10(summ[a]))
            cc = []
            for b in have:
                lb = np.abs(np.log10(summ[b]))
                cc.append("–" if a == b else f"{int(np.sum(la < lb))}/{ncell} ({np.mean(la - lb):+.3f})")
            w(f"| {LEVELS[a]} | " + " | ".join(cc) + " |")
        w("\nElement-wise APT, relative Frobenius difference between levels, mean over the molecules "
          "(B3 reported 0.7–4.5 % among the four levels on its 8 molecules):\n")  # fmt: skip
        w("| | " + " | ".join(LAB[b] for b in [REF, *have]) + " |")
        w("|---" * (len(have) + 2) + "|")
        for a in [REF, *have, *DIP]:
            w(f"| {LAB[a]} | " + " | ".join(f"{np.mean([rel(apts_all[m][a], apts_all[m][b]) for m in mols]):.3f}"
                                             for b in [REF, *have]) + " |")  # fmt: skip

    # ---------------- item 3: models against levels on these molecules
    sel = [b for b in bands if b[REF] > 1]
    w(f"\n## 2. Each dipole model against each level on these {len(mols)} molecules\n")
    w("MAE (km/mol) / median abs(log10 ratio), bands with B3LYP/6-31G(d,p) harmonic > 1 km/mol, on the "
      "6-31G(d,p) modes at its geometry (as B3). Training level in bold: MACE-POLAR-1 ωB97M-V/def2-TZVPD; "
      "MACE-MDP ωB97M-D3(BJ)/def2-TZVPPD; MACE4IR PBE0/FHI-aims tier-2, stand-in PBE0/def2-TZVP.\n")  # fmt: skip

    def mae_table(bs, title):
        w(f"\n### {title}\n")
        cols = [REF, *have]
        w("| dipole model | " + " | ".join(LAB[c] for c in cols) + " |")
        w("|---" * (len(cols) + 1) + "|")
        for m in DIP:
            cc = []
            for c in cols:
                d = np.mean([abs(b[m] - b[c]) for b in bs])
                lr = np.median([abs(np.log10(max(b[m], 1e-6) / b[c])) for b in bs])
                s = f"{d:.1f} / {lr:.3f}"
                cc.append(f"**{s}**" if TRAINED[m] == c else s)
            w(f"| {DIP[m]} | " + " | ".join(cc) + " |")
        if have:
            w("\nLevels against the reference, same metric: " + "; ".join(
                f"{LEVELS[c]} {np.mean([abs(b[c] - b[REF]) for b in bs]):.1f} / "
                f"{np.median([abs(np.log10(b[c] / b[REF])) for b in bs]):.3f}" for c in have))  # fmt: skip

    mae_table(sel, f"all bands ({len(sel)})")
    mae_table([b for b in sel if b["freq"] > 2700], f"X–H stretches > 2700 cm⁻¹ ({sum(b['freq'] > 2700 for b in sel)})")
    mae_table([b for b in sel if b["freq"] <= 2700], f"≤ 2700 cm⁻¹ ({sum(b['freq'] <= 2700 for b in sel)})")

    # ---------------- item 4: cost
    w("\n## 3. Cost\n")
    w("DFT on rune03 (Xeon Gold 6154), 8 cores each (Slurm wall from start/end epoch; program elapsed "
      "in brackets). B3LYP/6-31G(d,p): harmonic freq at the same geometry (this follow-up). Methanol's "
      "four levels are the B3 runs. ωB97M-V is ORCA NumFreq (numerical Hessian, 6N gradients), not "
      "representative of an analytic implementation; the APT needs only that, or a field-derivative "
      "run. ML: campaign VPT2 runs (local i7-6800K + RTX 2070 SUPER, Gaussian with 2 cores). One "
      "external call computes energy, forces, Hessian, dipole and APT at one geometry, i.e. one "
      "harmonic frequency + intensity calculation. Its wall time is printed per call in the batch "
      "logs (ml_alcohols.log, ml_inorganic.log; prep_bundle.call_times). 'ML call' = median over the "
      "runs of each run's median call (the first call of a run is listed apart; it includes warm-up); "
      "'runtime/calls' = whole run / number of calls, for comparison.\n")  # fmt: skip
    cost_rows = []
    hdr = ["molecule", "atoms", "B3LYP/6-31G(d,p) harm", *[LEVELS[v] + (" (NumFreq)" if v in NUMERICAL else "") for v in LEVELS],
           "ML call (median; min–max over runs)", "ML first call (median)", "ML runtime/calls (median)",
           "ML runs", "ML VPT2 run (median)"]  # fmt: skip
    w("| " + " | ".join(hdr) + " |")
    w("|---" * len(hdr) + "|")
    for mol in mols:
        n = len(bz[f"{mol}/m"])
        d0 = ROOT / "g_b3lyp_631gdp_harm" / mol
        t_ref, ok_ref = program_time(d0, mol)
        cc = [f"{hms(wall(d0))} ({hms(t_ref)}{'' if ok_ref or t_ref is None else ', not normal'})"]
        crow = dict(mol=mol, atoms=n, b3lyp_631gdp_wall_s=wall(d0), b3lyp_631gdp_prog_s=t_ref)
        for lev in LEVELS:
            d = run_dir(lev, mol)
            t, ok = program_time(d, mol)
            wl = wall(d)
            cc.append(f"{hms(wl)} ({hms(t)}{'' if ok or t is None else ', not normal'})")
            crow[f"{lev}_wall_s"], crow[f"{lev}_prog_s"] = wl, t
        runs = meta["ml_runs"][mol]
        call = np.array([np.median(r["calls_s"][1:]) for r in runs])
        first = np.array([r["calls_s"][0] for r in runs])
        per = np.array([r["runtime_s"] / r["n_calls"] for r in runs])
        tot = np.array([r["runtime_s"] for r in runs])
        cc += [f"{np.median(call):.2f} s ({call.min():.2f}–{call.max():.2f})", f"{np.median(first):.2f} s",
               f"{np.median(per):.2f} s", f"{len(runs)}", hms(np.median(tot))]  # fmt: skip
        crow.update(ml_call_median_s=float(np.median(call)), ml_call_min_s=float(call.min()),
                    ml_call_max_s=float(call.max()), ml_first_call_median_s=float(np.median(first)),
                    ml_runtime_per_call_median_s=float(np.median(per)), ml_runs=len(runs),
                    ml_vpt2_median_s=float(np.median(tot)))  # fmt: skip
        cost_rows.append(crow)
        w(f"| {mol} | {n} | " + " | ".join(cc) + " |")
    combos = [(e, d) for e in ("mace_anicc", "mace_omol", "mace_off", "mace_polar")
              for d in ("mace_polar1", "mace_mdp", "mace_ml")]  # fmt: skip
    w("\nML call by energy × dipole model (median call of the run, s; '–' = no run, MACE-ANI-CC has no S):\n")
    w("| energy model | dipole model | " + " | ".join(mols) + " |")
    w("|---" * (2 + len(mols)) + "|")
    for e, d in combos:
        cc = []
        for mol in mols:
            r = [r for r in meta["ml_runs"][mol] if r["energy"] == e and r["dipole"] == d]
            cc.append(f"{np.median(r[0]['calls_s'][1:]):.2f}" if r else "–")
        w(f"| {e} | {d} | " + " | ".join(cc) + " |")
    w("\nCampaign B3LYP/6-31G(d,p) opt + VPT2 (4 cores), elapsed per link (s): " + "; ".join(
        f"{m}: {' + '.join(f'{v:.0f}' for v in meta['campaign_b3lyp'][m]['elapsed_s'])}" for m in mols))  # fmt: skip
    tot = 0.0
    for lev_dir in ["g_b3lyp_631gdp_harm", *DIRS.values()]:
        for mol in mols:
            d = ROOT / lev_dir / mol
            if (d / "start_epoch.txt").exists() and wall(d):
                tot += wall(d)
    w(f"\nTotal wall time of this follow-up's jobs: {tot / 3600:.2f} h (sum over jobs, 8 cores each; "
      f"{tot * 8 / 3600:.1f} core-hours).")

    # ---------------- files
    keys = list(cells[0])
    with (OUT / "b7_cells.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(cells)
    keys = sorted({k for b in bands for k in b}, key=lambda k: (k not in ("mol", "mode", "freq"), k))
    with (OUT / "b7_bands.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(bands)
    with (OUT / "b7_cost.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(cost_rows[0]))
        wr.writeheader()
        wr.writerows(cost_rows)
    (OUT / "b7_tables.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))
    try:
        figure(cells, rows)
    except Exception as e:  # the figure is a convenience; the tables are the result
        print("figure skipped:", e)


def figure(cells, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5.0))
    xs = np.arange(len(cells))
    mk = {REF: ("k", "s"), "b3lyp_avtz": ("0.35", "o"), "pbe0_tzvp": ("0.5", "^"), "wb97mv_tzvpd": ("0.65", "v"),
          "wb97md3_tzvppd": ("0.8", "D"), "mace_polar1": ("#9467bd", "o"), "mace_mdp": ("#2ca02c", "o"),
          "mace_ml": ("#ff7f0e", "o"), "polar1_field": ("#9467bd", "x")}  # fmt: skip
    for j, r in enumerate(rows):
        q = np.array([c[r] / c["exp"] for c in cells])
        col, m = mk[r]
        ax.scatter(xs + (j - len(rows) / 2) * 0.07, q, color=col, marker=m, s=18, label=LAB[r],
                   facecolors="none" if r in LEVELS else None)  # fmt: skip
    ax.axhline(1, color="k", lw=0.6)
    ax.axhspan(0.75, 1.25, color="0.92", zorder=0)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{c['mol']}\n{c['region'].split('(')[-1].rstrip(')')}" for c in cells], rotation=90, fontsize=7)
    ax.set_ylabel("calc / exp")
    fig.legend(*ax.get_legend_handles_labels(), fontsize=6, ncol=5, loc="upper center", frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(OUT / "b7_calc_exp.png", dpi=150)


if __name__ == "__main__":
    main()
