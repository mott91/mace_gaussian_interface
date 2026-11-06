"""Follow-up 2, add-on: harmonic T1u intensities of C60 with all three dipole models.

MACE-MDP was never run on C60 (it does not fit on the 8 GB GPU there). Here every dipole model is
evaluated on the CPU in float64 with its unchanged campaign calculator (APT by autograd, units e):
- MACE-POLAR-1: calculators/mace_polar1.py, device="cpu";
- MACE-MDP:     calculators/mace_mdp.py (get_dielectric_derivatives), device="cpu";
- MACE4IR:      calculators/mace_loader.MACEDipoleCalculator (what mace_ml delegates to), cpu.

Geometries and harmonic modes (fchk "Current cartesian coordinates", "Cartesian Force Constants",
"Real atomic weights"; modes as in b1_b2_analyse.normal_modes):
- the three ML minima of the C60 VPT2 runs (ANI-cc, OMOL, OFF; each run's reference geometry);
- the B3LYP/6-31G(d,p) minimum (campaigns/c60/b3lyp_harmonic).
T1u identification without any model's intensity: at B3LYP the four T1u triplets are the four
bright triply degenerate levels (B3LYP APT). On an ML geometry, after Kabsch alignment onto the
B3LYP geometry (same atom order), each B3LYP T1u triplet picks the three ML modes with the largest
overlap with its subspace; intensities are summed over the three components.

Sanity check: POLAR-1 at the OMOL geometry must reproduce the APT stored in that Gaussian run's
fchk and its harmonic T1u intensities (75.0 / 70.3 / 0.1 / 33.7 km/mol, RESULTS_NOTES §12).
Writes c60_tables.md.
"""

from __future__ import annotations

import logging
import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")
import torch  # noqa: E402
from ase import Atoms  # noqa: E402
from b1_b2_analyse import KMMOL, normal_modes  # noqa: E402

from mace_gaussian.analysis.mode_matching import kabsch_rotation  # noqa: E402
from mace_gaussian.calculators.mace_loader import MACEDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.mace_mdp import MACEMDPDipoleCalculator  # noqa: E402
from mace_gaussian.calculators.mace_ml import DEFAULT_MACE_DIPOLE_MODEL  # noqa: E402
from mace_gaussian.calculators.mace_polar1 import MACEPolar1DipoleCalculator  # noqa: E402

C60 = ROOT / "campaigns/c60"
BOHR = 0.529177210903
GEOMS = {
    "B3LYP": C60 / "b3lyp_harmonic/c60_b3lyp_harm.fchk",
    "ANI-cc": C60 / "comparison_results/c60/mace_anicc_mace_polar1/gaussian_freq.fchk",
    "OMOL": C60 / "comparison_results/c60/mace_omol_mace_polar1/gaussian_freq.fchk",
    "OFF": C60 / "comparison_results/c60/mace_off_mace_polar1/gaussian_freq.fchk",
}
B3LYP_T1U = (85.1, 33.1, 27.1, 35.4)


def fchk(path: Path) -> dict:
    t = path.read_text()

    def arr(key):
        i = t.index(key + " ")
        n = int(t[i : t.index("\n", i)].split("N=")[1])
        return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])

    n = len(arr("Atomic numbers"))
    h = np.zeros((3 * n, 3 * n))
    h[np.tril_indices(3 * n)] = arr("Cartesian Force Constants")
    h = h + h.T - np.diag(np.diag(h))
    out = dict(z=arr("Atomic numbers").astype(int), x=arr("Current cartesian coordinates").reshape(-1, 3),
               m=arr("Real atomic weights"), hess=h)  # fmt: skip
    try:
        out["apt"] = arr("Dipole Derivatives").reshape(3 * n, 3)
    except ValueError:
        out["apt"] = None
    return out


def t1u_sets(freqs, disp, apt_b3lyp):
    """Indices of the four T1u triplets at B3LYP: the four brightest triply degenerate levels."""
    order = np.argsort(freqs)
    groups, cur = [], [order[0]]
    for i in order[1:]:
        if freqs[i] - freqs[cur[-1]] < 1.0:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    inten = KMMOL * np.sum((apt_b3lyp.T @ disp) ** 2, 0)
    trip = [g for g in groups if len(g) == 3]
    bright = sorted(trip, key=lambda g: -inten[g].sum())[:4]
    return sorted(bright, key=lambda g: freqs[g].mean())


def unit(v):
    return v / np.linalg.norm(v, axis=0)


def main():
    t0 = time.time()
    torch.set_num_threads(6)
    calcs = {"POLAR-1": MACEPolar1DipoleCalculator(device="cpu"), "MDP": MACEMDPDipoleCalculator(device="cpu"),
             "MACE4IR": MACEDipoleCalculator(DEFAULT_MACE_DIPOLE_MODEL, device="cpu")}  # fmt: skip
    f = {k: fchk(p) for k, p in GEOMS.items()}
    modes = {k: normal_modes(v["x"], v["m"], v["hess"]) for k, v in f.items()}
    fb, db = modes["B3LYP"]
    ref_sets = t1u_sets(fb, db, f["B3LYP"]["apt"])
    out = ["# C60: harmonic T1u intensities with all three dipole models (c60_dipoles.py)\n",
           "B3LYP T1u levels (cm-1): " + ", ".join(f"{fb[s].mean():.1f}" for s in ref_sets) + "\n"]  # fmt: skip
    w = out.append
    results, timing = {}, {}
    for geo, v in f.items():
        freqs, disp = modes[geo]
        if geo == "B3LYP":
            sets = ref_sets
            ovl = [1.0] * 4
        else:  # pair each B3LYP T1u triplet with the 3 ML modes of largest subspace overlap
            r = kabsch_rotation(v["x"] - v["x"].mean(0), f["B3LYP"]["x"] - f["B3LYP"]["x"].mean(0))
            n = len(v["m"])
            mw_ml = unit((disp.reshape(n, 3, -1).transpose(0, 2, 1) @ r).transpose(0, 2, 1).reshape(3 * n, -1))
            mw_b = unit(db)
            sets, ovl = [], []
            for s in ref_sets:
                p = np.sum((mw_b[:, s].T @ mw_ml) ** 2, axis=0)  # overlap of each ML mode with the triplet
                pick = list(np.argsort(-p)[:3])
                sets.append(pick)
                ovl.append(float(p[pick].sum() / 3))
        atoms = Atoms(numbers=v["z"], positions=v["x"] * BOHR)
        atoms.info.update(charge=0.0, spin=1.0)
        for name, calc in calcs.items():
            t1 = time.time()
            apt = np.asarray(calc.calculate_dipole_derivatives(atoms))
            timing[(geo, name)] = time.time() - t1
            inten = KMMOL * np.sum((apt.T @ disp) ** 2, 0)
            results[(geo, name)] = [float(inten[s].sum()) for s in sets]
            if geo == "B3LYP":  # element-wise APT error against B3LYP's own APT, same geometry
                ref = f["B3LYP"]["apt"]
                results[("apt_rel", name)] = float(np.linalg.norm(apt - ref) / np.linalg.norm(ref))
                results[("apt_max", name)] = float(np.abs(apt).max())
            if geo == "OMOL" and name == "POLAR-1":
                stored = v["apt"]
                results["check_apt"] = float(np.abs(apt - stored).max()) if stored is not None else None
                results["check_int"] = [float((KMMOL * np.sum((stored.T @ disp) ** 2, 0))[s].sum()) for s in sets]
        results[(geo, "freqs")] = [float(freqs[s].mean()) for s in sets]
        results[(geo, "overlap")] = ovl

    w("## Sanity check (POLAR-1 at the OMOL geometry)\n")
    w(f"- max |APT(now, CPU) - APT(fchk of the Gaussian run)| = {results['check_apt']:.1e} e")
    w("- T1u from the stored APT: " + " / ".join(f"{x:.1f}" for x in results["check_int"])
      + " km/mol; recomputed now: " + " / ".join(f"{x:.1f}" for x in results[("OMOL", "POLAR-1")])
      + " (RESULTS_NOTES §12: 75.0 / 70.3 / 0.1 / 33.7)\n")  # fmt: skip
    w("## T1u intensities (km/mol, sum of the three components) on each set of harmonic modes\n")
    w("| geometry and modes | T1u freq (cm-1) | mode overlap with B3LYP T1u | MDP | POLAR-1 | MACE4IR |")
    w("|---|---|---|---|---|---|")
    for geo in GEOMS:
        fr = " / ".join(f"{x:.0f}" for x in results[(geo, "freqs")])
        ov = " / ".join(f"{x:.3f}" for x in results[(geo, "overlap")])
        cells = [" / ".join(f"{x:.1f}" for x in results[(geo, n)]) for n in ("MDP", "POLAR-1", "MACE4IR")]
        w(f"| {geo} | {fr} | {ov} | " + " | ".join(cells) + " |")
    w(f"| B3LYP's own APT (reference) | | | {' / '.join(str(x) for x in B3LYP_T1U)} | | |")
    w("\nElement-wise APT at the B3LYP geometry, relative Frobenius error against B3LYP (largest "
      "element in e; B3LYP " + f"{np.abs(f['B3LYP']['apt']).max():.2f}" + "): " + ", ".join(
        f"{n} {results[('apt_rel', n)]:.2f} ({results[('apt_max', n)]:.2f})" for n in calcs))  # fmt: skip
    w("\nCPU time per APT evaluation (s): " + ", ".join(
        f"{n} {np.mean([timing[(g, n)] for g in GEOMS]):.1f}" for n in calcs) + f"; total {time.time() - t0:.0f} s")  # fmt: skip
    (HERE / "c60_tables.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
