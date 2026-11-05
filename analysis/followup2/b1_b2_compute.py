"""Follow-up 2, B1 + B2: dipole-model APTs at the B3LYP geometry (compute step).

B1: for every campaign molecule, the dipole and APT of one dipole model at the B3LYP/6-31G(d,p)
    minimum, read from the campaign fchk ("Current cartesian coordinates", Bohr). Computed with
    the pipeline's own routine (workflow.calculate_dipole_properties, as in every campaign call),
    so units and layout are the campaign's: dipole e*Bohr, APT (3N, 3) in e, row = Cartesian
    coordinate (atom-major), column = dipole component.
B2 (only for mace_polar1): the same quantities from the field response of the POLAR-1 energy,
    by central differences in a uniform field atoms.info["external_field"] (mace-torch 0.3.16):
      mu_field  = -[E(+D) - E(-D)] / 2D,   APT_field[i, b] = [F_i(+D e_b) - F_i(-D e_b)] / 2D,
    for D = 1e-4, 1e-3, 1e-2 (field in the model's unit, checked against the charge dipole).

Usage: python b1_b2_compute.py <dipole model>   (mace_ml | mace_mdp | mace_polar1)
Writes b1_apt_<model>.npz next to this script.
"""

from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")
from ase import Atoms  # noqa: E402

from mace_gaussian import workflow as wf  # noqa: E402
from mace_gaussian.calculators import dipole_factory  # noqa: E402

CAMP = ROOT / "campaigns/2026"
BOHR = 0.529177210903
FIELDS = (1e-4, 1e-3, 1e-2)


def fchk_arrays(mol: str) -> dict:
    t = (CAMP / mol / f"{mol}_freq_anharm.fchk").read_text()

    def arr(key):
        i = t.index(key)
        n = int(t[i : t.index("\n", i)].split("N=")[1])
        return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])

    return dict(z=arr("Atomic numbers").astype(int),
                x_bohr=arr("Current cartesian coordinates").reshape(-1, 3))  # fmt: skip


def molecules() -> list[str]:
    return sorted(p.name for p in CAMP.iterdir() if (p / f"{p.name}_freq_anharm.fchk").exists())


def field_response(calc, atoms: Atoms):
    """POLAR-1 energy/forces under a uniform field; returns {D: (mu_field, apt_field)}."""

    def ef(field):
        a = atoms.copy()
        a.info.update(atoms.info)
        a.info["external_field"] = list(map(float, field))
        calc.calculate(a, properties=["energy", "forces"])  # explicit: info is not in ASE's cache key
        return float(calc.results["energy"]), np.array(calc.results["forces"]), np.array(
            calc.results["dipole"]).reshape(-1)[:3]  # fmt: skip

    e0, f0, mu_charge = ef(np.zeros(3))
    out = {}
    for d in FIELDS:
        mu = np.zeros(3)
        apt = np.zeros((3 * len(atoms), 3))
        for b in range(3):
            fv = np.zeros(3)
            fv[b] = d
            ep, fp, _ = ef(fv)
            em, fm, _ = ef(-fv)
            mu[b] = -(ep - em) / (2 * d)
            apt[:, b] = ((fp - fm) / (2 * d)).reshape(-1)
        out[d] = (mu, apt)
    return mu_charge, out


def main():
    model = sys.argv[1]
    dcalc = dipole_factory.get_calculator(model)
    res = {}
    for mol in molecules():
        f = fchk_arrays(mol)
        atoms = Atoms(numbers=f["z"], positions=f["x_bohr"] * BOHR)
        atoms.info.update(charge=0.0, spin=1.0)
        mu, apt, _q, _pol = wf.calculate_dipole_properties(atoms, dcalc, 2, True)
        res[f"{mol}/mu"] = np.asarray(mu)  # e*Bohr
        res[f"{mol}/apt"] = np.asarray(apt)  # e, (3N, 3)
        msg = f"{mol:20s} |mu| {np.linalg.norm(mu) * 2.541746:.3f} D"
        if model == "mace_polar1":
            dcalc._ensure_calculator()
            mu_charge, fr = field_response(dcalc.calc, atoms)
            res[f"{mol}/mu_charge_eA"] = mu_charge
            for d, (muf, aptf) in fr.items():
                res[f"{mol}/mu_field_{d:g}"] = muf
                res[f"{mol}/apt_field_{d:g}"] = aptf
            mf = fr[1e-3][0]
            msg += (f"  charge-dipole {mu_charge.round(4)} e*A  field-dipole(1e-3) {mf.round(4)}  "
                    f"|APT_field - APT_charge| / |APT| = "
                    f"{np.linalg.norm(fr[1e-3][1] - apt) / np.linalg.norm(apt):.3f}")  # fmt: skip
        print(msg, flush=True)
    np.savez(HERE / f"b1_apt_{model}.npz", **res)


if __name__ == "__main__":
    main()
