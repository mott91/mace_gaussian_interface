"""Shared definitions for the surface scans along X-H stretches (Section 4.2).

Each method (B3LYP and the four MACE energy models) is scanned along its own normal mode from its
own optimized minimum, i.e. the geometry and mode its Gaussian VPT2 run of the 2026 campaign used:

    x(Q) = x0 + Q d,   d = L / sqrt(m)   (L: mass-weighted, normalized mode vector)

Q is the mass-weighted normal coordinate in amu^(1/2) Angstrom. Along this straight line the
curvature is exactly k(Q) = d^T H(x(Q)) d (attoJ amu^-1 Angstrom^-2, Gaussian's "k" units), so
dk/dQ = d3E/dQ3 and d2k/dQ2 = d4E/dQ4 are the cubic and quartic constants k(i,i,i), k(i,i,i,i)
that Gaussian prints for the VPT2 run (attoJ amu^-3/2 A^-3, attoJ amu^-2 A^-4).

Gaussian's VPT2 displacement h (default 0.010 Angstrom) is a Cartesian displacement along the unit
vector e = d / |d|, i.e. Q = h / |d| (checked against the B3LYP log in analyse.py).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAMP = ROOT / "campaigns" / "2026"
RUNS = CAMP / "comparison_results"
DATA = HERE / "data"

MODELS = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]  # fixed chapter order
METHODS = ["b3lyp", *MODELS]
# target modes: (key, molecule, B3LYP mode, 1-based in ascending frequency)
TARGETS = [
    ("methanol_oh", "methanol", 12),  # O-H stretch
    ("methanol_ch", "methanol", 11),  # C-H stretch that MACE-POLAR-1 places worst after VPT2
    ("glycine_nh", "glycine", 23),  # N-H stretch with the largest MACE-OMOL-0 error
]
STEPS = (0.005, 0.010, 0.020, 0.040)  # Gaussian displacement steps of Methods 3.7, Angstrom
DEFAULT_STEP = 0.010
FINE = np.round(np.arange(-0.050, 0.0501, 0.002), 4)  # target bond-length changes, Angstrom
WIDE = np.round(np.arange(-0.25, 0.4001, 0.01), 3)

HARTREE_AJ = 4.3597447222071  # attoJ
EV_AJ = 0.1602176634  # attoJ
BOHR_A = 0.529177210903
AMU_KG = 1.66053906660e-27
C_CM = 2.99792458e10
SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O"}


def k_to_wavenumber(k):
    """Curvature k (attoJ amu^-1 A^-2) -> harmonic wavenumber (cm-1); negative k -> negative."""
    k = np.asarray(k, float)
    w = np.sqrt(np.abs(k) * 1e-18 / (AMU_KG * 1e-20)) / (2 * np.pi * C_CM)
    return np.sign(k) * w


def fchk_arrays(path: Path, keys: list[str]) -> dict[str, np.ndarray]:
    """First occurrence of each array key in a formatted checkpoint."""
    t = path.read_text()
    out = {}
    for key in keys:
        i = t.index("\n" + key) + 1
        head = t[i : t.index("\n", i)]
        if "N=" not in head:  # scalar
            out[key] = np.array([float(head.split()[-1])])
            continue
        n = int(head.split("N=")[1])
        out[key] = np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])
    return out


def run_dir(mol: str, method: str) -> Path:
    if method == "b3lyp":
        return next(p for p in sorted((RUNS / mol).glob("b3lyp*")) if p.is_dir())
    return RUNS / mol / f"{method}_mace_ml"


def run_fchk(mol: str, method: str) -> Path:
    return next(iter(sorted(run_dir(mol, method).glob("*.fchk"))))


def run_log(mol: str, method: str) -> Path:
    return next(iter(sorted(run_dir(mol, method).glob("*.log"))))


@dataclass
class Mode:
    key: str
    mol: str
    method: str
    b3lyp_mode: int  # 1-based, ascending
    own_mode: int  # 1-based, ascending, in this method's own run
    z: np.ndarray  # atomic numbers
    x0: np.ndarray  # (N, 3) Angstrom, the run's optimized minimum
    m: np.ndarray  # atomic masses, amu
    d: np.ndarray  # (N, 3) Cartesian displacement per unit Q (amu^-1/2), +Q stretches X-H
    gauss_phase: int  # +1 if Gaussian's Vib-Mode vector stretches X-H, else -1
    h_atom: int
    x_atom: int
    freq: float  # harmonic wavenumber of the run, cm-1

    @property
    def dnorm(self) -> float:
        return float(np.linalg.norm(self.d))

    def geometry(self, q: float) -> np.ndarray:
        return self.x0 + q * self.d

    def bond(self, q) -> np.ndarray:
        q = np.atleast_1d(q)
        x = self.x0[None] + q[:, None, None] * self.d[None]
        r = np.linalg.norm(x[:, self.h_atom] - x[:, self.x_atom], axis=1)
        return r - np.linalg.norm(self.x0[self.h_atom] - self.x0[self.x_atom])

    @property
    def dr_dq(self) -> float:
        u = self.x0[self.h_atom] - self.x0[self.x_atom]
        u = u / np.linalg.norm(u)
        return float((self.d[self.h_atom] - self.d[self.x_atom]) @ u)

    @property
    def purity(self) -> dict:
        """How pure an X-H stretch: kinetic-energy share on the X-H pair (|L|^2 on H and X) and
        dr/dQ relative to an isolated X-H oscillator, 1/sqrt(mu)."""
        L2 = (self.d**2 * self.m[:, None]).sum(1)
        mu = self.m[self.h_atom] * self.m[self.x_atom] / (self.m[self.h_atom] + self.m[self.x_atom])
        on_x = [a for a in range(len(self.z)) if self.z[a] == 1
                and np.linalg.norm(self.x0[a] - self.x0[self.x_atom]) < 1.2]  # all H on the same X
        return {"share_XH": float(L2[self.h_atom] + L2[self.x_atom]), "share_H": float(L2[self.h_atom]),
                "share_XHn": float(L2[on_x].sum() + L2[self.x_atom]),
                "dr_dq_over_pure": self.dr_dq * float(np.sqrt(mu))}  # fmt: skip

    def q_for_bond(self, dr_targets) -> np.ndarray:
        """Q values whose bond-length change equals the targets (Newton on the exact geometry)."""
        q = np.asarray(dr_targets, float) / self.dr_dq
        for _ in range(20):
            eps = 1e-6
            f = self.bond(q) - dr_targets
            fp = (self.bond(q + eps) - self.bond(q - eps)) / (2 * eps)
            q = q - f / fp
        return q

    def q_for_step(self, h: float) -> float:
        """Gaussian displacement h (Cartesian Angstrom along d/|d|) as Q."""
        return h / self.dnorm


def b3lyp_band(mol: str, i: int, model: str) -> pd.Series:
    b = pd.read_csv(CAMP / "figures" / "supervisor" / "data" / "bands.csv")
    b = b[(b.molecule == mol) & (b.band_type == "fundamental") & (b.dipole_model == "mace_ml")
          & (b.energy_model == model) & (b.i == i)]  # fmt: skip
    assert len(b) == 1, (mol, i, model, len(b))
    return b.iloc[0]


def load_mode(key: str, method: str, mol: str | None = None, i: int | None = None) -> Mode:
    """The target mode `key` (TARGETS), or any B3LYP mode i of mol (key is then only a label)."""
    if mol is None:
        mol, i = next((m, i) for k, m, i in TARGETS if k == key)
    a = fchk_arrays(run_fchk(mol, method), ["Atomic numbers", "Current cartesian coordinates",
                                            "Real atomic weights", "Vib-E2", "Vib-Modes"])  # fmt: skip
    z = a["Atomic numbers"].astype(int)
    n = len(z)
    x0 = a["Current cartesian coordinates"].reshape(n, 3) * BOHR_A
    m = a["Real atomic weights"]
    nm = 3 * n - 6
    freqs = a["Vib-E2"][:nm]
    modes = a["Vib-Modes"].reshape(nm, n, 3)
    assert np.all(np.diff(freqs) > 0)
    # which of this run's modes is the target: the band table pairs B3LYP mode i with a MACE mode
    if method == "b3lyp":
        own = i
    else:
        f_ml = b3lyp_band(mol, i, method).ml_harmonic
        own = int(np.argmin(np.abs(freqs - f_ml))) + 1
        assert abs(freqs[own - 1] - f_ml) < 0.05, (key, method, freqs[own - 1], f_ml)
    v = modes[own - 1]  # Gaussian's Cartesian displacement vector
    L = v * np.sqrt(m)[:, None]
    L = L / np.linalg.norm(L)
    d = L / np.sqrt(m)[:, None]
    # the X-H bond: the hydrogen that moves most in the B3LYP mode (the same atom for every method;
    # in the antisymmetric NH2 stretch of glycine both hydrogens move equally) and its heavy atom
    if method == "b3lyp":
        h_atom = int(np.argmax(np.where(z == 1, (L**2).sum(1), -1)))
    else:
        h_atom = load_mode(key, "b3lyp", mol, i).h_atom
    dist = np.linalg.norm(x0 - x0[h_atom], axis=1)
    dist[h_atom] = np.inf
    dist[z == 1] = np.inf
    x_atom = int(np.argmin(dist))
    u = x0[h_atom] - x0[x_atom]
    phase = 1 if (v[h_atom] - v[x_atom]) @ u > 0 else -1
    if phase < 0:
        d = -d
    return Mode(key, mol, method, i, own, z, x0, m, d, phase, h_atom, x_atom, float(freqs[own - 1]))


def grid_points(mode: Mode) -> dict[str, np.ndarray]:
    """Q values of the three grids: fine (curvature), wide (energy), steps (+/- h, Gaussian)."""
    steps = np.array([s * mode.q_for_step(h) for h in STEPS for s in (-1, 1)])
    return {"fine": mode.q_for_bond(FINE), "wide": mode.q_for_bond(WIDE), "steps": steps}


def save(name: str, obj: dict) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj))


def load(name: str) -> dict:
    return json.loads((DATA / name).read_text())
