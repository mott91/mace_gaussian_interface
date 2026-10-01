"""Did an ML model end up in the same conformer as the B3LYP reference?

Every energy model re-optimizes on its own surface. For flexible molecules (acid and alcohol
chains, glycine, glucose) a model can slide into a different conformer; its spectrum is then
not comparable to the B3LYP one, whatever the model's accuracy. This compares the final ML
geometry with the final B3LYP geometry (both from the .fchk, same atom order):

- **RMSD** after optimal superposition (Kabsch), all atoms and heavy atoms. Same conformer:
  typically < 0.05 A; a different conformer usually > 0.2 A.
- **Torsions** around every rotatable bond, the decisive test: a single OH flip barely moves
  the RMSD but changes a dihedral by ~180 deg. Rotatable = acyclic single bond between two
  atoms that both carry another neighbour. Symmetric rotors (CH3, NH3, CF3: three identical
  terminal atoms) are skipped, because a 120 deg turn gives the same conformer.

``same_conformer`` is False when any torsion differs by more than ``TORSION_TOL_DEG``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..gaussian.fchk import parse_fchk_section

BOHR_TO_ANGSTROM = 0.529177210903
TORSION_TOL_DEG = 30.0


@dataclass
class ConformerResult:
    rmsd_all: float
    rmsd_heavy: float
    n_torsions: int
    max_torsion_dev: float  # deg, 0 when the molecule has no rotatable bond
    worst_torsion: tuple[int, int, int, int] | None  # 0-based atom indices
    same_conformer: bool


def read_fchk_geometry(fchk: Path) -> tuple[np.ndarray, np.ndarray]:
    """Atomic numbers and final Cartesian coordinates (Angstrom) of a .fchk."""
    content = Path(fchk).read_text()
    z = parse_fchk_section(content, "Atomic numbers", "I").astype(int)
    xyz = parse_fchk_section(content, "Current cartesian coordinates", "R").reshape(-1, 3)
    return z, xyz * BOHR_TO_ANGSTROM


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD between two point sets after optimal translation and rotation."""
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1))))


_COVALENT = {1: 0.31, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 14: 1.11, 15: 1.07,
             16: 1.05, 17: 1.02, 35: 1.20}  # fmt: skip


def bonds(z: np.ndarray, xyz: np.ndarray, scale: float = 1.25) -> list[tuple[int, int]]:
    """Covalent bonds from interatomic distances (sum of covalent radii x ``scale``)."""
    r = np.array([_COVALENT.get(int(k), 1.0) for k in z])
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    i, j = np.where(np.triu(d < scale * (r[:, None] + r[None]), k=1))
    return list(zip(i.tolist(), j.tolist()))


def _in_ring(a: int, b: int, nbrs: dict[int, set[int]]) -> bool:
    """Is bond a-b part of a ring (is b reachable from a without using a-b)?"""
    seen, stack = {a}, [n for n in nbrs[a] if n != b]
    while stack:
        n = stack.pop()
        if n == b:
            return True
        if n not in seen:
            seen.add(n)
            stack.extend(nbrs[n] - seen)
    return False


def _symmetric_rotor(center: int, other: int, z: np.ndarray, nbrs: dict[int, set[int]]) -> bool:
    ends = nbrs[center] - {other}
    return (
        len(ends) == 3
        and len({int(z[k]) for k in ends}) == 1
        and all(len(nbrs[k]) == 1 for k in ends)
    )


def rotatable_torsions(z: np.ndarray, xyz: np.ndarray) -> list[tuple[int, int, int, int]]:
    """One dihedral (i, j, k, l) per rotatable bond j-k; i, l = heaviest other neighbour."""
    nbrs: dict[int, set[int]] = {k: set() for k in range(len(z))}
    for a, b in bonds(z, xyz):
        nbrs[a].add(b)
        nbrs[b].add(a)
    out = []
    for j, k in bonds(z, xyz):
        if len(nbrs[j]) < 2 or len(nbrs[k]) < 2 or _in_ring(j, k, nbrs):
            continue
        if _symmetric_rotor(j, k, z, nbrs) or _symmetric_rotor(k, j, z, nbrs):
            continue
        i = max(nbrs[j] - {k}, key=lambda n: (z[n], -n))
        lt = max(nbrs[k] - {j}, key=lambda n: (z[n], -n))
        out.append((i, j, k, lt))
    return out


def dihedral(xyz: np.ndarray, i: int, j: int, k: int, lt: int) -> float:
    b0, b1, b2 = xyz[i] - xyz[j], xyz[k] - xyz[j], xyz[lt] - xyz[k]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def compare(z: np.ndarray, ref: np.ndarray, other: np.ndarray) -> ConformerResult:
    """Compare ``other`` (ML) with ``ref`` (B3LYP); same atoms in the same order."""
    heavy = z > 1
    torsions = rotatable_torsions(z, ref)
    devs = []
    for t in torsions:
        d = dihedral(other, *t) - dihedral(ref, *t)
        devs.append(abs((d + 180.0) % 360.0 - 180.0))
    worst = int(np.argmax(devs)) if devs else None
    max_dev = float(devs[worst]) if devs else 0.0
    return ConformerResult(
        rmsd_all=kabsch_rmsd(ref, other),
        rmsd_heavy=kabsch_rmsd(ref[heavy], other[heavy]) if heavy.sum() >= 3 else 0.0,
        n_torsions=len(torsions),
        max_torsion_dev=max_dev,
        worst_torsion=torsions[worst] if worst is not None else None,
        same_conformer=max_dev <= TORSION_TOL_DEG,
    )


def check_run(dft_fchk: Path, ml_fchk: Path) -> ConformerResult:
    z_ref, ref = read_fchk_geometry(dft_fchk)
    z_ml, ml = read_fchk_geometry(ml_fchk)
    if not np.array_equal(z_ref, z_ml):
        raise ValueError(f"atom order differs between {dft_fchk} and {ml_fchk}")
    return compare(z_ref, ref, ml)
