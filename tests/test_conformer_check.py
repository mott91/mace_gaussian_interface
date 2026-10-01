"""Conformer check: ML geometry vs B3LYP geometry (mace_gaussian.analysis.conformer_check)."""

import numpy as np
import pytest

from mace_gaussian.analysis.conformer_check import (
    compare,
    dihedral,
    kabsch_rmsd,
    rotatable_torsions,
)

# trans (Z) formic acid, planar: C, O(carbonyl), O(hydroxyl), H(C), H(O)
Z_HCOOH = np.array([6, 8, 8, 1, 1])
HCOOH = np.array(
    [
        [0.000, 0.000, 0.0],
        [0.660, 1.020, 0.0],
        [0.660, -1.130, 0.0],
        [-1.090, 0.000, 0.0],
        [1.610, -0.940, 0.0],
    ]
)
# ethanol: C(H3) C(H2) O H, then the five C-H hydrogens
Z_ETOH = np.array([6, 6, 8, 1, 1, 1, 1, 1, 1])
ETOH = np.array(
    [
        [-1.25, 0.20, 0.00],
        [0.00, -0.60, 0.00],
        [1.15, 0.25, 0.00],
        [1.95, -0.30, 0.00],  # O-H anti to C-C
        [-2.14, -0.43, 0.00],
        [-1.27, 0.84, 0.89],
        [-1.27, 0.84, -0.89],
        [0.03, -1.24, 0.89],
        [0.03, -1.24, -0.89],
    ]
)


def _rotate_about(xyz, axis_from, axis_to, atoms, deg):
    """Rotate ``atoms`` about the axis_from -> axis_to bond (Rodrigues)."""
    out = xyz.copy()
    k = xyz[axis_to] - xyz[axis_from]
    k /= np.linalg.norm(k)
    t = np.radians(deg)
    for a in atoms:
        v = xyz[a] - xyz[axis_to]
        v = v * np.cos(t) + np.cross(k, v) * np.sin(t) + k * np.dot(k, v) * (1 - np.cos(t))
        out[a] = xyz[axis_to] + v
    return out


def _random_rotation(seed=0):
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_rmsd_ignores_rotation_and_translation():
    moved = HCOOH @ _random_rotation().T + np.array([3.0, -1.0, 2.0])
    assert kabsch_rmsd(HCOOH, moved) == pytest.approx(0.0, abs=1e-9)


def test_same_conformer_after_rigid_motion():
    moved = ETOH @ _random_rotation(1).T + 5.0
    r = compare(Z_ETOH, ETOH, moved)
    assert r.same_conformer
    assert r.max_torsion_dev == pytest.approx(0.0, abs=1e-6)


def test_cooh_flip_is_flagged():
    cis = _rotate_about(HCOOH, 0, 2, [4], 180.0)  # H(O) around the C-O bond
    assert abs(dihedral(cis, 1, 0, 2, 4)) == pytest.approx(180.0, abs=1e-6)
    r = compare(Z_HCOOH, HCOOH, cis)
    assert not r.same_conformer
    assert r.max_torsion_dev == pytest.approx(180.0, abs=1e-6)
    assert r.worst_torsion == (1, 0, 2, 4)


def test_ethanol_gauche_is_flagged_methyl_turn_is_not():
    gauche = _rotate_about(ETOH, 1, 2, [3], 120.0)  # C-C-O-H anti -> gauche
    assert not compare(Z_ETOH, ETOH, gauche).same_conformer
    methyl = _rotate_about(ETOH, 1, 0, [4, 5, 6], 120.0)  # CH3 turn: same conformer
    r = compare(Z_ETOH, ETOH, methyl)
    assert r.same_conformer
    assert r.max_torsion_dev == pytest.approx(0.0, abs=1e-6)


def test_methyl_rotor_bond_is_not_a_torsion():
    tors = rotatable_torsions(Z_ETOH, ETOH)
    assert tors == [(0, 1, 2, 3)]  # only C-C-O-H; the CH3 side of C-C is a symmetric rotor
