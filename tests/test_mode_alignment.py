"""Mode overlaps must not depend on how Gaussian happened to orient each job.

Regression for the 2026-09-21 finding: ML and DFT checkpoints of the same molecule
differ by a rotation (each job is oriented independently, and every ML model
re-optimizes on its own surface), which scrambled the eigenvector assignment.
"""

from __future__ import annotations

import numpy as np

from mace_gaussian.analysis.mode_matching import (
    align_modes_to_reference,
    compute_mode_overlap,
    kabsch_rotation,
    match_modes,
)


def _rotation(angle: float, axis: int = 1) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    r = np.eye(3)
    a, b = [i for i in range(3) if i != axis]
    r[a, a] = r[b, b] = c
    r[a, b], r[b, a] = -s, s
    return r


def _fake_molecule(seed: int = 0):
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(4, 3))
    modes = rng.normal(size=(6, 4, 3))
    return coords, modes


def test_kabsch_recovers_a_known_rotation():
    coords, _ = _fake_molecule()
    r = _rotation(np.pi / 3)
    rotated = (coords - coords.mean(0)) @ r
    assert np.allclose(kabsch_rotation(rotated, coords), r.T, atol=1e-10)


def test_alignment_restores_overlaps_and_reports_rmsd():
    coords, modes = _fake_molecule()
    r = _rotation(np.pi)  # 180 degrees, as seen for methanol
    coords_rot, modes_rot = coords @ r, modes @ r

    # Without alignment a mode no longer matches its own counterpart...
    assert compute_mode_overlap(modes_rot[0], modes[0]) < 0.99
    aligned, rmsd = align_modes_to_reference(modes_rot, coords_rot, coords)
    assert rmsd < 1e-8
    for i in range(len(modes)):
        assert compute_mode_overlap(aligned[i], modes[i]) > 0.999


def test_assignment_is_identity_after_alignment():
    coords, modes = _fake_molecule(seed=3)
    r = _rotation(2.1, axis=0)
    aligned, _ = align_modes_to_reference(modes @ r, coords @ r, coords)
    matches = match_modes(aligned, modes, threshold=0.5)
    assert all(ref == i for i, (ref, _ov) in matches.items())
