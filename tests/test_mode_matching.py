"""Unit tests for mode_matching.py.

Tests mode overlap computation, mode matching/pairing, alignment matrix construction,
reduced mass calculation, and normalization. Includes both synthetic and real-data tests.

Does NOT import extract_mode_data_from_checkpoint (requires formchk).

Requirement: TEST-04
"""

import numpy as np
import pytest

from mace_gaussian.analysis.mode_matching import (
    compute_mode_overlap,
    compute_reduced_masses,
    create_alignment_matrix,
    match_modes,
    normalize_mode,
)

# ---------------------------------------------------------------------------
# Synthetic data tests (no fixtures needed)
# ---------------------------------------------------------------------------


class TestComputeModeOverlap:
    """Test the scalar product overlap between mode vectors."""

    def test_identical_modes_perfect_overlap(self):
        """Overlap of a mode with itself should be exactly 1.0."""
        mode = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        assert compute_mode_overlap(mode, mode) == pytest.approx(1.0, abs=1e-10)

    def test_opposite_modes_perfect_overlap(self):
        """Sign-flipped modes should still have overlap 1.0 (absolute dot product)."""
        mode = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        assert compute_mode_overlap(mode, -mode) == pytest.approx(1.0, abs=1e-10)

    def test_orthogonal_modes_zero_overlap(self):
        """Orthogonal modes should have zero overlap."""
        mode1 = np.array([[1.0, 0.0, 0.0]])
        mode2 = np.array([[0.0, 1.0, 0.0]])
        assert compute_mode_overlap(mode1, mode2) == pytest.approx(0.0, abs=1e-10)


class TestMatchModes:
    """Test mode-to-mode matching via maximum overlap."""

    def test_match_modes_identity(self):
        """Matching modes against themselves gives 1-to-1 identity mapping."""
        rng = np.random.RandomState(42)
        modes = rng.randn(3, 2, 3)  # 3 modes, 2 atoms, 3 coords

        matches = match_modes(modes, modes, threshold=0.5)

        for i in range(3):
            ref_idx, overlap = matches[i]
            assert ref_idx == i, f"Mode {i} should map to itself, got {ref_idx}"
            assert overlap == pytest.approx(1.0, abs=1e-10)

    def test_match_modes_permuted(self):
        """Matching permuted modes correctly recovers the original mapping."""
        # Create 3 clearly distinct modes (displacement on different atoms/axes)
        modes = np.zeros((3, 3, 3))
        modes[0, 0, 0] = 1.0  # Mode 0: atom 0, x-axis
        modes[1, 1, 1] = 1.0  # Mode 1: atom 1, y-axis
        modes[2, 2, 2] = 1.0  # Mode 2: atom 2, z-axis

        # Permute: [2, 0, 1]
        modes_permuted = modes[[2, 0, 1]]

        matches = match_modes(modes_permuted, modes, threshold=0.5)

        # Permuted mode 0 was originally mode 2, so should map to ref mode 2
        assert matches[0][0] == 2
        assert matches[1][0] == 0
        assert matches[2][0] == 1

    def test_match_modes_threshold(self):
        """Modes with moderate overlap are still recorded (with warning for low overlap)."""
        # Create two modes at 45 degrees in a 2D subspace
        mode1 = np.array([[1.0, 0.0, 0.0]])
        mode2 = np.array([[1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0]])

        modes_calc = mode1.reshape(1, 1, 3)
        modes_ref = mode2.reshape(1, 1, 3)

        matches = match_modes(modes_calc, modes_ref, threshold=0.9)

        # match_modes records best match regardless of threshold
        ref_idx, overlap = matches[0]
        assert ref_idx == 0
        expected_overlap = 1.0 / np.sqrt(2)  # cos(45 deg) ~ 0.707
        assert overlap == pytest.approx(expected_overlap, abs=1e-10)


class TestCreateAlignmentMatrix:
    """Test alignment matrix construction."""

    def test_create_alignment_matrix_shape(self):
        """Alignment matrix has shape (n_modes_calc, n_modes_ref)."""
        rng = np.random.RandomState(42)
        modes_a = rng.randn(3, 2, 3)
        modes_b = rng.randn(4, 2, 3)

        matrix = create_alignment_matrix(modes_a, modes_b)
        assert matrix.shape == (3, 4)

    def test_create_alignment_matrix_diagonal_dominance(self):
        """Self-alignment matrix should have diagonal values ~1.0."""
        rng = np.random.RandomState(42)
        modes = rng.randn(3, 2, 3)

        matrix = create_alignment_matrix(modes, modes)

        # Diagonal should be ~1.0 (self-overlap)
        assert np.all(np.diag(matrix) > 0.99)

        # Off-diagonal should be strictly less than diagonal
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert matrix[i, j] < matrix[i, i]


class TestComputeReducedMasses:
    """Test reduced mass computation from modes and atomic masses."""

    def test_compute_reduced_masses(self):
        """Reduced mass = sum of mass_atom * |displacement_atom|^2 for each mode."""
        # 1 mode, 2 atoms, each displaced along x by 1.0
        modes = np.array([[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]])
        masses = np.array([16.0, 1.0])

        reduced = compute_reduced_masses(modes, masses)

        # 16.0 * (1^2) + 1.0 * (1^2) = 17.0
        assert reduced[0] == pytest.approx(17.0)


class TestNormalizeMode:
    """Test mode norm computation."""

    def test_normalize_mode_pythagorean(self):
        """normalize_mode returns Euclidean norm of flattened mode vector."""
        result = normalize_mode(np.array([[3.0, 4.0, 0.0]]))
        assert result == pytest.approx(5.0)

    def test_normalize_mode_unit_vector(self):
        """Unit displacement on one atom gives norm 1.0."""
        result = normalize_mode(np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
        assert result == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Hungarian algorithm tests (bijective matching)
# ---------------------------------------------------------------------------


class TestHungarianMatching:
    """Test that match_modes() uses Hungarian algorithm for globally optimal bijective pairing."""

    def test_hungarian_prevents_duplicate_assignment(self):
        """Greedy fails here: calc modes 0 and 1 are both closest to ref mode 0.

        Hungarian should produce a bijective (1-to-1) assignment where no two
        calc modes map to the same ref mode.
        """
        # 3 atoms, 3 modes each
        # Ref mode 0: strong signal on atom 0 x-axis
        # Ref mode 1: strong signal on atom 1 y-axis
        # Ref mode 2: strong signal on atom 2 z-axis
        modes_ref = np.zeros((3, 3, 3))
        modes_ref[0, 0, 0] = 1.0
        modes_ref[1, 1, 1] = 1.0
        modes_ref[2, 2, 2] = 1.0

        # Calc mode 0: very similar to ref mode 0 (overlap ~0.98)
        modes_calc = np.zeros((3, 3, 3))
        modes_calc[0, 0, 0] = 0.98
        modes_calc[0, 0, 1] = 0.2  # small perturbation

        # Calc mode 1: also most similar to ref mode 0 under greedy (overlap ~0.95)
        # but has secondary overlap with ref mode 1
        modes_calc[1, 0, 0] = 0.7
        modes_calc[1, 1, 1] = 0.7  # significant component along ref mode 1

        # Calc mode 2: clearly matches ref mode 2
        modes_calc[2, 2, 2] = 1.0

        matches = match_modes(modes_calc, modes_ref, threshold=0.1)

        # All 3 calc modes should be matched
        assert len(matches) == 3

        # Bijection: all assigned ref indices must be unique
        assigned_refs = [matches[i][0] for i in range(3)]
        assert len(set(assigned_refs)) == 3, (
            f"Expected bijective assignment (3 unique refs), got {assigned_refs}"
        )

        # Calc mode 2 should match ref mode 2 (clear match)
        assert matches[2][0] == 2

    def test_unmatched_modes_more_calc_than_ref(self):
        """When n_calc > n_ref, extra calc modes get (None, 0.0)."""
        # 5 calc modes, 3 ref modes, 1 atom each
        modes_ref = np.zeros((3, 1, 3))
        modes_ref[0, 0, 0] = 1.0
        modes_ref[1, 0, 1] = 1.0
        modes_ref[2, 0, 2] = 1.0

        modes_calc = np.zeros((5, 1, 3))
        modes_calc[0, 0, 0] = 1.0  # matches ref 0
        modes_calc[1, 0, 1] = 1.0  # matches ref 1
        modes_calc[2, 0, 2] = 1.0  # matches ref 2
        # Calc modes 3 and 4: no good ref available
        modes_calc[3, 0, :] = [0.5, 0.5, 0.0]
        modes_calc[4, 0, :] = [0.0, 0.5, 0.5]

        matches = match_modes(modes_calc, modes_ref, threshold=0.1)

        # All 5 calc modes should have entries
        assert len(matches) == 5

        # 3 matched pairs should have valid ref indices
        matched = [(i, matches[i]) for i in range(5) if matches[i][0] is not None]
        unmatched = [(i, matches[i]) for i in range(5) if matches[i][0] is None]

        assert len(matched) == 3, f"Expected 3 matched pairs, got {len(matched)}"
        assert len(unmatched) == 2, f"Expected 2 unmatched pairs, got {len(unmatched)}"

        # Unmatched should have overlap 0.0
        for _idx, (ref_idx, overlap) in unmatched:
            assert ref_idx is None
            assert overlap == 0.0

    def test_unmatched_modes_more_ref_than_calc(self):
        """When n_calc < n_ref, all calc modes get valid matches (no None)."""
        # 3 calc modes, 5 ref modes, 1 atom each
        modes_ref = np.zeros((5, 1, 3))
        modes_ref[0, 0, 0] = 1.0
        modes_ref[1, 0, 1] = 1.0
        modes_ref[2, 0, 2] = 1.0
        modes_ref[3, 0, :] = [0.5, 0.5, 0.0]
        modes_ref[4, 0, :] = [0.0, 0.5, 0.5]

        modes_calc = np.zeros((3, 1, 3))
        modes_calc[0, 0, 0] = 1.0
        modes_calc[1, 0, 1] = 1.0
        modes_calc[2, 0, 2] = 1.0

        matches = match_modes(modes_calc, modes_ref, threshold=0.1)

        # All 3 calc modes matched, no None
        assert len(matches) == 3
        for i in range(3):
            ref_idx, overlap = matches[i]
            assert ref_idx is not None, f"Calc mode {i} should have a valid match"
            assert overlap > 0.5

    def test_low_overlap_pairs_kept(self):
        """Low-overlap pairs should be kept in dict with their overlap value preserved."""
        # Create nearly orthogonal modes
        modes_calc = np.zeros((2, 1, 3))
        modes_calc[0, 0, 0] = 1.0  # x-axis
        modes_calc[1, 0, 1] = 1.0  # y-axis

        modes_ref = np.zeros((2, 1, 3))
        modes_ref[0, 0, 0] = 1.0  # x-axis (good match for calc 0)
        # ref 1: mostly z-axis with tiny y component -> low overlap with calc 1
        modes_ref[1, 0, :] = [0.0, 0.1, 0.995]

        matches = match_modes(modes_calc, modes_ref, threshold=0.5)

        # Both calc modes should be in dict
        assert len(matches) == 2

        # Calc mode 0 -> ref 0 with high overlap
        assert matches[0][0] == 0
        assert matches[0][1] > 0.9

        # Calc mode 1 -> ref 1 with LOW overlap (kept, not dropped)
        assert matches[1][0] == 1
        assert matches[1][1] < 0.5  # below threshold but still present
        assert matches[1][1] > 0.0  # has some overlap value


# ---------------------------------------------------------------------------
# Real data tests (using .fchk fixtures)
# ---------------------------------------------------------------------------


class TestRealFchkModeMatching:
    """Test mode matching with actual water .fchk data."""

    def test_real_fchk_mode_overlap_self_match(self, water_dft_fchk):
        """Water DFT modes matched against themselves give overlap > 0.99 for all 3 modes."""
        from mace_gaussian.gaussian.fchk import extract_modes_from_fchk

        modes, _, _, _, _ = extract_modes_from_fchk(water_dft_fchk, force_harmonic=True)

        matches = match_modes(modes, modes)

        assert len(matches) == 3, "Water should have 3 vibrational modes"
        for i in range(3):
            ref_idx, overlap = matches[i]
            assert ref_idx == i, f"Mode {i} should self-match to index {i}"
            assert overlap > 0.99, f"Self-match overlap for mode {i} should be > 0.99"

    def test_real_fchk_dft_vs_ml_overlap(self, water_dft_fchk, water_ml_fchk):
        """Water DFT vs ML modes: at least 2 of 3 modes should have overlap > 0.5."""
        from mace_gaussian.gaussian.fchk import extract_modes_from_fchk

        modes_dft, _, _, _, _ = extract_modes_from_fchk(water_dft_fchk, force_harmonic=True)
        modes_ml, _, _, _, _ = extract_modes_from_fchk(water_ml_fchk, force_harmonic=True)

        matches = match_modes(modes_ml, modes_dft, threshold=0.1)

        good_matches = sum(1 for _, (_, overlap) in matches.items() if overlap > 0.5)
        assert good_matches >= 2, (
            f"Expected at least 2 modes with overlap > 0.5, got {good_matches}"
        )
