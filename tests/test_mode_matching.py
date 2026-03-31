"""Unit tests for mode_matching.py.

Tests mode overlap computation, mode matching/pairing, alignment matrix construction,
reduced mass calculation, and normalization. Includes both synthetic and real-data tests.

Does NOT import extract_mode_data_from_checkpoint (requires formchk).

Requirement: TEST-04
"""

import numpy as np
import pytest

from mace_gaussian.analysis.mode_matching import (
    DegenerateGroup,
    DegenerateGroupResult,
    compute_mode_overlap,
    compute_reduced_masses,
    compute_subspace_overlap,
    create_alignment_matrix,
    detect_degenerate_groups,
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


# ---------------------------------------------------------------------------
# Degenerate mode detection tests (Phase 19)
# ---------------------------------------------------------------------------


class TestDegenerateDetection:
    """Test detection of degenerate mode groups from DFT reference frequencies."""

    def test_no_degenerate_groups(self):
        """Frequencies well-separated (>5 cm-1) produce no degenerate groups."""
        freqs = np.array([500.0, 1000.0, 1500.0, 2000.0])
        groups = detect_degenerate_groups(freqs, threshold=5.0)
        assert len(groups) == 0

    def test_doubly_degenerate(self):
        """Two frequencies within 5 cm-1 form a doubly degenerate E group."""
        freqs = np.array([1000.0, 1003.0, 2000.0])
        groups = detect_degenerate_groups(freqs, threshold=5.0)
        assert len(groups) == 1
        g = groups[0]
        assert g.multiplicity == 2
        assert g.symmetry_label == "E"
        assert g.center_frequency == pytest.approx(1001.5, abs=0.1)
        assert set(g.ref_indices) == {0, 1}

    def test_triply_degenerate(self):
        """Three frequencies within 5 cm-1 form a triply degenerate T group."""
        freqs = np.array([1356.0, 1358.0, 1360.0, 3000.0])
        groups = detect_degenerate_groups(freqs, threshold=5.0)
        assert len(groups) == 1
        g = groups[0]
        assert g.multiplicity == 3
        assert g.symmetry_label == "T"
        assert g.center_frequency == pytest.approx(1358.0, abs=0.1)

    def test_multiple_groups(self):
        """Multiple degenerate groups detected correctly."""
        freqs = np.array([500.0, 502.0, 1000.0, 1002.0, 1004.0, 2000.0])
        groups = detect_degenerate_groups(freqs, threshold=5.0)
        assert len(groups) == 2
        # Sort by center frequency
        groups_sorted = sorted(groups, key=lambda g: g.center_frequency)
        assert groups_sorted[0].multiplicity == 2  # E at ~501
        assert groups_sorted[0].symmetry_label == "E"
        assert groups_sorted[1].multiplicity == 3  # T at ~1002
        assert groups_sorted[1].symmetry_label == "T"

    def test_threshold_boundary(self):
        """Frequencies exactly at threshold boundary are included; just beyond are not."""
        # Exactly 5 cm-1 apart -> should form group
        freqs_in = np.array([1000.0, 1005.0])
        groups_in = detect_degenerate_groups(freqs_in, threshold=5.0)
        assert len(groups_in) == 1

        # 6 cm-1 apart -> should NOT form group
        freqs_out = np.array([1000.0, 1006.0])
        groups_out = detect_degenerate_groups(freqs_out, threshold=5.0)
        assert len(groups_out) == 0

        # Exactly 0.5 cm-1 apart -> should form group
        freqs_in = np.array([1000.0, 1000.5])
        groups_in = detect_degenerate_groups(freqs_in, threshold=0.5)
        assert len(groups_in) == 1

        # 0.6 cm-1 apart -> should NOT form group
        freqs_out = np.array([1000.0, 1000.6])
        groups_out = detect_degenerate_groups(freqs_out, threshold=0.5)
        assert len(groups_out) == 0

class TestSubspaceOverlap:
    """Test subspace overlap computation via trace(M^T M) / k."""

    def test_perfect_subspace(self):
        """Identity-like sub-block gives overlap ~1.0."""
        # 3x3 alignment matrix with perfect 2x2 sub-block for a doubly degenerate group
        alignment = np.eye(3)
        overlap = compute_subspace_overlap(alignment, calc_indices=[0, 1], ref_indices=[0, 1])
        assert overlap == pytest.approx(1.0, abs=1e-10)

    def test_orthogonal_subspace(self):
        """Zero sub-block gives overlap 0.0."""
        alignment = np.zeros((3, 3))
        overlap = compute_subspace_overlap(alignment, calc_indices=[0, 1], ref_indices=[0, 1])
        assert overlap == pytest.approx(0.0, abs=1e-10)

    def test_partial_overlap(self):
        """Known sub-block gives trace(M^T M)/k within expected range."""
        # Create a 2x2 sub-block with known values
        alignment = np.array(
            [
                [0.8, 0.3, 0.0],
                [0.2, 0.9, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        # Sub-block for indices [0,1] x [0,1]:
        # M = [[0.8, 0.3], [0.2, 0.9]]
        # M^T M = [[0.8*0.8+0.2*0.2, 0.8*0.3+0.2*0.9], [0.3*0.8+0.9*0.2, 0.3*0.3+0.9*0.9]]
        #       = [[0.68, 0.42], [0.42, 0.90]]
        # trace = 0.68 + 0.90 = 1.58
        # overlap = 1.58 / 2 = 0.79
        overlap = compute_subspace_overlap(alignment, calc_indices=[0, 1], ref_indices=[0, 1])
        assert overlap == pytest.approx(0.79, abs=0.01)

    def test_rotated_degenerate_pair(self):
        """Two 2D modes rotated 45 deg within degenerate subspace give overlap ~1.0.

        Individual dot products would be ~0.707, but subspace overlap should be ~1.0
        because the subspace spanned is identical.
        """
        # DFT modes: [1,0] and [0,1] (standard basis)
        # ML modes: [cos45, sin45] and [-sin45, cos45] (rotated 45 deg)
        # The alignment matrix sub-block is the rotation matrix itself
        c = np.cos(np.pi / 4)
        s = np.sin(np.pi / 4)
        alignment = np.array(
            [
                [c, s, 0.0],
                [-s, c, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        # Sub-block M = [[c, s], [-s, c]] which is orthogonal
        # M^T M = I, trace = 2, overlap = 2/2 = 1.0
        overlap = compute_subspace_overlap(alignment, calc_indices=[0, 1], ref_indices=[0, 1])
        assert overlap == pytest.approx(1.0, abs=1e-10)

    def test_empty_indices(self):
        """Empty calc_indices returns 0.0."""
        alignment = np.eye(3)
        overlap = compute_subspace_overlap(alignment, calc_indices=[], ref_indices=[0, 1])
        assert overlap == pytest.approx(0.0, abs=1e-10)


class TestGroupAwareStatistics:
    """Test DegenerateGroupResult.statistics() with group-aware counting."""

    def test_statistics_count(self):
        """5 non-degenerate + 1 triple group = 6 total units (not 8)."""
        group = DegenerateGroup(
            ref_indices=[5, 6, 7],
            calc_indices=[5, 6, 7],
            center_frequency=1356.0,
            multiplicity=3,
            symmetry_label="T",
            subspace_overlap=0.95,
        )
        # Individual matches: 8 total modes (indices 0-7)
        matches = {i: (i, 0.9) for i in range(8)}
        result = DegenerateGroupResult(
            individual_matches=matches,
            groups=[group],
            non_degenerate_indices=[0, 1, 2, 3, 4],
        )
        stats = result.statistics()
        assert stats["total_units"] == 6  # 5 non-deg + 1 group

    def test_average_overlap(self):
        """Average overlap weighted correctly across groups and non-degenerate modes."""
        group = DegenerateGroup(
            ref_indices=[2, 3],
            calc_indices=[2, 3],
            center_frequency=1000.0,
            multiplicity=2,
            symmetry_label="E",
            subspace_overlap=0.8,
        )
        # 2 non-degenerate modes with overlap 1.0 each, 1 group with overlap 0.8
        matches = {0: (0, 1.0), 1: (1, 1.0), 2: (2, 0.7), 3: (3, 0.7)}
        result = DegenerateGroupResult(
            individual_matches=matches,
            groups=[group],
            non_degenerate_indices=[0, 1],
        )
        stats = result.statistics()
        # avg = (1.0 + 1.0 + 0.8) / 3 = 0.933...
        assert stats["avg_overlap"] == pytest.approx(2.8 / 3, abs=0.01)
        # confident: 1.0 >= 0.5, 1.0 >= 0.5, 0.8 >= 0.5 -> 3
        assert stats["confident_count"] == 3


class TestGroupRegressionData:
    """Test DegenerateGroupResult.regression_data() for group-averaged data points."""

    def test_group_averaged(self):
        """Triply degenerate group produces 1 data point with averaged freq/intensity."""
        group = DegenerateGroup(
            ref_indices=[0, 1, 2],
            calc_indices=[0, 1, 2],
            center_frequency=1358.0,
            multiplicity=3,
            symmetry_label="T",
            subspace_overlap=0.95,
        )
        matches = {0: (0, 0.7), 1: (1, 0.7), 2: (2, 0.7), 3: (3, 0.9)}
        result = DegenerateGroupResult(
            individual_matches=matches,
            groups=[group],
            non_degenerate_indices=[3],
        )
        freqs_calc = np.array([1350.0, 1355.0, 1360.0, 3000.0])
        freqs_ref = np.array([1356.0, 1358.0, 1360.0, 3050.0])
        ints_calc = np.array([10.0, 12.0, 14.0, 50.0])
        ints_ref = np.array([11.0, 13.0, 15.0, 55.0])

        dft_f, ml_f, dft_i, ml_i, mask = result.regression_data(
            freqs_calc, freqs_ref, ints_calc, ints_ref
        )
        # 1 group + 1 non-degenerate = 2 data points
        assert len(dft_f) == 2
        # The group-averaged point (could be first or second depending on order)
        # Find which index is the degenerate one
        deg_idx = np.where(mask)[0][0]
        assert dft_f[deg_idx] == pytest.approx(np.mean([1356.0, 1358.0, 1360.0]), abs=0.1)
        assert ml_f[deg_idx] == pytest.approx(np.mean([1350.0, 1355.0, 1360.0]), abs=0.1)
        assert dft_i[deg_idx] == pytest.approx(np.mean([11.0, 13.0, 15.0]), abs=0.1)
        assert ml_i[deg_idx] == pytest.approx(np.mean([10.0, 12.0, 14.0]), abs=0.1)

    def test_non_degenerate_unchanged(self):
        """Non-degenerate modes produce individual data points unchanged."""
        matches = {0: (0, 0.9), 1: (1, 0.85)}
        result = DegenerateGroupResult(
            individual_matches=matches,
            groups=[],
            non_degenerate_indices=[0, 1],
        )
        freqs_calc = np.array([1000.0, 2000.0])
        freqs_ref = np.array([1010.0, 2020.0])
        ints_calc = np.array([30.0, 60.0])
        ints_ref = np.array([35.0, 65.0])

        dft_f, ml_f, dft_i, ml_i, mask = result.regression_data(
            freqs_calc, freqs_ref, ints_calc, ints_ref
        )
        assert len(dft_f) == 2
        assert not np.any(mask)  # No degenerate points
        assert dft_f[0] == pytest.approx(1010.0)
        assert ml_f[0] == pytest.approx(1000.0)
        assert dft_i[0] == pytest.approx(35.0)
        assert ml_i[0] == pytest.approx(30.0)
