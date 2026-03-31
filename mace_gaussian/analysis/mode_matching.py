#!/usr/bin/env python3
"""
Mode Matching via Normal Mode Overlap

Translates build_modes.gk awk script to Python.
Matches vibrational modes between two calculations using scalar product
(dot product) rather than just frequency ordering.

This handles cases where mode ordering changes between DFT and ML calculations.

Uses .fchk (formatted checkpoint) files from Gaussian for clean parsing.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.optimize import linear_sum_assignment

from ..gaussian.fchk import extract_modes_from_fchk, get_fchk_from_chk

logger = logging.getLogger(__name__)


def extract_mode_data_from_checkpoint(
    chk_or_fchk_file: str, force_harmonic: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Extract vibrational modes, frequencies, coordinates, and masses from Gaussian checkpoint.

    Can accept either .chk or .fchk file. If .chk is provided, converts to .fchk.

    Parameters
    ----------
    chk_or_fchk_file : str
        Path to Gaussian .chk or .fchk file
    force_harmonic : bool, optional
        If True, force extraction of harmonic frequencies/modes only.
        Default: False. Set to True for mode overlap heatmaps.

    Returns
    -------
    modes : np.ndarray
        Vibrational mode vectors, shape (n_modes, n_atoms, 3)
    frequencies : np.ndarray
        Vibrational frequencies in cm^-1, shape (n_modes,)
    coords : np.ndarray
        Atomic coordinates in Angstroms, shape (n_atoms, 3)
    masses : np.ndarray
        Atomic masses in amu, shape (n_atoms,)
    n_atoms : int
        Number of atoms
    """
    file_path = Path(chk_or_fchk_file)

    # Convert .chk to .fchk if necessary
    if file_path.suffix == ".chk":
        fchk_file = get_fchk_from_chk(str(file_path))
    elif file_path.suffix == ".fchk":
        fchk_file = str(file_path)
    else:
        raise ValueError(f"Expected .chk or .fchk file, got: {file_path.suffix}")

    # Extract data from .fchk (force harmonic if requested)
    modes, frequencies, coords, masses, n_atoms = extract_modes_from_fchk(
        fchk_file, force_harmonic=force_harmonic
    )

    logger.info(f"Extracted {modes.shape[0]} modes for {n_atoms} atoms from {fchk_file}")

    return modes, frequencies, coords, masses, n_atoms


def compute_reduced_masses(modes: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Compute reduced masses for each vibrational mode.

    μ_i = Σ_atoms m_atom * (L_i,atom · L_i,atom)

    where L_i,atom is the displacement vector for atom in mode i.

    Parameters
    ----------
    modes : np.ndarray
        Shape (n_modes, n_atoms, 3)
    masses : np.ndarray
        Shape (n_atoms,)

    Returns
    -------
    reduced_masses : np.ndarray
        Shape (n_modes,)
    """
    n_modes, n_atoms, _ = modes.shape
    reduced_masses = np.zeros(n_modes)

    for i in range(n_modes):
        mode_i = modes[i]  # (n_atoms, 3)
        # Compute sum of mass * (L·L) for each atom
        for atom_idx in range(n_atoms):
            displacement = mode_i[atom_idx]  # (3,)
            reduced_masses[i] += masses[atom_idx] * np.dot(displacement, displacement)

    return reduced_masses


def normalize_mode(mode: np.ndarray) -> float:
    """
    Compute the norm of a mode vector.

    Parameters
    ----------
    mode : np.ndarray
        Shape (n_atoms, 3)

    Returns
    -------
    norm : float
    """
    return np.sqrt(np.sum(mode**2))


def compute_mode_overlap(mode1: np.ndarray, mode2: np.ndarray) -> float:
    """
    Compute normalized scalar product (overlap) between two modes.

    overlap = |mode1 . mode2| / (||mode1|| x ||mode2||)

    Returns value between 0 and 1:
    - 1: modes are perfectly aligned (same or opposite direction)
    - 0: modes are orthogonal (completely different)

    Parameters
    ----------
    mode1, mode2 : np.ndarray
        Shape (n_atoms, 3)

    Returns
    -------
    overlap : float
        Between 0 and 1
    """
    # Flatten modes
    mode1_flat = mode1.flatten()
    mode2_flat = mode2.flatten()

    # Compute norms
    norm1 = np.linalg.norm(mode1_flat)
    norm2 = np.linalg.norm(mode2_flat)

    # Compute dot product
    dot_product = np.dot(mode1_flat, mode2_flat)

    # Normalized overlap (take absolute value to handle sign flip)
    overlap = np.abs(dot_product) / (norm1 * norm2)

    return overlap


def match_modes(
    modes_calc: np.ndarray, modes_ref: np.ndarray, threshold: float = 0.5
) -> dict[int, tuple[int | None, float]]:
    """
    Match modes between calculation and reference via Hungarian algorithm.

    Uses globally optimal bijective assignment (linear_sum_assignment) to find
    the 1-to-1 mode pairing that maximizes total overlap. This prevents the
    greedy failure where two calc modes claim the same ref mode.

    Parameters
    ----------
    modes_calc : np.ndarray
        Modes from current calculation, shape (n_modes, n_atoms, 3)
    modes_ref : np.ndarray
        Modes from reference calculation, shape (n_modes_ref, n_atoms, 3)
    threshold : float
        Minimum overlap threshold for logging warnings (default: 0.5).
        Low-overlap pairs are kept in the dict regardless.

    Returns
    -------
    matches : dict
        Maps calc_mode_idx -> (ref_mode_idx | None, overlap).
        When n_calc > n_ref, unmatched calc modes get (None, 0.0).
    """
    n_modes_calc = modes_calc.shape[0]
    n_modes_ref = modes_ref.shape[0]

    logger.info(f"Matching {n_modes_calc} calculation modes to {n_modes_ref} reference modes")

    # Build full overlap matrix using existing function
    overlap_matrix = create_alignment_matrix(modes_calc, modes_ref)

    # Hungarian algorithm: negate overlap for maximization (scipy minimizes)
    row_ind, col_ind = linear_sum_assignment(-overlap_matrix)

    # Build matches dict from assignment
    matches: dict[int, tuple[int | None, float]] = {}
    assigned_calc = set(row_ind)

    for r, c in zip(row_ind, col_ind):
        overlap = overlap_matrix[r, c]
        matches[r] = (int(c), float(overlap))

        if overlap >= threshold:
            logger.debug(f"Mode {r:3d} -> Ref mode {c:3d} (overlap: {overlap:.4f})")
        else:
            logger.warning(f"Mode {r:3d} -> Ref mode {c:3d} (low overlap: {overlap:.4f})")

    # Unmatched calc modes (when n_calc > n_ref)
    for i in range(n_modes_calc):
        if i not in assigned_calc:
            matches[i] = (None, 0.0)
            logger.warning(f"Mode {i:3d} -> No match available (n_calc > n_ref)")

    return matches


def create_alignment_matrix(modes_calc: np.ndarray, modes_ref: np.ndarray) -> np.ndarray:
    """
    Create full alignment matrix showing overlap between all mode pairs.

    Parameters
    ----------
    modes_calc : np.ndarray
        Shape (n_modes_calc, n_atoms, 3)
    modes_ref : np.ndarray
        Shape (n_modes_ref, n_atoms, 3)

    Returns
    -------
    alignment_matrix : np.ndarray
        Shape (n_modes_calc, n_modes_ref)
        Element [i,j] is overlap between calc mode i and ref mode j
    """
    n_modes_calc = modes_calc.shape[0]
    n_modes_ref = modes_ref.shape[0]

    alignment_matrix = np.zeros((n_modes_calc, n_modes_ref))

    for i in range(n_modes_calc):
        for j in range(n_modes_ref):
            alignment_matrix[i, j] = compute_mode_overlap(modes_calc[i], modes_ref[j])

    return alignment_matrix


def plot_mode_overlap_heatmap(
    alignment_matrix: np.ndarray,
    output_file: Optional[str] = None,
    calc_label: str = "ML Calculation",
    ref_label: str = "DFT Reference",
    freqs_calc: Optional[np.ndarray] = None,
    freqs_ref: Optional[np.ndarray] = None,
    matches: Optional[dict[int, tuple[int, float]]] = None,
    x_labels: Optional[list[str]] = None,
    y_labels: Optional[list[str]] = None,
) -> None:
    """
    Plot heatmap of mode overlap matrix with elegant pastel styling.

    Parameters
    ----------
    alignment_matrix : np.ndarray
        Shape (n_modes_calc, n_modes_ref), overlap values between 0 and 1
    output_file : str, optional
        If provided, save plot to this file
    calc_label : str
        Label for calculation modes (y-axis)
    ref_label : str
        Label for reference modes (x-axis)
    freqs_calc : np.ndarray, optional
        Frequencies for calculation modes (for labels)
    freqs_ref : np.ndarray, optional
        Frequencies for reference modes (for labels)
    matches : dict, optional
        If provided, highlight the best matches with markers
    """
    n_modes_calc, n_modes_ref = alignment_matrix.shape

    # Set elegant style with clean background
    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(10, 8), facecolor="white")
    ax.set_facecolor("#fafafa")

    # Create elegant pastel red colormap: ivory -> soft red -> medium red
    from matplotlib.colors import LinearSegmentedColormap

    colors = [
        "#f8f8f5",  # Warm ivory (no overlap)
        "#fad4d4",  # Soft pastel red (weak overlap)
        "#f5a5a5",  # Medium pastel red (medium overlap)
        "#e87d7d",  # Stronger pastel red (strong overlap)
    ]
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list("pastel_elegant", colors, N=n_bins)

    # Create heatmap with subtle styling
    im = ax.imshow(
        alignment_matrix,
        cmap=cmap,
        vmin=0,
        vmax=1,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )

    # Add elegant colorbar with refined styling
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.ax.tick_params(labelsize=9, length=3, width=0.5, colors="#5a5a5a")
    cbar.set_label(
        "Mode Overlap",
        rotation=270,
        labelpad=18,
        fontsize=10,
        color="#4a4a4a",
        family="sans-serif",
        weight="normal",
    )
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor("#d0d0d0")

    # Set ticks
    ax.set_xticks(np.arange(n_modes_ref))
    ax.set_yticks(np.arange(n_modes_calc))

    # Create refined labels with frequencies (use overrides if provided)
    if x_labels is None:
        if freqs_ref is not None:
            x_labels = [f"{i}\n{freqs_ref[i]:.0f}" for i in range(n_modes_ref)]
        else:
            x_labels = [f"{i}" for i in range(n_modes_ref)]

    if y_labels is None:
        if freqs_calc is not None:
            y_labels = [f"{i}: {freqs_calc[i]:.0f}" for i in range(n_modes_calc)]
        else:
            y_labels = [f"{i}" for i in range(n_modes_calc)]

    ax.set_xticklabels(x_labels, fontsize=8.5, color="#4a4a4a", family="sans-serif")
    ax.set_yticklabels(y_labels, fontsize=8.5, color="#4a4a4a", family="sans-serif")

    # Elegant axis labels with proper LaTeX formatting
    ax.set_xlabel(
        f"{ref_label} Mode (cm$^{{-1}}$)",
        fontsize=14,
        color="#3a3a3a",
        family="sans-serif",
        weight="normal",
        labelpad=8,
    )
    ax.set_ylabel(
        f"{calc_label} Mode (cm$^{{-1}}$)",
        fontsize=14,
        color="#3a3a3a",
        family="sans-serif",
        weight="normal",
        labelpad=8,
    )
    ax.set_title(
        "Vibrational Mode Overlap Matrix",
        fontsize=16,
        color="#2a2a2a",
        family="sans-serif",
        weight="normal",
        pad=15,
    )

    # Subtle tick styling
    ax.tick_params(axis="both", which="major", length=4, width=0.5, colors="#6a6a6a")

    # Grid lines removed for cleaner appearance
    ax.grid(False)

    # Remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["left"].set_color("#c0c0c0")
    ax.spines["bottom"].set_color("#c0c0c0")

    # Add overlap values as text with refined color logic
    for i in range(n_modes_calc):
        for j in range(n_modes_ref):
            value = alignment_matrix[i, j]
            # Softer text color choices for better readability
            if value > 0.7:
                text_color = "#2a2a2a"  # Dark gray for high overlap
                weight = "semibold"
            elif value > 0.3:
                text_color = "#4a4a4a"  # Medium gray
                weight = "normal"
            else:
                text_color = "#7a7a7a"  # Light gray for low overlap
                weight = "normal"

            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=7.5,
                family="monospace",
                weight=weight,
            )

    # Draw borders around Hungarian-matched cells
    if matches is not None:
        for calc_idx, (ref_idx, overlap) in matches.items():
            if ref_idx is None:
                continue
            linestyle = "-" if overlap >= 0.5 else "--"
            border_color = "#4a4a4a"
            rect = Rectangle(
                (ref_idx - 0.5, calc_idx - 0.5),
                1,
                1,
                linewidth=2,
                edgecolor=border_color,
                facecolor="none",
                linestyle=linestyle,
                zorder=5,
            )
            ax.add_patch(rect)

        # Legend for border styles
        legend_elements = [
            Line2D(
                [0],
                [0],
                color="#4a4a4a",
                linewidth=2,
                linestyle="-",
                label="Confident match ($\\geq$ 0.5)",
            ),
            Line2D(
                [0],
                [0],
                color="#4a4a4a",
                linewidth=2,
                linestyle="--",
                label="Uncertain match (< 0.5)",
            ),
        ]
        ax.legend(handles=legend_elements, loc="upper left", fontsize=8, framealpha=0.9)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved mode overlap heatmap to {output_file}")
        print(f"Saved heatmap: {output_file}")
    else:
        plt.show()

    plt.close()
    plt.style.use("default")  # Reset style


# ---------------------------------------------------------------------------
# Degenerate mode handling (Phase 19)
# ---------------------------------------------------------------------------


@dataclass
class DegenerateGroup:
    """A group of nearly-degenerate modes in the DFT reference.

    Modes within a frequency threshold are grouped together. The subspace
    overlap (trace(M^T M)/k) replaces individual dot-product overlaps for
    these modes, avoiding misleadingly low scores when eigenvectors are
    rotated within the degenerate subspace.
    """

    ref_indices: list[int]
    calc_indices: list[int] = field(default_factory=list)
    center_frequency: float = 0.0
    multiplicity: int = 0
    symmetry_label: str = ""
    subspace_overlap: float = 0.0


@dataclass
class DegenerateGroupResult:
    """Wraps Hungarian match results with degenerate group awareness.

    Provides group-aware statistics and regression data where degenerate
    groups count as single units (D-06) and produce averaged data points (D-07).
    """

    individual_matches: dict[int, tuple[int | None, float]]
    groups: list[DegenerateGroup]
    non_degenerate_indices: list[int]

    def effective_overlap(self, ref_idx: int) -> float:
        """Get overlap for a ref mode.

        Returns group subspace overlap if degenerate, individual overlap otherwise.
        """
        for group in self.groups:
            if ref_idx in group.ref_indices:
                return group.subspace_overlap
        # Non-degenerate: find the calc mode matched to this ref mode
        for _calc_idx, (matched_ref, overlap) in self.individual_matches.items():
            if matched_ref == ref_idx:
                return overlap
        return 0.0

    def statistics(self) -> dict:
        """Group-aware match statistics.

        Returns dict with:
        - total_units: number of non-degenerate modes + number of groups
        - avg_overlap: mean of effective overlaps (one per unit)
        - confident_count: units with effective overlap >= 0.5
        """
        overlaps = []
        # One overlap per non-degenerate mode
        for ref_idx in self.non_degenerate_indices:
            overlaps.append(self.effective_overlap(ref_idx))
        # One overlap per group
        for group in self.groups:
            overlaps.append(group.subspace_overlap)

        total_units = len(overlaps)
        avg_overlap = float(np.mean(overlaps)) if overlaps else 0.0
        confident_count = sum(1 for o in overlaps if o >= 0.5)

        return {
            "total_units": total_units,
            "avg_overlap": avg_overlap,
            "confident_count": confident_count,
        }

    def regression_data(
        self,
        freqs_calc: np.ndarray,
        freqs_ref: np.ndarray,
        ints_calc: np.ndarray,
        ints_ref: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """One data point per non-degenerate mode + one per group (averaged).

        Returns (dft_freqs, ml_freqs, dft_ints, ml_ints, is_degenerate_mask).
        For degenerate groups: mean of DFT/ML freqs and intensities in group.
        """
        dft_freqs_list = []
        ml_freqs_list = []
        dft_ints_list = []
        ml_ints_list = []
        is_deg_list = []

        # Non-degenerate modes: individual data points
        for ref_idx in self.non_degenerate_indices:
            # Find the calc mode matched to this ref
            calc_idx = None
            for ci, (ri, _overlap) in self.individual_matches.items():
                if ri == ref_idx:
                    calc_idx = ci
                    break
            if calc_idx is not None:
                dft_freqs_list.append(freqs_ref[ref_idx])
                ml_freqs_list.append(freqs_calc[calc_idx])
                dft_ints_list.append(ints_ref[ref_idx])
                ml_ints_list.append(ints_calc[calc_idx])
                is_deg_list.append(False)

        # Degenerate groups: averaged data points
        for group in self.groups:
            if group.calc_indices:
                dft_freqs_list.append(float(np.mean(freqs_ref[group.ref_indices])))
                ml_freqs_list.append(float(np.mean(freqs_calc[group.calc_indices])))
                dft_ints_list.append(float(np.mean(ints_ref[group.ref_indices])))
                ml_ints_list.append(float(np.mean(ints_calc[group.calc_indices])))
                is_deg_list.append(True)

        return (
            np.array(dft_freqs_list),
            np.array(ml_freqs_list),
            np.array(dft_ints_list),
            np.array(ml_ints_list),
            np.array(is_deg_list, dtype=bool),
        )


def detect_degenerate_groups(
    ref_frequencies: np.ndarray, threshold: float = 5.0
) -> list[DegenerateGroup]:
    """Identify groups of nearly-degenerate modes in DFT reference frequencies.

    Modes are grouped if consecutive frequencies (sorted) differ by <= threshold.
    Groups of size 1 are non-degenerate and excluded from the result.

    Parameters
    ----------
    ref_frequencies : np.ndarray
        DFT reference frequencies in cm-1.
    threshold : float
        Maximum frequency difference (cm-1) to consider modes degenerate.

    Returns
    -------
    List of DegenerateGroup with ref_indices, center_frequency, multiplicity,
    and symmetry_label populated. calc_indices and subspace_overlap are left
    at defaults (to be filled after Hungarian matching).
    """
    if len(ref_frequencies) == 0:
        return []

    sorted_indices = np.argsort(ref_frequencies)
    sorted_freqs = ref_frequencies[sorted_indices]

    groups: list[DegenerateGroup] = []
    current_group = [int(sorted_indices[0])]

    for i in range(1, len(sorted_freqs)):
        if sorted_freqs[i] - sorted_freqs[i - 1] <= threshold:
            current_group.append(int(sorted_indices[i]))
        else:
            if len(current_group) >= 2:
                _add_degenerate_group(groups, current_group, ref_frequencies)
            current_group = [int(sorted_indices[i])]

    # Don't forget the last group
    if len(current_group) >= 2:
        _add_degenerate_group(groups, current_group, ref_frequencies)

    return groups


def _add_degenerate_group(
    groups: list[DegenerateGroup],
    indices: list[int],
    freqs: np.ndarray,
) -> None:
    """Helper to create and append a DegenerateGroup."""
    k = len(indices)
    label = {2: "E", 3: "T"}.get(k, f"deg-{k}")
    groups.append(
        DegenerateGroup(
            ref_indices=indices,
            calc_indices=[],
            center_frequency=float(np.mean(freqs[indices])),
            multiplicity=k,
            symmetry_label=label,
            subspace_overlap=0.0,
        )
    )


def compute_subspace_overlap(
    alignment_matrix: np.ndarray,
    calc_indices: list[int],
    ref_indices: list[int],
) -> float:
    """Compute subspace overlap for a degenerate group.

    Uses trace(M^T M) / k where M is the sub-block of the alignment matrix
    and k = min(len(calc_indices), len(ref_indices)).

    Parameters
    ----------
    alignment_matrix : np.ndarray
        Full overlap matrix, shape (n_calc, n_ref).
    calc_indices : list[int]
        ML mode indices matched to this group.
    ref_indices : list[int]
        DFT mode indices in this group.

    Returns
    -------
    Subspace overlap in [0, 1]. Value of 1 means perfect subspace alignment.
    """
    k = min(len(calc_indices), len(ref_indices))
    if k == 0:
        return 0.0

    M = alignment_matrix[np.ix_(calc_indices, ref_indices)]
    return float(np.trace(M.T @ M) / k)


def build_degenerate_result(
    matches: dict[int, tuple[int | None, float]],
    alignment_matrix: np.ndarray,
    freqs_ref: np.ndarray,
    threshold: float = 5.0,
) -> DegenerateGroupResult:
    """Build a DegenerateGroupResult from Hungarian matches and DFT frequencies.

    Detects degenerate groups, assigns calc indices from matches, computes
    subspace overlap for each group, and identifies non-degenerate ref indices.

    Parameters
    ----------
    matches : dict
        From match_modes(): calc_idx -> (ref_idx | None, overlap).
    alignment_matrix : np.ndarray
        Full overlap matrix, shape (n_calc, n_ref).
    freqs_ref : np.ndarray
        DFT reference frequencies in cm-1.
    threshold : float
        Frequency threshold for degenerate grouping.

    Returns
    -------
    DegenerateGroupResult with groups populated with calc_indices and
    subspace_overlap values.
    """
    groups = detect_degenerate_groups(freqs_ref, threshold=threshold)

    # Build reverse mapping: ref_idx -> calc_idx
    ref_to_calc: dict[int, int] = {}
    for calc_idx, (ref_idx, _overlap) in matches.items():
        if ref_idx is not None:
            ref_to_calc[ref_idx] = calc_idx

    # Assign calc_indices and compute subspace overlap for each group
    grouped_ref_indices: set[int] = set()
    for group in groups:
        group.calc_indices = [ref_to_calc[ri] for ri in group.ref_indices if ri in ref_to_calc]
        group.subspace_overlap = compute_subspace_overlap(
            alignment_matrix, group.calc_indices, group.ref_indices
        )
        grouped_ref_indices.update(group.ref_indices)

    # Non-degenerate ref indices
    all_ref_indices = set(range(len(freqs_ref)))
    non_degenerate = sorted(all_ref_indices - grouped_ref_indices)

    return DegenerateGroupResult(
        individual_matches=matches,
        groups=groups,
        non_degenerate_indices=non_degenerate,
    )


def collapse_alignment_matrix(
    alignment_matrix: np.ndarray,
    groups: list[DegenerateGroup],
    freqs_calc: np.ndarray,
    freqs_ref: np.ndarray,
    matches: dict[int, tuple[int | None, float]],
) -> tuple[np.ndarray, list[str], list[str]]:
    """Collapse alignment matrix by merging degenerate groups into single entries.

    Non-degenerate modes become individual rows/columns labeled by frequency.
    Degenerate groups become single rows/columns labeled with symmetry and
    center frequency (e.g. "T2 @1356"). Group cells use subspace overlap;
    cross cells (one side group, other not) use max overlap from the sub-block.

    Parameters
    ----------
    alignment_matrix : np.ndarray
        Full overlap matrix, shape (n_calc, n_ref).
    groups : list[DegenerateGroup]
        Degenerate groups from build_degenerate_result().
    freqs_calc : np.ndarray
        ML calculation frequencies.
    freqs_ref : np.ndarray
        DFT reference frequencies.
    matches : dict
        Hungarian matches: calc_idx -> (ref_idx | None, overlap).

    Returns
    -------
    (collapsed_matrix, calc_labels, ref_labels) where collapsed_matrix has
    one row per calc entry and one column per ref entry.
    """
    n_calc = alignment_matrix.shape[0]
    n_ref = alignment_matrix.shape[1]

    # Build ref entries: (indices_list, label)
    ref_entries: list[tuple[list[int], str]] = []
    used_ref: set[int] = set()
    for g in groups:
        ref_entries.append((g.ref_indices, f"{g.symmetry_label} @{g.center_frequency:.0f}"))
        used_ref.update(g.ref_indices)
    for j in range(n_ref):
        if j not in used_ref:
            ref_entries.append(([j], f"{freqs_ref[j]:.0f}"))

    # Build calc entries: (indices_list, label)
    calc_entries: list[tuple[list[int], str]] = []
    used_calc: set[int] = set()
    for g in groups:
        if g.calc_indices:
            calc_entries.append(
                (
                    g.calc_indices,
                    f"{g.symmetry_label} @{g.center_frequency:.0f}",
                )
            )
            used_calc.update(g.calc_indices)
    for i in range(n_calc):
        if i not in used_calc:
            calc_entries.append(([i], f"{freqs_calc[i]:.0f}"))

    # Build collapsed matrix
    n_rows = len(calc_entries)
    n_cols = len(ref_entries)
    collapsed = np.zeros((n_rows, n_cols))

    for ri, (ref_idxs, _) in enumerate(ref_entries):
        for ci, (calc_idxs, _) in enumerate(calc_entries):
            is_ref_group = len(ref_idxs) > 1
            is_calc_group = len(calc_idxs) > 1

            if is_ref_group and is_calc_group:
                # Both grouped: use precomputed subspace overlap only for
                # matching groups (same ref indices); cross-group cells
                # use max of sub-block
                ref_set = set(ref_idxs)
                calc_set = set(calc_idxs)
                matched_group = None
                for g in groups:
                    if set(g.ref_indices) == ref_set and set(g.calc_indices) == calc_set:
                        matched_group = g
                        break
                if matched_group is not None:
                    collapsed[ci, ri] = matched_group.subspace_overlap
                else:
                    sub = alignment_matrix[np.ix_(calc_idxs, ref_idxs)]
                    collapsed[ci, ri] = float(np.max(sub))
            elif not is_ref_group and not is_calc_group:
                # Both non-degenerate: copy from original
                collapsed[ci, ri] = alignment_matrix[calc_idxs[0], ref_idxs[0]]
            else:
                # Cross: one side group, other not -- use max of sub-block
                sub = alignment_matrix[np.ix_(calc_idxs, ref_idxs)]
                collapsed[ci, ri] = float(np.max(sub))

    calc_labels = [label for _, label in calc_entries]
    ref_labels = [label for _, label in ref_entries]

    return collapsed, calc_labels, ref_labels


# Example usage
if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) < 3:
        print("Usage: python mode_matching.py <calc.chk/.fchk> <ref.chk/.fchk> [output.png]")
        print("\nOptional: Provide output filename to save heatmap")
        sys.exit(1)

    calc_file = sys.argv[1]
    ref_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    print("Extracting modes from calculation...")
    modes_calc, freqs_calc, coords_calc, masses_calc, n_atoms_calc = (
        extract_mode_data_from_checkpoint(calc_file)
    )

    print("Extracting modes from reference...")
    modes_ref, freqs_ref, coords_ref, masses_ref, n_atoms_ref = extract_mode_data_from_checkpoint(
        ref_file
    )

    if n_atoms_calc != n_atoms_ref:
        print(f"ERROR: Different number of atoms ({n_atoms_calc} vs {n_atoms_ref})")
        sys.exit(1)

    print(f"\nMatching {modes_calc.shape[0]} modes...")
    matches = match_modes(modes_calc, modes_ref)

    print("\n" + "=" * 60)
    print("MODE MATCHING RESULTS")
    print("=" * 60)
    print(f"{'Calc Mode':<12} {'Ref Mode':<12} {'Overlap':<12}")
    print("-" * 60)

    for calc_idx in sorted(matches.keys()):
        ref_idx, overlap = matches[calc_idx]
        print(f"{calc_idx:<12} {ref_idx:<12} {overlap:>11.4f}")

    print("=" * 60)

    # Create alignment matrix and plot heatmap
    print("\nCreating overlap heatmap...")
    alignment_matrix = create_alignment_matrix(modes_calc, modes_ref)

    # Determine output filename
    if output_file is None:
        calc_name = Path(calc_file).stem
        ref_name = Path(ref_file).stem
        output_file = f"mode_overlap_{calc_name}_vs_{ref_name}.png"

    plot_mode_overlap_heatmap(
        alignment_matrix,
        output_file=output_file,
        calc_label="ML Calculation",
        ref_label="DFT Reference",
        freqs_calc=freqs_calc,
        freqs_ref=freqs_ref,
        matches=matches,
    )

    print(f"\n✓ Heatmap saved to: {output_file}")
