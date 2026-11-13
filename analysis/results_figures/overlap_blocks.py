"""Off-diagonal blocks of the mode-overlap matrices (Section 4.1, "Normal Modes"; Fig. 4.2).

Builds the same matrices as make_thesis_figures.fig_overlap_matrices (runs with MACE4IR dipoles,
MACE modes Kabsch-rotated into the B3LYP frame, |mass-weighted scalar product|), checks that the
Hungarian pairs and their overlaps equal the band table, and then describes the blocks:

- block: B3LYP modes that one MACE mode shares, i.e. connected through MACE modes with overlap
  >= THRESH to more than one B3LYP mode (connected components, size >= 2);
- mode character from the B3LYP eigenvector: share of the mass-weighted displacement (kinetic
  energy) on hydroxyl H, on C-H hydrogens and on the heavy atoms (C, O).

Also: subspace overlap of the degenerate B3LYP pairs (CO2 bend, NH3/PH3 E modes), where the
overlap of a single mode depends on the arbitrary basis inside the pair:
sqrt(sum of squared overlaps between the pair and its assigned MACE modes / 2).

python analysis/results_figures/overlap_blocks.py [molecule ...]   (default beta_glucose)
Run from a scratch directory (importing mace_gaussian writes .model.pt into the cwd).
"""

import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
CAMP = ROOT / "campaigns" / "2026"
sys.path.insert(0, str(ROOT))
logging.disable(logging.WARNING)
from mace_gaussian.analysis.mode_matching import (  # noqa: E402
    align_modes_to_reference,
    create_alignment_matrix,
    extract_mode_data_from_checkpoint,
)

spec = importlib.util.spec_from_file_location("thesis_stats", ROOT / "analysis" / "followup3" / "thesis_stats.py")
ts = importlib.util.module_from_spec(spec)
sys.modules["thesis_stats"] = ts
spec.loader.exec_module(ts)

ORDER = ["mace_anicc", "mace_omol", "mace_polar", "mace_off"]
THRESH = 0.25


def character(modes, coords, masses):
    """(n_modes, 3) shares of the mass-weighted displacement on hydroxyl H, C-H H, heavy atoms."""
    sym = np.round(masses).astype(int)
    d = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
    cls = np.full(len(sym), 2)  # heavy atom
    for a in np.where(sym == 1)[0]:
        partner = min((b for b in range(len(sym)) if b != a), key=lambda b: d[a, b])
        cls[a] = 0 if sym[partner] == 16 else 1
    q2 = (modes.reshape(len(modes), -1, 3) ** 2).sum(-1)
    q2 /= q2.sum(1, keepdims=True)
    return np.stack([q2[:, cls == k].sum(1) for k in range(3)], 1)


def components(edges, n):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return [sorted(g) for g in groups.values() if len(g) > 1]


def main(mols):
    _, _, fund, _, _ = ts.load()
    f = fund[fund.dipole_model == "mace_ml"]
    for mol in mols:
        base = CAMP / "comparison_results" / mol
        dft_fchk = next(iter(sorted(base.glob("b3lyp*/*.fchk"))))
        md, fd, xd, masses, _ = extract_mode_data_from_checkpoint(str(dft_fchk), force_harmonic=True)
        ch = character(np.asarray(md), np.asarray(xd), np.asarray(masses))
        lab = pd.read_csv(ts.HERE / "rotor_labels.csv")
        rot = {(x.molecule, int(x.i)) for x in lab.itertuples() if x.rotor}
        print(f"\n### {mol}: {len(fd)} modes, blocks = B3LYP modes linked by a MACE mode with overlap >= {THRESH}\n")
        all_blocks = {}
        for m in ORDER:
            fchk = base / f"{m}_mace_ml" / "gaussian_freq.fchk"
            if not fchk.exists():
                continue
            mm, _, xm, _, _ = extract_mode_data_from_checkpoint(str(fchk), force_harmonic=True)
            mm, _ = align_modes_to_reference(mm, xm, xd)
            ov = np.abs(create_alignment_matrix(mm, md))  # [mace, b3lyp]
            r, c = linear_sum_assignment(-ov)
            assign = dict(zip(r, c))  # mace -> b3lyp
            pair_ov = {int(b) + 1: ov[a, b] for a, b in zip(r, c)}
            # check against the band table (degenerate sets are averaged there: compare the minimum)
            bt = f[(f.molecule == mol) & (f.energy_model == m)]
            if bt.ndeg.max() == 1:
                diff = max(abs(pair_ov[int(i)] - o) for i, o in zip(bt.i, bt.overlap))
                assert diff < 0.005, (mol, m, diff)
            edges = [(b1, b2) for a in range(ov.shape[0]) for b1 in np.where(ov[a] >= THRESH)[0]
                     for b2 in np.where(ov[a] >= THRESH)[0] if b1 < b2]  # fmt: skip
            blocks = components(edges, len(fd))
            all_blocks[m] = blocks
            print(f"{m}: {len(blocks)} blocks (pairs equal the band table)")
            for g in blocks:
                off = max(ov[a, b] for a in range(ov.shape[0]) for b in g if assign[a] != b)
                share = ch[g].mean(0)
                nrot = sum((mol, i + 1) in rot for i in g)
                print(f"   modes {g[0] + 1}-{g[-1] + 1} ({len(g)}), {fd[g[0]]:.0f}-{fd[g[-1]]:.0f} cm-1, "
                      f"min pair overlap {min(pair_ov[i + 1] for i in g):.2f}, largest off-pair {off:.2f}; "
                      f"share OH-H {share[0]:.2f}, CH-H {share[1]:.2f}, heavy {share[2]:.2f}; rotors {nrot}")  # fmt: skip
        # the same blocks in every model?
        ms = list(all_blocks)
        sets = {m: {frozenset(g) for g in all_blocks[m]} for m in ms}
        common = set.intersection(*sets.values()) if sets else set()
        print(f"\nidentical blocks in all {len(ms)} models: {sorted((min(g) + 1, max(g) + 1) for g in common)}")
        cover = {m: set().union(*all_blocks[m]) if all_blocks[m] else set() for m in ms}
        both = set.intersection(*cover.values())
        print(f"modes inside a block in every model: {len(both)}; per model: "
              + ", ".join(f"{m} {len(cover[m])}" for m in ms))
        print("\nB3LYP mode character (OH-H / CH-H / heavy), modes 1-66 by frequency:")
        for i in range(len(fd)):
            print(f"   {i + 1:3d} {fd[i]:7.1f}  {ch[i, 0]:.2f} {ch[i, 1]:.2f} {ch[i, 2]:.2f}"
                  + ("  rotor" if (mol, i + 1) in rot else ""))  # fmt: skip


def degenerate(mols=("carbon_dioxide", "ammonia", "phosphine")):
    print("\n### degenerate pairs: single-mode overlap (band table) and subspace overlap\n")
    for mol in mols:
        base = CAMP / "comparison_results" / mol
        md, fd, xd, _, _ = extract_mode_data_from_checkpoint(str(sorted(base.glob("b3lyp*/*.fchk"))[0]), force_harmonic=True)
        fd = np.asarray(fd)
        groups = []
        for i in range(len(fd)):
            if i and fd[i] - fd[i - 1] <= ts.DEGENERATE_CM:
                groups[-1].append(i)
            else:
                groups.append([i])
        for m in ORDER:
            fchk = base / f"{m}_mace_ml" / "gaussian_freq.fchk"
            if not fchk.exists():
                continue
            mm, _, xm, _, _ = extract_mode_data_from_checkpoint(str(fchk), force_harmonic=True)
            mm, _ = align_modes_to_reference(mm, xm, xd)
            ov = np.abs(create_alignment_matrix(mm, md))
            r, c = linear_sum_assignment(-ov)
            inv = {b: a for a, b in zip(r, c)}
            for g in (g for g in groups if len(g) > 1):
                sub = np.sqrt((ov[np.ix_([inv[b] for b in g], g)] ** 2).sum() / len(g))
                single = np.mean([ov[inv[b], b] for b in g])
                print(f"   {mol:15s} {m:11s} modes {g[0] + 1}-{g[-1] + 1} ({fd[g[0]]:.0f} cm-1): "
                      f"single {single:.3f}, subspace {sub:.4f}")  # fmt: skip


if __name__ == "__main__":
    main(sys.argv[1:] or ["beta_glucose"])
    degenerate()
