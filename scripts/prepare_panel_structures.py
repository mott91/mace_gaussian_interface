#!/usr/bin/env python3
"""Starting structures for the 2026 supervisor panel (thesis/MOLECULES.md).

Fetches each molecule's 3D record from PubChem, sets the conformer-defining torsions with
RDKit, relaxes the rest with MMFF94 while those torsions are held, and writes
``molecules/panel_2026/<name>.xyz`` plus ``molecules/panel_2026.txt`` (batch list).

Conformer conventions (Methods: B3LYP is the reference, so what matters is that ML and DFT
start from the *same* conformer, not that it is the global minimum):

- chains all-anti (C-C-C-C and C-C-C-O 180 deg, C-C-O-H 180 deg),
- carboxyl groups syn (O=C-O-H 0 deg, the Z form), the alpha C-C bond eclipsing C=O,
- named conformers: cis-HCOOH (O=C-O-H 180), carbonic acid cis-cis, glycine Ip,
  H2O2 skew (H-O-O-H 115, never the planar saddle),
- beta-D-glucopyranose as PubChem gives it (ring pucker printed for a manual check).

Usage::

    python scripts/prepare_panel_structures.py            # all, skip existing files
    python scripts/prepare_panel_structures.py --force    # refetch everything
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "molecules" / "panel_2026"
SDF_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF?record_type=3d"

ANTI_CHAIN = [
    ("[#6]-[#6]-[#6]-[#6]", 180.0),
    ("[#6]-[#6]-[#6]-[#8;X2;H1]", 180.0),
    ("[#6]-[#6]-[#8]-[#1]", 180.0),
]
CARBOXYL_SYN = [
    ("[#6]-[#6]-[#6]=[#8]", 0.0),  # alpha C-C eclipses C=O
    ("[#8]=[#6]-[#8]-[#1]", 0.0),  # Z (syn) COOH
]


@dataclass
class Entry:
    name: str  # file stem
    query: str | None  # PubChem name; None = built in code
    torsions: list[tuple[str, float]] = field(default_factory=list)
    first_only: set[str] = field(default_factory=set)  # SMARTS applied to one match only
    group: str = ""


PANEL = [
    # carboxylic acids C1-C6
    Entry("formic_acid_trans", "formic acid", [("[#8]=[#6]-[#8]-[#1]", 0.0)], group="acids"),
    Entry("formic_acid_cis", "formic acid", [("[#8]=[#6]-[#8]-[#1]", 180.0)], group="acids"),
    Entry("acetic_acid", "acetic acid", CARBOXYL_SYN, group="acids"),
    Entry("propionic_acid", "propionic acid", ANTI_CHAIN + CARBOXYL_SYN, group="acids"),
    Entry("butyric_acid", "butyric acid", ANTI_CHAIN + CARBOXYL_SYN, group="acids"),
    Entry("valeric_acid", "valeric acid", ANTI_CHAIN + CARBOXYL_SYN, group="acids"),
    Entry("caproic_acid", "hexanoic acid", ANTI_CHAIN + CARBOXYL_SYN, group="acids"),
    # alcohols
    Entry("methanol", "methanol", group="alcohols"),
    Entry("ethanol", "ethanol", ANTI_CHAIN, group="alcohols"),
    Entry("1-propanol", "1-propanol", ANTI_CHAIN, group="alcohols"),
    Entry("1-butanol", "1-butanol", ANTI_CHAIN, group="alcohols"),
    Entry(
        "isopropanol",
        "2-propanol",
        [("[#6]-[#6]-[#8]-[#1]", 180.0)],
        first_only={"[#6]-[#6]-[#8]-[#1]"},
        group="alcohols",
    ),
    Entry(
        "tert-butanol",
        "tert-butanol",
        [("[#6]-[#6]-[#8]-[#1]", 180.0)],
        first_only={"[#6]-[#6]-[#8]-[#1]"},
        group="alcohols",
    ),
    # in-house VCI set (methanol, formic acid above)
    Entry("carbonic_acid", "carbonic acid", [("[#8]=[#6]-[#8]-[#1]", 0.0)], group="VCI"),
    Entry("carbon_dioxide", "carbon dioxide", group="VCI"),
    # biomolecules
    Entry(
        "glycine",
        "glycine",
        [("[#7]-[#6]-[#6]=[#8]", 0.0), ("[#8]=[#6]-[#8]-[#1]", 0.0), ("[#1]-[#7]-[#6]-[#6]", 60.0)],
        first_only={"[#1]-[#7]-[#6]-[#6]"},  # rotating one N-H turns the whole NH2 group
        group="bio",
    ),
    Entry("beta_glucose", "beta-D-glucopyranose", group="bio"),
    # inorganic
    Entry("sulfur_dioxide", "sulfur dioxide", group="inorganic"),
    Entry("phosphine", "phosphine", group="inorganic"),
    Entry("ammonia", "ammonia", group="inorganic"),
    Entry(
        "hydrogen_peroxide",
        "hydrogen peroxide",
        [("[#1]-[#8]-[#8]-[#1]", 115.0)],
        group="inorganic",
    ),
    Entry("nitric_acid", "nitric acid", group="inorganic"),
    # ML only
    Entry("c60", None, group="ML only"),
]


def fetch_mol(query: str) -> Chem.Mol:
    r = requests.get(SDF_URL.format(name=query), timeout=30)
    r.raise_for_status()
    mol = Chem.MolFromMolBlock(r.text, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read the PubChem record for {query!r}")
    return mol


def set_torsions(mol: Chem.Mol, entry: Entry) -> list[tuple[tuple[int, ...], float]]:
    """Set every matching torsion; return the (atoms, target) list for constraints."""
    conf = mol.GetConformer()
    held = []
    for smarts, angle in entry.torsions:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        if smarts in entry.first_only:
            matches = matches[:1]
        for m in matches:
            rdMolTransforms.SetDihedralDeg(conf, *m, angle)
            held.append((m, angle))
    return held


def relax(mol: Chem.Mol, held) -> None:
    """MMFF94 relaxation with the conformer-defining torsions held fixed."""
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return  # no MMFF parameters (e.g. SO2, PH3): keep the PubChem geometry
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    for (i, j, k, l), angle in held:  # noqa: E741
        ff.MMFFAddTorsionConstraint(i, j, k, l, False, angle, angle, 1e4)
    ff.Minimize(maxIts=2000)


def write_xyz(mol: Chem.Mol, path: Path, comment: str) -> None:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    path.write_text("\n".join(lines) + "\n")


def report(mol: Chem.Mol, held) -> str:
    conf = mol.GetConformer()
    seen, parts = set(), []
    for m, target in held:
        key = tuple(sorted((m[1], m[2])))
        if key in seen:
            continue
        seen.add(key)
        sym = "-".join(mol.GetAtomWithIdx(i).GetSymbol() for i in m)
        parts.append(f"{sym} {rdMolTransforms.GetDihedralDeg(conf, *m):6.1f} (target {target:g})")
    return "; ".join(parts) if parts else "as fetched"


def glucose_pucker(mol: Chem.Mol) -> str:
    """Ring O-C1-C2-C3 / C3-C4-C5-O signs: 4C1 has the alternating chair pattern with C4 up,
    C1 down. Printed for a manual check rather than decided here."""
    ring = mol.GetRingInfo().AtomRings()
    six = [r for r in ring if len(r) == 6]
    if not six:
        return "no six-ring"
    r = six[0]
    conf = mol.GetConformer()
    dih = [
        rdMolTransforms.GetDihedralDeg(conf, r[i], r[(i + 1) % 6], r[(i + 2) % 6], r[(i + 3) % 6])
        for i in range(6)
    ]
    return "ring torsions " + " ".join(f"{d:+.0f}" for d in dih)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    written, failed = [], []
    for e in PANEL:
        if args.only and e.name not in args.only:
            continue
        path = OUT / f"{e.name}.xyz"
        if path.exists() and not args.force:
            print(f"  skip  {e.name:20s} (exists)")
            written.append(path)
            continue
        try:
            if e.query is None:
                from ase.build import molecule
                from ase.io import write

                write(path, molecule("C60"), format="xyz", comment="C60 (ASE g2 geometry)")
                print(f"  ok    {e.name:20s} ASE built-in")
            else:
                mol = fetch_mol(e.query)
                held = set_torsions(mol, e)
                relax(mol, held)
                write_xyz(mol, path, f"{e.name}: PubChem '{e.query}', {report(mol, held)}")
                extra = f"  [{glucose_pucker(mol)}]" if e.name == "beta_glucose" else ""
                print(f"  ok    {e.name:20s} {report(mol, held)}{extra}")
            written.append(path)
        except Exception as exc:  # keep going; list failures at the end
            print(f"  FAIL  {e.name:20s} {exc!r}")
            failed.append(e.name)

    (REPO / "molecules" / "panel_2026.txt").write_text(
        "\n".join(str(p.relative_to(REPO)) for p in written) + "\n"
    )
    print(
        f"\n{len(written)} structures -> {OUT.relative_to(REPO)}; list: molecules/panel_2026.txt"
    )
    if failed:
        print("failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
