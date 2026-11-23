"""LaTeX for the appendix table of the reference fundamentals (figures-ch4-tweaks, item 6;
figures-ch4-last-fixes, item 4: footnote letters in order of appearance, values as the sources print
them, mode descriptions in one style, Jacox cited).

  Table 1  measured gas-phase fundamentals of the seven experiment molecules: every band of
           mace_gaussian/analysis/band_origins.json for them; footnote marks from
           analysis/followup3/out/experiment_bands.csv (the bands and pairings of the experiment
           comparison) and from the CO2 nu1 rule of thesis_stats.section_experiment
  Table 2  VCI fundamentals on CCSD(T)-F12 surfaces (thesis/data/vci_fundamentals.csv); the six bands
           flagged in thesis/data/vci_comparison_campaign2026.csv are marked, the other 30 enter the
           coupled-cluster comparison

The blocks are two tables; Table 1 (58 rows) is longer than a page and comes in two parts (\\ContinuedFloat).

    python analysis/results_figures/reference_fundamentals_table.py [out.tex]
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
THESIS = Path.home() / "thesis"
ORIGINS = ROOT / "mace_gaussian" / "analysis" / "band_origins.json"
EXP_BANDS = ROOT / "analysis" / "followup3" / "out" / "experiment_bands.csv"
VCI = THESIS / "data" / "vci_fundamentals.csv"
VCI_CMP = THESIS / "data" / "vci_comparison_campaign2026.csv"

NAME = {"formic_acid_trans": r"\textit{trans}-Formic acid", "acetic_acid": "Acetic acid", "methanol": "Methanol",
        "carbonic_acid": "Carbonic acid", "carbon_dioxide": "CO$_2$", "ammonia": "NH$_3$", "phosphine": "PH$_3$",
        "nitric_acid": "Nitric acid"}  # fmt: skip
# order of the appendix tables: acids, alcohols, small molecules in the order of Table 3.3
EXP_ORDER = ["formic_acid_trans", "acetic_acid", "methanol", "carbon_dioxide", "ammonia", "phosphine", "nitric_acid"]
VCI_ORDER = ["formic_acid_trans", "methanol", "carbon_dioxide", "carbonic_acid"]
# short source keys; the caption resolves them
SOURCE_KEY = [("Shimanouchi", "Sh"), ("Jacox", "Ja"), ("PCCP 28, 19068", "Di26"), ("PCCP 22, 17932", "Di20")]
VCI_KEY = {"Dinu 2024": "Di24", "Dinu 2026": "Di26", "Dinu 2020": "Di20", "Schlagin 2024": "Sc24"}

# Values as the sources print them. band_origins.json stores numbers, so trailing zeros are lost there
# (580.30 -> 580.3). Shimanouchi prints integers (WebBook "Selected Freq."), Dinu one decimal (the json
# keeps that form; checked below), the Jacox compilation as below: the nine gas-phase values of the WebBook
# HNO3 page (ID=C7697372, Mask=800), copied 2026-11-23 and checked against the json numbers.
JACOX_PRINTED = {"ν1": "3550.0", "ν2": "1709.57", "ν3": "1325.74", "ν4": "1303.52", "ν5": "879.11",
                 "ν6": "646.83", "ν7": "580.30", "ν8": "763.15", "ν9": "458.23"}  # fmt: skip

# Mode descriptions in one style (figures-ch4-last-fixes 4d): full words in lower case, bonds with an en
# dash (O--H, C--H, C--O, C--C, N--O), C=O and the groups (CH3, CH2, NO2, ONO2, COH, OCO, CCO) as they
# are; the content unchanged. Abbreviations: Shimanouchi str = stretch, deform = deformation, d- =
# degenerate, s- = symmetric, op- = out-of-plane; Sym/Deg = symmetric/degenerate; asym/anti =
# antisymmetric; Jacox "a-stretch" = antisymmetric stretch [I: no key on the WebBook page; usual sense];
# Jacox "OPLA" has no key on the WebBook page and stays. Every description must be in this map.
DESC = {
    # Shimanouchi (acetic acid)
    "OH str": "O--H stretch", "CH3 d-str": "CH3 degenerate stretch", "CH3 s-str": "CH3 symmetric stretch",
    "C=O str": "C=O stretch", "CH3 d-deform": "CH3 degenerate deformation",
    "CH3 s-deform": "CH3 symmetric deformation", "OH bend": "O--H bend", "C-O str": "C--O stretch",
    "CH3 rock": "CH3 rock", "CC str": "C--C stretch", "OCO deform": "OCO deformation",
    "CCO deform": "CCO deformation", "C=O op-bend": "C=O out-of-plane bend", "C-O torsion": "C--O torsion",
    # Shimanouchi (methanol, NH3, PH3)
    "OH stretch": "O--H stretch", "CH3 d-stretch": "CH3 degenerate stretch",
    "CH3 s-stretch": "CH3 symmetric stretch", "CO stretch": "C--O stretch", "torsion": "torsion",
    "symmetric stretch": "symmetric stretch", "symmetric deformation": "symmetric deformation",
    "degenerate stretch": "degenerate stretch", "degenerate deformation": "degenerate deformation",
    "Sym str": "symmetric stretch", "Sym deform": "symmetric deformation", "Deg str": "degenerate stretch",
    "Deg deform": "degenerate deformation",
    # Dinu (formic acid, CO2; experiment and VCI)
    "CH stretch": "C--H stretch", "C=O stretch": "C=O stretch", "CH in-plane bend": "C--H in-plane bend",
    "COH in-plane bend": "COH in-plane bend", "C-O stretch": "C--O stretch",
    "OCO in-plane bend": "OCO in-plane bend", "CH out-of-plane bend": "C--H out-of-plane bend",
    "COH out-of-plane bend": "COH out-of-plane bend", "sym stretch": "symmetric stretch", "bend": "bend",
    "asym stretch": "antisymmetric stretch",
    # Jacox (HNO3)
    "NO2 a-stretch": "NO2 antisymmetric stretch", "Mixed": "mixed", "ON stretch": "N--O stretch",
    "NO2 scissors": "NO2 scissors", "NO2 rock": "NO2 rock", "ONO2 OPLA": "ONO2 OPLA", "Torsion": "torsion",
    # Dinu 2024 (methanol VCI)
    "CH3 stretch": "CH3 stretch", "CH2 asym stretch": "CH2 antisymmetric stretch",
    "CH3 sym stretch": "CH3 symmetric stretch", "CH3 deformation": "CH3 deformation",
    "CH2 asym deformation": "CH2 antisymmetric deformation", "CH3 sym deformation": "CH3 symmetric deformation",
    "OH in-plane bend": "O--H in-plane bend", "CH3 asym rock": "CH3 antisymmetric rock",
    "CH3 rock / CO": "CH3 rock / C--O", "OH torsion": "O--H torsion",
    # Schlagin 2024 (carbonic acid VCI)
    "OH sym stretch": "O--H symmetric stretch", "OH anti stretch": "O--H antisymmetric stretch",
    "C-O anti stretch": "C--O antisymmetric stretch", "OH in-plane bend sym": "O--H symmetric in-plane bend",
    "OH in-plane bend anti": "O--H antisymmetric in-plane bend", "C-O sym stretch": "C--O symmetric stretch",
    "OCO out-of-plane": "OCO out-of-plane", "C=O in-plane bend": "C=O in-plane bend",
    "OH out-of-plane sym": "O--H symmetric out-of-plane", "C-O in-plane bend": "C--O in-plane bend",
    "OH out-of-plane anti": "O--H antisymmetric out-of-plane",
}  # fmt: skip
DESC_USED: dict[str, str] = {}


def nu(label: str) -> str:
    return f"$\\nu_{{{int(re.sub(r'[^0-9]', '', label))}}}$"


def num(label: str) -> int:
    return int(re.sub(r"[^0-9]", "", label))


def tex_desc(s: str) -> str:
    new = DESC[s]
    DESC_USED[s] = new
    new = re.sub(r"\b(CH|NO|ONO)(\d)", lambda m: f"{m.group(1)}$_{m.group(2)}$", new)
    return new.replace("&", r"\&")


def printed(stored: str, label: str, src: str) -> str:
    """The value as the source prints it; stored = the json number as written in the file."""
    v = float(stored)
    if src == "Sh":
        assert v == int(v), stored
        return f"{v:.0f}"
    if src == "Ja":
        assert float(JACOX_PRINTED[label]) == v, (label, stored)
        return JACOX_PRINTED[label]
    assert re.fullmatch(r"\d+\.\d", stored), stored  # Dinu: one decimal
    return stored


def source_key(text: str) -> str:
    return next(k for pat, k in SOURCE_KEY if pat in text)


def experiment_block() -> tuple[dict[str, list[tuple[str, str]]], dict[str, int]]:
    """Rows per molecule as (footnote mark, LaTeX row)."""
    origins = json.loads(ORIGINS.read_text(), parse_float=str, parse_int=str)  # numbers as written
    e = pd.read_csv(EXP_BANDS).set_index(["molecule", "label"])
    rows, counts = {}, {"bands": 0, "compared": 0}
    for mol in EXP_ORDER:
        rows[mol] = []
        src = source_key(origins[mol]["source"])
        bands = sorted(origins[mol]["bands"], key=lambda b: num(b["label"]))
        for i, b in enumerate(bands):
            r = e.loc[(mol, b["label"])]
            assert float(r.exp) == float(b["freq_cm"]), (mol, b["label"])
            a = str(r.assignment)
            if mol in ("ammonia", "phosphine"):
                mark = "symtop"  # reported separately (thesis_stats.SYMTOP)
            elif a.startswith("excluded"):
                mark = next(r for k, r in (("gas-phase value missing", "liquid"), ("NIST: value undefined", "undefined"),
                                           ("a' band without", "nopartner")) if k in a)  # fmt: skip
            elif mol == "carbon_dioxide" and b["label"] == "ν1":
                mark = "fermi"  # Fermi dyad, dropped as co2_nu1 in section_experiment
            elif mol == "acetic_acid" and b["label"] == "ν7":
                mark = "pairing"  # in the comparison; uncertain pairing
            else:
                mark = ""
            counts["bands"] += 1
            counts["compared"] += mark in ("", "pairing")
            sup = f"\\textsuperscript{{@{mark}@}}" if mark else ""
            rows[mol].append((mark, f"    {NAME[mol] if i == 0 else '':28s} & {nu(b['label'])}{sup} & "
                              f"{tex_desc(b['description'])} & {printed(b['freq_cm'], b['label'], src)} & "
                              f"{src if i == 0 else ''} \\\\"))  # fmt: skip
    return rows, counts


def vci_block() -> tuple[dict[str, list[tuple[str, str]]], dict[str, int]]:
    v = pd.read_csv(VCI, dtype={"vci": str})  # values as written (one decimal, as in the sources)
    c = pd.read_csv(VCI_CMP).set_index(["molecule", "label"])
    rows, counts = {}, {"bands": 0, "compared": 0}
    for mol in VCI_ORDER:
        rows[mol] = []
        d = v[v.molecule == mol].assign(n=lambda x: x.label.map(num)).sort_values("n")
        for i, (_, b) in enumerate(d.iterrows()):
            r = c.loc[(mol, b.label)]
            assert float(r.vci) == float(b.vci), (mol, b.label)
            flag = "" if pd.isna(r.flag) else str(r.flag)
            mark = "" if not flag else ("torsion" if flag.startswith("torsion") else "resonance")
            assert re.fullmatch(r"\d+\.\d", b.vci), b.vci
            counts["bands"] += 1
            counts["compared"] += not mark
            src = VCI_KEY[next(k for k in VCI_KEY if b.source.startswith(k))]
            sup = f"\\textsuperscript{{@{mark}@}}" if mark else ""
            rows[mol].append((mark, f"    {NAME[mol] if i == 0 else '':28s} & {nu(b.label)}{sup} & "
                              f"{tex_desc(b.description)} & {b.vci} & {src if i == 0 else ''} \\\\"))
    return rows, counts


HEAD = r"""  \footnotesize
  \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllrl@{}}
    \toprule
    Molecule & Band & Mode & %s / cm$^{-1}$ & Source \\
    \midrule"""
FOOT = r"""    \bottomrule
  \end{tabular*}"""

SOURCES = {
    "Sh": r"Sh, Shimanouchi\cite{shimanouchi1972} via the NIST WebBook\cite{nist_webbook}",
    "Ja": r"Ja, Jacox\cite{jacox1998} via the NIST WebBook\cite{nist_webbook}",
    "Di20": r"Di20, Dinu \textit{et al.}\cite{dinu2020co2ch4}",
    "Di24": r"Di24, Dinu \textit{et al.}\cite{dinu2024methanol}",
    "Di26": r"Di26, Dinu \textit{et al.}\cite{dinu2026formicI}",
    "Sc24": r"Sc24, Schlagin \textit{et al.}\cite{schlagin2024carbonic}",
}
# footnote reasons; the letters are assigned in the order of first appearance (letters())
EXP_NOTES = {
    "symtop": r"NH$_3$ and PH$_3$: reported separately, not in the 46 bands.",
    "liquid": r"Not compared: no gas-phase value in the source; 1165\,cm$^{-1}$ is the liquid value.",
    "undefined": r"Not compared: the source gives no defined gas-phase value (strong coupling of internal and "
                 r"overall rotation).",
    "nopartner": r"Not compared: no B3LYP a$'$ mode to pair it with.",
    "fermi": r"Not compared: Fermi dyad with 2$\nu_2$.",
    "pairing": r"In the comparison, but the pairing is uncertain: assigned to the B3LYP band at 1338\,cm$^{-1}$.",
}  # fmt: skip
VCI_NOTES = {
    "torsion": "Not compared: large-amplitude torsion.",
    "resonance": "Not compared: in a resonance (polyad, Fermi dyad or doublet) in the source.",
}


def letters(rows: dict, mols: list[str]) -> dict[str, str]:
    """Footnote letter of every reason, in the order in which the reasons first appear."""
    order = list(dict.fromkeys(m for mol in mols for m, _ in rows[mol] if m))
    return {m: "abcdefghij"[i] for i, m in enumerate(order)}


# Table 1 (58 rows) is about 200 pt taller than a page: two parts of one table (\ContinuedFloat)
EXP_PARTS = [["formic_acid_trans", "acetic_acid", "methanol"], ["carbon_dioxide", "ammonia", "phosphine", "nitric_acid"]]


def table(intro: str, label: str | None, col: str, rows: dict, mols: list[str], notes: dict, cont: bool,
          let: dict[str, str]) -> str:
    marks = list(dict.fromkeys(m for mol in mols for m, _ in rows[mol] if m))
    marks.sort(key=lambda m: let[m])
    # sources in the order of their first row
    srcs = list(dict.fromkeys(k for mol in mols for _, r in rows[mol] for k in SOURCES if r.endswith(f"& {k} \\\\")))
    cap = (intro + " Sources: " + "; ".join(SOURCES[k] for k in srcs) + ". "
           + " ".join(f"\\textsuperscript{{{let[m]}}}~{notes[m]}" for m in marks))  # fmt: skip
    body = "\n    \\addlinespace\n".join("\n".join(r for _, r in rows[mol]) for mol in mols)
    for m, x in let.items():
        body = body.replace(f"@{m}@", x)
    assert "@" not in body
    return "\n".join([r"\begin{table}[htbp]", r"  \centering"] + ([r"  \ContinuedFloat"] if cont else [])
                     + [f"  \\caption{{{cap.strip()}}}"] + ([f"  \\label{{{label}}}"] if label else [])
                     + [HEAD % col, body, FOOT, r"\end{table}", ""])  # fmt: skip


def main(out: Path) -> None:
    exp_rows, ec = experiment_block()
    vci_rows, vc = vci_block()
    assert ec == {"bands": 58, "compared": 46}, ec
    assert vc == {"bands": 36, "compared": 30}, vc
    exp_intro = (r"Measured gas-phase fundamentals of the seven molecules of the experiment comparison, with the "
                 r"band labels of the dotted lines in the spectra and the mode descriptions of the sources. "
                 r"Bands without a mark enter the 46-band comparison.")  # fmt: skip
    head = (f"""% Generated by mace_gaussian/analysis/results_figures/reference_fundamentals_table.py
% (figures-ch4-tweaks, item 6) from mace_gaussian/analysis/band_origins.json,
% analysis/followup3/out/experiment_bands.csv, thesis/data/vci_fundamentals.csv and
% thesis/data/vci_comparison_campaign2026.csv. Do not edit by hand.
% Table 1: {ec['bands']} bands, {ec['compared']} of them in the experiment comparison; two parts
% (\\ContinuedFloat of the caption package, same table number), since 58 rows do not fit one page.
% Table 2: {vc['bands']} bands, {vc['compared']} of them in the coupled-cluster comparison.
% Footnote letters in the order of first appearance (table 1: part 1, then part 2); the VCI table
% starts again at a. Values as the sources print them; mode descriptions in one style (DESC).
""")  # fmt: skip
    let_exp = letters(exp_rows, EXP_PARTS[0] + EXP_PARTS[1])
    let_vci = letters(vci_rows, VCI_ORDER)
    print("footnotes, table 1:", let_exp, " table 2:", let_vci)
    tex = head + "\n".join([
        table(exp_intro, "tab:app_exp_fundamentals", "Exp.", exp_rows, EXP_PARTS[0], EXP_NOTES, False, let_exp),
        table(r"Measured gas-phase fundamentals (continued).", None, "Exp.", exp_rows, EXP_PARTS[1], EXP_NOTES, True,
              let_exp),
        table(r"VCI fundamentals on CCSD(T)-F12 potential energy surfaces, the reference of the coupled-cluster "
              r"comparison, with the mode descriptions of the sources. Bands without a mark enter the 30-band "
              r"comparison.", "tab:app_vci_fundamentals", "VCI", vci_rows, VCI_ORDER, VCI_NOTES, False, let_vci),
    ])  # fmt: skip
    for old_desc, new_desc in DESC_USED.items():
        if old_desc != new_desc:
            print(f"desc  {old_desc!r} -> {new_desc!r}")
    out.write_text(tex)
    print(f"wrote {out}: table 1 {ec}, table 2 {vc}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else THESIS / "docs" / "handoff" / "reference_fundamentals_table.tex")
