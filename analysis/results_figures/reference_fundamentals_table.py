"""LaTeX for the appendix table of the reference fundamentals (figures-ch4-tweaks, item 6).

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


def nu(label: str) -> str:
    return f"$\\nu_{{{int(re.sub(r'[^0-9]', '', label))}}}$"


def num(label: str) -> int:
    return int(re.sub(r"[^0-9]", "", label))


def tex_desc(s: str) -> str:
    s = re.sub(r"\b(CH|NO|ONO)(\d)", lambda m: f"{m.group(1)}$_{m.group(2)}$", s)
    return s.replace("&", r"\&")


def fmt(v: float, src: str) -> str:
    # Shimanouchi gives integers; the other sources as stored
    return f"{v:.0f}" if src == "Sh" else f"{v}"


def source_key(text: str) -> str:
    return next(k for pat, k in SOURCE_KEY if pat in text)


def experiment_block() -> tuple[dict[str, list[tuple[str, str]]], dict[str, int]]:
    """Rows per molecule as (footnote mark, LaTeX row)."""
    origins = json.loads(ORIGINS.read_text())
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
                mark = "a"  # reported separately (thesis_stats.SYMTOP)
            elif a.startswith("excluded"):
                mark = {"gas-phase value missing": "b", "NIST: value undefined": "c", "a' band without": "d"}[
                    next(k for k in ("gas-phase value missing", "NIST: value undefined", "a' band without") if k in a)
                ]  # fmt: skip
            elif mol == "carbon_dioxide" and b["label"] == "ν1":
                mark = "e"  # Fermi dyad, dropped as co2_nu1 in section_experiment
            elif mol == "acetic_acid" and b["label"] == "ν7":
                mark = "f"  # in the comparison; uncertain pairing
            else:
                mark = ""
            counts["bands"] += 1
            counts["compared"] += mark in ("", "f")
            sup = f"\\textsuperscript{{{mark}}}" if mark else ""
            rows[mol].append((mark, f"    {NAME[mol] if i == 0 else '':28s} & {nu(b['label'])}{sup} & "
                              f"{tex_desc(b['description'])} & {fmt(b['freq_cm'], src)} & {src if i == 0 else ''} \\\\"))
    return rows, counts


def vci_block() -> tuple[dict[str, list[tuple[str, str]]], dict[str, int]]:
    v = pd.read_csv(VCI)
    c = pd.read_csv(VCI_CMP).set_index(["molecule", "label"])
    rows, counts = {}, {"bands": 0, "compared": 0}
    for mol in VCI_ORDER:
        rows[mol] = []
        d = v[v.molecule == mol].assign(n=lambda x: x.label.map(num)).sort_values("n")
        for i, (_, b) in enumerate(d.iterrows()):
            r = c.loc[(mol, b.label)]
            assert float(r.vci) == float(b.vci), (mol, b.label)
            flag = "" if pd.isna(r.flag) else str(r.flag)
            mark = "" if not flag else ("a" if flag.startswith("torsion") else "b")
            counts["bands"] += 1
            counts["compared"] += not mark
            src = VCI_KEY[next(k for k in VCI_KEY if b.source.startswith(k))]
            sup = f"\\textsuperscript{{{mark}}}" if mark else ""
            rows[mol].append((mark, f"    {NAME[mol] if i == 0 else '':28s} & {nu(b.label)}{sup} & "
                              f"{tex_desc(b.description)} & {b.vci:.1f} & {src if i == 0 else ''} \\\\"))
    return rows, counts


HEAD = r"""  \footnotesize
  \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllrl@{}}
    \toprule
    Molecule & Band & Mode (source) & %s / cm$^{-1}$ & Source \\
    \midrule"""
FOOT = r"""    \bottomrule
  \end{tabular*}"""

SOURCES = {
    "Sh": r"Sh, Shimanouchi\cite{shimanouchi1972} via the NIST WebBook\cite{nist_webbook}",
    "Ja": r"Ja, the Jacox compilation via the NIST WebBook\cite{nist_webbook}",
    "Di20": r"Di20, Dinu \textit{et al.}\cite{dinu2020co2ch4}",
    "Di24": r"Di24, Dinu \textit{et al.}\cite{dinu2024methanol}",
    "Di26": r"Di26, Dinu \textit{et al.}\cite{dinu2026formicI}",
    "Sc24": r"Sc24, Schlagin \textit{et al.}\cite{schlagin2024carbonic}",
}
EXP_NOTES = {
    "a": r"NH$_3$ and PH$_3$: reported separately, not in the 46 bands.",
    "b": r"Not compared: no gas-phase value in the source; 1165\,cm$^{-1}$ is the liquid value.",
    "c": r"Not compared: the source gives no defined gas-phase value (strong coupling of internal and "
         r"overall rotation).",
    "d": r"Not compared: no B3LYP a$'$ mode to pair it with.",
    "e": r"Not compared: Fermi dyad with 2$\nu_2$.",
    "f": r"In the comparison, but the pairing is uncertain: assigned to the B3LYP band at 1338\,cm$^{-1}$.",
}  # fmt: skip
VCI_NOTES = {
    "a": "Not compared: large-amplitude torsion.",
    "b": "Not compared: in a resonance (polyad, Fermi dyad or doublet) in the source.",
}
# Table 1 (58 rows) is about 200 pt taller than a page: two parts of one table (\ContinuedFloat)
EXP_PARTS = [["formic_acid_trans", "acetic_acid", "methanol"], ["carbon_dioxide", "ammonia", "phosphine", "nitric_acid"]]


def table(intro: str, label: str | None, col: str, rows: dict, mols: list[str], notes: dict, cont: bool) -> str:
    marks = sorted({m for mol in mols for m, _ in rows[mol] if m})
    srcs = [k for k in SOURCES if any(f"& {k} \\\\" in r for mol in mols for _, r in rows[mol])]
    cap = (intro + " Sources: " + "; ".join(SOURCES[k] for k in srcs) + ". "
           + " ".join(f"\\textsuperscript{{{m}}}~{notes[m]}" for m in marks))  # fmt: skip
    body = "\n    \\addlinespace\n".join("\n".join(r for _, r in rows[mol]) for mol in mols)
    return "\n".join([r"\begin{table}[htbp]", r"  \centering"] + ([r"  \ContinuedFloat"] if cont else [])
                     + [f"  \\caption{{{cap}}}"] + ([f"  \\label{{{label}}}"] if label else [])
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
% The Jacox compilation has no bib entry yet (see the report); its key Ja cites the WebBook.
""")  # fmt: skip
    tex = head + "\n".join([
        table(exp_intro, "tab:app_exp_fundamentals", "Exp.", exp_rows, EXP_PARTS[0], EXP_NOTES, False),
        table(r"Measured gas-phase fundamentals (continued).", None, "Exp.", exp_rows, EXP_PARTS[1], EXP_NOTES, True),
        table(r"VCI fundamentals on CCSD(T)-F12 potential energy surfaces, the reference of the coupled-cluster "
              r"comparison, with the mode descriptions of the sources. Bands without a mark enter the 30-band "
              r"comparison.", "tab:app_vci_fundamentals", "VCI", vci_rows, VCI_ORDER, VCI_NOTES, False),
    ])  # fmt: skip
    out.write_text(tex)
    print(f"wrote {out}: table 1 {ec}, table 2 {vc}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else THESIS / "docs" / "handoff" / "reference_fundamentals_table.tex")
