"""Per-mode master table: every method's number for every DFT normal mode in one row.

The rest of the report is organized per method (one section per ML model). The
thesis needs the transpose: one row per DFT fundamental, with the DFT harmonic and
VPT2 frequencies, each ML model's harmonic and VPT2 frequency for the *same* normal
mode (via the eigenvector mapping the analysis already computed), and, where a
curated value exists, the experimental band origin.

Inputs come from the ``analysis_results`` dict produced by
``ComparisonWorkflow.run_full_analysis``. Each comparison carries

- ``dft_spectrum`` / ``ml_spectrum``: ``SpectrumData`` whose fundamental mode ids
  ``F{k}`` live in the checkpoint (ascending-frequency) index space of their own
  calculation,
- ``mode_mapping``: ML checkpoint index -> DFT checkpoint index (0-based), from
  ``extract_mode_mapping``; ``mode_overlaps`` gives the eigenvector overlap,
- ``_dft_results`` / ``_ml_results``: the raw results.json dicts, whose
  ``frequencies.harmonic`` list is in checkpoint order and supplies the harmonic
  column in anharmonic mode.

Experimental band origins are read from ``band_origins.json`` next to this module
(Shimanouchi 1972 via the NIST WebBook, see the file's ``_meta``). Bands are
assigned to DFT modes by nearest frequency (VPT2 when available, else harmonic),
one DFT mode per unit of degeneracy, unless the band lists explicit ``modes``.
"""

from __future__ import annotations

import csv
import html
import json
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

BAND_ORIGINS_PATH = Path(__file__).parent / "band_origins.json"

# Largest |nu_exp - nu_DFT| accepted by the automatic band assignment.
ASSIGNMENT_WINDOW_CM = 200.0

# Short labels for the LaTeX export; unknown names fall back to the escaped raw name.
_ENERGY_LABELS = {
    "mace_omol": "OMOL",
    "mace_off": "OFF",
    "mace_mp": "MP",
    "mace_anicc": "ANI-cc",
    "mace_polar": "POLAR",
}
_DIPOLE_LABELS = {
    "mace_ml": "",
    "mace_polar1": "/P1",
    "mace_mdp": "/MDP",
    "espaloma": "/esp",
    "xtb": "/xtb",
}


# Dipole calculators, most preferred first, used to pick one representative run per
# energy model for the frequency views (frequencies do not depend on the dipole model).
_DIPOLE_PREFERENCE = ("mace_ml", "mace_polar1", "mace_mdp", "espaloma", "xtb")


def split_method(method: str) -> tuple[str, str | None]:
    """``mace_off_espaloma`` -> ``("mace_off", "espaloma")``; unknown suffix -> ``(name, None)``."""
    for dip in _DIPOLE_PREFERENCE:
        if method.endswith("_" + dip):
            return method[: -len(dip) - 1], dip
    return method, None


def frequency_methods(methods: list[dict]) -> list[dict]:
    """One entry per energy model, keeping the run with the preferred dipole calculator.

    Harmonic and VPT2 frequencies come from the energy model alone, so the report's
    5 x 3 method grid collapses to five frequency columns. The label drops the
    dipole suffix.
    """
    best: dict[str, tuple[int, dict]] = {}
    order: list[str] = []
    for m in methods:
        energy, dip = split_method(m["name"])
        rank = len(_DIPOLE_PREFERENCE)
        if dip in _DIPOLE_PREFERENCE:
            rank = _DIPOLE_PREFERENCE.index(dip)
        if energy not in best:
            order.append(energy)
        if energy not in best or rank < best[energy][0]:
            best[energy] = (rank, {**m, "label": display_name(energy)})
    return [best[e][1] for e in order]


def display_name(method: str) -> str:
    """``mace_omol_mace_ml`` -> ``OMOL``, ``mace_off_espaloma`` -> ``OFF/esp``."""
    for energy, label in _ENERGY_LABELS.items():
        if method.startswith(energy + "_"):
            return label + _DIPOLE_LABELS.get(
                method[len(energy) + 1 :], "/" + method[len(energy) + 1 :]
            )
        if method == energy:
            return label
    return method


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------


def load_band_origins(molecule: str, path: Path = BAND_ORIGINS_PATH) -> dict | None:
    """Curated band-origin entry for ``molecule`` or None."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        logger.warning("Could not read band origins from %s: %s", path, e)
        return None
    entry = data.get(molecule)
    return entry if isinstance(entry, dict) and entry.get("bands") else None


def _fundamentals(spectrum) -> dict[int, tuple[float, float]]:
    """{checkpoint index (1-based): (frequency, intensity)} for the fundamentals."""
    out: dict[int, tuple[float, float]] = {}
    for i, (label, mode_id) in enumerate(zip(spectrum.labels, spectrum.mode_ids)):
        if label == "fundamental" and mode_id.startswith("F"):
            try:
                k = int(mode_id[1:])
            except ValueError:
                continue
            out[k] = (float(spectrum.frequencies[i]), float(spectrum.intensities[i]))
    return out


def _harmonic_list(results: dict | None) -> list[dict]:
    if not results:
        return []
    return list(results.get("frequencies", {}).get("harmonic", []) or [])


def _method_cells(comp: dict, mode: str, dft_modes: list[int]) -> tuple[dict[int, dict], str]:
    """Per-DFT-mode cells for one ML method, plus how ML modes were paired."""
    ml_fund = _fundamentals(comp["ml_spectrum"])
    ml_harm = _harmonic_list(comp.get("_ml_results"))
    overlaps = comp.get("mode_overlaps") or {}
    mapping = comp.get("mode_mapping")

    # dft index (1-based) -> ml index (1-based). With several ML modes claiming one
    # DFT mode, keep the best overlap.
    if mapping:
        pairing = "eigenvector"
        dft_to_ml: dict[int, int] = {}
        for ml_idx, dft_idx in mapping.items():
            k, j = dft_idx + 1, ml_idx + 1
            if k not in dft_to_ml or overlaps.get(ml_idx, 0) > overlaps.get(dft_to_ml[k] - 1, 0):
                dft_to_ml[k] = j
    else:
        pairing = "index"
        dft_to_ml = {k: k for k in dft_modes}

    cells: dict[int, dict] = {}
    for k in dft_modes:
        j = dft_to_ml.get(k)
        cell: dict[str, Any] = {
            "ml_mode": j,
            "harmonic": None,
            "vpt2": None,
            "intensity": None,
            "overlap": None,
        }
        if j is not None:
            if mode == "anharmonic":
                if j in ml_fund:
                    cell["vpt2"], cell["intensity"] = ml_fund[j]
                if 0 < j <= len(ml_harm):
                    cell["harmonic"] = float(ml_harm[j - 1]["freq_cm"])
            elif j in ml_fund:
                cell["harmonic"], cell["intensity"] = ml_fund[j]
            if mapping:
                ov = overlaps.get(j - 1)
                cell["overlap"] = float(ov) if ov is not None else None
        cells[k] = cell
    return cells, pairing


def assign_band_origins(
    bands: list[dict], reference: dict[int, float], window: float = ASSIGNMENT_WINDOW_CM
) -> dict[int, dict]:
    """Map curated bands onto DFT modes.

    Parameters
    ----------
    bands:
        Entries from ``band_origins.json``; ``degeneracy`` (default 1) is how many
        DFT modes a band may take, ``modes`` (1-based) forces the assignment.
    reference:
        {DFT mode index: reference frequency} used for nearest matching.

    Returns
    -------
    {DFT mode index: band dict + ``assignment`` ("explicit" | "nearest")}
    """
    taken: dict[int, dict] = {}
    candidates: list[tuple[float, int, int]] = []  # (distance, band index, dft mode)
    for b_i, band in enumerate(bands):
        if band.get("modes"):
            for k in band["modes"]:
                if k in reference and k not in taken:
                    taken[k] = {**band, "assignment": "explicit"}
            continue
        for k, nu in reference.items():
            d = abs(nu - float(band["freq_cm"]))
            if d <= window:
                candidates.append((d, b_i, k))

    remaining = {
        b_i: int(band.get("degeneracy", 1))
        for b_i, band in enumerate(bands)
        if not band.get("modes")
    }
    for d, b_i, k in sorted(candidates):
        if k in taken or remaining.get(b_i, 0) <= 0:
            continue
        taken[k] = {**bands[b_i], "assignment": "nearest", "distance_cm": d}
        remaining[b_i] -= 1
    return taken


def build_master_table(
    analysis_results: dict, band_origins: dict | None = None, band_origins_path: Path | None = None
) -> dict:
    """Pivot ``analysis_results`` into one row per DFT normal mode.

    Returns a JSON-serializable dict::

        {
          "molecule": str, "mode": "anharmonic" | "harmonic",
          "methods": [{"name", "label", "pairing"}, ...],
          "experimental": {"source": str, ...} | None,
          "rows": [
            {"dft_mode": k, "dft_harmonic": float | None, "dft_vpt2": float | None,
             "dft_intensity": float | None,
             "experimental": {"label", "symmetry", "description", "freq_cm", "rating",
                              "assignment"} | None,
             "methods": {name: {"ml_mode", "harmonic", "vpt2", "intensity", "overlap"}}},
            ...
          ]
        }
    """
    comparisons = analysis_results.get("comparisons") or []
    molecule = analysis_results.get("molecule", "")
    mode = analysis_results.get("mode", "anharmonic")
    table: dict[str, Any] = {
        "molecule": molecule,
        "mode": mode,
        "methods": [],
        "experimental": None,
        "rows": [],
    }
    if not comparisons:
        return table

    first = comparisons[0]
    dft_fund = _fundamentals(first["dft_spectrum"])
    dft_harm = _harmonic_list(first.get("_dft_results"))
    dft_modes = sorted(set(dft_fund) | set(range(1, len(dft_harm) + 1)))
    if not dft_modes:
        return table

    rows: dict[int, dict] = {}
    for k in dft_modes:
        harmonic = vpt2 = intensity = None
        if mode == "anharmonic":
            if k in dft_fund:
                vpt2, intensity = dft_fund[k]
            if 0 < k <= len(dft_harm):
                harmonic = float(dft_harm[k - 1]["freq_cm"])
        elif k in dft_fund:
            harmonic, intensity = dft_fund[k]
        rows[k] = {
            "dft_mode": k,
            "dft_harmonic": harmonic,
            "dft_vpt2": vpt2,
            "dft_intensity": intensity,
            "experimental": None,
            "methods": {},
        }

    for comp in comparisons:
        name = comp["name"]
        cells, pairing = _method_cells(comp, mode, dft_modes)
        table["methods"].append({"name": name, "label": display_name(name), "pairing": pairing})
        for k in dft_modes:
            rows[k]["methods"][name] = cells[k]

    if band_origins is None:
        band_origins = load_band_origins(molecule, band_origins_path or BAND_ORIGINS_PATH)
    if band_origins:
        reference = {
            k: (r["dft_vpt2"] if r["dft_vpt2"] is not None else r["dft_harmonic"])
            for k, r in rows.items()
        }
        reference = {k: v for k, v in reference.items() if v is not None}
        assigned = assign_band_origins(band_origins["bands"], reference)
        for k, band in assigned.items():
            rows[k]["experimental"] = {
                key: band.get(key)
                for key in (
                    "label",
                    "symmetry",
                    "description",
                    "freq_cm",
                    "rating",
                    "note",
                    "assignment",
                    "distance_cm",
                )
            }
        table["experimental"] = {
            "source": band_origins.get("source"),
            "cas": band_origins.get("cas"),
            "n_bands": len(band_origins["bands"]),
            "n_assigned": len(assigned),
            "assignment_window_cm": ASSIGNMENT_WINDOW_CM,
        }

    table["rows"] = [rows[k] for k in dft_modes]
    return table


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt(v: float | None, digits: int = 1) -> str:
    return "" if v is None else f"{v:.{digits}f}"


def _delta(v: float | None, ref: float | None, digits: int = 1) -> str:
    if v is None or ref is None:
        return ""
    return f"{v - ref:+.{digits}f}"


def render_master_table_html(table: dict) -> str:
    """HTML ``<section id="master-table">`` for the report."""
    rows = table.get("rows") or []
    methods = frequency_methods(table.get("methods") or [])
    anharm = table.get("mode") == "anharmonic"
    if not rows:
        return (
            '<section id="master-table"><h2>Per-mode master table</h2>'
            '<div class="warning-box">No modes available.</div></section>'
        )

    esc = html.escape
    freq_cols = ["harm", "VPT2"] if anharm else ["harm"]
    # Header row 1: groups; row 2: sub-columns.
    # Experiment columns only where a band origin exists for at least one mode
    has_exp = any(r.get("experimental") for r in rows)
    head1 = ['<th rowspan="2">Mode</th>']
    if has_exp:
        head1.append('<th colspan="2">Experiment</th>')
    head1.append(f'<th colspan="{len(freq_cols)}">DFT</th>')
    for m in methods:
        title = f"paired by {m['pairing']}"
        head1.append(
            f'<th colspan="{len(freq_cols) + 1}" title="{esc(title)}">{esc(m["label"])}</th>'
        )
    head2 = (["<th>band</th>", "<th>ν₀</th>"] if has_exp else []) + [
        f"<th>{c}</th>" for c in freq_cols
    ]
    for _ in methods:
        head2 += [f"<th>{c}</th>" for c in freq_cols] + [
            f"<th title='{esc(('VPT2' if anharm else 'harmonic') + ' minus DFT')}'>Δ</th>"
        ]

    body = []
    for r in rows:
        exp = r.get("experimental") or {}
        exp_label = exp.get("label") or ""
        if exp:
            tip = " ".join(
                str(x)
                for x in (
                    exp.get("symmetry"),
                    exp.get("description"),
                    f"rating {exp.get('rating')}",
                )
                if x
            )
            exp_label = f'<span title="{esc(tip)}">{esc(exp_label)}</span>'
        cells = [f"<td>{r['dft_mode']}</td>"]
        if has_exp:
            cells += [f"<td>{exp_label}</td>", f"<td>{_fmt(exp.get('freq_cm'))}</td>"]
        cells.append(f"<td>{_fmt(r['dft_harmonic'])}</td>")
        dft_ref = r["dft_harmonic"]
        if anharm:
            cells.append(f"<td>{_fmt(r['dft_vpt2'])}</td>")
            dft_ref = r["dft_vpt2"]
        for m in methods:
            c = r["methods"].get(m["name"]) or {}
            ov = c.get("overlap")
            style = ' class="metric-warning"' if ov is not None and ov < 0.7 else ""
            tip = f"ML mode {c.get('ml_mode')}" + (f", overlap {ov:.2f}" if ov is not None else "")
            cells.append(f'<td title="{esc(tip)}"{style}>{_fmt(c.get("harmonic"))}</td>')
            val = c.get("harmonic")
            if anharm:
                cells.append(f'<td title="{esc(tip)}"{style}>{_fmt(c.get("vpt2"))}</td>')
                val = c.get("vpt2")
            cells.append(f"<td{style}>{_delta(val, dft_ref)}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")

    exp_info = table.get("experimental")
    note = ""
    if exp_info and exp_info.get("source"):
        note = (
            f"<p>Experimental band origins: {esc(exp_info['source'])} "
            f"({exp_info['n_assigned']} of {exp_info['n_bands']} bands assigned, "
            f"nearest DFT {'VPT2' if anharm else 'harmonic'} frequency within "
            f"{exp_info['assignment_window_cm']:.0f} cm⁻¹).</p>"
        )
    else:
        note = "<p>No curated experimental band origins for this molecule.</p>"
    caption = (
        "One row per DFT normal mode (checkpoint index). One column group per energy model "
        "(frequencies do not depend on the dipole model; intensities per run are in "
        "master_table.csv). ML columns are paired to the DFT mode "
        "by eigenvector overlap; a highlighted cell has overlap below 0.7. Δ is "
        f"{'VPT2' if anharm else 'harmonic'} minus DFT {'VPT2' if anharm else 'harmonic'}, in cm⁻¹."
    )
    return (
        '<section id="master-table"><h2>Per-mode master table</h2>'
        f"<p>{caption}</p>{note}"
        '<div class="table-scroll" style="overflow:auto"><table class="data-table">'
        f"<thead><tr>{''.join(head1)}</tr><tr>{''.join(head2)}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div></section>"
    )


def write_master_table_csv(table: dict, out_path: Path) -> None:
    """Wide CSV: one row per DFT mode, one column group per method."""
    methods = [m["name"] for m in table.get("methods") or []]
    header = [
        "dft_mode",
        "exp_label",
        "exp_freq_cm",
        "exp_symmetry",
        "exp_description",
        "exp_rating",
        "dft_harmonic",
        "dft_vpt2",
        "dft_intensity",
    ]
    for m in methods:
        header += [f"{m}_ml_mode", f"{m}_harmonic", f"{m}_vpt2", f"{m}_intensity", f"{m}_overlap"]
    with Path(out_path).open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in table.get("rows") or []:
            exp = r.get("experimental") or {}
            row = [
                r["dft_mode"],
                exp.get("label", ""),
                exp.get("freq_cm", ""),
                exp.get("symmetry", ""),
                exp.get("description", ""),
                exp.get("rating", ""),
                _fmt(r["dft_harmonic"], 4),
                _fmt(r["dft_vpt2"], 4),
                _fmt(r["dft_intensity"], 4),
            ]
            for m in methods:
                c = r["methods"].get(m) or {}
                row += [
                    c.get("ml_mode", ""),
                    _fmt(c.get("harmonic"), 4),
                    _fmt(c.get("vpt2"), 4),
                    _fmt(c.get("intensity"), 4),
                    _fmt(c.get("overlap"), 3),
                ]
            w.writerow(row)


def _tex(s: Any) -> str:
    return (
        str(s)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("#", r"\#")
    )


def _tex_label(label: str) -> str:
    """``\u03bd3`` (Greek nu + digits) -> ``$\\nu_{3}$``; anything else is escaped verbatim."""
    m = re.fullmatch("\u03bd(\\d+)", label or "")
    return f"$\\nu_{{{m.group(1)}}}$" if m else _tex(label)


def write_master_table_latex(table: dict, out_path: Path, digits: int = 0) -> None:
    """booktabs ``tabular`` (no floating environment) for ``\\input`` in the thesis.

    Columns: mode, experimental band and origin, DFT harmonic and VPT2, then per
    method the VPT2 (or harmonic) frequency and its deviation from DFT. Cells with
    eigenvector overlap below 0.7 are marked with a dagger.
    """
    rows = table.get("rows") or []
    methods = frequency_methods(table.get("methods") or [])
    anharm = table.get("mode") == "anharmonic"
    n_dft = 2 if anharm else 1
    col_spec = "l l r " + "r " * n_dft + "".join("r r " for _ in methods)
    lines = [
        "% generated by mace_gaussian.analysis.master_table; do not edit by hand",
        f"% molecule: {table.get('molecule', '')}, mode: {table.get('mode', '')}",
        f"\\begin{{tabular}}{{{col_spec.strip()}}}",
        "\\toprule",
    ]
    head1 = ["", "\\multicolumn{2}{c}{Experiment}", f"\\multicolumn{{{n_dft}}}{{c}}{{DFT}}"]
    head1 += [f"\\multicolumn{{2}}{{c}}{{{_tex(m['label'])}}}" for m in methods]
    lines.append(" & ".join(head1) + " \\\\")
    cmid = ["\\cmidrule(lr){2-3}", f"\\cmidrule(lr){{4-{3 + n_dft}}}"]
    col = 4 + n_dft
    for _ in methods:
        cmid.append(f"\\cmidrule(lr){{{col}-{col + 1}}}")
        col += 2
    lines.append(" ".join(cmid))
    head2 = ["Mode", "band", "$\\nu_0$", "harm."] + (["VPT2"] if anharm else [])
    head2 += ["VPT2" if anharm else "harm.", "$\\Delta$"] * len(methods)
    lines.append(" & ".join(head2) + " \\\\")
    lines.append("\\midrule")
    for r in rows:
        exp = r.get("experimental") or {}
        cells = [
            str(r["dft_mode"]),
            _tex_label(exp.get("label", "")),
            _fmt(exp.get("freq_cm"), digits),
            _fmt(r["dft_harmonic"], digits),
        ]
        dft_ref = r["dft_harmonic"]
        if anharm:
            cells.append(_fmt(r["dft_vpt2"], digits))
            dft_ref = r["dft_vpt2"]
        for m in methods:
            c = r["methods"].get(m["name"]) or {}
            val = c.get("vpt2") if anharm else c.get("harmonic")
            mark = "$^\\dagger$" if (c.get("overlap") is not None and c["overlap"] < 0.7) else ""
            cells += [_fmt(val, digits) + mark, _delta(val, dft_ref, digits)]
        lines.append(" & ".join(cells) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    Path(out_path).write_text("\n".join(lines), encoding="utf-8")
