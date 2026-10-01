"""Band-level data for the supervisor figures: every VPT2 band and χ, ML next to DFT.

Two tables, both in the DFT checkpoint (ascending-frequency) mode numbering:

- **bands**: one row per (molecule, ML run, band). Band types are fundamentals ``i``,
  overtones ``2i`` and combinations ``i+j``. ML bands are put into the DFT numbering
  through the eigenvector mode mapping, so a row always compares the same motion.
- **chi**: one row per (molecule, energy model, anharmonic constant). The supervisor's
  χ is computed from band positions, the way experiment would see it:
  ``x_ii = (overtone - 2 fundamental) / 2`` and ``x_ij = combination - v_i - v_j``.
  Gaussian's printed (deperturbed) X matrix sits next to it; where the two disagree by
  more than ``RESONANCE_TOL_CM`` the band positions were shifted by a resonance
  treatment, and the row is flagged ``resonant``.

Frequencies and χ depend only on the energy model; the dipole model changes intensities
only. The chi table therefore uses one run per energy model.

Degenerate modes (linear molecules, symmetric tops: entries carrying ``l``) are left out
of the chi table, because their overtone splits by l and the simple formula does not hold.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ..gaussian.parser import GaussianLogParser
from .analyze_spectra import gaussian_mode_to_checkpoint_index

DIPOLE_MODELS = ("mace_polar1", "mace_ml", "espaloma")
RESONANCE_TOL_CM = 1.0


@dataclass
class BandRow:
    molecule: str
    energy_model: str
    dipole_model: str
    band_type: str  # fundamental | overtone | combination
    i: int  # DFT mode, 1-based checkpoint order
    j: int | None  # second mode of a combination band
    l: int | None  # noqa: E741 -- vibrational angular momentum of a degenerate overtone (HCN 2v2: 0, +-2)
    dft_harmonic: float
    dft_anharmonic: float
    ml_harmonic: float
    ml_anharmonic: float
    dft_intensity: float
    ml_intensity: float
    overlap: float | None  # smallest eigenvector overlap of the modes involved


@dataclass
class ChiRow:
    molecule: str
    energy_model: str
    kind: str  # diagonal | offdiagonal
    i: int
    j: int
    dft_bands: float
    ml_bands: float
    dft_xmatrix: float | None
    ml_xmatrix: float | None
    dft_resonant: bool
    ml_resonant: bool
    overlap: float | None


def split_run_name(name: str) -> tuple[str, str]:
    """``mace_omol_mace_ml`` -> (``mace_omol``, ``mace_ml``)."""
    for dip in DIPOLE_MODELS:
        if name.endswith("_" + dip):
            return name[: -len(dip) - 1], dip
    return name, ""


def _results_log(run_dir: Path, results: dict) -> Path | None:
    """The log whose fundamentals match results.json (DFT folders can hold several)."""
    want = [round(e["freq_cm"], 2) for e in results["frequencies"].get("anharmonic", [])]
    logs = sorted(run_dir.glob("*.log"))
    for log in logs:
        got = [
            round(e["freq_cm"], 2)
            for e in GaussianLogParser(str(log)).parse_anharmonic_frequencies()
        ]
        if got == want:
            return log
    return logs[0] if logs else None


class _Run:
    """One Gaussian run: bands keyed in its own checkpoint numbering (1-based)."""

    def __init__(self, run_dir: Path):
        self.results = json.loads((run_dir / "results.json").read_text())
        freqs = self.results["frequencies"]
        anharm = freqs.get("anharmonic", [])
        g2c = gaussian_mode_to_checkpoint_index(anharm)
        self.g2c = g2c
        self.fund = {g2c.get(e["mode"], e["mode"]): e for e in anharm}
        # Symmetric tops label every mode with l; only modes with two components are degenerate
        components: dict[int, list[int]] = {}
        for e in anharm:
            if "mode_gaussian" in e:
                components.setdefault(e["mode_gaussian"], []).append(g2c.get(e["mode"], e["mode"]))
        self.degenerate = {k for ks in components.values() if len(ks) > 1 for k in ks}
        self.over = {
            (g2c.get(e["mode"], e["mode"]), e.get("l")): e
            for e in freqs.get("overtones", [])
            if e.get("overtone_level", 2) == 2
        }
        self.comb = {}
        for e in freqs.get("combination_bands", []):
            a, b = g2c.get(e["mode1"], e["mode1"]), g2c.get(e["mode2"], e["mode2"])
            self.comb[(min(a, b), max(a, b))] = e
        log = _results_log(run_dir, self.results)
        raw = GaussianLogParser(str(log)).parse_x_matrix() if log else {}
        self.x = {}
        for (a, b), v in raw.items():
            ca, cb = g2c.get(a, a), g2c.get(b, b)
            self.x[(min(ca, cb), max(ca, cb))] = v


def _translate(mapping: dict[int, int], k: int) -> int | None:
    """ML checkpoint mode (1-based) -> DFT checkpoint mode (1-based)."""
    d = mapping.get(k - 1)
    return None if d is None else d + 1


def build_tables(
    molecule: str,
    results_base: Path = Path("comparison_results"),
    analysis_base: Path = Path("analysis_results"),
) -> tuple[list[BandRow], list[ChiRow]]:
    report = json.loads((analysis_base / molecule / "report_data.json").read_text())
    dft_dirs = sorted(p for p in (results_base / molecule).glob("b3lyp*") if p.is_dir())
    if not dft_dirs:
        return [], []
    dft = _Run(dft_dirs[0])

    bands: list[BandRow] = []
    chi: list[ChiRow] = []
    chi_done: set[str] = set()
    for comp in report["comparisons"]:
        name = comp["name"]
        run_dir = results_base / molecule / name
        if not (run_dir / "results.json").exists():
            continue
        energy, dipole = split_run_name(name)
        ml = _Run(run_dir)
        mapping = {int(k): int(v) for k, v in (comp.get("mode_mapping") or {}).items()}
        overlaps = {int(k): float(v) for k, v in (comp.get("mode_overlaps") or {}).items()}
        ovl = _overlap_lookup(mapping, overlaps)
        bands += _band_rows(molecule, energy, dipole, dft, ml, mapping, ovl)
        if energy not in chi_done:
            chi_done.add(energy)
            chi += _chi_rows(molecule, energy, dft, ml, mapping, ovl)
    return bands, chi


def _overlap_lookup(mapping: dict[int, int], overlaps: dict[int, float]):
    """Smallest eigenvector overlap of the given DFT modes (1-based), None if any unmapped."""
    to_ml = {d + 1: m + 1 for m, d in mapping.items()}

    def ovl(*dft_modes: int) -> float | None:
        vals = [overlaps.get(to_ml[d] - 1) for d in dft_modes if d in to_ml]
        vals = [v for v in vals if v is not None]
        return min(vals) if len(vals) == len(dft_modes) else None

    return ovl


def _band_rows(molecule, energy, dipole, dft: _Run, ml: _Run, mapping, ovl) -> list[BandRow]:
    def row(kind, i, j, d, m):
        return BandRow(
            molecule,
            energy,
            dipole,
            kind,
            i,
            j,
            d.get("l"),
            d["freq_harmonic"],
            d.get("freq_cm", d.get("freq_anharmonic")),
            m["freq_harmonic"],
            m.get("freq_cm", m.get("freq_anharmonic")),
            d["ir_intensity"],
            m["ir_intensity"],
            ovl(i) if j is None else ovl(i, j),
        )

    rows = []
    for k, m in ml.fund.items():
        i = _translate(mapping, k)
        if i in dft.fund:
            rows.append(row("fundamental", i, None, dft.fund[i], m))
    for (k, ell), m in ml.over.items():
        i = _translate(mapping, k)
        if (i, ell) in dft.over:
            rows.append(row("overtone", i, None, dft.over[(i, ell)], m))
    for (a, b), m in ml.comb.items():
        ia, ib = _translate(mapping, a), _translate(mapping, b)
        if ia is None or ib is None:
            continue
        key = (min(ia, ib), max(ia, ib))
        if key in dft.comb:
            rows.append(row("combination", key[0], key[1], dft.comb[key], m))
    return rows


def _chi_rows(molecule, energy, dft: _Run, ml: _Run, mapping, ovl) -> list[ChiRow]:
    def from_bands(run: _Run, i: int, j: int) -> float | None:
        fi, fj = run.fund.get(i), run.fund.get(j)
        if i == j:
            o = run.over.get((i, None))
            if fi is None or o is None:
                return None
            return (o["freq_anharmonic"] - 2 * fi["freq_cm"]) / 2
        c = run.comb.get((min(i, j), max(i, j)))
        if fi is None or fj is None or c is None:
            return None
        return c["freq_anharmonic"] - fi["freq_cm"] - fj["freq_cm"]

    def resonant(bands_val: float, x: float | None) -> bool:
        return x is not None and abs(bands_val - x) > RESONANCE_TOL_CM

    rows = []
    ml_to_dft = {k: _translate(mapping, k) for k in ml.fund}
    dft_to_ml = {d: k for k, d in ml_to_dft.items() if d is not None}
    modes = sorted(i for i in dft.fund if i in dft_to_ml)
    for a in modes:
        for b in modes:
            if b < a:
                continue
            if a in dft.degenerate or b in dft.degenerate:
                continue
            ma, mb = dft_to_ml[a], dft_to_ml[b]
            if ma in ml.degenerate or mb in ml.degenerate:
                continue
            d_val, m_val = from_bands(dft, a, b), from_bands(ml, ma, mb)
            if d_val is None or m_val is None:
                continue
            d_x = dft.x.get((a, b))
            m_x = ml.x.get((min(ma, mb), max(ma, mb)))
            rows.append(
                ChiRow(
                    molecule,
                    energy,
                    "diagonal" if a == b else "offdiagonal",
                    a,
                    b,
                    d_val,
                    m_val,
                    d_x,
                    m_x,
                    resonant(d_val, d_x),
                    resonant(m_val, m_x),
                    ovl(a) if a == b else ovl(a, b),
                )
            )
    return rows


def write_tables(molecules: list[str], out_dir: Path, **kw) -> tuple[Path, Path]:
    import csv

    out_dir.mkdir(parents=True, exist_ok=True)
    all_bands, all_chi = [], []
    for mol in molecules:
        b, c = build_tables(mol, **kw)
        all_bands += b
        all_chi += c
    paths = []
    for name, rows, cls in (("bands.csv", all_bands, BandRow), ("chi.csv", all_chi, ChiRow)):
        path = out_dir / name
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(cls.__dataclass_fields__))
            w.writeheader()
            for r in rows:
                w.writerow(asdict(r))
        paths.append(path)
    return paths[0], paths[1]
