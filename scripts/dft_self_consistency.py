#!/usr/bin/env python3
"""Harness self-consistency check: B3LYP through the external interface vs. native.

Runs the production ``run_frequency_calculation`` with Gaussian itself as the
"ML model" (``gaussian_b3lyp`` energy and dipole calculators), then runs a native
``freq(anharm) B3LYP/6-31G(d,p)`` job at exactly the geometry the harness wrote
into its own .gjf, and compares the two logs mode by mode. Same physics, same
geometry, same VPT2 settings: any difference is harness error.

Usage::

    python scripts/dft_self_consistency.py molecules/water.xyz
    python scripts/dft_self_consistency.py molecules/water.xyz --fd-dipole-derivs

Outputs land in ``validation_results/dft_self_consistency/<variant>/<molecule>/``:
``gaussian_b3lyp/`` (harness run), ``b3lyp_native/`` (reference run),
``summary.md`` and ``summary.json``.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from ase.io import read

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from mace_gaussian.calculators import dipole_factory  # noqa: E402
from mace_gaussian.calculators.gaussian_reference import (  # noqa: E402
    CALCULATOR_NAME,
    REFERENCE_BASIS,
    REFERENCE_METHOD,
    GaussianReferenceJob,
    set_shared_job,
)
from mace_gaussian.dft_baseline import (  # noqa: E402
    DFT_BASELINES,
    create_gaussian_dft_input,
    run_gaussian_dft,
)
from mace_gaussian.gaussian.fchk import convert_chk_to_fchk  # noqa: E402
from mace_gaussian.gaussian.parser import parse_gaussian_log  # noqa: E402
from mace_gaussian.utils.results import ResultsManager  # noqa: E402
from mace_gaussian.utils.units import HARTREE_TO_EV  # noqa: E402
from mace_gaussian.workflow import run_frequency_calculation  # noqa: E402

logger = logging.getLogger("dft_self_consistency")

NATIVE_NAME = "b3lyp_native"


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def run_harness(atoms, name, results_mgr, charge, multiplicity, nproc, fd_dipole) -> Path:
    """Stage 1: production external-interface run with Gaussian as the model."""
    set_shared_job(GaussianReferenceJob(nproc=nproc))
    dip = dipole_factory.get_calculator(CALCULATOR_NAME)
    dip.analytic_derivatives = not fd_dipole

    ok = run_frequency_calculation(
        atoms,
        name,
        CALCULATOR_NAME,
        CALCULATOR_NAME,
        results_mgr,
        charge=charge,
        multiplicity=multiplicity,
        keep_scratch=True,
    )
    if not ok:
        raise SystemExit("harness run failed; see log above")
    return results_mgr.create_frequency_directory(name, CALCULATOR_NAME, CALCULATOR_NAME)


def read_gjf_geometry(gjf: Path):
    """Return an Atoms object with the coordinates exactly as written in the .gjf."""
    from ase import Atoms

    lines = gjf.read_text().splitlines()
    # Route, blank, title, blank, "charge mult", coordinates..., blank
    idx = next(i for i, ln in enumerate(lines) if ln.startswith("#"))
    idx = next(i for i in range(idx + 1, len(lines)) if not lines[i].strip())  # after route
    idx = next(i for i in range(idx + 1, len(lines)) if not lines[i].strip())  # after title
    charge, mult = (int(x) for x in lines[idx + 1].split())
    symbols, positions = [], []
    for ln in lines[idx + 2 :]:
        if not ln.strip():
            break
        s, x, y, z = ln.split()
        symbols.append(s)
        positions.append([float(x), float(y), float(z)])
    return Atoms(symbols, positions=positions), charge, mult


def run_native(atoms, name, results_mgr, charge, multiplicity, nproc) -> Path:
    """Stage 2: native freq(anharm) at the harness geometry, no optimization."""
    native_dir = results_mgr.create_frequency_directory(name, NATIVE_NAME, NATIVE_NAME)
    print("=" * 60)
    print(f"NATIVE REFERENCE: {REFERENCE_METHOD}/{REFERENCE_BASIS} freq(anharm), no opt")
    print("=" * 60)
    t0 = time.time()
    create_gaussian_dft_input(
        atoms,
        "gaussian_dft.gjf",
        REFERENCE_METHOD,
        REFERENCE_BASIS,
        charge,
        multiplicity,
        title="self-consistency reference at the harness geometry",
        output_dir=native_dir,
        nproc=nproc,
        optimize=False,
    )
    ok, log = run_gaussian_dft("gaussian_dft.gjf", cwd=str(native_dir))
    if not ok:
        raise SystemExit(f"native run failed; see {log}")
    chk = native_dir / "gaussian_dft.chk"
    if chk.exists():
        convert_chk_to_fchk(str(chk), str(native_dir / "gaussian_dft.fchk"))
    parsed = parse_gaussian_log(log)
    results_mgr.save_frequency_results(
        molecule_name=name,
        energy_calculator=NATIVE_NAME,
        dipole_calculator=NATIVE_NAME,
        calculator_type="dft",
        frequencies_data={k: parsed.get(k, []) for k in _BAND_KINDS},
        energy=(parsed.get("final_energy_hartree") or 0.0) * HARTREE_TO_EV,
        dipole=parsed.get("dipole_moment"),
        runtime=time.time() - t0,
        gaussian_log=str(log),
        gaussian_gjf=str(native_dir / "gaussian_dft.gjf"),
        calculation_parameters={"vpt2_diagnostics": parsed.get("vpt2_diagnostics")},
        gaussian_timing=parsed.get("timing"),
    )
    print(f"  Completed in {time.time() - t0:.1f} seconds\n")
    return native_dir


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

_BAND_KINDS = ("harmonic", "anharmonic", "overtones", "combination_bands")


def _key(kind: str, entry: dict, index: int) -> tuple:
    if kind == "harmonic":
        return (index + 1,)
    if kind == "anharmonic":
        return (entry["mode"],)
    if kind == "overtones":
        return (entry["mode"], entry["overtone_level"])
    return (entry["mode1"], entry["mode2"])


def _freq(kind: str, entry: dict) -> float:
    return entry["freq_cm"] if kind in ("harmonic", "anharmonic") else entry["freq_anharmonic"]


def compare_logs(harness_log: Path, native_log: Path, tol_freq: float, tol_int_rel: float) -> dict:
    """Pair every band of the two logs and report the differences."""
    h = parse_gaussian_log(str(harness_log))
    n = parse_gaussian_log(str(native_log))
    int_abs_floor = 0.05  # km/mol: below this, relative tolerance is meaningless

    kinds = {}
    worst_freq = 0.0
    worst_int = 0.0
    failures = []
    for kind in _BAND_KINDS:
        rows = []
        h_map = {_key(kind, e, i): e for i, e in enumerate(h.get(kind, []))}
        n_map = {_key(kind, e, i): e for i, e in enumerate(n.get(kind, []))}
        for key in sorted(set(h_map) | set(n_map)):
            he, ne = h_map.get(key), n_map.get(key)
            if he is None or ne is None:
                failures.append(f"{kind} {key}: present on one side only")
                rows.append({"key": list(key), "missing": "harness" if he is None else "native"})
                continue
            fh, fn = _freq(kind, he), _freq(kind, ne)
            ih, in_ = he["ir_intensity"], ne["ir_intensity"]
            d_freq = fh - fn
            d_int = ih - in_
            int_tol = max(int_abs_floor, tol_int_rel * max(abs(ih), abs(in_)))
            ok_f, ok_i = abs(d_freq) <= tol_freq, abs(d_int) <= int_tol
            worst_freq = max(worst_freq, abs(d_freq))
            worst_int = max(worst_int, abs(d_int) / max(abs(in_), int_abs_floor))
            if not ok_f:
                failures.append(f"{kind} {key}: freq differs by {d_freq:+.4f} cm^-1")
            if not ok_i:
                failures.append(f"{kind} {key}: intensity differs by {d_int:+.4f} km/mol")
            rows.append(
                {
                    "key": list(key),
                    "harness_freq": fh,
                    "native_freq": fn,
                    "delta_freq": d_freq,
                    "harness_int": ih,
                    "native_int": in_,
                    "delta_int": d_int,
                    "pass": ok_f and ok_i,
                }
            )
        kinds[kind] = rows

    e_h, e_n = h.get("final_energy_hartree"), n.get("final_energy_hartree")
    d_h, d_n = h.get("dipole_moment") or {}, n.get("dipole_moment") or {}
    return {
        "bands": kinds,
        "energy_hartree": {"harness": e_h, "native": e_n, "delta": (e_h or 0) - (e_n or 0)},
        "dipole_debye": {
            "harness": d_h.get("magnitude"),
            "native": d_n.get("magnitude"),
            "delta": (d_h.get("magnitude") or 0) - (d_n.get("magnitude") or 0),
        },
        "vpt2_diagnostics": {
            "harness": h.get("vpt2_diagnostics"),
            "native": n.get("vpt2_diagnostics"),
        },
        "max_abs_delta_freq_cm": worst_freq,
        "max_rel_delta_intensity": worst_int,
        "tolerances": {"freq_cm": tol_freq, "intensity_rel": tol_int_rel},
        "failures": failures,
        "passed": not failures,
    }


def write_summary(result: dict, meta: dict, out_dir: Path) -> None:
    (out_dir / "summary.json").write_text(json.dumps({**meta, **result}, indent=2, default=float))
    lines = [
        f"# DFT self-consistency check: {meta['molecule']} ({meta['variant']})",
        "",
        f"Harness: `{CALCULATOR_NAME}` through `run_frequency_calculation` "
        f"({meta['harness_reference_jobs']} inner Gaussian jobs). "
        f"Native: `{REFERENCE_METHOD}/{REFERENCE_BASIS} freq(anharm)` at the harness geometry.",
        "",
        f"**Verdict: {'PASS' if result['passed'] else 'FAIL'}** "
        f"(max |Δν| = {result['max_abs_delta_freq_cm']:.4f} cm⁻¹ against "
        f"{result['tolerances']['freq_cm']} cm⁻¹; "
        f"max relative |ΔI| = {result['max_rel_delta_intensity']:.2%} against "
        f"{result['tolerances']['intensity_rel']:.0%})",
        "",
        f"Final energy: harness {result['energy_hartree']['harness']:.10f} Ha, "
        f"native {result['energy_hartree']['native']:.10f} Ha, "
        f"Δ = {result['energy_hartree']['delta']:.2e} Ha",
        f"Dipole: harness {result['dipole_debye']['harness']} D, native "
        f"{result['dipole_debye']['native']} D",
        "",
    ]
    titles = {
        "harmonic": "Harmonic fundamentals",
        "anharmonic": "VPT2 fundamentals",
        "overtones": "Overtones",
        "combination_bands": "Combination bands",
    }
    for kind, rows in result["bands"].items():
        if not rows:
            continue
        lines += [
            f"## {titles[kind]}",
            "",
            "| band | freq harness | freq native | Δ (cm⁻¹) | I harness | I native "
            "| ΔI (km/mol) | ok |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            key = "/".join(str(k) for k in r["key"])
            if "missing" in r:
                lines.append(f"| {key} | missing on {r['missing']} side | | | | | | ✗ |")
                continue
            lines.append(
                f"| {key} | {r['harness_freq']:.4f} | {r['native_freq']:.4f} | "
                f"{r['delta_freq']:+.4f} | {r['harness_int']:.4f} | {r['native_int']:.4f} | "
                f"{r['delta_int']:+.4f} | {'✓' if r['pass'] else '✗'} |"
            )
        lines.append("")
    if result["failures"]:
        lines += ["## Failures", ""] + [f"- {f}" for f in result["failures"]] + [""]
    (out_dir / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("xyz", type=Path)
    ap.add_argument("--name", help="molecule name (default: xyz stem)")
    ap.add_argument(
        "--out", type=Path, default=REPO / "validation_results" / "dft_self_consistency"
    )
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--multiplicity", type=int, default=1)
    ap.add_argument("--nproc", type=int, default=4, help="cores for the inner and native jobs")
    ap.add_argument(
        "--fd-dipole-derivs",
        action="store_true",
        help="use the base-class finite-difference dipole derivatives (the ML path) "
        "instead of Gaussian's analytic ones",
    )
    ap.add_argument("--tol-freq", type=float, default=0.1, help="cm^-1")
    ap.add_argument("--tol-int-rel", type=float, default=0.01, help="relative")
    ap.add_argument(
        "--compare-only", action="store_true", help="skip both runs, re-compare existing logs"
    )
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    cfg = DFT_BASELINES["b3lyp"]
    assert (cfg["method"], cfg["basis"]) == (REFERENCE_METHOD, REFERENCE_BASIS), (
        "reference calculator and DFT_BASELINES['b3lyp'] disagree"
    )

    name = args.name or args.xyz.stem
    variant = "fd_dipole" if args.fd_dipole_derivs else "analytic"
    (args.out / variant).mkdir(parents=True, exist_ok=True)
    results_mgr = ResultsManager(str(args.out / variant))
    mol_dir = results_mgr.create_molecule_directory(name)
    harness_dir = mol_dir / CALCULATOR_NAME
    native_dir = mol_dir / NATIVE_NAME

    t0 = time.time()
    n_ref_jobs = None
    if not args.compare_only:
        atoms = read(str(args.xyz))
        harness_dir = run_harness(
            atoms,
            name,
            results_mgr,
            args.charge,
            args.multiplicity,
            args.nproc,
            args.fd_dipole_derivs,
        )
        from mace_gaussian.calculators.gaussian_reference import get_shared_job

        n_ref_jobs = dict(get_shared_job().n_jobs)
        ref_atoms, charge, mult = read_gjf_geometry(harness_dir / "gaussian_freq.gjf")
        native_dir = run_native(ref_atoms, name, results_mgr, charge, mult, args.nproc)

    result = compare_logs(
        harness_dir / "gaussian_freq.log",
        native_dir / "gaussian_dft.log",
        args.tol_freq,
        args.tol_int_rel,
    )
    meta = {
        "molecule": name,
        "variant": variant,
        "reference": f"{REFERENCE_METHOD}/{REFERENCE_BASIS}",
        "harness_reference_jobs": n_ref_jobs,
        "wall_time_s": time.time() - t0,
        "harness_dir": str(harness_dir),
        "native_dir": str(native_dir),
    }
    write_summary(result, meta, mol_dir)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
