#!/usr/bin/env python3
"""Re-parse Gaussian logs and refresh the frequency data in results.json files.

results.json stores parsed frequency/intensity data captured at calculation
time. When the log parser is fixed (e.g. the DS(anharm) intensity overwrite or
the imaginary-frequency regex), the stored data goes stale. This script
re-parses each calculation's Gaussian log in place — no recomputation — and
rewrites only the parser-derived fields: frequencies (harmonic, anharmonic,
overtones, combination_bands) and dipole.

Usage:
    python scripts/refresh_results_from_logs.py [molecule ...]

With no arguments, refreshes every molecule under comparison_results/.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mace_gaussian.gaussian.parser import parse_gaussian_log

RESULTS_BASE = Path("comparison_results")
SKIP_DIRS = {"geometry_opt", "experimental"}
LOG_NAMES = ["gaussian_freq.log", "gaussian_dft.log"]


def find_log(calc_dir: Path) -> Path | None:
    for name in LOG_NAMES:
        candidate = calc_dir / name
        if candidate.exists():
            return candidate
    # Older DFT baselines are named after the molecule (hydrogen_cyanide_freq_anharm.log)
    logs = sorted(calc_dir.glob("*.log"))
    return logs[0] if len(logs) == 1 else None


def refresh_one(calc_dir: Path) -> str:
    json_path = calc_dir / "results.json"
    if not json_path.exists():
        return "no results.json"
    log_path = find_log(calc_dir)
    if log_path is None:
        return "no gaussian log"

    with json_path.open() as f:
        results = json.load(f)

    try:
        parsed = parse_gaussian_log(str(log_path))
    except Exception as e:  # keep going over the batch; report per-dir
        return f"parse failed: {e}"

    results["frequencies"] = {
        "harmonic": parsed.get("harmonic", []),
        "anharmonic": parsed.get("anharmonic", []),
        "overtones": parsed.get("overtones", []),
        "combination_bands": parsed.get("combination_bands", []),
    }
    if parsed.get("dipole_moment") is not None:
        results["dipole"] = parsed["dipole_moment"]

    with json_path.open("w") as f:
        json.dump(results, f, indent=2)
    return "refreshed"


def main() -> None:
    global RESULTS_BASE
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("molecules", nargs="*")
    ap.add_argument("--campaign", default=None, help="refresh campaigns/<NAME>/ instead")
    args = ap.parse_args()
    if args.campaign is not None:
        from mace_gaussian.campaign import campaign_paths

        RESULTS_BASE = campaign_paths(args.campaign).comparison
    molecules = args.molecules or sorted(p.name for p in RESULTS_BASE.iterdir() if p.is_dir())
    counts: dict[str, int] = {}
    for mol in molecules:
        mol_dir = RESULTS_BASE / mol
        if not mol_dir.is_dir():
            print(f"{mol}: not found, skipping")
            continue
        for calc_dir in sorted(mol_dir.iterdir()):
            if not calc_dir.is_dir() or calc_dir.name in SKIP_DIRS:
                continue
            status = refresh_one(calc_dir)
            counts[status] = counts.get(status, 0) + 1
            if status not in ("refreshed", "no results.json"):
                print(f"  {mol}/{calc_dir.name}: {status}")
    print("Summary:", ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
