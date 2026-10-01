#!/usr/bin/env python3
"""Protocol test before the 2026 campaign: does symmetry change the B3LYP VPT2 baseline?

The existing baselines run ``opt freq(anharm)`` with Gaussian's default symmetry handling,
so ethane ends up exactly D3d and gets symmetric-top VPT2, while the ML runs (ASE geometry,
slightly broken symmetry) get asymmetric-top VPT2. This reruns three baselines from the same
starting geometry with ``nosymm`` (and ``nosymm opt=tight``) on the cluster and compares
the fundamentals with the existing baseline.

Results go to ``validation_results/symmetry_test/``; ``comparison_results`` is not touched.

Usage::

    python scripts/symmetry_test.py submit      # write inputs, sbatch on rune03
    python scripts/symmetry_test.py status
    python scripts/symmetry_test.py retrieve    # copy back, parse, print the comparison
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from mace_gaussian import slurm  # noqa: E402
from mace_gaussian.gaussian.parser import GaussianLogParser  # noqa: E402

HOST = "mot@tci5"
OUT = REPO / "validation_results" / "symmetry_test"
MOLECULES = ["ethane", "methane", "ammonia"]
VARIANTS = {
    "nosymm": "# opt freq(anharm) b3lyp/6-31G(d,p) nosymm",
    "nosymm_tight": "# opt=tight freq(anharm) b3lyp/6-31G(d,p) nosymm",
}
JOBS = OUT / "jobs.json"


def _start_geometry(mol: str) -> list[str]:
    """Charge/multiplicity and atom lines of the existing baseline input (same start)."""
    gjf = REPO / "comparison_results" / mol / "b3lyp_6-31Gdp" / "gaussian_dft.gjf"
    blocks = gjf.read_text().strip().split("\n\n")
    return blocks[2].strip().splitlines()  # link0+route / title / charge+atoms


def submit() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    for mol in MOLECULES:
        geom = _start_geometry(mol)
        for variant, route in VARIANTS.items():
            name = f"symtest_{mol}_{variant}"
            d = OUT / name
            d.mkdir(exist_ok=True)
            gjf = d / f"{name}.gjf"
            gjf.write_text(
                f"%chk={name}_freq_anharm.chk\n%mem=4GB\n%NProcShared=4\n{route}\n\n"
                f"symmetry test: {mol} {variant}\n\n" + "\n".join(geom) + "\n\n"
            )
            jobs.append({"name": name, "gjf_path": str(gjf)})
    ids = slurm.submit_dft_jobs(jobs, HOST, REPO / "templates" / "slurm_dft.sh", str(OUT))
    JOBS.write_text(json.dumps(ids, indent=1))
    for name, jid in ids.items():
        print(f"  {name:32s} job {jid}")
    missing = {j["name"] for j in jobs} - set(ids)
    if missing:
        print("NOT submitted:", ", ".join(sorted(missing)))


def status() -> None:
    """One non-blocking sacct query (slurm.poll_jobs would wait for completion)."""
    ids = json.loads(JOBS.read_text())
    res = slurm._ssh_run(HOST, f"sacct -n -X -P -o JobID,State,Elapsed -j {','.join(ids.values())}")
    by_id = {line.split("|")[0]: line.split("|")[1:] for line in res.stdout.split()}
    for name, jid in ids.items():
        state, elapsed = by_id.get(jid, ["?", ""])
        print(f"  {name:32s} {jid:>9s}  {state:10s} {elapsed}")


def _fundamentals(log: Path) -> list[tuple[float, float]]:
    rows = GaussianLogParser(str(log)).parse_anharmonic_frequencies()
    return sorted((r["freq_harmonic"], r["freq_cm"]) for r in rows)


def retrieve() -> None:
    ids = json.loads(JOBS.read_text())
    ok = slurm.retrieve_results(HOST, list(ids), results_dir=str(OUT))
    print("retrieved:", {k: v for k, v in ok.items()})
    for mol in MOLECULES:
        old = _fundamentals(
            REPO / "comparison_results" / mol / "b3lyp_6-31Gdp" / "gaussian_dft.log"
        )
        cols = {"symm (old)": old}
        for variant in VARIANTS:
            log = (
                OUT
                / f"symtest_{mol}_{variant}"
                / "b3lyp_6-31Gdp"
                / f"symtest_{mol}_{variant}_freq_anharm.log"
            )
            if log.exists():
                cols[variant] = _fundamentals(log)
        print(f"\n{mol}: harmonic -> VPT2 fundamental (cm-1), sorted by harmonic")
        print("  " + "".join(f"{c:>24s}" for c in cols))
        n = max(len(v) for v in cols.values())
        for k in range(n):
            cells = []
            for v in cols.values():
                cells.append(f"{v[k][0]:9.1f} -> {v[k][1]:9.1f}  " if k < len(v) else " " * 24)
            print("  " + "".join(cells))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    {"submit": submit, "status": status, "retrieve": retrieve}[cmd]()
