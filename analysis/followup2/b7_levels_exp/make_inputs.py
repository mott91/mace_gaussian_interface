"""Follow-up 2, B7: higher levels for the B4 molecules (closing the gap between B3 and B4).

B3 computed harmonic APTs at four better levels for 8 molecules, B4 compared region-integrated
intensities with NIST gas-phase spectra for 5 molecules; only methanol is in both. This writes the
four B3 levels for the other four B4 molecules, with exactly the B3 settings
(b3_reference/make_inputs.py): the B3LYP/6-31G(d,p) campaign minimum (fchk, Bohr -> Angstrom)
written unchanged, nosymm / no reorientation, 8 cores on rune03.

- g_b3lyp_avtz_fix  B3LYP/aug-cc-pVTZ (Gaussian, analytic freq)
- g_pbe0_tzvp_fix   PBE0/def2-TZVP (Gaussian, analytic freq): MACE4IR stand-in
- o_wb97mv_fix      wB97M-V/def2-TZVPD (ORCA 6.1.1, NumFreq): POLAR-1 / OMol25
- o_wb97md3bj_fix   wB97M-D3(BJ)/def2-TZVPPD (ORCA 6.1.1, analytic Freq): MDP / SPICE
- g_b3lyp_631gdp_harm  B3LYP/6-31G(d,p) harmonic freq at the same geometry and 8 cores, for the
                    cost table only (the campaign runs are opt + VPT2 on 4 cores). Also for
                    methanol, whose four higher levels are reused from B3.

    python analysis/followup2/b7_levels_exp/make_inputs.py
writes <level>/<mol>/{input, job.sh} and jobs.txt; submit.sh copies to the cluster and submits.
"""

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CAMP = Path("/home/mot/mace_gaussian/campaigns/2026")
REMOTE = "/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b7"
NEW = ["ethanol", "isopropanol", "1-butanol", "sulfur_dioxide"]
COST_ONLY = ["methanol"]  # four higher levels from B3; only the 6-31G(d,p) timing job is new
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
BOHR = 0.529177210903
NCPU = 8


def fchk_geometry(mol):
    t = (CAMP / mol / f"{mol}_freq_anharm.fchk").read_text()

    def arr(key):
        i = t.index(key)
        n = int(t[i : t.index("\n", i)].split("N=")[1])
        return np.array([float(v) for v in t[t.index("\n", i) + 1 :].split()[:n]])

    z = arr("Atomic numbers").astype(int)
    x = arr("Current cartesian coordinates").reshape(-1, 3) * BOHR
    return [f"{SYM[a]:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}" for a, p in zip(z, x)]


# identical to B3 (b3_reference/make_inputs.py)
GAUSS = {
    "g_b3lyp_avtz_fix": "# freq b3lyp/aug-cc-pVTZ nosymm int=ultrafine",
    "g_pbe0_tzvp_fix": "# freq pbe1pbe/def2TZVP nosymm int=ultrafine",
}
ORCA = {
    "o_wb97mv_fix": "! wB97M-V def2-TZVPD NumFreq TightSCF DefGrid3",
    "o_wb97md3bj_fix": "! wB97M-D3BJ def2-TZVPPD Freq TightSCF DefGrid3",
}
# the campaign route without opt and anharm
COST = {"g_b3lyp_631gdp_harm": "# freq b3lyp/6-31G(d,p) nosymm int=ultrafine"}

SLURM = """#!/bin/bash
#SBATCH --job-name=f2b7_{tag}
#SBATCH --output={rdir}/slurm_%j.out
#SBATCH --error={rdir}/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={ncpu}
#SBATCH --mem=40G
#SBATCH --partition=CPU_rune
#SBATCH --nodelist=rune03
cd {rdir}
date +%s > start_epoch.txt
{body}
date +%s > end_epoch.txt
"""


def gauss_job(mol, lev, route, geo, jobs):
    d = HERE / lev / mol
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{mol}.gjf").write_text(
        f"%chk={mol}.chk\n%mem=32GB\n%NProcShared={NCPU}\n{route}\n\n"
        f"{mol} {lev} (follow-up 2, B7)\n\n0 1\n{geo}\n\n"
    )
    body = (f"module load rune/gaussian/g16-c.02\nexport OMP_NUM_THREADS={NCPU}\n"
            f"g16 {mol}.gjf\nformchk {mol}.chk {mol}.fchk\nrm -f {mol}.chk")  # fmt: skip
    rdir = f"{REMOTE}/{lev}/{mol}"
    (d / "job.sh").write_text(SLURM.format(tag=f"{lev}_{mol}", rdir=rdir, ncpu=NCPU, body=body))
    jobs.append(f"{lev}/{mol}")


def orca_job(mol, lev, route, geo, jobs):
    d = HERE / lev / mol
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{mol}.inp").write_text(
        f"{route}\n%pal nprocs {NCPU} end\n%maxcore 4000\n"
        f"%freq CentralDiff true end\n* xyz 0 1\n{geo}\n*\n"
    )
    body = ("module load chaos/orca/6.1.1\n"
            f"$(which orca) {mol}.inp > {mol}.out\nrm -f {mol}.gbw {mol}*.tmp*")  # fmt: skip
    rdir = f"{REMOTE}/{lev}/{mol}"
    js = SLURM.format(tag=f"{lev}_{mol}", rdir=rdir, ncpu=NCPU, body=body)
    # ORCA starts 8 MPI ranks: 8 tasks on ONE node (as B3 after its failed first submission)
    js = js.replace("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8",
                    "#SBATCH --nodes=1\n#SBATCH --ntasks=8\n#SBATCH --cpus-per-task=1")
    (d / "job.sh").write_text(js)
    jobs.append(f"{lev}/{mol}")


def main():
    jobs = []
    for mol in NEW + COST_ONLY:
        geo = "\n".join(fchk_geometry(mol))
        if mol in NEW:
            for lev, route in GAUSS.items():
                gauss_job(mol, lev, route, geo, jobs)
            for lev, route in ORCA.items():
                orca_job(mol, lev, route, geo, jobs)
        for lev, route in COST.items():
            gauss_job(mol, lev, route, geo, jobs)
    (HERE / "jobs.txt").write_text("\n".join(jobs) + "\n")
    print(len(jobs), "jobs")


if __name__ == "__main__":
    main()
