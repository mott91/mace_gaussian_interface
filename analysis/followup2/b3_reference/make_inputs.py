"""Follow-up 2, B3: inputs for the reference-level check (harmonic only).

Eight small molecules that carry the intensity differences of A1. Five levels:

- g_b3lyp_avtz_opt   B3LYP/aug-cc-pVTZ, own geometry (opt=verytight + freq): how far the
                     reference itself moves.
- g_b3lyp_avtz_fix   B3LYP/aug-cc-pVTZ at the B3LYP/6-31G(d,p) geometry (freq only): APT
                     with geometry and normal modes held fixed.
- g_pbe0_tzvp_fix    PBE0/def2-TZVP at that geometry: stand-in for MACE4IR's PBE0 /
                     FHI-aims tier-2 (QCML; arXiv 2508.19118, sec. 4.3).
- o_wb97mv_fix       wB97M-V/def2-TZVPD (OMol25: POLAR-1's training level), ORCA NumFreq
                     (no analytic Hessian with VV10 in ORCA 6.1.1).
- o_wb97md3bj_fix    wB97M-D3(BJ)/def2-TZVPPD (SPICE: MACE-MDP's training level), ORCA Freq.

The fixed geometry is the B3LYP/6-31G(d,p) minimum of the campaign (fchk, Bohr -> Angstrom),
written unchanged, so with nosymm (Gaussian) / no reorientation (ORCA) every APT is in the
frame of the campaign fchk.
"""

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CAMP = Path("/home/mot/mace_gaussian/campaigns/2026")
REMOTE = "/scratch_rune03a/mot/calculations/mace_gaussian/followup2_b3"
MOLS = ["methanol", "formic_acid_trans", "acetic_acid", "carbon_dioxide", "nitric_acid",
        "carbonic_acid", "hydrogen_peroxide", "glycine"]  # fmt: skip
SYM = {1: "H", 6: "C", 7: "N", 8: "O"}
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


GAUSS = {
    "g_b3lyp_avtz_opt": "# opt=verytight freq b3lyp/aug-cc-pVTZ nosymm int=ultrafine",
    "g_b3lyp_avtz_fix": "# freq b3lyp/aug-cc-pVTZ nosymm int=ultrafine",
    "g_pbe0_tzvp_fix": "# freq pbe1pbe/def2TZVP nosymm int=ultrafine",
}
ORCA = {
    "o_wb97mv_fix": "! wB97M-V def2-TZVPD NumFreq TightSCF DefGrid3",
    "o_wb97md3bj_fix": "! wB97M-D3BJ def2-TZVPPD Freq TightSCF DefGrid3",
}

SLURM = """#!/bin/bash
#SBATCH --job-name=f2b3_{tag}
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


def main():
    jobs = []
    for mol in MOLS:
        geo = "\n".join(fchk_geometry(mol))
        for lev, route in GAUSS.items():
            d = HERE / lev / mol
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{mol}.gjf").write_text(
                f"%chk={mol}.chk\n%mem=32GB\n%NProcShared={NCPU}\n{route}\n\n"
                f"{mol} {lev} (follow-up 2, B3)\n\n0 1\n{geo}\n\n"
            )
            body = (f"module load rune/gaussian/g16-c.02\nexport OMP_NUM_THREADS={NCPU}\n"
                    f"g16 {mol}.gjf\nformchk {mol}.chk {mol}.fchk\nrm -f {mol}.chk")  # fmt: skip
            rdir = f"{REMOTE}/{lev}/{mol}"
            (d / "job.sh").write_text(SLURM.format(tag=f"{lev}_{mol}", rdir=rdir, ncpu=NCPU, body=body))
            jobs.append(f"{lev}/{mol}")
        for lev, route in ORCA.items():
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
            # ORCA starts 8 MPI ranks: 8 tasks on ONE node (rune03 scratch is local; a first
            # submission without --nodes=1 put ranks on rune04 and failed)
            js = js.replace("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8",
                            "#SBATCH --nodes=1\n#SBATCH --ntasks=8\n#SBATCH --cpus-per-task=1")
            (d / "job.sh").write_text(js)
            jobs.append(f"{lev}/{mol}")
    (HERE / "jobs.txt").write_text("\n".join(jobs) + "\n")
    print(len(jobs), "jobs")


if __name__ == "__main__":
    main()
