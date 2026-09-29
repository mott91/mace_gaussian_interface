"""Gaussian itself as the "ML model": the harness self-consistency check.

Every ML result in this project passes through the external interface: Gaussian
hands a geometry to ``workflow.run_next_calculation``, the calculators return
energy, gradient, Hessian, dipole and dipole derivatives in ASE-side units, and
``gaussian.io.write_gaussian_output`` converts them back to atomic units. A unit
factor, a sign or an index-ordering error anywhere in that pipe would shift every
ML spectrum by the same amount, indistinguishable from a real model error.

This module plugs a Gaussian single-point (B3LYP/6-31G(d,p) by default) into that
same pipe. With Gaussian on both ends, the physics is identical and any
difference to a native ``freq(anharm)`` run is harness error. The expected
agreement is finite-difference noise, well below 0.1 cm^-1.

Two adapters share one cached ``GaussianReferenceJob`` so each geometry is
computed once:

- ``GaussianReferenceCalculator``: an ASE calculator (energy, forces,
  ``get_hessian``), registered as ``gaussian_b3lyp`` in ``workflow.calculator``.
- ``GaussianReferenceDipoleCalculator``: a ``DipoleCalculatorBase`` returning the
  analytic dipole, dipole derivatives and polarizability from the same job,
  registered as ``gaussian_b3lyp`` in the dipole factory.

The inner job always runs with ``nosymm``. Without it Gaussian rotates the
molecule into its standard orientation and the gradient, Hessian and dipole would
come back in a different frame from the coordinates the outer job passed in.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from ..gaussian.fchk import convert_chk_to_fchk, parse_fchk_section
from ..utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV
from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)

REFERENCE_METHOD = "b3lyp"
REFERENCE_BASIS = "6-31G(d,p)"
CALCULATOR_NAME = "gaussian_b3lyp"

# Unit factors from Gaussian's fchk (atomic units) to the ASE-side contract.
_HARTREE_PER_BOHR_TO_EV_PER_A = HARTREE_TO_EV / BOHR_TO_ANGSTROM
_HARTREE_PER_BOHR2_TO_EV_PER_A2 = HARTREE_TO_EV / BOHR_TO_ANGSTROM**2
_BOHR3_TO_A3 = BOHR_TO_ANGSTROM**3

# Tolerance for the frame check: fchk coordinates must equal the input geometry.
_FRAME_TOL_A = 1e-6


@dataclass
class ReferenceResult:
    """One Gaussian single point, in the units the harness expects from an ML model.

    ``hessian``, ``dipole_derivatives`` and ``polarizability`` are None for a
    gradient-only job (``level="force"``).
    """

    energy: float  # eV
    forces: np.ndarray  # eV/Angstrom, shape (natoms, 3)
    dipole: np.ndarray  # e*Bohr, shape (3,)
    hessian: np.ndarray | None  # eV/Angstrom^2, shape (3*natoms, 3*natoms)
    dipole_derivatives: np.ndarray | None  # e, shape (3*natoms, 3): row = atom coordinate
    polarizability: np.ndarray | None  # Angstrom^3, shape (3, 3)
    level: str  # "force" or "freq"


def parse_reference_fchk(
    fchk_path: str | Path, expected_positions_A: np.ndarray | None = None
) -> ReferenceResult:
    """Read a ``force nosymm`` or ``freq nosymm`` fchk into harness units.

    Parameters
    ----------
    fchk_path:
        Formatted checkpoint written by the inner job.
    expected_positions_A:
        If given, the fchk coordinates (Bohr) must match these Angstrom positions
        to ``_FRAME_TOL_A``. This is the guard against a reoriented inner job.

    Raises
    ------
    ValueError
        Frame mismatch, or a required section is missing.
    """
    content = Path(fchk_path).read_text()
    natoms = int(_parse_fchk_scalar(content, "Number of atoms"))
    n3 = 3 * natoms

    positions = parse_fchk_section(content, "Current cartesian coordinates").reshape(natoms, 3)
    positions_A = positions * BOHR_TO_ANGSTROM
    if expected_positions_A is not None:
        shift = np.abs(positions_A - np.asarray(expected_positions_A)).max()
        if shift > _FRAME_TOL_A:
            raise ValueError(
                f"fchk coordinates differ from the input geometry by up to {shift:.2e} A: "
                "the inner job reoriented the molecule (missing nosymm?)"
            )

    # fchk stores the scalar as "Total Energy   R   -7.64E+01" without an N= count.
    energy_hartree = _parse_fchk_scalar(content, "Total Energy")
    gradient = parse_fchk_section(content, "Cartesian Gradient").reshape(natoms, 3)
    dipole = parse_fchk_section(content, "Dipole Moment")

    hessian = dipole_derivatives = polarizability = None
    level = "force"
    if "Cartesian Force Constants" in content:
        level = "freq"
        tri = parse_fchk_section(content, "Cartesian Force Constants")
        hessian_au = _lower_triangle_to_full(tri, n3)
        hessian = hessian_au * _HARTREE_PER_BOHR2_TO_EV_PER_A2
        # Verified layout (water probe, 2026-09-18): the 9*natoms values are
        # d(mu_x, mu_y, mu_z)/dR for R = x1, y1, z1, x2, ... , i.e. C-order (3N, 3),
        # which is exactly the (3*natoms, 3) array write_gaussian_output consumes.
        dipole_derivatives = parse_fchk_section(content, "Dipole Derivatives").reshape(n3, 3)
        if "Polarizability" in content:
            p = parse_fchk_section(content, "Polarizability")  # xx, xy, yy, xz, yz, zz in Bohr^3
            polarizability = (
                np.array([[p[0], p[1], p[3]], [p[1], p[2], p[4]], [p[3], p[4], p[5]]])
                * _BOHR3_TO_A3
            )

    return ReferenceResult(
        energy=float(energy_hartree * HARTREE_TO_EV),
        forces=-gradient * _HARTREE_PER_BOHR_TO_EV_PER_A,
        dipole=dipole,
        hessian=hessian,
        dipole_derivatives=dipole_derivatives,
        polarizability=polarizability,
        level=level,
    )


def _parse_fchk_scalar(content: str, name: str) -> float:
    for line in content.splitlines():
        if line.startswith(name) and " N=" not in line:
            return float(line.split()[-1])
    raise ValueError(f"Section '{name}' not found in .fchk file")


def _lower_triangle_to_full(tri: np.ndarray, n: int) -> np.ndarray:
    """Expand Gaussian's row-major lower triangle ((1,1),(2,1),(2,2),(3,1),...) to (n, n)."""
    if tri.size != n * (n + 1) // 2:
        raise ValueError(f"Expected {n * (n + 1) // 2} lower-triangle values, got {tri.size}")
    full = np.zeros((n, n))
    idx = 0
    for i in range(n):
        full[i, : i + 1] = tri[idx : idx + i + 1]
        idx += i + 1
    return full + np.tril(full, -1).T


class GaussianReferenceJob:
    """Run one Gaussian single point per geometry and cache the result.

    The cache holds the most recent geometry only; the external loop asks for the
    energy, then the Hessian, then the dipole of the same geometry, and LBFGS asks
    for forces at a new geometry each step. A ``force`` job is upgraded to ``freq``
    in place when the Hessian or dipole derivatives are requested.
    """

    def __init__(
        self,
        method: str = REFERENCE_METHOD,
        basis: str = REFERENCE_BASIS,
        nproc: int = 4,
        mem: str = "2GB",
        keep_files: bool = False,
        workdir: str | Path | None = None,
    ) -> None:
        self.method = method
        self.basis = basis
        self.nproc = nproc
        self.mem = mem
        self.keep_files = keep_files
        self.workdir = Path(workdir) if workdir is not None else None
        self.n_jobs = {"force": 0, "freq": 0}
        self._cache_key: tuple | None = None
        self._cache: ReferenceResult | None = None

    @staticmethod
    def available() -> bool:
        return shutil.which("g16") is not None and shutil.which("formchk") is not None

    @property
    def route(self) -> str:
        return f"{self.method}/{self.basis}"

    def compute(self, atoms, need_hessian: bool) -> ReferenceResult:
        key = self._key(atoms)
        cached = self._cache is not None and self._cache_key == key
        if cached and (self._cache.level == "freq" or not need_hessian):
            return self._cache
        level = "freq" if need_hessian else "force"
        result = self._run(atoms, level)
        self._cache_key, self._cache = key, result
        return result

    @staticmethod
    def _key(atoms) -> tuple:
        return (
            tuple(atoms.get_chemical_symbols()),
            np.round(atoms.get_positions(), 10).tobytes(),
            int(atoms.info.get("charge", 0)),
            int(atoms.info.get("spin", 1)),
        )

    def _run(self, atoms, level: str) -> ReferenceResult:
        charge = int(atoms.info.get("charge", 0))
        multiplicity = int(atoms.info.get("spin", 1))
        positions = atoms.get_positions()
        keyword = "freq" if level == "freq" else "force"

        if self.workdir is not None:
            self.workdir.mkdir(parents=True, exist_ok=True)
        tmp = tempfile.mkdtemp(prefix=f"gref_{level}_", dir=self.workdir)
        tmp_path = Path(tmp)
        try:
            gjf = tmp_path / "ref.gjf"
            with gjf.open("w") as f:
                f.write(f"%chk=ref.chk\n%mem={self.mem}\n%NProcShared={self.nproc}\n")
                f.write(f"# {self.route} {keyword} nosymm\n\n")
                f.write(f"reference {level} job\n\n{charge} {multiplicity}\n")
                for s, p in zip(atoms.get_chemical_symbols(), positions):
                    f.write(f"{s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")
                f.write("\n")

            with (tmp_path / "console.txt").open("wb") as console:
                proc = subprocess.run(
                    ["g16", "ref.gjf"], cwd=tmp, stdout=console, stderr=subprocess.STDOUT
                )
            log = tmp_path / "ref.log"
            log_text = log.read_text(errors="replace") if log.exists() else ""
            if proc.returncode != 0 or "Normal termination" not in log_text:
                raise RuntimeError(
                    f"Reference Gaussian job failed (exit {proc.returncode}) in {tmp}: "
                    f"{log_text[-2000:]}"
                )
            fchk = convert_chk_to_fchk(str(tmp_path / "ref.chk"), str(tmp_path / "ref.fchk"))
            result = parse_reference_fchk(fchk, expected_positions_A=positions)
        except Exception:
            self.keep_files = True  # leave the evidence in place
            raise
        finally:
            if not self.keep_files:
                shutil.rmtree(tmp, ignore_errors=True)

        self.n_jobs[level] += 1
        logger.debug(
            "Reference %s job done (E=%.8f eV, |F|max=%.2e eV/A)",
            level,
            result.energy,
            np.abs(result.forces).max(),
        )
        return result


_shared_job: GaussianReferenceJob | None = None


def get_shared_job() -> GaussianReferenceJob:
    """The one job instance both adapters use, so a geometry is computed once."""
    global _shared_job
    if _shared_job is None:
        _shared_job = GaussianReferenceJob()
    return _shared_job


def set_shared_job(job: GaussianReferenceJob | None) -> None:
    global _shared_job
    _shared_job = job


class GaussianReferenceCalculator(Calculator):
    """ASE calculator: energy and forces from the reference job, plus ``get_hessian``."""

    implemented_properties: ClassVar[list[str]] = ["energy", "forces"]

    def __init__(self, job: GaussianReferenceJob | None = None, **kwargs) -> None:
        super().__init__(**kwargs)
        self.job = job if job is not None else get_shared_job()

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        result = self.job.compute(self.atoms, need_hessian=False)
        self.results = {"energy": result.energy, "forces": result.forces}

    def get_hessian(self, atoms=None) -> np.ndarray:
        """Analytic Hessian in eV/Angstrom^2, shape (3N, 3N), as ``calculate_hessian`` expects."""
        if atoms is None:
            atoms = self.atoms
        return self.job.compute(atoms, need_hessian=True).hessian


class GaussianReferenceDipoleCalculator(DipoleCalculatorBase):
    """Dipole, dipole derivatives and polarizability from the reference job.

    ``analytic_derivatives=True`` (default) returns Gaussian's analytic dipole
    derivatives; ``False`` falls back to the base-class central differences, which
    exercises the same finite-difference code the ML dipole models go through.
    """

    def __init__(self, job: GaussianReferenceJob | None = None) -> None:
        self._job = job
        self.analytic_derivatives = True
        super().__init__(CALCULATOR_NAME)

    @property
    def job(self) -> GaussianReferenceJob:
        if self._job is None:
            self._job = get_shared_job()
        return self._job

    def _check_availability(self) -> bool:
        self.available = GaussianReferenceJob.available()
        if self.available:
            logger.info("✓ Gaussian reference dipole calculator available")
        else:
            logger.info("✗ Gaussian reference dipole calculator unavailable (no g16/formchk)")
        return self.available

    def calculate_dipole(self, atoms, **kwargs):
        result = self.job.compute(atoms, need_hessian=False)
        return result.dipole.copy(), None

    def calculate_dipole_derivatives(self, atoms, displacement=0.01, **kwargs) -> np.ndarray:
        if not self.analytic_derivatives:
            return super().calculate_dipole_derivatives(atoms, displacement=displacement, **kwargs)
        return self.job.compute(atoms, need_hessian=True).dipole_derivatives.copy()

    def calculate_polarizability(self, atoms) -> np.ndarray:
        """Static polarizability in Angstrom^3, shape (3, 3)."""
        return self.job.compute(atoms, need_hessian=True).polarizability.copy()
