"""Tests for the Gaussian reference calculator (harness self-consistency check).

Everything here runs without g16: the fixture ``b3lyp_freq_nosymm.fchk`` is a
B3LYP/6-31G(d,p) ``freq nosymm`` water job, and a fake job loads it instead of
running Gaussian. The round-trip test is the offline half of the check: fchk
atomic units -> ASE-side units -> ``write_gaussian_output`` -> atomic units must
reproduce the fchk numbers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from mace_gaussian.calculators.gaussian_reference import (
    CALCULATOR_NAME,
    GaussianReferenceCalculator,
    GaussianReferenceDipoleCalculator,
    GaussianReferenceJob,
    ReferenceResult,
    _lower_triangle_to_full,
    parse_reference_fchk,
)
from mace_gaussian.gaussian.io import write_gaussian_output
from mace_gaussian.utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV

FCHK = Path(__file__).parent / "fixtures" / "water" / "b3lyp_freq_nosymm.fchk"

# Input geometry of the fixture job (Angstrom), as written in the .gjf.
WATER_POSITIONS = np.array(
    [
        [0.0, 0.0, 0.11779],
        [0.0, 0.75545, -0.47116],
        [0.0, -0.75545, -0.47116],
    ]
)
FIXTURE_ENERGY_HARTREE = -76.41963434002696


def _water() -> Atoms:
    atoms = Atoms("OHH", positions=WATER_POSITIONS)
    atoms.info["charge"] = 0.0
    atoms.info["spin"] = 1.0
    return atoms


class FakeJob(GaussianReferenceJob):
    """Loads the fixture instead of running g16; counts jobs like the real one."""

    def _run(self, atoms, level):
        result = parse_reference_fchk(FCHK, expected_positions_A=atoms.get_positions())
        if level == "force":
            result = ReferenceResult(
                energy=result.energy,
                forces=result.forces,
                dipole=result.dipole,
                hessian=None,
                dipole_derivatives=None,
                polarizability=None,
                level="force",
            )
        self.n_jobs[level] += 1
        return result


class TestParseReferenceFchk:
    def test_units_shapes_and_frame(self):
        r = parse_reference_fchk(FCHK, expected_positions_A=WATER_POSITIONS)
        assert r.level == "freq"
        assert r.energy == pytest.approx(FIXTURE_ENERGY_HARTREE * HARTREE_TO_EV, rel=1e-12)
        assert r.forces.shape == (3, 3)
        assert r.hessian.shape == (9, 9)
        assert np.allclose(r.hessian, r.hessian.T)
        # Translational invariance: each row sums to zero over the atoms.
        row_sums = r.hessian.reshape(9, 3, 3).sum(axis=1)
        assert np.allclose(row_sums, 0.0, atol=1e-6)
        # Stretch curvature (H y and z entries) is large and positive.
        assert np.diag(r.hessian)[4] > 30
        assert r.dipole.shape == (3,)
        assert r.dipole_derivatives.shape == (9, 3)
        assert r.polarizability.shape == (3, 3)
        assert np.allclose(r.polarizability, r.polarizability.T)
        # Water in the yz plane, dipole along -z, about 2 Debye = 0.8 e*Bohr.
        assert r.dipole[0] == pytest.approx(0.0, abs=1e-8)
        assert r.dipole[2] == pytest.approx(-0.80, abs=0.02)

    def test_dipole_derivative_layout_is_row_per_coordinate(self):
        """Row 3*i + k holds d(mu_x, mu_y, mu_z)/dR_{i,k}.

        Checked against central differences of the Gaussian dipole with H1
        displaced +-0.005 A along y (2026-09-18 probe): d mu_y/dy = 0.208 e and
        d mu_z/dy = 0.087 e. The other layout would put 0.119 in the z slot.
        """
        r = parse_reference_fchk(FCHK)
        h1_y = r.dipole_derivatives[3 * 1 + 1]
        assert h1_y[0] == pytest.approx(0.0, abs=1e-8)
        assert h1_y[1] == pytest.approx(0.2077, abs=2e-3)
        assert h1_y[2] == pytest.approx(0.0865, abs=2e-3)

    def test_dipole_derivative_sum_rule(self):
        """Neutral molecule: summing d mu/dR over atoms gives the zero charge tensor."""
        r = parse_reference_fchk(FCHK)
        per_atom = r.dipole_derivatives.reshape(3, 3, 3)  # atom, coordinate, dipole component
        assert np.allclose(per_atom.sum(axis=0), 0.0, atol=1e-6)

    def test_frame_check_rejects_reoriented_job(self):
        shifted = WATER_POSITIONS + np.array([0.1, 0.0, 0.0])
        with pytest.raises(ValueError, match="reoriented"):
            parse_reference_fchk(FCHK, expected_positions_A=shifted)

    def test_lower_triangle_roundtrip(self):
        rng = np.random.default_rng(0)
        a = rng.normal(size=(6, 6))
        sym = a + a.T
        tri = np.concatenate([sym[i, : i + 1] for i in range(6)])
        assert np.allclose(_lower_triangle_to_full(tri, 6), sym)
        with pytest.raises(ValueError):
            _lower_triangle_to_full(tri[:-1], 6)


class TestHarnessRoundTrip:
    def test_write_gaussian_output_reproduces_fchk_atomic_units(self, tmp_path):
        """ASE units -> write_gaussian_output -> file must equal the fchk values.

        This applies exactly the conversions run_next_calculation applies to an
        ML model's output (energy, gradient, Hessian) and checks the file Gaussian
        would read against the original atomic-unit numbers.
        """
        from mace_gaussian.gaussian.fchk import parse_fchk_section

        r = parse_reference_fchk(FCHK)
        natoms = 3
        # Same conversion as workflow.calculate_hessian applies to calculator.get_hessian()
        hessian_for_gaussian = r.hessian * (BOHR_TO_ANGSTROM**2) / HARTREE_TO_EV
        polar_bohr3 = r.polarizability / BOHR_TO_ANGSTROM**3
        polar6 = np.array(
            [
                polar_bohr3[0, 0],
                polar_bohr3[0, 1],
                polar_bohr3[1, 1],
                polar_bohr3[0, 2],
                polar_bohr3[1, 2],
                polar_bohr3[2, 2],
            ]
        )
        out = tmp_path / "gau.out"
        write_gaussian_output(
            str(out),
            natoms,
            r.energy,
            -r.forces,
            r.dipole,
            r.dipole_derivatives,
            hessian_for_gaussian,
            deriv=2,
            polarizability=polar6,
        )

        vals = [float(x) for x in out.read_text().replace("D", "E").split()]
        content = FCHK.read_text()
        expected = np.concatenate(
            [
                [FIXTURE_ENERGY_HARTREE],
                parse_fchk_section(content, "Dipole Moment"),
                parse_fchk_section(content, "Cartesian Gradient"),
                parse_fchk_section(content, "Polarizability"),
                parse_fchk_section(content, "Dipole Derivatives"),
                parse_fchk_section(content, "Cartesian Force Constants"),
            ]
        )
        assert len(vals) == len(expected)
        np.testing.assert_allclose(vals, expected, rtol=1e-9, atol=1e-12)


class TestAdapters:
    def test_energy_calculator_and_dipole_share_one_job(self):
        job = FakeJob()
        atoms = _water()
        atoms.calc = GaussianReferenceCalculator(job=job)
        forces = atoms.get_forces()
        assert forces.shape == (3, 3)
        assert job.n_jobs == {"force": 1, "freq": 0}

        hessian = atoms.calc.get_hessian(atoms)
        assert hessian.shape == (9, 9)
        assert job.n_jobs == {"force": 1, "freq": 1}

        dip = GaussianReferenceDipoleCalculator(job=job)
        dipole, charges = dip.calculate_dipole(atoms)
        derivs = dip.calculate_dipole_derivatives(atoms, displacement=0.005)
        polar = dip.calculate_polarizability(atoms)
        assert charges is None
        assert dipole.shape == (3,) and derivs.shape == (9, 3) and polar.shape == (3, 3)
        # All served from the cached freq job.
        assert job.n_jobs == {"force": 1, "freq": 1}

    def test_new_geometry_invalidates_cache(self):
        job = FakeJob()
        atoms = _water()
        atoms.calc = GaussianReferenceCalculator(job=job)
        atoms.get_forces()
        atoms.positions[1, 1] += 0.01
        with pytest.raises(ValueError, match="reoriented"):
            # FakeJob checks the fixture frame, so a moved atom is detected: proves
            # the cache key changed and _run was called again.
            atoms.get_forces()
        assert job.n_jobs["force"] == 1  # second call raised before counting

    def test_finite_difference_derivatives_use_base_class(self):
        job = FakeJob()
        dip = GaussianReferenceDipoleCalculator(job=job)
        dip.analytic_derivatives = False
        atoms = _water()
        # Base-class central differences displace atoms, which FakeJob rejects
        # (frame check): the point is only that the base-class path is taken.
        with pytest.raises(ValueError, match="reoriented"):
            dip.calculate_dipole_derivatives(atoms, displacement=0.005)

    def test_registered_in_factory_and_workflow(self):
        from mace_gaussian.calculators import dipole_factory

        assert CALCULATOR_NAME in dipole_factory.calculators
        assert isinstance(
            dipole_factory.calculators[CALCULATOR_NAME], GaussianReferenceDipoleCalculator
        )

        from mace_gaussian.workflow import calculator

        assert isinstance(calculator(CALCULATOR_NAME), GaussianReferenceCalculator)
