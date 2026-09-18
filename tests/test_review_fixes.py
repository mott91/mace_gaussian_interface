"""Regression tests for the 2026-09 review fixes (docs/explained/REVIEW_FINDINGS.md).

H2: anharmonic mode IDs live in .fchk (ascending-frequency) index space.
H3: mode vectors are mass-weighted, so a calculation overlapped with itself is the identity.
H4: parse_final_energy returns the energy, not a thermal correction.
M1: dipole failures are counted and flagged, not silently zeroed.
M2: the Gaussian deadline is enforced inside the wait loop.
M3: the wait loop returns as soon as a request arrives.
M4: geometry_optimisation reports the optimizer's own convergence verdict.
M5: batch leaderboard metrics come from the eigenvector-matched CSVs.
"""

import json
import time
from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).parent / "fixtures"


# --- H2 -----------------------------------------------------------------------


class TestGaussianModeToCheckpointIndex:
    def test_water_symmetry_order_is_translated(self):
        from mace_gaussian.analysis.analyze_spectra import gaussian_mode_to_checkpoint_index

        # Gaussian's numbering for water: Mode(1)=3799, Mode(2)=1665, Mode(3)=3912
        anh = [
            {"mode": 1, "freq_harmonic": 3799.22, "freq_cm": 3624.0},
            {"mode": 2, "freq_harmonic": 1665.30, "freq_cm": 1615.2},
            {"mode": 3, "freq_harmonic": 3912.42, "freq_cm": 3722.3},
        ]
        assert gaussian_mode_to_checkpoint_index(anh) == {2: 1, 1: 2, 3: 3}

    def test_missing_freq_harmonic_falls_back_to_identity(self):
        from mace_gaussian.analysis.analyze_spectra import gaussian_mode_to_checkpoint_index

        assert gaussian_mode_to_checkpoint_index([{"mode": 1, "freq_cm": 1.0}]) == {}
        assert gaussian_mode_to_checkpoint_index([]) == {}

    def test_spectrum_ids_match_fchk_order(self):
        """Anharmonic F-labels must name the same modes as the .fchk harmonic order."""
        from mace_gaussian.analysis.analyze_spectra import SpectrumAnalyzer

        results = json.loads((FIXTURES / "water" / "results.json").read_text())
        spec = SpectrumAnalyzer().extract_spectrum_data(
            results, include_overtones=True, include_combinations=True, use_harmonic=False
        )
        fund = {
            mid: f
            for mid, f, lab in zip(spec.mode_ids, spec.frequencies, spec.labels)
            if lab == "fundamental"
        }
        # F1 must be the lowest-harmonic mode (the bend), F3 the highest
        harm_by_mode = {e["mode"]: e["freq_harmonic"] for e in results["frequencies"]["anharmonic"]}
        lowest = min(harm_by_mode, key=harm_by_mode.get)
        highest = max(harm_by_mode, key=harm_by_mode.get)
        anh_by_mode = {e["mode"]: e["freq_cm"] for e in results["frequencies"]["anharmonic"]}
        assert fund["F1"] == pytest.approx(anh_by_mode[lowest])
        assert fund["F3"] == pytest.approx(anh_by_mode[highest])


# --- H3 -----------------------------------------------------------------------


class TestMassWeightedOverlap:
    @pytest.mark.parametrize("fchk", ["water/dft_b3lyp.fchk", "CH4_ase/dft_b3lyp.fchk"])
    def test_self_overlap_is_identity(self, fchk):
        from mace_gaussian.analysis.mode_matching import (
            create_alignment_matrix,
            extract_mode_data_from_checkpoint,
        )

        modes, *_ = extract_mode_data_from_checkpoint(str(FIXTURES / fchk), force_harmonic=True)
        S = create_alignment_matrix(modes, modes)
        assert np.allclose(S, np.eye(len(modes)), atol=1e-6)


# --- H4 -----------------------------------------------------------------------


class TestParseFinalEnergyAnchored:
    def test_external_log_returns_energy_not_thermal_correction(self):
        from mace_gaussian.gaussian.parser import GaussianLogParser

        parser = GaussianLogParser(str(FIXTURES / "water" / "ml_mace_mp_esp.log"))
        # Fixture line 6: " Energy=   -0.520130339     NIter=   0."
        assert parser.parse_final_energy() == pytest.approx(-0.520130339)

    def test_dft_scf_done_line(self, tmp_path):
        from mace_gaussian.gaussian.parser import GaussianLogParser

        log = tmp_path / "dft.log"
        log.write_text(
            " SCF Done:  E(RB3LYP) =  -76.4196339639     A.U. after    7 cycles\n"
            " Thermal correction to Energy=                    0.024197\n"
            " Thermal correction to Gibbs Free Energy=         0.003705\n"
        )
        assert GaussianLogParser(str(log)).parse_final_energy() == pytest.approx(-76.4196339639)

    def test_truncated_dft_fixture_has_no_energy(self):
        """The DFT fixture is cut before SCF Done; the old regex 'found' 0.003705 there."""
        from mace_gaussian.gaussian.parser import GaussianLogParser

        parser = GaussianLogParser(str(FIXTURES / "water" / "dft_b3lyp.log"))
        assert parser.parse_final_energy() is None


# --- M1 -----------------------------------------------------------------------


class _BrokenDipole:
    name = "broken"

    def calculate_dipole(self, atoms, **kw):
        raise RuntimeError("CUDA out of memory")


class _CountingDipole:
    """Base-class finite differences with a dipole that fails on the third call."""

    name = "counting"

    def __init__(self):
        self.calls = 0

    def calculate_dipole(self, atoms, **kw):
        self.calls += 1
        if self.calls >= 3:
            raise RuntimeError("boom")
        return np.zeros(3), None


class TestDipoleFallbackIsFlagged:
    def test_failure_is_counted_on_atoms(self):
        from ase import Atoms

        from mace_gaussian.workflow import calculate_dipole_properties

        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
        dip, derivs, _charges, _pol = calculate_dipole_properties(
            atoms, _BrokenDipole(), deriv=2, calculate_derivatives=True
        )
        assert np.all(dip == 0) and np.all(derivs == 0)
        assert atoms.info["dipole_fallback_count"] == 1
        assert "CUDA out of memory" in atoms.info["dipole_fallback_last_error"]

    def test_finite_difference_failure_propagates(self):
        """base.py used to swallow the exception and return a half-filled zero array."""
        from ase import Atoms

        from mace_gaussian.calculators.base import DipoleCalculatorBase

        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
        calc = _CountingDipole()
        with pytest.raises(RuntimeError, match="boom"):
            DipoleCalculatorBase.calculate_dipole_derivatives(calc, atoms)
        # positions restored by the finally block
        assert atoms.positions[1, 2] == pytest.approx(0.74)


# --- M2 / M3 ------------------------------------------------------------------


class _FakeProc:
    def __init__(self, returncode=None):
        self.returncode = returncode

    def poll(self):
        return self.returncode


@pytest.fixture
def rep_socket():
    import zmq

    ctx = zmq.Context()
    sock = ctx.socket(zmq.REP)
    sock.setsockopt(zmq.LINGER, 0)
    sock.bind("inproc://review-test")
    yield ctx, sock
    sock.close()
    ctx.term()


class TestWaitLoop:
    def test_deadline_raises_while_idle(self, rep_socket):
        from mace_gaussian.gaussian.zmq_server import is_calc_finished
        from mace_gaussian.utils.exceptions import GaussianTimeoutError

        _, sock = rep_socket
        t0 = time.time()
        with pytest.raises(GaussianTimeoutError):
            is_calc_finished(_FakeProc(None), sock, deadline=time.time() - 1)
        assert time.time() - t0 < 3  # one poll interval, not forever

    def test_exited_process_returns_true(self, rep_socket):
        from mace_gaussian.gaussian.zmq_server import is_calc_finished

        _, sock = rep_socket
        assert is_calc_finished(_FakeProc(0), sock, deadline=time.time() + 60) is True

    def test_message_returns_false_immediately(self, rep_socket):
        import zmq

        from mace_gaussian.gaussian.zmq_server import is_calc_finished

        ctx, sock = rep_socket
        req = ctx.socket(zmq.REQ)
        req.setsockopt(zmq.LINGER, 0)
        req.connect("inproc://review-test")
        req.send_string("in|out")
        t0 = time.time()
        assert is_calc_finished(_FakeProc(None), sock, deadline=time.time() + 60) is False
        assert time.time() - t0 < 0.5  # M3: no unconditional 1 s sleep
        assert sock.recv_string() == "in|out"
        req.close()


# --- M4 -----------------------------------------------------------------------


class TestGeometryOptimisationReportsConvergence:
    def test_returns_optimizer_verdict(self):
        from ase import Atoms
        from ase.calculators.emt import EMT

        from mace_gaussian.workflow import OPT_FMAX, geometry_optimisation

        assert pytest.approx(1e-4) == OPT_FMAX
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.9]])
        atoms.calc = EMT()
        _mol, steps, converged = geometry_optimisation(atoms, fmax=1e-2)
        assert isinstance(converged, bool) and converged
        assert steps > 0


# --- M5 -----------------------------------------------------------------------


class TestBatchMetricsFromMatchedCsv:
    def test_imaginary_pairs_excluded_and_pearson_r2(self, tmp_path):
        import pandas as pd

        from mace_gaussian.analysis.batch_report import _metrics_from_matched_csv

        df = pd.DataFrame(
            {
                "DFT_Frequency_cm": [-274.0, 1000.0, 2000.0, 3000.0],
                "ML_Frequency_cm": [1356.0, 1010.0, 1990.0, 3020.0],
            }
        )
        p = tmp_path / "comparison_x.csv"
        df.to_csv(p, index=False)
        r2, rmse, n = _metrics_from_matched_csv(p)
        assert n == 3  # the -274 row is dropped
        assert rmse == pytest.approx(np.sqrt((100 + 100 + 400) / 3))
        assert 0.99 < r2 <= 1.0

    def test_missing_csv_falls_back_to_sorted(self, tmp_path):
        from mace_gaussian.analysis.batch_report import _compute_combo_metrics

        combo = tmp_path / "mace_x_y"
        combo.mkdir()
        (combo / "results.json").write_text(
            json.dumps({"frequencies": {"harmonic": [{"freq_cm": f} for f in (1010, 1990)]}})
        )
        row = _compute_combo_metrics(
            "mol", "mace_x_y", combo, [1000.0, 2000.0], 3, matched_csv=tmp_path / "nope.csv"
        )
        assert row["pairing"] == "sorted" and row["n_freqs"] == 2
