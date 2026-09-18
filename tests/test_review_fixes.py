"""Regression tests for the 2026-09 review fixes (docs/explained/REVIEW_FINDINGS.md).

H2: anharmonic mode IDs live in .fchk (ascending-frequency) index space.
H3: mode vectors are mass-weighted, so a calculation overlapped with itself is the identity.
H4: parse_final_energy returns the energy, not a thermal correction.
"""

import json
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
