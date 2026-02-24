"""Reference output regression tests for water and CH4 results.json files.

Validates that committed ML pipeline outputs (mace_mp + espaloma) maintain expected
structure and physically reasonable frequency values. Also verifies test marker
infrastructure and coverage reporting.

Requirements covered: TEST-05, TEST-06, TEST-07, TEST-09.
"""

import json

import pytest

# ---------------------------------------------------------------------------
# Water reference output tests (TEST-05)
# ---------------------------------------------------------------------------


class TestWaterReferenceOutput:
    """Validate committed water ML results.json structure and values."""

    def test_water_reference_output_structure(self, water_results_json):
        """Water results.json has expected top-level keys and frequency sections."""
        data = json.loads(water_results_json.read_text())

        assert isinstance(data, dict)

        # Top-level keys
        expected_top_keys = {
            "molecule",
            "energy_calculator",
            "dipole_calculator",
            "calculator_type",
            "timestamp",
            "frequencies",
            "energy_eV",
            "dipole",
            "runtime_s",
            "files",
        }
        assert expected_top_keys.issubset(data.keys()), (
            f"Missing top-level keys: {expected_top_keys - data.keys()}"
        )

        # Frequency sections
        freq = data["frequencies"]
        expected_freq_keys = {"harmonic", "anharmonic", "overtones", "combination_bands"}
        assert expected_freq_keys.issubset(freq.keys()), (
            f"Missing frequency keys: {expected_freq_keys - freq.keys()}"
        )

        # Water has exactly 3 normal modes (3N-6 = 3*3-6 = 3)
        assert len(freq["harmonic"]) == 3, (
            f"Water should have 3 harmonic modes, got {len(freq['harmonic'])}"
        )
        assert len(freq["anharmonic"]) == 3, (
            f"Water should have 3 anharmonic fundamentals, got {len(freq['anharmonic'])}"
        )
        assert len(freq["overtones"]) == 3, (
            f"Water should have 3 first overtones, got {len(freq['overtones'])}"
        )

    def test_water_reference_output_values(self, water_results_json):
        """Water ML frequencies are in physically reasonable ranges."""
        data = json.loads(water_results_json.read_text())

        # Metadata checks
        assert data["molecule"] == "water"
        assert data["energy_calculator"] == "mace_mp"
        assert data["dipole_calculator"] == "espaloma"
        assert data["calculator_type"] == "ml"

        # Harmonic frequencies: ML values are ~1462, ~3852, ~3975 cm-1
        # All should be in 400-4500 cm-1 range for water fundamentals
        harmonic_freqs = [f["freq_cm"] for f in data["frequencies"]["harmonic"]]
        for freq in harmonic_freqs:
            assert 400 <= freq <= 4500, (
                f"Water harmonic frequency {freq} cm-1 outside reasonable range [400, 4500]"
            )

        # Verify specific known values (ML mace_mp/espaloma results)
        # Bending mode ~1462 cm-1, symmetric stretch ~3852, asymmetric stretch ~3975
        assert any(1400 <= f <= 1550 for f in harmonic_freqs), (
            "No water bending mode found in 1400-1550 cm-1 range"
        )
        assert any(3700 <= f <= 4000 for f in harmonic_freqs), (
            "No water stretching mode found in 3700-4000 cm-1 range"
        )

        # IR intensities should be non-negative
        for entry in data["frequencies"]["harmonic"]:
            assert entry["ir_intensity"] >= 0, "IR intensity cannot be negative"

        # Dipole should be non-zero for water (polar molecule)
        assert data["dipole"]["magnitude"] > 0, "Water dipole magnitude should be > 0"


# ---------------------------------------------------------------------------
# CH4 reference output tests (TEST-06)
# ---------------------------------------------------------------------------


class TestCH4ReferenceOutput:
    """Validate committed CH4 ML results.json structure and values."""

    def test_ch4_reference_output_structure(self, ch4_results_json):
        """CH4 results.json has expected top-level keys and frequency sections."""
        data = json.loads(ch4_results_json.read_text())

        assert isinstance(data, dict)

        # Top-level keys
        expected_top_keys = {
            "molecule",
            "energy_calculator",
            "dipole_calculator",
            "calculator_type",
            "frequencies",
            "energy_eV",
            "dipole",
            "runtime_s",
            "files",
        }
        assert expected_top_keys.issubset(data.keys())

        # Frequency sections
        freq = data["frequencies"]
        expected_freq_keys = {"harmonic", "anharmonic", "overtones", "combination_bands"}
        assert expected_freq_keys.issubset(freq.keys())

        # CH4: 3*5-6 = 9 normal modes, but degenerate modes are collapsed in harmonic
        # The JSON shows 4 unique harmonic frequencies (t2 bending x3, e bending x2,
        # a1 stretch x1, t2 stretch x3 -- but stored as 4 unique entries)
        assert len(freq["harmonic"]) == 4, (
            f"CH4 should have 4 unique harmonic frequency entries, got {len(freq['harmonic'])}"
        )

        # 9 anharmonic fundamentals (one per mode)
        assert len(freq["anharmonic"]) == 9, (
            f"CH4 should have 9 anharmonic fundamentals, got {len(freq['anharmonic'])}"
        )

        # 9 first overtones (one per mode)
        assert len(freq["overtones"]) == 9, (
            f"CH4 should have 9 first overtones, got {len(freq['overtones'])}"
        )

    def test_ch4_reference_output_values(self, ch4_results_json):
        """CH4 ML frequencies are in physically reasonable ranges."""
        data = json.loads(ch4_results_json.read_text())

        # Metadata checks
        assert data["molecule"] == "CH4_ase"
        assert data["energy_calculator"] == "mace_mp"
        assert data["dipole_calculator"] == "espaloma"
        assert data["calculator_type"] == "ml"

        # Harmonic frequencies: ML values are ~1163, ~1410, ~3076, ~3163 cm-1
        harmonic_freqs = [f["freq_cm"] for f in data["frequencies"]["harmonic"]]
        for freq in harmonic_freqs:
            assert 800 <= freq <= 3500, (
                f"CH4 harmonic frequency {freq} cm-1 outside reasonable range [800, 3500]"
            )

        # CH4 has fundamentals around ~1300 (bending) and ~3000 (stretching)
        assert any(1100 <= f <= 1500 for f in harmonic_freqs), (
            "No CH4 bending mode found in 1100-1500 cm-1 range"
        )
        assert any(2900 <= f <= 3300 for f in harmonic_freqs), (
            "No CH4 stretching mode found in 2900-3300 cm-1 range"
        )

        # IR intensities should be non-negative
        for entry in data["frequencies"]["harmonic"]:
            assert entry["ir_intensity"] >= 0, "IR intensity cannot be negative"

        # CH4 is nonpolar -- dipole magnitude should be 0
        assert data["dipole"]["magnitude"] == 0.0, (
            f"CH4 dipole should be 0 (nonpolar), got {data['dipole']['magnitude']}"
        )


# ---------------------------------------------------------------------------
# Fixture existence meta-tests
# ---------------------------------------------------------------------------


class TestReferenceOutputsCommitted:
    """Verify that reference fixture files are actually present on disk."""

    def test_reference_outputs_are_committed(self, fixtures_dir):
        """Both water and CH4 results.json fixtures exist."""
        water_json = fixtures_dir / "water" / "results.json"
        ch4_json = fixtures_dir / "CH4_ase" / "results.json"

        assert water_json.exists(), f"Water results.json not found at {water_json}"
        assert ch4_json.exists(), f"CH4 results.json not found at {ch4_json}"

    def test_reference_outputs_are_valid_json(self, fixtures_dir):
        """Both results.json files are parseable JSON."""
        for molecule_dir in ["water", "CH4_ase"]:
            json_path = fixtures_dir / molecule_dir / "results.json"
            data = json.loads(json_path.read_text())
            assert isinstance(data, dict), f"{molecule_dir}/results.json is not a JSON object"


# ---------------------------------------------------------------------------
# Test marker infrastructure (TEST-07)
# ---------------------------------------------------------------------------

# Marker filtering is verified by running:
#   pytest -m "not gpu and not gaussian and not slow" --co -q
# and confirming all collected tests are unit tests (no special markers).


def test_unit_test_baseline():
    """This test should always be collected when running with -m 'not gpu and not gaussian'."""
    assert True


@pytest.mark.gpu
def test_gpu_placeholder():
    """Placeholder GPU test -- verifies marker infrastructure works."""
    pytest.skip("No GPU test implemented yet -- marker infrastructure verification only")


@pytest.mark.gaussian
def test_gaussian_placeholder():
    """Placeholder Gaussian test -- verifies marker infrastructure works."""
    pytest.skip("No Gaussian test implemented yet -- marker infrastructure verification only")


# ---------------------------------------------------------------------------
# Results manager version metadata tests (Phase 2)
# ---------------------------------------------------------------------------


class TestResultsManagerMetadata:
    """Verify that ResultsManager embeds version_info in saved JSON."""

    def test_frequency_results_contain_version_info(self, tmp_path):
        """save_frequency_results() writes version_info to results.json."""
        from mace_gaussian.utils.results import ResultsManager

        mgr = ResultsManager(base_output_dir=str(tmp_path))
        mgr.save_frequency_results(
            molecule_name="test_mol",
            energy_calculator="mace_mp",
            dipole_calculator="espaloma",
            calculator_type="ml",
            frequencies_data={
                "harmonic": [{"freq_cm": 1000.0, "ir_intensity": 10.0}],
                "anharmonic": [],
                "overtones": [],
                "combination_bands": [],
            },
            energy=0.0,
            dipole=None,
            runtime=1.0,
        )

        results_path = tmp_path / "test_mol" / "mace_mp_espaloma" / "results.json"
        assert results_path.exists()

        data = json.loads(results_path.read_text())
        assert "version_info" in data
        assert "tool_version" in data["version_info"]
        assert "python_version" in data["version_info"]
        assert "platform" in data["version_info"]

    def test_frequency_results_contain_calculation_parameters(self, tmp_path):
        """save_frequency_results() writes calculation_parameters to results.json."""
        from mace_gaussian.utils.results import ResultsManager

        mgr = ResultsManager(base_output_dir=str(tmp_path))
        params = {"basis": "6-31G(d,p)", "method": "B3LYP"}
        mgr.save_frequency_results(
            molecule_name="test_mol",
            energy_calculator="dft",
            dipole_calculator="dft",
            calculator_type="dft",
            frequencies_data={
                "harmonic": [],
                "anharmonic": [],
                "overtones": [],
                "combination_bands": [],
            },
            energy=0.0,
            dipole=None,
            runtime=1.0,
            calculation_parameters=params,
        )

        results_path = tmp_path / "test_mol" / "dft" / "results.json"
        data = json.loads(results_path.read_text())
        assert data["calculation_parameters"] == params

    def test_optimization_results_contain_version_info(self, tmp_path):
        """save_optimization_results() writes version_info to results.json."""
        import numpy as np
        from ase import Atoms

        from mace_gaussian.utils.results import ResultsManager

        # Create real ASE Atoms objects (mock fails ase.io.write)
        atoms = Atoms(
            "OH2",
            positions=np.array(
                [[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]]
            ),
        )

        mgr = ResultsManager(base_output_dir=str(tmp_path))
        mgr.save_optimization_results(
            molecule_name="test_mol",
            initial_atoms=atoms,
            final_atoms=atoms,
            calculator_name="mace_mp",
            initial_energy=-76.0,
            final_energy=-76.1,
            converged=True,
            num_steps=42,
            runtime=5.0,
        )

        results_path = tmp_path / "test_mol" / "geometry_opt" / "results.json"
        assert results_path.exists()

        data = json.loads(results_path.read_text())
        assert "version_info" in data
        assert "tool_version" in data["version_info"]
        assert data["num_steps"] == 42
