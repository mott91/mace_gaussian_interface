"""Tests for GaussianLogParser covering harmonic, anharmonic, overtone,
combination band, energy parsing, and edge cases.

Covers requirements: TEST-01 (harmonic parser), TEST-02 (anharmonic/overtone/combo),
TEST-08 (acoh bug documentation).
"""

import pytest

from mace_gaussian.gaussian.parser import GaussianLogParser
from mace_gaussian.utils.exceptions import GaussianParseError

# --- Expected values (verified against actual parser output from trimmed fixtures) ---

# Water B3LYP/6-31G(d,p) harmonic frequencies
WATER_HARMONIC_EXPECTED = [
    {"freq_cm": 1665.3001, "ir_intensity": 70.3304},
    {"freq_cm": 3799.2203, "ir_intensity": 1.6374},
    {"freq_cm": 3912.4174, "ir_intensity": 20.2278},
]

# Water B3LYP/6-31G(d,p) anharmonic fundamentals
WATER_ANHARMONIC_EXPECTED = [
    {"mode": 1, "freq_cm": 3624.001, "ir_intensity": 0.62829225, "freq_harmonic": 3799.22},
    {"mode": 2, "freq_cm": 1615.153, "ir_intensity": 69.29467704, "freq_harmonic": 1665.3},
    {"mode": 3, "freq_cm": 3722.273, "ir_intensity": 16.09760189, "freq_harmonic": 3912.417},
]

# Water B3LYP/6-31G(d,p) overtones
WATER_OVERTONES_EXPECTED = [
    {
        "mode": 1,
        "overtone_level": 2,
        "freq_harmonic": 7598.441,
        "freq_anharmonic": 7161.57,
        "ir_intensity": 0.77569767,
    },
    {
        "mode": 2,
        "overtone_level": 2,
        "freq_harmonic": 3330.6,
        "freq_anharmonic": 3194.746,
        "ir_intensity": 0.80602582,
    },
    {
        "mode": 3,
        "overtone_level": 2,
        "freq_harmonic": 7824.835,
        "freq_anharmonic": 7346.675,
        "ir_intensity": 0.04914658,
    },
]

# Water energy from archive section
WATER_ENERGY_EXPECTED = 0.003705  # Hartrees


# --- Harmonic frequency tests ---


class TestParseHarmonicFrequencies:
    """Tests for parse_harmonic_frequencies() method."""

    def test_water_harmonic_count_and_values(self, water_dft_log):
        """Water DFT should have exactly 3 harmonic modes with correct freq and IR values."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_harmonic_frequencies()

        assert len(result) == 3

        for actual, expected in zip(result, WATER_HARMONIC_EXPECTED):
            assert actual["freq_cm"] == pytest.approx(expected["freq_cm"], abs=0.01)
            assert actual["ir_intensity"] == pytest.approx(expected["ir_intensity"], abs=0.01)

    def test_ch4_harmonic_degenerate_preserved(self, ch4_dft_log):
        """CH4 DFT has 9 vibrational modes (3*5-6=9). Degenerate modes are
        preserved -- the parser deduplicates repeated frequency blocks (Gaussian
        prints them twice for anharmonic calc) but keeps all modes within a block.

        CH4 symmetry (Td) produces:
          - 1 A1 mode at 3046.54 cm-1
          - 1 E mode (2-fold) at 1578.57 cm-1 (IR-inactive)
          - 1 T2 mode (3-fold) at 1356.19 cm-1
          - 1 T2 mode (3-fold) at 3162.41 cm-1
          Total: 1 + 2 + 3 + 3 = 9 modes
        """
        parser = GaussianLogParser(ch4_dft_log)
        result = parser.parse_harmonic_frequencies()

        # Parser preserves all 9 vibrational modes including degenerates
        assert len(result) == 9

        # All frequencies must be positive (no imaginary modes)
        for entry in result:
            assert entry["freq_cm"] > 0

    @pytest.mark.parametrize(
        "molecule,log_fixture,expected_unique_count",
        [
            ("water", "water_dft_log", 3),
            ("CH4", "ch4_dft_log", 9),  # 9 modes preserved (degenerate modes kept)
        ],
    )
    def test_harmonic_mode_count_parametrized(
        self, molecule, log_fixture, expected_unique_count, request
    ):
        """Verify harmonic mode count per molecule.

        Gaussian logs repeat frequency blocks (initial + anharmonic pre-analysis).
        The parser deduplicates repeated blocks but preserves degenerate modes
        within each block.
        """
        log_path = request.getfixturevalue(log_fixture)
        parser = GaussianLogParser(log_path)
        result = parser.parse_harmonic_frequencies()
        assert len(result) == expected_unique_count

    def test_ml_log_harmonic_parsing(self, water_ml_log):
        """Parser works on ML-generated Gaussian logs (external calculation with ML dipoles).

        This verifies the parser is not DFT-specific and correctly handles logs
        produced via the ZMQ external helper script.
        """
        parser = GaussianLogParser(water_ml_log)
        result = parser.parse_harmonic_frequencies()

        assert len(result) == 3

        # ML frequencies differ from DFT but should all be positive
        for entry in result:
            assert entry["freq_cm"] > 0
            assert "ir_intensity" in entry


# --- Anharmonic frequency tests ---


class TestParseAnharmonicFrequencies:
    """Tests for parse_anharmonic_frequencies() method."""

    def test_water_anharmonic_count_and_values(self, water_dft_log):
        """Water DFT should have 3 anharmonic fundamentals matching expected values."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_anharmonic_frequencies()

        assert len(result) == 3

        for actual, expected in zip(result, WATER_ANHARMONIC_EXPECTED):
            assert actual["mode"] == expected["mode"]
            assert actual["freq_cm"] == pytest.approx(expected["freq_cm"], abs=0.01)
            assert actual["ir_intensity"] == pytest.approx(expected["ir_intensity"], abs=0.001)
            assert actual["freq_harmonic"] == pytest.approx(expected["freq_harmonic"], abs=0.01)

    def test_anharmonic_keys_present(self, water_dft_log):
        """Each anharmonic entry must have mode, freq_cm, ir_intensity, freq_harmonic keys."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_anharmonic_frequencies()

        required_keys = {"mode", "freq_cm", "ir_intensity", "freq_harmonic"}
        for entry in result:
            assert required_keys.issubset(entry.keys())


# --- Overtone tests ---


class TestParseOvertones:
    """Tests for parse_overtones() method."""

    def test_water_overtones_count_and_values(self, water_dft_log):
        """Water DFT should have 3 overtones matching expected values."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_overtones()

        assert len(result) == 3

        for actual, expected in zip(result, WATER_OVERTONES_EXPECTED):
            assert actual["mode"] == expected["mode"]
            assert actual["overtone_level"] == expected["overtone_level"]
            assert actual["freq_harmonic"] == pytest.approx(expected["freq_harmonic"], abs=0.01)
            assert actual["freq_anharmonic"] == pytest.approx(expected["freq_anharmonic"], abs=0.01)
            assert actual["ir_intensity"] == pytest.approx(expected["ir_intensity"], abs=0.001)

    def test_overtone_keys_present(self, water_dft_log):
        """Each overtone entry must have expected keys."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_overtones()

        required_keys = {
            "mode",
            "overtone_level",
            "freq_harmonic",
            "freq_anharmonic",
            "ir_intensity",
        }
        for entry in result:
            assert required_keys.issubset(entry.keys())


# --- Combination band tests ---


class TestParseCombinationBands:
    """Tests for parse_combination_bands() method."""

    def test_water_combination_bands_documents_duplication_bug(self, water_dft_log):
        """Water DFT combination bands return 6 entries (NOT 3).

        Known bug: parser captures combination bands from both I(anharm) and
        DS(anharm) sections, yielding duplicates. True unique count is 3.
        See Phase 1 research.
        """
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_combination_bands()

        # Documents current (buggy) behavior: 6 instead of 3
        assert len(result) == 6

        # Verify all entries have expected keys
        required_keys = {"mode1", "mode2", "freq_harmonic", "freq_anharmonic", "ir_intensity"}
        for entry in result:
            assert required_keys.issubset(entry.keys())

    def test_combination_band_frequencies_positive(self, water_dft_log):
        """All combination band frequencies should be positive."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_combination_bands()

        for entry in result:
            assert entry["freq_harmonic"] > 0
            assert entry["freq_anharmonic"] > 0
            assert entry["ir_intensity"] >= 0


# --- Energy parsing tests ---


class TestParseFinalEnergy:
    """Tests for parse_final_energy() method."""

    def test_water_energy_extraction(self, water_dft_log):
        """Energy parsing extracts Hartree energy from 'Energy=' pattern in archive section."""
        parser = GaussianLogParser(water_dft_log)
        result = parser.parse_final_energy()

        assert result is not None
        assert result == pytest.approx(WATER_ENERGY_EXPECTED, abs=0.000001)


# --- Edge case and error handling tests ---


class TestEdgeCases:
    """Tests for error handling and edge cases."""

    def test_parser_file_not_found(self):
        """GaussianLogParser should raise FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            GaussianLogParser("nonexistent_file_that_does_not_exist.log")

    def test_acoh_anharmonic_parsing(self, acoh_ml_log):
        """Acetic acid ML log uses 'Vibrational Energies at Anharmonic Level' format
        (Format B) with H/L overlap prefix indicators. Parser now handles both formats.

        Expected: 18 modes (8 atoms -> 3*8-6 = 18)
        """
        parser = GaussianLogParser(acoh_ml_log)
        result = parser.parse_anharmonic_frequencies()
        assert len(result) == 18  # 8 atoms -> 3*8-6 = 18 modes
        for entry in result:
            assert "mode" in entry
            assert "freq_cm" in entry
            assert "freq_harmonic" in entry
            assert "ir_intensity" in entry


# --- Parser error handling tests (Phase 2) ---


MINIMAL_HARMONIC_LOG = """\
 Frequencies --   1595.1234  3657.0567  3755.9876
 IR Inten    --     75.1234   2.5678    44.3210
"""

EMPTY_LOG = """\
Gaussian 16 output
Normal termination of Gaussian 16.
"""


class TestParserErrorHandling:
    """Tests for GaussianParseError behavior in parser methods."""

    def test_harmonic_raises_on_empty_log(self, tmp_path):
        """parse_harmonic_frequencies() raises GaussianParseError when no frequencies found."""
        log_file = tmp_path / "empty.log"
        log_file.write_text(EMPTY_LOG)

        parser = GaussianLogParser(str(log_file))
        with pytest.raises(GaussianParseError, match="No harmonic frequencies found"):
            parser.parse_harmonic_frequencies()

    def test_anharmonic_default_returns_empty(self, tmp_path):
        """parse_anharmonic_frequencies() (default strict=False) returns empty list
        when no Fundamental Bands section exists."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        result = parser.parse_anharmonic_frequencies()
        assert result == []

    def test_anharmonic_strict_raises(self, tmp_path):
        """parse_anharmonic_frequencies(strict=True) raises GaussianParseError
        when no Fundamental Bands section exists."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        with pytest.raises(GaussianParseError, match="No Fundamental Bands"):
            parser.parse_anharmonic_frequencies(strict=True)

    def test_overtones_default_returns_empty(self, tmp_path):
        """parse_overtones() (default strict=False) returns empty list
        when no Overtones section exists."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        result = parser.parse_overtones()
        assert result == []

    def test_overtones_strict_raises(self, tmp_path):
        """parse_overtones(strict=True) raises GaussianParseError."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        with pytest.raises(GaussianParseError, match="No Overtones section"):
            parser.parse_overtones(strict=True)

    def test_combination_bands_default_returns_empty(self, tmp_path):
        """parse_combination_bands() (default strict=False) returns empty list
        when no Combination Bands section exists."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        result = parser.parse_combination_bands()
        assert result == []

    def test_combination_bands_strict_raises(self, tmp_path):
        """parse_combination_bands(strict=True) raises GaussianParseError."""
        log_file = tmp_path / "harmonic_only.log"
        log_file.write_text(MINIMAL_HARMONIC_LOG)

        parser = GaussianLogParser(str(log_file))
        with pytest.raises(GaussianParseError, match="No Combination Bands"):
            parser.parse_combination_bands(strict=True)
