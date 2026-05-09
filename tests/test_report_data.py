"""Tests for mace_gaussian.analysis.report_data (Phase 23 Wave 0 RED tests)."""

import csv
import json

import numpy as np
import pytest

from mace_gaussian.analysis.report_data import export_report_data


class TestExportReportData:
    def test_writes_json_file(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        assert out.exists()

    def test_json_has_schema_version(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["schema_version"] == 1

    def test_json_has_molecule_and_mode(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["molecule"] == "water"
        assert data["mode"] == "anharmonic"

    def test_json_comparisons_shape(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert len(data["comparisons"]) == 2
        comp = data["comparisons"][0]
        assert comp["name"] == "mace_off_espaloma"
        assert "metrics" in comp
        assert "runtime" in comp
        assert "hardware" in comp
        assert "spectrum_ml" in comp
        assert "spectrum_dft" in comp

    def test_json_metrics_keys(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        metrics = data["comparisons"][0]["metrics"]
        expected = {
            "r2_freq",
            "r2_intensity",
            "rmse_freq",
            "mae_freq",
            "rmse_intensity",
            "mae_intensity",
            "max_error_freq",
            "num_matched",
            "num_dft_only",
            "num_ml_only",
            "match_rate",
        }
        assert set(metrics.keys()) == expected

    def test_json_spectrum_shape(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        spec = data["comparisons"][0]["spectrum_ml"]
        assert "frequencies_cm" in spec
        assert "intensities" in spec
        assert "labels" in spec
        assert "mode_ids" in spec
        assert isinstance(spec["frequencies_cm"], list)

    def test_json_experimental_present(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["experimental"] is not None
        assert "source" in data["experimental"]
        assert "wavenumbers_cm" in data["experimental"]

    def test_json_experimental_none(self, fake_analysis_results, tmp_path):
        fake_analysis_results["experimental"] = None
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["experimental"] is None

    def test_numpy_arrays_serialized_to_lists(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        freqs = data["comparisons"][0]["spectrum_ml"]["frequencies_cm"]
        assert isinstance(freqs, list)
        assert all(isinstance(f, float) for f in freqs)

    def test_experimental_agreement_in_comparison(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["comparisons"][0]["experimental_agreement"] == pytest.approx(0.92)

    def test_nan_experimental_agreement_becomes_null(self, fake_analysis_results, tmp_path):
        fake_analysis_results["comparisons"][0]["experimental_agreement"] = float("nan")
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        assert data["comparisons"][0]["experimental_agreement"] is None


class TestSummaryMetricsCsv:
    def test_writes_csv(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        csv_path = tmp_path / "summary_metrics.csv"
        assert csv_path.exists()

    def test_csv_has_header_and_rows(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        csv_path = tmp_path / "summary_metrics.csv"
        with csv_path.open() as f:
            reader = csv.reader(f)
            rows = list(reader)
        assert rows[0][0] == "method"
        assert len(rows) == 3  # header + 2 comparisons

    def test_csv_method_names(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        csv_path = tmp_path / "summary_metrics.csv"
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            names = [row["method"] for row in reader]
        assert "mace_off_espaloma" in names
        assert "mace_anicc_mace_ml" in names


class TestRoundTrip:
    def test_json_is_valid(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        # Should not raise
        data = json.loads(out.read_text())
        assert isinstance(data, dict)

    def test_runtime_keys(self, fake_analysis_results, tmp_path):
        out = tmp_path / "report_data.json"
        export_report_data(fake_analysis_results, out)
        data = json.loads(out.read_text())
        rt = data["comparisons"][0]["runtime"]
        assert "ml_pipeline_s" in rt
        assert "dft_pipeline_s" in rt
        assert "ml_gaussian_elapsed_s" in rt
        assert "dft_gaussian_elapsed_s" in rt
        assert "speedup" in rt
