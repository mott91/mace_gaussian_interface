"""RED tests for Phase 23 report_data module."""

import json

import pytest

report_data = pytest.importorskip("mace_gaussian.analysis.report_data")


def test_export_report_data_writes_json_file(fake_analysis_results, tmp_path):
    out = tmp_path / "report_data.json"
    report_data.export_report_data(fake_analysis_results, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_export_report_data_roundtrip(fake_analysis_results, tmp_path):
    out = tmp_path / "report_data.json"
    report_data.export_report_data(fake_analysis_results, out)
    with out.open() as f:
        payload = json.load(f)
    assert payload["schema_version"] == 1
    assert payload["molecule"] == "water"
    assert payload["mode"] == "anharmonic"
    assert "comparisons" in payload
    assert len(payload["comparisons"]) == 1
    comp = payload["comparisons"][0]
    assert comp["name"] == "mace_off_espaloma"
    assert "metrics" in comp
    assert comp["metrics"]["r2_freq"] == pytest.approx(0.987)


def test_export_report_data_handles_numpy_arrays(fake_analysis_results, tmp_path):
    out = tmp_path / "report_data.json"
    report_data.export_report_data(fake_analysis_results, out)
    with out.open() as f:
        payload = json.load(f)
    comp = payload["comparisons"][0]
    assert isinstance(comp["spectrum_ml"]["frequencies_cm"], list)
    assert isinstance(comp["spectrum_ml"]["intensities"], list)


def test_export_report_data_without_experimental(fake_analysis_results, tmp_path):
    fake_analysis_results["experimental"] = None
    out = tmp_path / "report_data.json"
    report_data.export_report_data(fake_analysis_results, out)
    with out.open() as f:
        payload = json.load(f)
    assert payload["experimental"] is None
