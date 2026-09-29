"""Smoke tests for the per-molecule report (report_v2.py, the default since 2026-09-18)."""

from __future__ import annotations

from mace_gaussian.analysis.report_v2 import ReportV2Generator


def test_v2_report_builds_from_fixture(fake_analysis_results, tmp_path):
    out = ReportV2Generator("water", tmp_path, mode="anharmonic").generate(fake_analysis_results)
    html = out.read_text(encoding="utf-8")
    assert out.name == "report.html"
    assert (tmp_path / "report_data.json").exists()
    assert (tmp_path / "master_table.tex").exists()
    # Two runs of two different energy models -> two model sections.
    assert html.count('class="comparison-section model-section"') == 2
    for anchor in (
        "executive-summary",
        "combined",
        "per-mode-errors",
        "master-table",
        "model-1",
        "summary-table",
    ):
        assert f'id="{anchor}"' in html
    assert html.count("cdn.plot.ly/plotly") == 1
    assert "OFF" in html and "ANI-cc" in html


def test_v2_report_harmonic_mode_has_no_overtones(fake_analysis_results, tmp_path):
    fake_analysis_results["mode"] = "harmonic"
    out = ReportV2Generator("water", tmp_path, mode="harmonic").generate(fake_analysis_results)
    html = out.read_text(encoding="utf-8")
    assert 'id="overtones"' not in html
    assert "Harmonic: ML" in html


def test_v2_report_empty(tmp_path):
    out = ReportV2Generator("x", tmp_path).generate({"molecule": "x", "comparisons": []})
    assert "No comparisons" in out.read_text()
