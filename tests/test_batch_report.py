"""Tests for batch report generation."""

from click.testing import CliRunner

from mace_gaussian.analysis.batch_report import (
    aggregate_results,
    generate_batch_report,
)
from mace_gaussian.cli import cli


def test_aggregate_results_with_real_data():
    """Aggregate from the actual comparison_results/ directory."""
    df = aggregate_results("comparison_results")
    assert not df.empty, "Expected non-empty DataFrame from comparison_results"
    assert set(df.columns) >= {"molecule", "combo", "r2", "rmse"}
    assert "water" in df["molecule"].values
    assert (df["r2"] >= -1).all() and (df["r2"] <= 1).all()
    assert (df["rmse"] >= 0).all()


def test_generate_batch_report_creates_html(tmp_path):
    """Generate report from real data into a temp directory."""
    output_dir = str(tmp_path / "report_out")
    generate_batch_report(
        results_dir="comparison_results",
        output_dir=output_dir,
    )
    html_file = tmp_path / "report_out" / "batch_report.html"
    assert html_file.exists()
    content = html_file.read_text()
    assert "Leaderboard" in content
    assert "data:image/png;base64," in content
    assert len(content) > 1000


def test_report_cli_help():
    """CLI report command shows help with expected options."""
    runner = CliRunner()
    result = runner.invoke(cli, ["report", "--help"])
    assert result.exit_code == 0
    assert "--results-dir" in result.output
    assert "--output-dir" in result.output


def test_aggregate_results_empty_dir(tmp_path):
    """Aggregating an empty directory returns empty DataFrame."""
    df = aggregate_results(str(tmp_path))
    assert df.empty
