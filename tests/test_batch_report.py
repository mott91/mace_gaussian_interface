"""Tests for batch report generation."""

import json
from pathlib import Path

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


# ---- Phase 23 T-23-01 batch-path XSS mitigation RED test ----
# This test is RED until Plan 05 Task 1 wraps batch_report.py's f-string
# interpolations of molecule/combo/hardware strings with html.escape.
# It builds a minimal fake comparison_results/ tree with a malicious
# molecule directory name, runs generate_batch_report, and asserts the
# raw <script> tag is NOT present in the rendered HTML.

def _write_min_results_json(path: Path, freqs: list[float]) -> None:
    """Write a minimal Phase 21 results.json schema with harmonic freqs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "frequencies": {
                    "harmonic": [{"freq_cm": f, "ir_intensity": 1.0} for f in freqs],
                },
                "runtime_s": 1.0,
                "gaussian_timing": {"total_elapsed_s": 1.0},
                "version_info": {"cpu_model": "TestCPU", "gpu_name": ""},
                "hardware": {"cpu_model": "TestCPU", "node": "testnode"},
            }
        )
    )


def test_batch_report_escapes_molecule_name(tmp_path):
    """T-23-01 mitigation: batch report must HTML-escape malicious molecule names.

    Payload uses `<img onerror=...>` rather than `<script>...</script>` because
    the forward slash in a closing `</script>` tag gets interpreted as a path
    separator when used as a directory name, fragmenting the tree.
    """
    malicious = "<img src=x onerror=alert(1)>"
    results_dir = tmp_path / "comparison_results"
    mol_dir = results_dir / malicious
    freqs = [1600.0, 3700.0, 3800.0]

    _write_min_results_json(mol_dir / "b3lyp_6-31Gdp" / "results.json", freqs)
    _write_min_results_json(mol_dir / "mace_off_espaloma" / "results.json", freqs)

    out_dir = tmp_path / "batch_report_out"
    report_file = generate_batch_report(
        results_dir=str(results_dir),
        output_dir=str(out_dir),
    )

    html = Path(report_file).read_text()

    assert malicious not in html, (
        "batch_report did not escape malicious molecule name -- "
        "T-23-01 mitigation missing. Expected html.escape wrapping in "
        "mace_gaussian/analysis/batch_report.py."
    )
    assert "&lt;img src=x onerror=alert(1)&gt;" in html, (
        "Expected HTML-escaped form of payload in batch report output"
    )
