"""Campaign isolation: one name decides every folder, legacy layout unchanged."""

from pathlib import Path
from unittest import mock

import pytest

from mace_gaussian import batch
from mace_gaussian.campaign import DEFAULT_REMOTE_BASE, campaign_paths


def test_legacy_layout_unchanged():
    p = campaign_paths(None)
    assert p.comparison == Path("comparison_results")
    assert p.analysis == Path("analysis_results")
    assert p.harmonic == Path("analysis_results_harmonic")
    assert p.figures.parts[-2:] == ("thesis", "figures")
    assert p.remote_base == DEFAULT_REMOTE_BASE


def test_campaign_owns_every_folder():
    p = campaign_paths("2026")
    for folder in (p.comparison, p.analysis, p.harmonic):
        assert folder.parts[:2] == ("campaigns", "2026")
    assert "campaigns" in p.figures.parts and "2026" in p.figures.parts
    assert p.remote_base == f"{DEFAULT_REMOTE_BASE}/campaigns/2026"
    # the harmonic analysis appends "_harmonic" to the analysis folder
    assert str(p.analysis) + "_harmonic" == str(p.harmonic)


@pytest.mark.parametrize("bad", ["", "../x", "a/b", " 2026", "-x"])
def test_campaign_name_validated(bad):
    with pytest.raises(ValueError):
        campaign_paths(bad)


def test_template_resources():
    assert batch.template_resources("#SBATCH --cpus-per-task=8\n#SBATCH --mem=16G\n") == (8, "15GB")
    assert batch.template_resources("#SBATCH --cpus-per-task=4\n#SBATCH --mem=4G\n") == (4, "4GB")
    assert batch.template_resources("no sbatch lines") == (4, "4GB")


def test_analyses_go_to_the_given_folder():
    with (
        mock.patch("mace_gaussian.analysis.analyze_molecule") as anh,
        mock.patch("mace_gaussian.analysis.analyze_molecule_harmonic") as har,
    ):
        batch._run_analyses(
            "h2o2", "campaigns/x/comparison_results", {}, "campaigns/x/analysis_results"
        )
    assert anh.call_args.kwargs == {
        "base_results_dir": "campaigns/x/comparison_results",
        "output_dir": "campaigns/x/analysis_results",
    }
    assert har.call_args.kwargs["output_dir"] == "campaigns/x/analysis_results"
    assert har.call_args.kwargs["base_results_dir"] == "campaigns/x/comparison_results"


def test_cli_rejects_campaign_with_output_dir(tmp_path):
    from click.testing import CliRunner

    from mace_gaussian.cli import cli

    listfile = tmp_path / "b.txt"
    listfile.write_text("")
    res = CliRunner().invoke(
        cli, ["batch", str(listfile), "--campaign", "x", "--output-dir", "elsewhere"]
    )
    assert res.exit_code == 2
    assert "either --campaign or --output-dir" in res.output
