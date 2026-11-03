"""The ML (external) frequency route: campaign default, and the step-size override."""

from ase.build import molecule

from mace_gaussian.gaussian import io
from mace_gaussian.gaussian.io import ase_to_gjf


def _route_lines(tmp_path):
    ase_to_gjf(molecule("NH3"), "ml.gjf", output_dir=tmp_path)
    return [ln for ln in (tmp_path / "ml.gjf").read_text().splitlines() if ln.startswith("#")]


def test_default_route_is_the_campaign_route(tmp_path):
    lines = _route_lines(tmp_path)
    assert lines[0] == "# freq (anharm)"
    assert lines[1].startswith("# external=")


def test_freq_keyword_can_be_replaced(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "ML_FREQ_KEYWORD", "freq=(anharm,step=200)")
    assert _route_lines(tmp_path)[0] == "# freq=(anharm,step=200)"
