"""The batch run refreshes the thesis figures through a subprocess that must never break it."""

from __future__ import annotations

import subprocess
from types import SimpleNamespace

from mace_gaussian import batch


def test_refresh_passes_directories_and_reports_success(monkeypatch, tmp_path):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return SimpleNamespace(returncode=0, stdout="gallery: x (23 figures)\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert batch.refresh_thesis_figures(str(tmp_path / "a"), str(tmp_path / "c")) is True
    cmd = calls[0]
    assert cmd[1].endswith("make_thesis_figures.py")
    assert cmd[cmd.index("--analysis-dir") + 1] == str((tmp_path / "a").resolve())
    assert cmd[cmd.index("--comparison-dir") + 1] == str((tmp_path / "c").resolve())


def test_refresh_failure_is_reported_not_raised(monkeypatch, tmp_path):
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(returncode=1, stdout="", stderr="LaTeX error"),
    )
    assert batch.refresh_thesis_figures(str(tmp_path), str(tmp_path)) is False


def test_refresh_skips_when_script_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(batch, "FIGURE_SCRIPT", tmp_path / "missing.py")
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))
    assert batch.refresh_thesis_figures(str(tmp_path), str(tmp_path)) is False
    assert called == []
