"""The NIST JDX parser must return wavenumbers whatever the file's x unit is."""

from __future__ import annotations

import numpy as np

from mace_gaussian.analysis import nist_fetcher


def _install_fake_jcamp(monkeypatch, payload):
    """Point nist_fetcher's ``import jcamp`` at a stub returning ``payload``."""
    import sys

    stub = type("J", (), {"jcamp_readfile": staticmethod(lambda _path: payload)})
    monkeypatch.setitem(sys.modules, "jcamp", stub)


def test_micrometer_axis_is_converted_to_wavenumbers(monkeypatch, tmp_path):
    payload = {
        "x": [2.5, 5.0, 10.0],
        "y": [0.9, 0.5, 0.95],  # fractional transmittance
        "xunits": "MICROMETERS",
        "yunits": "TRANSMITTANCE",
        "title": "methanol",
    }
    _install_fake_jcamp(monkeypatch, payload)
    jdx = tmp_path / "x.jdx"
    jdx.write_text("dummy")
    spec = nist_fetcher._parse_jdx_file(jdx, "methanol")
    assert spec is not None
    np.testing.assert_allclose(spec.wavenumbers, [1000.0, 2000.0, 4000.0])
    # sorted ascending in wavenumber, strongest absorption (T=0.5 at 5 um) at 2000 cm-1
    assert spec.absorbance[1] == 1.0


def test_wavenumber_axis_is_left_alone(monkeypatch, tmp_path):
    payload = {
        "x": [4000.0, 3000.0, 1000.0],
        "y": [0.1, 0.9, 0.4],
        "xunits": "1/CM",
        "yunits": "ABSORBANCE",
    }
    _install_fake_jcamp(monkeypatch, payload)
    jdx = tmp_path / "x.jdx"
    jdx.write_text("dummy")
    spec = nist_fetcher._parse_jdx_file(jdx, "water")
    np.testing.assert_allclose(spec.wavenumbers, [1000.0, 3000.0, 4000.0])
    np.testing.assert_allclose(spec.absorbance, [0.4 / 0.9, 1.0, 0.1 / 0.9])
