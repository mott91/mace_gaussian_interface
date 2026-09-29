"""Tests for the per-mode master table (analysis/master_table.py)."""

# ruff: noqa: RUF001  (Greek nu in band labels is intentional)

from __future__ import annotations

import csv
import json

import numpy as np
import pytest

from mace_gaussian.analysis.analyze_spectra import SpectrumData
from mace_gaussian.analysis.master_table import (
    assign_band_origins,
    build_master_table,
    display_name,
    frequency_methods,
    render_master_table_html,
    write_master_table_csv,
    write_master_table_latex,
)


def _spectrum(freqs, ints, ids=None):
    n = len(freqs)
    return SpectrumData(
        frequencies=np.array(freqs, dtype=float),
        intensities=np.array(ints, dtype=float),
        labels=["fundamental"] * n,
        mode_ids=ids or [f"F{i + 1}" for i in range(n)],
    )


@pytest.fixture
def water_like_results():
    """Two ML methods; the second has its modes 2 and 3 swapped relative to DFT.

    DFT checkpoint order (ascending harmonic): bend 1665, sym str 3799, asym str 3912.
    ML method B's checkpoint order puts the asym stretch before the sym stretch,
    so the eigenvector mapping is {0: 0, 1: 2, 2: 1}.
    """
    dft_results = {
        "frequencies": {
            "harmonic": [{"freq_cm": 1665.3}, {"freq_cm": 3799.2}, {"freq_cm": 3912.4}],
        }
    }
    dft_spec = _spectrum([1615.1, 3624.5, 3722.7], [69.3, 0.6, 16.1])
    ml_a_results = {
        "frequencies": {
            "harmonic": [{"freq_cm": 1621.9}, {"freq_cm": 3818.1}, {"freq_cm": 3929.0}],
        }
    }
    ml_a_spec = _spectrum([1580.0, 3643.8, 3740.0], [80.5, 2.1, 20.0])
    ml_b_results = {
        "frequencies": {
            "harmonic": [{"freq_cm": 1600.0}, {"freq_cm": 3950.0}, {"freq_cm": 3800.0}],
        }
    }
    ml_b_spec = _spectrum([1560.0, 3760.0, 3620.0], [70.0, 20.0, 1.0])
    return {
        "molecule": "water",
        "mode": "anharmonic",
        "comparisons": [
            {
                "name": "mace_omol_mace_ml",
                "dft_spectrum": dft_spec,
                "ml_spectrum": ml_a_spec,
                "mode_mapping": {0: 0, 1: 1, 2: 2},
                "mode_overlaps": {0: 0.99, 1: 0.98, 2: 0.97},
                "_dft_results": dft_results,
                "_ml_results": ml_a_results,
            },
            {
                "name": "mace_off_espaloma",
                "dft_spectrum": dft_spec,
                "ml_spectrum": ml_b_spec,
                "mode_mapping": {0: 0, 1: 2, 2: 1},
                "mode_overlaps": {0: 0.95, 1: 0.60, 2: 0.94},
                "_dft_results": dft_results,
                "_ml_results": ml_b_results,
            },
        ],
    }


class TestBuild:
    def test_one_row_per_dft_mode_with_dft_columns(self, water_like_results):
        t = build_master_table(water_like_results)
        assert t["mode"] == "anharmonic"
        assert [r["dft_mode"] for r in t["rows"]] == [1, 2, 3]
        assert [r["dft_harmonic"] for r in t["rows"]] == [1665.3, 3799.2, 3912.4]
        assert [r["dft_vpt2"] for r in t["rows"]] == [1615.1, 3624.5, 3722.7]
        assert [m["name"] for m in t["methods"]] == ["mace_omol_mace_ml", "mace_off_espaloma"]
        assert all(m["pairing"] == "eigenvector" for m in t["methods"])

    def test_ml_columns_follow_the_eigenvector_mapping(self, water_like_results):
        t = build_master_table(water_like_results)
        row_sym = t["rows"][1]  # DFT mode 2 = symmetric stretch
        a = row_sym["methods"]["mace_omol_mace_ml"]
        b = row_sym["methods"]["mace_off_espaloma"]
        assert a["ml_mode"] == 2 and a["harmonic"] == 3818.1 and a["vpt2"] == 3643.8
        # Method B stores the symmetric stretch as its mode 3.
        assert b["ml_mode"] == 3 and b["harmonic"] == 3800.0 and b["vpt2"] == 3620.0
        assert b["overlap"] == pytest.approx(0.94)
        row_asym = t["rows"][2]
        assert row_asym["methods"]["mace_off_espaloma"]["ml_mode"] == 2
        assert row_asym["methods"]["mace_off_espaloma"]["overlap"] == pytest.approx(0.60)

    def test_experimental_band_origins_assigned_by_nearest_vpt2(self, water_like_results):
        t = build_master_table(water_like_results)  # uses the shipped band_origins.json
        labels = [(r["experimental"] or {}).get("label") for r in t["rows"]]
        assert labels == ["ν2", "ν1", "ν3"]
        assert t["rows"][0]["experimental"]["freq_cm"] == 1595
        assert t["rows"][0]["experimental"]["assignment"] == "nearest"
        assert t["experimental"]["n_assigned"] == 3

    def test_index_pairing_without_mapping_and_without_raw_results(self, fake_analysis_results):
        """The report fixture has no mode_mapping and no _ml_results: must not crash."""
        t = build_master_table(fake_analysis_results)
        assert len(t["rows"]) == 3
        assert all(m["pairing"] == "index" for m in t["methods"])
        cell = t["rows"][0]["methods"]["mace_off_espaloma"]
        assert cell["vpt2"] == 1600.0 and cell["harmonic"] is None

    def test_harmonic_mode_has_no_vpt2_columns(self, water_like_results):
        water_like_results["mode"] = "harmonic"
        t = build_master_table(water_like_results)
        assert all(r["dft_vpt2"] is None for r in t["rows"])
        # In harmonic mode the spectrum *is* the harmonic spectrum.
        assert t["rows"][0]["dft_harmonic"] == 1615.1
        assert t["rows"][0]["methods"]["mace_omol_mace_ml"]["harmonic"] == 1580.0

    def test_empty_comparisons(self):
        t = build_master_table({"molecule": "x", "mode": "anharmonic", "comparisons": []})
        assert t["rows"] == [] and t["methods"] == []


class TestAssignment:
    def test_degenerate_band_takes_several_modes(self):
        bands = [
            {"label": "ν4", "freq_cm": 1306, "degeneracy": 3},
            {"label": "ν2", "freq_cm": 1534, "degeneracy": 2},
            {"label": "ν1", "freq_cm": 2917, "degeneracy": 1},
            {"label": "ν3", "freq_cm": 3019, "degeneracy": 3},
        ]
        ref = {1: 1327, 2: 1327, 3: 1327, 4: 1545, 5: 1545, 6: 2927, 7: 3025, 8: 3025, 9: 3025}
        out = assign_band_origins(bands, ref)
        assert [out[k]["label"] for k in range(1, 10)] == [
            "ν4",
            "ν4",
            "ν4",
            "ν2",
            "ν2",
            "ν1",
            "ν3",
            "ν3",
            "ν3",
        ]

    def test_explicit_modes_override_and_window_limits(self):
        bands = [
            {"label": "a", "freq_cm": 1000, "modes": [2]},
            {"label": "b", "freq_cm": 5000},  # nothing within the window
        ]
        out = assign_band_origins(bands, {1: 1000, 2: 3000})
        assert out == {2: {**bands[0], "assignment": "explicit"}}


class TestOutputs:
    def test_html_section_contains_values_and_anchor(self, water_like_results):
        t = build_master_table(water_like_results)
        h = render_master_table_html(t)
        assert 'id="master-table"' in h
        assert "3818.1" in h and "3624.5" in h and "1595" in h
        assert "metric-warning" in h  # the 0.60 overlap cell is flagged
        assert "+19.3" in h  # 3643.8 - 3624.5

    def test_csv_and_latex_files(self, water_like_results, tmp_path):
        t = build_master_table(water_like_results)
        write_master_table_csv(t, tmp_path / "m.csv")
        with (tmp_path / "m.csv").open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 3
        assert rows[1]["mace_off_espaloma_ml_mode"] == "3"
        assert rows[1]["mace_off_espaloma_vpt2"] == "3620.0000"
        assert rows[0]["exp_label"] == "ν2"

        write_master_table_latex(t, tmp_path / "m.tex")
        tex = (tmp_path / "m.tex").read_text()
        assert "\\begin{tabular}" in tex and "\\toprule" in tex
        assert "OMOL" in tex and "{OFF}" in tex  # collapsed to one column per energy model
        assert "$\\nu_{2}$" in tex
        assert "$^\\dagger$" in tex  # low-overlap marker
        assert "mace_off" not in tex  # raw names never leak unescaped

    def test_json_serializable(self, water_like_results):
        json.dumps(build_master_table(water_like_results))


def test_frequency_methods_collapse_to_one_per_energy_model():
    methods = [
        {"name": "mace_off_espaloma", "label": "x", "pairing": "eigenvector"},
        {"name": "mace_off_mace_ml", "label": "x", "pairing": "eigenvector"},
        {"name": "mace_omol_mace_polar1", "label": "x", "pairing": "eigenvector"},
        {"name": "custom", "label": "custom", "pairing": "index"},
    ]
    out = frequency_methods(methods)
    assert [(m["name"], m["label"]) for m in out] == [
        ("mace_off_mace_ml", "OFF"),
        ("mace_omol_mace_polar1", "OMOL"),
        ("custom", "custom"),
    ]


def test_display_name():
    assert display_name("mace_omol_mace_ml") == "OMOL"
    assert display_name("mace_off_espaloma") == "OFF/esp"
    assert display_name("mace_anicc_mace_polar1") == "ANI-cc/P1"
    assert display_name("something_else") == "something_else"
