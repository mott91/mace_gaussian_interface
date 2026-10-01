"""Band and chi tables for the supervisor figures (mace_gaussian.analysis.band_table)."""

import json

import pytest

from mace_gaussian.analysis.band_table import build_tables, split_run_name

X_BLOCK = """ Total Anharmonic X Matrix (in cm^-1)
 ------------------------------------
                1             2
      1 {x11}
      2 {x21}  {x22}

"""


def _run(path, fund, over, comb, x):
    """A 2-mode run. fund: {mode: (harm, anharm)} in Gaussian numbering."""
    path.mkdir(parents=True)
    freqs = {
        "anharmonic": [
            {"mode": m, "freq_harmonic": h, "freq_cm": a, "ir_intensity": 10.0 * m}
            for m, (h, a) in fund.items()
        ],
        "overtones": [
            {
                "mode": m,
                "overtone_level": 2,
                "freq_harmonic": 2 * fund[m][0],
                "freq_anharmonic": v,
                "ir_intensity": 1.0,
            }
            for m, v in over.items()
        ],
        "combination_bands": [
            {
                "mode1": 2,
                "mode2": 1,
                "freq_harmonic": fund[1][0] + fund[2][0],
                "freq_anharmonic": comb,
                "ir_intensity": 0.5,
            }
        ],
    }
    (path / "results.json").write_text(json.dumps({"frequencies": freqs}))
    (path / "gaussian_freq.log").write_text(
        X_BLOCK.format(**{k: f"{v:.6E}".replace("E", "D") for k, v in x.items()})
    )


@pytest.fixture
def swapped(tmp_path):
    """DFT: mode 1 = 1000, mode 2 = 2000. ML lists the same motions the other way round
    (its 2000 band is its checkpoint mode 2, but the eigenvector mapping says ML ckpt 0 is
    DFT ckpt 1). Mode 1 of the DFT has a resonance: its overtone is pushed off 2v-2x."""
    base, ana = tmp_path / "comparison_results", tmp_path / "analysis_results"
    _run(
        base / "mol" / "b3lyp_6-31Gdp",
        fund={1: (1050.0, 1000.0), 2: (2100.0, 2000.0)},
        over={1: 1990.0, 2: 3960.0},  # x11 = -5, x22 = -20
        comb=2990.0,  # x12 = -10
        x={"x11": -8.0, "x21": -10.0, "x22": -20.0},  # x11 deperturbed differs -> resonant
    )
    _run(
        base / "mol" / "mace_omol_mace_ml",
        fund={1: (1040.0, 990.0), 2: (2090.0, 1985.0)},
        over={1: 1972.0, 2: 3930.0},  # x11 = -4, x22 = -20
        comb=2965.0,  # x12 = -10
        x={"x11": -4.0, "x21": -10.0, "x22": -20.0},
    )
    (ana / "mol").mkdir(parents=True)
    report = {
        "comparisons": [
            {
                "name": "mace_omol_mace_ml",
                "mode_mapping": {"0": 1, "1": 0},
                "mode_overlaps": {"0": 0.9, "1": 0.8},
            }
        ]
    }
    (ana / "mol" / "report_data.json").write_text(json.dumps(report))
    return build_tables("mol", results_base=base, analysis_base=ana)


def test_split_run_name():
    assert split_run_name("mace_omol_mace_polar1") == ("mace_omol", "mace_polar1")
    assert split_run_name("mace_anicc_espaloma") == ("mace_anicc", "espaloma")


def test_bands_follow_the_mode_mapping(swapped):
    bands, _ = swapped
    by = {(r.band_type, r.i, r.j): r for r in bands}
    assert len(bands) == 5
    # DFT mode 1 (1000) pairs with ML ckpt mode 2 (1985), not ML ckpt mode 1 (990)
    assert by[("fundamental", 1, None)].ml_anharmonic == 1985.0
    assert by[("overtone", 1, None)].ml_anharmonic == 3930.0
    assert by[("combination", 1, 2)].ml_anharmonic == 2965.0
    assert by[("fundamental", 1, None)].overlap == pytest.approx(0.8)
    assert by[("combination", 1, 2)].overlap == pytest.approx(0.8)


def test_chi_from_band_positions_and_resonance_flag(swapped):
    _, chi = swapped
    by = {(r.i, r.j): r for r in chi}
    assert set(by) == {(1, 1), (1, 2), (2, 2)}
    d11 = by[(1, 1)]
    assert d11.kind == "diagonal"
    assert d11.dft_bands == pytest.approx(-5.0)
    assert d11.ml_bands == pytest.approx(-20.0)  # ML mode paired with DFT mode 1
    assert d11.dft_resonant and not d11.ml_resonant
    assert by[(1, 2)].dft_bands == pytest.approx(-10.0)
    assert by[(1, 2)].ml_bands == pytest.approx(-10.0)
    assert not by[(1, 2)].dft_resonant
