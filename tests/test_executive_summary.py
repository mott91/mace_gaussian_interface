"""RED tests for Phase 23 executive_summary module."""

import math

import numpy as np
import pytest

executive_summary = pytest.importorskip("mace_gaussian.analysis.executive_summary")


def test_compute_experimental_agreement_identical_spectra():
    x = np.linspace(400, 4000, 100)
    spec = np.exp(-((x - 2000) ** 2) / 1000)
    r = executive_summary.compute_experimental_agreement(spec, spec)
    assert r == pytest.approx(1.0, abs=1e-6)


def test_compute_experimental_agreement_none_or_zero_returns_nan():
    r = executive_summary.compute_experimental_agreement(np.zeros(100), np.zeros(100))
    assert math.isnan(r)


def test_rank_methods_deterministic(fake_metrics):
    comparisons = [
        {
            "name": "method_a",
            "metrics": fake_metrics,
            "speedup": 5.0,
            "experimental_agreement": 0.9,
        },
        {
            "name": "method_b",
            "metrics": fake_metrics,
            "speedup": 3.0,
            "experimental_agreement": 0.7,
        },
    ]
    ranked1 = executive_summary.rank_methods(comparisons)
    ranked2 = executive_summary.rank_methods(comparisons)
    assert [r["name"] for r in ranked1] == [r["name"] for r in ranked2]
    assert len(ranked1) == 2
    assert "composite_score" in ranked1[0]


def test_build_verdict_with_experimental(fake_metrics):
    ranked = [
        {
            "name": "mace_off_espaloma",
            "composite_score": 0.1,
            "r2_freq": 0.99,
            "rmse_freq": 15.0,
            "r2_intensity": 0.9,
            "speedup": 4.5,
            "experimental_agreement": 0.88,
        }
    ]
    verdict = executive_summary.build_verdict(ranked, has_experimental=True)
    assert "mace_off_espaloma" in verdict
    assert "experiment" in verdict.lower()


def test_build_verdict_no_experimental():
    ranked = [
        {
            "name": "mace_off_espaloma",
            "composite_score": 0.1,
            "r2_freq": 0.99,
            "rmse_freq": 15.0,
            "r2_intensity": 0.9,
            "speedup": 4.5,
            "experimental_agreement": None,
        }
    ]
    verdict = executive_summary.build_verdict(ranked, has_experimental=False)
    assert "mace_off_espaloma" in verdict
    assert "DFT" in verdict or "speedup" in verdict.lower()
