"""Tests for mace_gaussian.analysis.executive_summary (Phase 23 Wave 0 RED tests)."""

import math

import numpy as np
import pytest

from mace_gaussian.analysis.executive_summary import (
    build_verdict,
    compute_experimental_agreement,
    rank_methods,
)


# ── compute_experimental_agreement ──────────────────────────────────────────


class TestComputeExperimentalAgreement:
    def test_identical_spectra_returns_one(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        assert compute_experimental_agreement(arr, arr) == pytest.approx(1.0)

    def test_none_ml_returns_nan(self):
        assert math.isnan(compute_experimental_agreement(None, np.array([1.0, 2.0])))

    def test_none_exp_returns_nan(self):
        assert math.isnan(compute_experimental_agreement(np.array([1.0, 2.0]), None))

    def test_zero_exp_returns_nan(self):
        ml = np.array([1.0, 2.0, 3.0])
        exp = np.zeros(3)
        assert math.isnan(compute_experimental_agreement(ml, exp))

    def test_zero_ml_returns_nan(self):
        ml = np.zeros(3)
        exp = np.array([1.0, 2.0, 3.0])
        assert math.isnan(compute_experimental_agreement(ml, exp))

    def test_negative_correlation_clipped_to_zero(self):
        ml = np.array([1.0, 2.0, 3.0, 4.0])
        exp = np.array([4.0, 3.0, 2.0, 1.0])  # perfectly anti-correlated
        result = compute_experimental_agreement(ml, exp)
        assert result == pytest.approx(0.0)

    def test_returns_python_float(self):
        arr = np.array([1.0, 2.0, 3.0])
        result = compute_experimental_agreement(arr, arr)
        assert type(result) is float


# ── rank_methods ────────────────────────────────────────────────────────────


class TestRankMethods:
    def test_empty_comparisons(self):
        assert rank_methods([]) == []

    def test_single_method_scores_zero_rmse_norm(self, fake_metrics):
        comps = [
            {
                "name": "method_a",
                "metrics": fake_metrics,
                "speedup": 5.0,
                "experimental_agreement": 0.9,
            }
        ]
        ranked = rank_methods(comps)
        assert len(ranked) == 1
        assert ranked[0]["name"] == "method_a"
        # Single method → rmse_norm = 0 (min == max)
        assert ranked[0]["composite_score"] >= 0.0

    def test_two_methods_better_first(self, fake_metrics, fake_metrics_b):
        comps = [
            {
                "name": "worse",
                "metrics": fake_metrics_b,
                "speedup": 3.0,
                "experimental_agreement": 0.7,
            },
            {
                "name": "better",
                "metrics": fake_metrics,
                "speedup": 10.0,
                "experimental_agreement": 0.95,
            },
        ]
        ranked = rank_methods(comps)
        assert ranked[0]["name"] == "better"
        assert ranked[1]["name"] == "worse"
        assert ranked[0]["composite_score"] < ranked[1]["composite_score"]

    def test_output_keys(self, fake_metrics):
        comps = [
            {
                "name": "m",
                "metrics": fake_metrics,
                "speedup": 5.0,
                "experimental_agreement": 0.9,
            }
        ]
        ranked = rank_methods(comps)
        expected_keys = {
            "name",
            "composite_score",
            "r2_freq",
            "r2_intensity",
            "rmse_freq",
            "mae_freq",
            "rmse_intensity",
            "mae_intensity",
            "speedup",
            "experimental_agreement",
        }
        assert set(ranked[0].keys()) == expected_keys

    def test_no_experimental_redistributes_weights(self, fake_metrics, fake_metrics_b):
        """When all experimental_agreement are NaN, exp weight goes to rmse."""
        comps = [
            {
                "name": "a",
                "metrics": fake_metrics,
                "speedup": 5.0,
                "experimental_agreement": float("nan"),
            },
            {
                "name": "b",
                "metrics": fake_metrics_b,
                "speedup": 3.0,
                "experimental_agreement": float("nan"),
            },
        ]
        ranked = rank_methods(comps)
        assert ranked[0]["experimental_agreement"] is None
        assert ranked[1]["experimental_agreement"] is None

    def test_deterministic(self, fake_metrics, fake_metrics_b):
        comps = [
            {
                "name": "a",
                "metrics": fake_metrics,
                "speedup": 5.0,
                "experimental_agreement": 0.9,
            },
            {
                "name": "b",
                "metrics": fake_metrics_b,
                "speedup": 3.0,
                "experimental_agreement": 0.7,
            },
        ]
        r1 = rank_methods(comps)
        r2 = rank_methods(comps)
        assert [x["name"] for x in r1] == [x["name"] for x in r2]
        assert [x["composite_score"] for x in r1] == [x["composite_score"] for x in r2]


# ── build_verdict ───────────────────────────────────────────────────────────


class TestBuildVerdict:
    def test_empty_list(self):
        assert build_verdict([], has_experimental=False) == "No ML comparisons available."

    def test_with_experimental(self):
        ranked = [
            {
                "name": "mace_off_espaloma",
                "composite_score": 0.1,
                "r2_freq": 0.995,
                "r2_intensity": 0.85,
                "rmse_freq": 20.0,
                "speedup": 9.6,
                "experimental_agreement": 0.92,
            }
        ]
        verdict = build_verdict(ranked, has_experimental=True)
        assert "mace_off_espaloma" in verdict
        assert "closest to experiment" in verdict
        assert "0.995" in verdict
        assert "20.0" in verdict
        # Agreement value intentionally absent: near-noise metric hidden from display
        assert "0.92" not in verdict

    def test_without_experimental(self):
        ranked = [
            {
                "name": "mace_anicc_mace_ml",
                "composite_score": 0.2,
                "r2_freq": 0.980,
                "r2_intensity": 0.70,
                "rmse_freq": 30.0,
                "speedup": 8.0,
                "experimental_agreement": None,
            }
        ]
        verdict = build_verdict(ranked, has_experimental=False)
        assert "mace_anicc_mace_ml" in verdict
        assert "closest to DFT" in verdict
        assert "speedup" in verdict
