"""Phase 23 executive summary ranking and verdict logic (D-01, D-02).

Computes the 'best method' composite score for per-molecule reports and
emits a one-line verdict for the executive summary card.

Scoring weights — symmetric across frequency and intensity so dipole-calculator
choice is reflected as strongly as energy-calculator choice, all against the DFT
baseline:
- 25% normalized RMSE_freq
- 25% (1 - R²_freq)
- 25% normalized RMSE_intensity
- 25% (1 - R²_intensity)

experimental_agreement is still computed and exported but carries zero weight: the
ranking measures how well each ML model reproduces DFT, not experiment.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.stats import pearsonr

_DEFAULT_WEIGHTS = {
    "rmse_freq": 0.25,
    "r2_freq": 0.25,
    "rmse_int": 0.25,
    "r2_int": 0.25,
    "exp": 0.0,
}


def _min_max_normalize(arr: np.ndarray) -> np.ndarray:
    """Min-max normalize to [0, 1]; returns zeros when all values are equal."""
    lo = float(arr.min())
    hi = float(arr.max())
    span = max(hi - lo, 1e-9)
    return (arr - lo) / span


def compute_experimental_agreement(
    ml_broadened: np.ndarray | None,
    exp_on_grid: np.ndarray | None,
) -> float:
    """Pearson correlation of broadened ML spectrum vs experimental on common grid.

    Returns value in [0, 1]. Negative correlations clipped to 0. NaN when either
    input is None or the experimental is all zero.
    """
    if ml_broadened is None or exp_on_grid is None:
        return float("nan")
    if np.allclose(exp_on_grid, 0) or np.allclose(ml_broadened, 0):
        return float("nan")
    r, _ = pearsonr(ml_broadened, exp_on_grid)
    if math.isnan(r):
        return float("nan")
    return float(max(0.0, r))


def rank_methods(
    comparisons: list[dict[str, Any]],
    weights: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    """Return comparisons sorted by composite score (lower is better).

    Each output dict carries name, composite_score, r2_freq, r2_intensity,
    rmse_freq, speedup, experimental_agreement (None when unavailable).
    """
    if not comparisons:
        return []

    # Detect experimental availability across the pool
    has_exp = any(
        not (
            c.get("experimental_agreement") is None
            or (
                isinstance(c.get("experimental_agreement"), float)
                and math.isnan(c["experimental_agreement"])
            )
        )
        for c in comparisons
    )
    if weights is None:
        weights = _DEFAULT_WEIGHTS

    rmse_freq_arr = np.array([c["metrics"].rmse_freq for c in comparisons], dtype=float)
    rmse_freq_norm = _min_max_normalize(rmse_freq_arr)
    rmse_int_arr = np.array([c["metrics"].rmse_intensity for c in comparisons], dtype=float)
    rmse_int_norm = _min_max_normalize(rmse_int_arr)

    scored: list[dict[str, Any]] = []
    for c, rn_freq, rn_int in zip(comparisons, rmse_freq_norm, rmse_int_norm):
        m = c["metrics"]
        exp_agree_raw = c.get("experimental_agreement")
        if exp_agree_raw is None or (
            isinstance(exp_agree_raw, float) and math.isnan(exp_agree_raw)
        ):
            exp_agree_val = 0.0
            exp_agree_out: float | None = None
        else:
            exp_agree_val = float(exp_agree_raw)
            exp_agree_out = exp_agree_val
        score = (
            weights["rmse_freq"] * float(rn_freq)
            + weights["r2_freq"] * (1.0 - float(m.r2_freq))
            + weights["rmse_int"] * float(rn_int)
            + weights["r2_int"] * (1.0 - float(m.r2_intensity))
            + weights["exp"] * (1.0 - exp_agree_val)
        )
        n_freq = int(getattr(m, "num_peaks", 0))
        n_int = max(0, n_freq - int(getattr(m, "num_intensity_filtered", 0)))
        scored.append(
            {
                "name": c["name"],
                "composite_score": float(score),
                "r2_freq": float(m.r2_freq),
                "r2_intensity": float(m.r2_intensity),
                "rmse_freq": float(m.rmse_freq),
                "mae_freq": float(m.mae_freq),
                "slope_freq": float(getattr(m, "slope_freq", float("nan"))),
                "rmse_intensity": float(m.rmse_intensity),
                "mae_intensity": float(m.mae_intensity),
                "n_freq": n_freq,
                "n_int": n_int,
                "speedup": float(c.get("speedup", 0.0)),
                "experimental_agreement": exp_agree_out if has_exp else None,
            }
        )
    scored.sort(key=lambda x: x["composite_score"])
    return scored


def build_verdict(ranked: list[dict[str, Any]], has_experimental: bool) -> str:
    """Human-readable one-liner for the executive summary card (D-01).

    Every number in the ranking is measured against the DFT reference, so the verdict
    says so regardless of whether an experimental spectrum was available. (The old
    wording "closest to experiment" quoted DFT metrics while pointing at a hidden,
    near-noise spectrum-correlation score; report review item 1, 2026-09-18.)
    The metric R\u00b2 is not quoted: on a handful of modes it is always ~1 and says nothing.
    """
    if not ranked:
        return "No ML comparisons available."
    best = ranked[0]
    slope = best.get("slope_freq")
    slope_txt = f", slope {slope:.3f}" if slope is not None and math.isfinite(slope) else ""
    return (
        f"{best['name']} is closest to DFT "
        f"(MAE {best.get('mae_freq', float('nan')):.1f} cm\u207b\u00b9, "
        f"RMSE {best['rmse_freq']:.1f} cm\u207b\u00b9{slope_txt}, "
        f"{best['speedup']:.1f}\u00d7 speedup)."
    )
