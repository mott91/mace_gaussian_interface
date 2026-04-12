"""Phase 23 executive summary ranking and verdict logic (D-01, D-02).

Computes the 'best method' composite score for per-molecule reports and
emits a one-line verdict for the executive summary card.

Scoring weights (documented in 23-RESEARCH.md A2):
- 40% normalized RMSE (primary accuracy metric)
- 30% (1 - R2_freq)  (linearity of frequency correlation)
- 20% (1 - R2_intensity) (intensity accuracy)
- 10% (1 - experimental_agreement) when experimental data available; else 0 and redistribute to RMSE
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.stats import pearsonr

_DEFAULT_WEIGHTS = {"rmse": 0.4, "r2_freq": 0.3, "r2_int": 0.2, "exp": 0.1}
_WEIGHTS_NO_EXP = {"rmse": 0.5, "r2_freq": 0.3, "r2_int": 0.2, "exp": 0.0}


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
        weights = _DEFAULT_WEIGHTS if has_exp else _WEIGHTS_NO_EXP

    rmses = np.array([c["metrics"].rmse_freq for c in comparisons], dtype=float)
    rmse_min = float(rmses.min())
    rmse_max = float(rmses.max())
    span = max(rmse_max - rmse_min, 1e-9)
    rmse_norm = (rmses - rmse_min) / span

    scored: list[dict[str, Any]] = []
    for c, rn in zip(comparisons, rmse_norm):
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
            weights["rmse"] * float(rn)
            + weights["r2_freq"] * (1.0 - float(m.r2_freq))
            + weights["r2_int"] * (1.0 - float(m.r2_intensity))
            + weights["exp"] * (1.0 - exp_agree_val)
        )
        scored.append(
            {
                "name": c["name"],
                "composite_score": float(score),
                "r2_freq": float(m.r2_freq),
                "r2_intensity": float(m.r2_intensity),
                "rmse_freq": float(m.rmse_freq),
                "speedup": float(c.get("speedup", 0.0)),
                "experimental_agreement": exp_agree_out if has_exp else None,
            }
        )
    scored.sort(key=lambda x: x["composite_score"])
    return scored


def build_verdict(ranked: list[dict[str, Any]], has_experimental: bool) -> str:
    """Human-readable one-liner for the executive summary card (D-01)."""
    if not ranked:
        return "No ML comparisons available."
    best = ranked[0]
    if has_experimental and best.get("experimental_agreement") is not None:
        return (
            f"{best['name']} is closest to experiment "
            f"(R\u00b2_freq={best['r2_freq']:.3f}, RMSE={best['rmse_freq']:.1f} cm\u207b\u00b9, "
            f"experimental agreement={best['experimental_agreement']:.2f})."
        )
    return (
        f"{best['name']} is closest to DFT "
        f"(R\u00b2_freq={best['r2_freq']:.3f}, RMSE={best['rmse_freq']:.1f} cm\u207b\u00b9, "
        f"{best['speedup']:.1f}\u00d7 speedup)."
    )
