"""Phase 23 structured data export (D-07).

Writes report_data.json + summary_metrics.csv (+ optional spectrum_grid.csv)
alongside each report.html for downstream thesis figure generation.

Downstream consumers read these files -- never the HTML -- so schema stability
is contractual. Schema version is tracked in ``schema_version``.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA_VERSION = 1


def export_report_data(analysis_results: dict[str, Any], output_path: Path) -> None:
    """Write report_data.json and companion CSV files.

    Parameters
    ----------
    analysis_results:
        The full analysis results dict from ``analysis_workflow``.
    output_path:
        Path for the JSON file; CSVs land next to it.
    """
    output_path = Path(output_path)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "molecule": analysis_results.get("molecule", ""),
        "mode": analysis_results.get("mode", "anharmonic"),
        "bandwidth_fwhm_cm": float(analysis_results.get("bandwidth_fwhm", 10.0)),
        "experimental": _serialize_experimental(analysis_results.get("experimental")),
        "comparisons": [_serialize_comparison(c) for c in analysis_results.get("comparisons", [])],
        "executive_summary": analysis_results.get("executive_summary") or {},
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        json.dump(payload, f, indent=2, default=_json_default)

    _write_summary_metrics_csv(payload, output_path.parent / "summary_metrics.csv")

    freq_grid = analysis_results.get("freq_grid")
    if freq_grid is not None:
        _write_spectrum_grid_csv(
            freq_grid,
            analysis_results,
            output_path.parent / "spectrum_grid.csv",
        )


# ---------------------------------------------------------------------------
# Internal serializers
# ---------------------------------------------------------------------------


def _serialize_comparison(c: dict[str, Any]) -> dict[str, Any]:
    m = c["metrics"]
    exp_agree = c.get("experimental_agreement")
    if isinstance(exp_agree, float) and math.isnan(exp_agree):
        exp_agree = None
    return {
        "name": c["name"],
        "metrics": {
            "r2_freq": float(m.r2_freq),
            "r2_intensity": float(m.r2_intensity),
            "rmse_freq": float(m.rmse_freq),
            "mae_freq": float(m.mae_freq),
            "rmse_intensity": float(m.rmse_intensity),
            "mae_intensity": float(m.mae_intensity),
            "max_error_freq": float(m.max_error_freq),
            "num_matched": int(m.num_matched),
            "num_dft_only": int(m.num_dft_only),
            "num_ml_only": int(m.num_ml_only),
            "match_rate": float(m.match_rate),
        },
        "runtime": {
            "ml_pipeline_s": float(c.get("ml_runtime", 0.0)),
            "dft_pipeline_s": float(c.get("dft_runtime", 0.0)),
            "ml_gaussian_elapsed_s": float(
                c.get("ml_gaussian_timing", {}).get("total_elapsed_s", 0.0)
            ),
            "dft_gaussian_elapsed_s": float(
                c.get("dft_gaussian_timing", {}).get("total_elapsed_s", 0.0)
            ),
            "speedup": float(c.get("speedup", 0.0)),
        },
        "hardware": {
            "ml": c.get("ml_hardware", {}),
            "dft": c.get("dft_hardware", {}),
        },
        "spectrum_ml": _serialize_spectrum(c["ml_spectrum"]),
        "spectrum_dft": _serialize_spectrum(c["dft_spectrum"]),
        "experimental_agreement": exp_agree,
    }


def _serialize_spectrum(spec: Any) -> dict[str, Any]:
    return {
        "frequencies_cm": _to_list(spec.frequencies),
        "intensities": _to_list(spec.intensities),
        "labels": list(spec.labels),
        "mode_ids": list(spec.mode_ids),
    }


def _serialize_experimental(exp: Any) -> dict[str, Any] | None:
    if exp is None:
        return None
    return {
        "source": exp.source,
        "molecule_name": exp.molecule_name,
        "cas_number": exp.cas_number,
        "wavenumbers_cm": _to_list(exp.wavenumbers),
        "absorbance_normalized": _to_list(exp.absorbance),
    }


def _to_list(arr: Any) -> list[float]:
    if arr is None:
        return []
    if isinstance(arr, np.ndarray):
        return arr.tolist()
    return list(arr)


def _json_default(o: Any) -> Any:
    """Fallback serializer for numpy types that slip through."""
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if hasattr(o, "tolist"):
        return o.tolist()
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


# ---------------------------------------------------------------------------
# CSV companions
# ---------------------------------------------------------------------------


def _write_summary_metrics_csv(payload: dict[str, Any], out_path: Path) -> None:
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "method",
                "r2_freq",
                "r2_intensity",
                "rmse_freq",
                "mae_freq",
                "rmse_intensity",
                "mae_intensity",
                "speedup",
                "experimental_agreement",
            ]
        )
        for c in payload["comparisons"]:
            writer.writerow(
                [
                    c["name"],
                    c["metrics"]["r2_freq"],
                    c["metrics"]["r2_intensity"],
                    c["metrics"]["rmse_freq"],
                    c["metrics"]["mae_freq"],
                    c["metrics"]["rmse_intensity"],
                    c["metrics"]["mae_intensity"],
                    c["runtime"]["speedup"],
                    c.get("experimental_agreement", ""),
                ]
            )


def _write_spectrum_grid_csv(
    freq_grid: np.ndarray,
    analysis_results: dict[str, Any],
    out_path: Path,
) -> None:
    """Write frequency grid CSV. Full broadened columns land in Plan 04 integration."""
    header = ["freq_cm"]
    rows = [[float(x)] for x in np.asarray(freq_grid).tolist()]
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
