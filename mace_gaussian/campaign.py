"""Where a set of results lives: the legacy folders, or one isolated campaign.

The 2026 benchmark campaign recomputes molecules that already exist in the legacy
``comparison_results/`` (methanol, ammonia) under a new DFT protocol. Old and new data must
never meet: an analysis that pairs a new ML run with an old B3LYP baseline, or a cluster
retrieval that picks up an old log after a failed job, gives a wrong reference without any
error. So every entry point resolves its folders here, from a single campaign name:

- ``None`` (default): exactly the legacy layout, unchanged behaviour;
- ``"2026"``: everything under ``campaigns/2026/`` locally and ``<remote>/campaigns/2026``
  on the cluster.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_REMOTE_BASE = "/scratch_rune03a/mot/calculations/mace_gaussian"
_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


@dataclass(frozen=True)
class CampaignPaths:
    name: str | None
    comparison: Path  # raw results: <molecule>/<run>/results.json, logs, fchk
    analysis: Path  # anharmonic analyses (reports, report_data.json)
    harmonic: Path  # harmonic analyses (analysis_workflow appends "_harmonic")
    figures: Path  # thesis-style figures and galleries
    remote_base: str  # cluster scratch root for DFT jobs


def campaign_paths(name: str | None = None, root: Path = REPO) -> CampaignPaths:
    """Folders for the legacy layout (``name=None``) or the campaign ``name``."""
    if name is None:
        return CampaignPaths(
            name=None,
            comparison=Path("comparison_results"),
            analysis=Path("analysis_results"),
            harmonic=Path("analysis_results_harmonic"),
            figures=root / "thesis" / "figures",
            remote_base=DEFAULT_REMOTE_BASE,
        )
    if not _NAME.match(name):
        raise ValueError(f"campaign name {name!r}: letters, digits, '_' and '-' only")
    base = Path("campaigns") / name
    return CampaignPaths(
        name=name,
        comparison=base / "comparison_results",
        analysis=base / "analysis_results",
        harmonic=base / "analysis_results_harmonic",
        figures=root / "campaigns" / name / "figures",
        remote_base=f"{DEFAULT_REMOTE_BASE}/campaigns/{name}",
    )
