"""One visual language for every figure: fixed colors per energy model.

Palette shared with thesis_style.py (colorblind-validated categorical slots). DFT is always black, experiment always grey.
Dipole models are distinguished by marker symbol, never by color, because they only
change intensities. Matplotlib/seaborn thesis figures should import from here too.
"""

from __future__ import annotations

DFT_COLOR = "#000000"
EXP_COLOR = "#7F7F7F"

# Same five slots as thesis_style.SERIES (the user's figstyle palette), so the
# HTML report and the thesis figures agree on every model's colour.
ENERGY_MODEL_COLORS: dict[str, str] = {
    "mace_omol": "#2a78d6",  # blue, slot 1
    "mace_mp": "#eb6834",  # orange, slot 2
    "mace_off": "#1baf7a",  # green, slot 3
    "mace_polar": "#eda100",  # yellow, slot 4
    "mace_anicc": "#e87ba4",  # pink, slot 5
}
_FALLBACK_COLORS = ["#56B4E9", "#F0E442", "#999999", "#8B4513", "#2F4F4F"]

# Dipole models get their own hues (Matplotlib's default cycle, distinct from the
# energy-model palette) for figures that compare dipole models within one energy model.
DIPOLE_COLORS: dict[str, str] = {  # Matplotlib tab10 order
    "mace_ml": "#1f77b4",  # blue
    "mace_polar1": "#ff7f0e",  # orange
    "espaloma": "#2ca02c",  # green
    "xtb": "#d62728",  # red
}

DIPOLE_SYMBOLS: dict[str, str] = {
    "mace_ml": "circle",
    "mace_polar1": "diamond",
    "espaloma": "square",
    "xtb": "cross",
}


def model_color(energy_model: str) -> str:
    """Color for an energy model; unknown names get a stable fallback."""
    if energy_model in ENERGY_MODEL_COLORS:
        return ENERGY_MODEL_COLORS[energy_model]
    return _FALLBACK_COLORS[hash(energy_model) % len(_FALLBACK_COLORS)]


def dipole_color(dipole_model: str | None) -> str:
    return DIPOLE_COLORS.get(dipole_model or "", "#333333")


def dipole_symbol(dipole_model: str | None) -> str:
    return DIPOLE_SYMBOLS.get(dipole_model or "", "circle")


def hex_to_rgba(color: str, alpha: float) -> str:
    c = color.lstrip("#")
    r, g, b = (int(c[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
