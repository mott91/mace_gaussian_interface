"""RED tests for Phase 23 _shared_css module."""
import pytest


def test_build_css_returns_style_block():
    from mace_gaussian.analysis._shared_css import build_css

    css = build_css()
    assert "<style>" in css
    assert "</style>" in css


def test_build_css_includes_required_class_names():
    from mace_gaussian.analysis._shared_css import build_css

    css = build_css()
    required = [
        "executive-summary",
        "method-card",
        "best-method",
        "data-table",
        "comparison-section",
        "verdict",
        "plotly-graph-div",
    ]
    for cls in required:
        assert cls in css, f"Missing class: {cls}"


def test_build_css_is_static_string():
    from mace_gaussian.analysis._shared_css import build_css

    # Calling twice returns identical content
    assert build_css() == build_css()
