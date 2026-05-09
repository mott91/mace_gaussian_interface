"""RED tests for Phase 23 plotly_builders module."""
import numpy as np
import pytest

plotly_builders = pytest.importorskip("mace_gaussian.analysis.plotly_builders")


def test_build_spectrum_figure_returns_plotly_figure():
    import plotly.graph_objects as go

    freq_grid = np.linspace(400, 4000, 200)
    dft_norm = np.zeros_like(freq_grid)
    ml_norm = np.zeros_like(freq_grid)
    fig = plotly_builders.build_spectrum_figure(
        freq_grid=freq_grid,
        dft_norm=dft_norm,
        ml_norm=ml_norm,
        ml_name="test",
    )
    assert isinstance(fig, go.Figure)
    # DFT + ML = 2 traces minimum
    assert len(fig.data) >= 2


def test_spectrum_includes_experimental():
    freq_grid = np.linspace(400, 4000, 200)
    dft_norm = np.zeros_like(freq_grid)
    ml_norm = np.zeros_like(freq_grid)
    exp_on_grid = np.zeros_like(freq_grid)
    fig = plotly_builders.build_spectrum_figure(
        freq_grid=freq_grid,
        dft_norm=dft_norm,
        ml_norm=ml_norm,
        ml_name="test",
        experimental_norm=exp_on_grid,
    )
    # DFT + ML + experimental = 3 traces
    assert len(fig.data) == 3
    names = [t.name for t in fig.data]
    assert any("xperimental" in n or "NIST" in n for n in names)


def test_spectrum_x_axis_reversed_for_spectroscopic_convention():
    freq_grid = np.linspace(400, 4000, 200)
    fig = plotly_builders.build_spectrum_figure(
        freq_grid=freq_grid,
        dft_norm=np.zeros_like(freq_grid),
        ml_norm=np.zeros_like(freq_grid),
        ml_name="test",
    )
    # High wavenumber on left = autorange reversed
    assert fig.layout.xaxis.autorange == "reversed"


def test_build_regression_figure_returns_plotly_figure():
    import plotly.graph_objects as go

    dft_freqs = np.array([1000.0, 2000.0, 3000.0])
    ml_freqs = np.array([1010.0, 2005.0, 2990.0])
    fig = plotly_builders.build_regression_figure(
        dft_freqs=dft_freqs,
        ml_freqs=ml_freqs,
        ml_name="test",
    )
    assert isinstance(fig, go.Figure)


def test_regression_figure_marks_low_overlap_with_open_marker():
    """A fundamental with eigenvector overlap below the threshold renders as
    a hollow ('-open' suffix) symbol in its own legend trace, separate from
    the high-overlap fundamentals.  Reviewers can then visually discount the
    flagged points without re-running the analysis.
    """
    dft_freqs = np.array([1000.0, 2000.0, 3000.0])
    ml_freqs = np.array([1010.0, 2050.0, 2990.0])
    mode_ids = ["F1", "F2", "F3"]
    # F2 has overlap 0.4 — below the 0.7 threshold — F1 and F3 are clean.
    mode_overlaps = [0.95, 0.4, 0.92]

    fig = plotly_builders.build_regression_figure(
        dft_freqs=dft_freqs,
        ml_freqs=ml_freqs,
        ml_name="test",
        mode_ids=mode_ids,
        mode_overlaps=mode_overlaps,
    )

    open_traces = [t for t in fig.data if getattr(t.marker, "symbol", None) == "circle-open"]
    closed_traces = [t for t in fig.data if getattr(t.marker, "symbol", None) == "circle"]
    assert len(open_traces) == 1, "low-overlap fundamental should produce one hollow trace"
    assert len(closed_traces) == 1, "high-overlap fundamentals should produce one filled trace"
    assert "low overlap" in open_traces[0].name
    # Two high-overlap fundamentals (F1, F3); one low-overlap (F2).
    assert len(closed_traces[0].x) == 2
    assert len(open_traces[0].x) == 1


def test_build_combined_spectrum_figure_multi_method():
    import plotly.graph_objects as go

    freq_grid = np.linspace(400, 4000, 200)
    ml_norms = {
        "method_a": np.zeros_like(freq_grid),
        "method_b": np.zeros_like(freq_grid),
    }
    fig = plotly_builders.build_combined_spectrum_figure(
        freq_grid=freq_grid,
        dft_norm=np.zeros_like(freq_grid),
        ml_norms=ml_norms,
    )
    assert isinstance(fig, go.Figure)
    # 1 DFT + 2 ML methods = at least 3 traces
    assert len(fig.data) >= 3


def test_combined_spectrum_x_axis_reversed():
    freq_grid = np.linspace(400, 4000, 200)
    ml_norms = {"method_a": np.zeros_like(freq_grid)}
    fig = plotly_builders.build_combined_spectrum_figure(
        freq_grid=freq_grid,
        dft_norm=np.zeros_like(freq_grid),
        ml_norms=ml_norms,
    )
    assert fig.layout.xaxis.autorange == "reversed"


def test_experimental_on_grid_returns_none_for_none():
    freq_grid = np.linspace(400, 4000, 200)
    result = plotly_builders.experimental_on_grid(None, freq_grid)
    assert result is None


def test_experimental_on_grid_normalizes():
    from types import SimpleNamespace

    freq_grid = np.linspace(400, 4000, 200)
    exp = SimpleNamespace(
        wavenumbers=np.linspace(400, 4000, 100),
        absorbance=np.random.rand(100) * 5.0,
    )
    result = plotly_builders.experimental_on_grid(exp, freq_grid)
    assert result is not None
    assert len(result) == len(freq_grid)
    assert np.max(result) <= 1.0 + 1e-10
