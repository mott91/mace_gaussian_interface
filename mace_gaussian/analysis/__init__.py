"""Analysis and reporting subpackage."""

# Note: analysis_workflow imports seaborn/matplotlib (heavy plotting deps).
# We do NOT eagerly import it here so that lightweight submodules like
# mode_matching can be imported in test environments without those deps.
# Callers that need analyze_molecule must import from .analysis_workflow directly.


def analyze_molecule(*args, **kwargs):
    """Lazy wrapper — imports analysis_workflow on first call."""
    from .analysis_workflow import analyze_molecule as _fn

    return _fn(*args, **kwargs)


def analyze_molecule_harmonic(*args, **kwargs):
    """Lazy wrapper — imports analysis_workflow on first call."""
    from .analysis_workflow import analyze_molecule_harmonic as _fn

    return _fn(*args, **kwargs)


def run_analysis_main(*args, **kwargs):
    """Lazy wrapper — imports analysis_workflow on first call."""
    from .analysis_workflow import run_analysis_main as _fn

    return _fn(*args, **kwargs)


def run_analysis_harmonic_main(*args, **kwargs):
    """Lazy wrapper — imports analysis_workflow on first call."""
    from .analysis_workflow import run_analysis_harmonic_main as _fn

    return _fn(*args, **kwargs)


__all__ = [
    "analyze_molecule",
    "analyze_molecule_harmonic",
    "run_analysis_main",
    "run_analysis_harmonic_main",
]
