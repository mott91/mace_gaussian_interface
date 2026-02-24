"""MACE-Gaussian interface: ML potentials + Gaussian for IR spectroscopy."""

__version__ = "0.2.0"

# run_pipeline is imported lazily here to avoid circular imports.
# workflow.py uses relative imports, so this is safe once workflow.py
# is inside the package (Plan 02). For now, leave this commented out
# — it will be activated in Plan 02 after workflow.py is moved.
# from .workflow import run_pipeline

__all__ = ["__version__"]
