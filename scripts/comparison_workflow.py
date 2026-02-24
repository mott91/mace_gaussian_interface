#!/usr/bin/env python3
"""Research script for comparison workflow — delegates to mace_gaussian.analysis."""
from mace_gaussian.analysis.analysis_workflow import ComparisonWorkflow, analyze_molecule, analyze_molecule_harmonic  # noqa: F401

if __name__ == "__main__":
    import sys

    print("Use run_analysis.py or run_analysis_harmonic.py for analysis.")
    print("Or import from mace_gaussian.analysis.analysis_workflow directly.")
    sys.exit(0)
