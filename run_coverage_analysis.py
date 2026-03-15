#!/usr/bin/env python3
"""Frequency range coverage analysis across molecules.

Usage: python run_coverage_analysis.py water aspirin [molecule3 ...]
"""
import sys

from mace_gaussian.analysis.coverage_analysis import CoverageAnalyzer


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_coverage_analysis.py molecule1 [molecule2 ...]")
        print("Example: python run_coverage_analysis.py water aspirin")
        sys.exit(1)

    molecules = sys.argv[1:]
    analyzer = CoverageAnalyzer(molecules=molecules)
    output_dir = analyzer.run()
    print(f"\nCoverage analysis complete. Results in: {output_dir}/")


if __name__ == "__main__":
    main()
