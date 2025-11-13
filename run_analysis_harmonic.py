#!/usr/bin/env python3
"""
Harmonic-only IR spectral analysis

This script analyzes ONLY harmonic fundamental frequencies with proper
eigenvector-based mode matching. Overtones and combinations are excluded
since they don't have meaningful mode matching via eigenvector dot products.

Usage:
    python run_analysis_harmonic.py <molecule_name>

Example:
    python run_analysis_harmonic.py water
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from comparison_workflow import analyze_molecule_harmonic


def main():
    """Main entry point"""
    if len(sys.argv) < 2:
        print("=" * 60)
        print("HARMONIC-ONLY IR SPECTRAL ANALYSIS")
        print("=" * 60)
        print("\nUsage:")
        print("  python run_analysis_harmonic.py <molecule_name>")
        print("\nExample:")
        print("  python run_analysis_harmonic.py water")
        print("\nThis will:")
        print("  1. Analyze ONLY harmonic fundamental frequencies")
        print("  2. Use eigenvector dot product mode matching")
        print("  3. Generate comparison plots (fundamentals only)")
        print("  4. Create HTML report")
        print("  5. Save to analysis_results_harmonic/")
        print("\nNote: Overtones and combinations are excluded")
        print("=" * 60)
        sys.exit(1)

    molecule_name = sys.argv[1]

    print("\n" + "=" * 60)
    print(f"HARMONIC ANALYSIS: {molecule_name.upper()}")
    print("=" * 60 + "\n")

    try:
        results = analyze_molecule_harmonic(molecule_name)

        print("\n" + "=" * 60)
        print("SUCCESS!")
        print("=" * 60)
        print(f"\nResults directory: {results['output_dir']}")
        print(f"HTML Report: {results['output_dir']}/report.html")
        print("\nOpen the HTML report in your browser to view the analysis!")
        print("=" * 60 + "\n")

    except FileNotFoundError as e:
        print("\n" + "=" * 60)
        print("ERROR: Files not found")
        print("=" * 60)
        print(f"\n{e}")
        print("\nMake sure you have:")
        print(f"  - comparison_results/{molecule_name}/<dft>/results.json (DFT baseline)")
        print(f"  - comparison_results/{molecule_name}/<calculator>/results.json (ML results)")
        print("=" * 60 + "\n")
        sys.exit(1)

    except Exception as e:
        print("\n" + "=" * 60)
        print("ERROR: Analysis failed")
        print("=" * 60)
        print(f"\n{e}")
        print("\nCheck the logs above for details.")
        print("=" * 60 + "\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
