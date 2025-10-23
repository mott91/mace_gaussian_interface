#!/usr/bin/env python3
"""
Simple convenience script to run IR spectral analysis

Usage:
    python run_analysis.py <molecule_name>
    
Example:
    python run_analysis.py acoh
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from comparison_workflow import analyze_molecule


def main():
    """Main entry point"""
    if len(sys.argv) < 2:
        print("=" * 60)
        print("IR SPECTRAL ANALYSIS FRAMEWORK")
        print("=" * 60)
        print("\nUsage:")
        print("  python run_analysis.py <molecule_name>")
        print("\nExample:")
        print("  python run_analysis.py acoh")
        print("\nThis will:")
        print("  1. Find DFT anharmonic baseline")
        print("  2. Find all ML calculation results")
        print("  3. Generate comparison plots")
        print("  4. Calculate statistical metrics")
        print("  5. Create HTML report")
        print("\n" + "=" * 60)
        sys.exit(1)
    
    molecule_name = sys.argv[1]
    
    print("\n" + "=" * 60)
    print(f"ANALYZING: {molecule_name.upper()}")
    print("=" * 60 + "\n")
    
    try:
        results = analyze_molecule(molecule_name)
        
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
        print(f"  - comparison_results/{molecule_name}/freq_anharm/results.json (DFT baseline)")
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
