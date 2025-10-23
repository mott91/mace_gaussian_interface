"""
Comprehensive Analysis Workflow

Scans results directories, compares all ML calculators against DFT baseline,
generates plots, tables, and comprehensive HTML report.
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
from datetime import datetime
import logging

from analyze_spectra import SpectrumAnalyzer, SpectrumData, ComparisonMetrics

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ComparisonWorkflow:
    """Orchestrates full comparison analysis workflow"""
    
    def __init__(self, 
                 molecule_name: str,
                 base_results_dir: str = "comparison_results",
                 output_dir: str = "analysis_results"):
        """
        Initialize workflow
        
        Parameters
        ----------
        molecule_name : str
            Name of molecule to analyze
        base_results_dir : str
            Base directory containing comparison_results
        output_dir : str
            Directory for analysis outputs
        """
        self.molecule_name = molecule_name
        self.base_results_dir = Path(base_results_dir)
        self.molecule_dir = self.base_results_dir / molecule_name
        
        # Create output directories
        self.output_dir = Path(output_dir) / molecule_name
        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"
        
        for dir_path in [self.output_dir, self.plots_dir, self.data_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize analyzer
        self.analyzer = SpectrumAnalyzer(
            freq_range=(400, 4000),
            bandwidth_fwhm=8.0,
            freq_step=0.5
        )
        
        logger.info(f"Initialized workflow for {molecule_name}")
    
    def find_dft_baseline(self):
        """Find DFT baseline by looking for calculator_type: 'dft'"""
        for item in self.molecule_dir.iterdir():
            if item.is_dir():
                json_path = item / "results.json"
                if json_path.exists():
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                        if data.get('calculator_type') == 'dft':
                            return json_path
        return None
    
    def find_ml_results(self) -> List[Tuple[str, Path]]:
        """
        Find all ML calculation results
        
        Returns
        -------
        list of (name, path) tuples
            List of ML calculator combinations and their result paths
        """
        ml_results = []
        
        # Scan molecule directory for ML results
        for item in self.molecule_dir.iterdir():
            if item.is_dir() and item.name != "geometry_opt" and item.name != "freq_anharm":
                json_path = item / "results.json"
                if json_path.exists():
                    # Check if it's an ML calculation
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                        if data.get('calculator_type') == 'ml':
                            name = f"{data['energy_calculator']}_{data['dipole_calculator']}"
                            ml_results.append((name, json_path))
                            logger.info(f"Found ML result: {name}")
        
        return ml_results
    
    def create_comparison_table(self,
                               dft_spectrum: SpectrumData,
                               ml_spectrum: SpectrumData,
                               ml_name: str) -> pd.DataFrame:
        """
        Create detailed comparison table
        
        Parameters
        ----------
        dft_spectrum, ml_spectrum : SpectrumData
            Spectra to compare
        ml_name : str
            Name of ML calculator
            
        Returns
        -------
        pd.DataFrame
            Comparison table
        """
        # Match peaks
        dft_freq, ml_freq, dft_int, ml_int = self.analyzer.match_peaks(
            dft_spectrum, ml_spectrum
        )
        
        # Create DataFrame
        df = pd.DataFrame({
            'DFT_Frequency_cm': dft_freq,
            'ML_Frequency_cm': ml_freq,
            'Freq_Difference_cm': ml_freq - dft_freq,
            'Percent_Error': 100 * (ml_freq - dft_freq) / dft_freq,
            'DFT_Intensity_km_mol': dft_int,
            'ML_Intensity_km_mol': ml_int,
            'Intensity_Difference': ml_int - dft_int
        })
        
        # Sort by DFT frequency
        df = df.sort_values('DFT_Frequency_cm').reset_index(drop=True)
        
        return df
    
    def run_single_comparison(self,
                            ml_name: str,
                            ml_path: Path,
                            dft_path: Path) -> Dict:
        """
        Run complete comparison for one ML calculator
        
        Parameters
        ----------
        ml_name : str
            Name of ML calculator
        ml_path : Path
            Path to ML results.json
        dft_path : Path
            Path to DFT results.json
            
        Returns
        -------
        dict
            Comparison results including metrics and file paths
        """
        logger.info(f"Comparing {ml_name} vs DFT...")
        
        # Load results
        ml_results = self.analyzer.load_results(ml_path)
        dft_results = self.analyzer.load_results(dft_path)
        
        # Extract spectra (with overtones and combinations)
        ml_spectrum = self.analyzer.extract_spectrum_data(
            ml_results, include_overtones=True, include_combinations=True
        )
        dft_spectrum = self.analyzer.extract_spectrum_data(
            dft_results, include_overtones=True, include_combinations=True
        )
        
        logger.info(f"  ML peaks: {len(ml_spectrum.frequencies)}, "
                   f"DFT peaks: {len(dft_spectrum.frequencies)}")
        
        # Calculate metrics
        metrics = self.analyzer.calculate_metrics(ml_spectrum, dft_spectrum)
        
        # Create plots
        spectrum_plot_path = self.plots_dir / f"spectrum_{ml_name}.pdf"
        regression_plot_path = self.plots_dir / f"regression_{ml_name}.pdf"
        
        self.analyzer.plot_spectra_comparison(
            ml_spectrum, dft_spectrum, ml_name, str(spectrum_plot_path)
        )
        
        self.analyzer.plot_regression(
            ml_spectrum, dft_spectrum, metrics, ml_name, str(regression_plot_path)
        )
        
        # Create comparison table
        comparison_df = self.create_comparison_table(
            dft_spectrum, ml_spectrum, ml_name
        )
        
        table_path = self.data_dir / f"comparison_{ml_name}.csv"
        comparison_df.to_csv(table_path, index=False, float_format='%.4f')
        
        # Get runtime info
        ml_runtime = ml_results.get('runtime_s', 0)
        dft_runtime = dft_results.get('runtime_s', 0)
        speedup = dft_runtime / ml_runtime if ml_runtime > 0 else 0
        
        return {
            'name': ml_name,
            'metrics': metrics,
            'ml_peaks': len(ml_spectrum.frequencies),
            'dft_peaks': len(dft_spectrum.frequencies),
            'ml_runtime': ml_runtime,
            'dft_runtime': dft_runtime,
            'speedup': speedup,
            'spectrum_plot': spectrum_plot_path.name,
            'regression_plot': regression_plot_path.name,
            'table_file': table_path.name,
            'comparison_df': comparison_df
        }
    
    def run_full_analysis(self) -> Dict:
        """
        Run complete analysis workflow
        
        Returns
        -------
        dict
            Complete analysis results for all comparisons
        """
        logger.info("=" * 60)
        logger.info(f"STARTING ANALYSIS FOR {self.molecule_name.upper()}")
        logger.info("=" * 60)
        
        # Find DFT baseline
        dft_path = self.find_dft_baseline()
        if dft_path is None:
            raise FileNotFoundError(
                f"No DFT baseline found for {self.molecule_name}"
            )
        
        # Find all ML results
        ml_results = self.find_ml_results()
        if not ml_results:
            raise FileNotFoundError(
                f"No ML results found for {self.molecule_name}"
            )
        
        logger.info(f"Found {len(ml_results)} ML calculations to compare")
        
        # Run comparisons
        comparisons = []
        for ml_name, ml_path in ml_results:
            try:
                result = self.run_single_comparison(ml_name, ml_path, dft_path)
                comparisons.append(result)
                logger.info(f"  [OK] {ml_name}: MAE={result['metrics'].mae_freq:.2f} cm^-1, "
                          f"R^2={result['metrics'].r2_freq:.4f}")
            except Exception as e:
                logger.error(f"  [FAIL] {ml_name} failed: {e}")
        
        # Save summary metrics
        metrics_summary = {
            'molecule': self.molecule_name,
            'analysis_date': datetime.now().isoformat(),
            'num_ml_calculators': len(comparisons),
            'comparisons': [
                {
                    'name': c['name'],
                    'mae_freq': c['metrics'].mae_freq,
                    'rmse_freq': c['metrics'].rmse_freq,
                    'r2_freq': c['metrics'].r2_freq,
                    'r2_intensity': c['metrics'].r2_intensity,
                    'speedup': c['speedup']
                }
                for c in comparisons
            ]
        }
        
        with open(self.data_dir / "metrics_summary.json", 'w') as f:
            json.dump(metrics_summary, f, indent=2)
        
        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info("=" * 60)
        
        return {
            'molecule': self.molecule_name,
            'comparisons': comparisons,
            'output_dir': self.output_dir
        }
    
    def generate_html_report(self, analysis_results: Dict):
        """
        Generate comprehensive HTML report
        
        Parameters
        ----------
        analysis_results : dict
            Results from run_full_analysis()
        """
        from html_report_generator import HTMLReportGenerator
        
        generator = HTMLReportGenerator(
            molecule_name=self.molecule_name,
            output_dir=self.output_dir
        )
        
        generator.generate_report(analysis_results)
        
        logger.info(f"HTML report generated: {self.output_dir}/report.html")


def analyze_molecule(molecule_name: str,
                    base_results_dir: str = "comparison_results",
                    output_dir: str = "analysis_results") -> Dict:
    """
    Convenience function to run complete analysis
    
    Parameters
    ----------
    molecule_name : str
        Name of molecule to analyze
    base_results_dir : str
        Base directory containing comparison_results
    output_dir : str
        Directory for analysis outputs
        
    Returns
    -------
    dict
        Analysis results
    """
    workflow = ComparisonWorkflow(
        molecule_name=molecule_name,
        base_results_dir=base_results_dir,
        output_dir=output_dir
    )
    
    results = workflow.run_full_analysis()
    workflow.generate_html_report(results)
    
    return results


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python comparison_workflow.py <molecule_name>")
        print("Example: python comparison_workflow.py acoh")
        sys.exit(1)
    
    molecule_name = sys.argv[1]
    analyze_molecule(molecule_name)
