"""
Comprehensive Analysis Workflow - FIXED VERSION

Changes:
- PNG instead of PDF for HTML display
- Combined plots with all ML methods vs wb97xd DFT baseline  
- Filters wb97xd DFT as the reference baseline
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
    
    def find_dft_baseline(self, prefer_wb97xd: bool = True) -> Optional[Path]:
        """
        Find DFT anharmonic baseline results
        
        Scans all subdirectories for results.json files with calculator_type='dft'
        Prefers wb97xd if multiple DFT results exist
        
        Parameters
        ----------
        prefer_wb97xd : bool
            If True, prefer wb97xd DFT method when multiple DFT results exist
        
        Returns
        -------
        Path or None
            Path to DFT baseline results.json
        """
        dft_results = []
        
        # Scan all subdirectories for DFT results
        for item in self.molecule_dir.iterdir():
            if item.is_dir() and item.name != "geometry_opt":
                json_path = item / "results.json"
                if json_path.exists():
                    try:
                        with open(json_path, 'r') as f:
                            data = json.load(f)
                            if data.get('calculator_type') == 'dft':
                                dft_results.append((item.name, json_path))
                    except (json.JSONDecodeError, KeyError) as e:
                        logger.warning(f"Error reading {json_path}: {e}")
                        continue
        
        if not dft_results:
            logger.warning("No DFT anharmonic baseline found!")
            return None
        
        # Prefer wb97xd if requested and available
        if prefer_wb97xd:
            for name, path in dft_results:
                if 'wb97' in name.lower():
                    logger.info(f"Found DFT baseline (wb97xd): {path}")
                    return path
        
        # Otherwise return first DFT result
        logger.info(f"Found DFT baseline: {dft_results[0][1]}")
        return dft_results[0][1]
    
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
            if item.is_dir() and item.name != "geometry_opt":
                json_path = item / "results.json"
                if json_path.exists():
                    try:
                        with open(json_path, 'r') as f:
                            data = json.load(f)
                            # Check if it's an ML calculation
                            if data.get('calculator_type') == 'ml':
                                name = f"{data['energy_calculator']}_{data['dipole_calculator']}"
                                ml_results.append((name, json_path))
                                logger.info(f"Found ML result: {name}")
                    except (json.JSONDecodeError, KeyError) as e:
                        logger.warning(f"Error reading {json_path}: {e}")
                        continue
        
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
        # Match modes
        dft_freq, ml_freq, dft_int, ml_int, _ = self.analyzer.match_by_mode(
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
        
        # Create plots - SAVE AS PNG NOT PDF!
        spectrum_plot_path = self.plots_dir / f"spectrum_{ml_name}.png"
        regression_plot_path = self.plots_dir / f"regression_{ml_name}.png"
        
        self.analyzer.plot_spectra_comparison(
            ml_spectrum, dft_spectrum, ml_name, 
            molecule_name=self.molecule_name,
            save_path=str(spectrum_plot_path)
        )
        
        self.analyzer.plot_regression(
            ml_spectrum, dft_spectrum, metrics, ml_name,
            molecule_name=self.molecule_name,
            save_path=str(regression_plot_path)
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
            'comparison_df': comparison_df,
            'ml_spectrum': ml_spectrum,  # Store for combined plots
            'dft_spectrum': dft_spectrum
        }
    
    def create_combined_plots(self, comparisons: List[Dict], dft_spectrum: SpectrumData):
        """
        Create combined plots with all ML methods vs DFT
        
        Parameters
        ----------
        comparisons : list
            List of comparison results from run_single_comparison
        dft_spectrum : SpectrumData
            DFT reference spectrum
        """
        if not comparisons:
            return
        
        # Create combined spectrum plot (400-4000 cm⁻¹)
        combined_spectrum_path = self.plots_dir / "spectrum_combined.png"
        self.analyzer.plot_combined_spectra(
            ml_spectra=[c['ml_spectrum'] for c in comparisons],
            ml_names=[c['name'] for c in comparisons],
            dft_spectrum=dft_spectrum,
            molecule_name=self.molecule_name,
            save_path=str(combined_spectrum_path)
        )

        # Create extended combined spectrum plot (400-8000 cm⁻¹, includes overtones)
        combined_spectrum_extended_path = self.plots_dir / "spectrum_combined_extended.png"
        self.analyzer.plot_combined_spectra_extended(
            ml_spectra=[c['ml_spectrum'] for c in comparisons],
            ml_names=[c['name'] for c in comparisons],
            dft_spectrum=dft_spectrum,
            molecule_name=self.molecule_name,
            save_path=str(combined_spectrum_extended_path)
        )

        # Create combined regression plot
        combined_regression_path = self.plots_dir / "regression_combined.png"
        self.analyzer.plot_combined_regression(
            ml_spectra=[c['ml_spectrum'] for c in comparisons],
            ml_names=[c['name'] for c in comparisons],
            dft_spectrum=dft_spectrum,
            metrics_list=[c['metrics'] for c in comparisons],
            molecule_name=self.molecule_name,
            save_path=str(combined_regression_path)
        )
        
        logger.info("Created combined plots")
    
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
        
        # Find DFT baseline (prefer wb97xd)
        dft_path = self.find_dft_baseline(prefer_wb97xd=True)
        if dft_path is None:
            logger.error("=" * 60)
            logger.error(f"NO DFT BASELINE FOUND FOR {self.molecule_name.upper()}")
            logger.error("=" * 60)
            logger.error(f"Searched in: {self.molecule_dir}")
            logger.error("Could not find any directory with calculator_type='dft' in results.json")
            logger.error("=" * 60)
            
            # Return empty result instead of raising exception
            return {
                'molecule': self.molecule_name,
                'comparisons': [],
                'output_dir': self.output_dir,
                'error': 'No DFT baseline found'
            }
        
        # Find all ML results
        ml_results = self.find_ml_results()
        if not ml_results:
            logger.warning("=" * 60)
            logger.warning(f"NO ML RESULTS FOUND FOR {self.molecule_name.upper()}")
            logger.warning("=" * 60)
            logger.warning(f"Searched in: {self.molecule_dir}")
            logger.warning("No directories with calculator_type='ml' found.")
            logger.warning("=" * 60)
            
            return {
                'molecule': self.molecule_name,
                'comparisons': [],
                'output_dir': self.output_dir,
                'error': 'No ML results found'
            }
        
        logger.info(f"Found DFT baseline: {dft_path.parent.name}")
        logger.info(f"Found {len(ml_results)} ML calculations to compare")
        
        # Load DFT spectrum once for combined plots
        dft_results = self.analyzer.load_results(dft_path)
        dft_spectrum = self.analyzer.extract_spectrum_data(
            dft_results, include_overtones=True, include_combinations=True
        )
        
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
                import traceback
                traceback.print_exc()
        
        # Create combined plots
        self.create_combined_plots(comparisons, dft_spectrum)
        
        # Only save summary if we have comparisons
        if comparisons:
            # Save summary metrics
            metrics_summary = {
                'molecule': self.molecule_name,
                'analysis_date': datetime.now().isoformat(),
                'num_ml_calculators': len(comparisons),
                'dft_method': dft_path.parent.name,
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
            'output_dir': self.output_dir,
            'has_combined_plots': True
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
        
        # Only show success message if we have comparisons
        if analysis_results.get('comparisons'):
            logger.info(f"HTML report generated: {self.output_dir}/report.html")
        else:
            logger.warning(f"Error report generated: {self.output_dir}/report.html")


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
        print("Example: python comparison_workflow.py water")
        sys.exit(1)
    
    molecule_name = sys.argv[1]
    analyze_molecule(molecule_name)