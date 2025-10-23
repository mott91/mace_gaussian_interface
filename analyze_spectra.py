"""
Comprehensive IR Spectral Analysis Script

Analyzes ML vs DFT frequency calculations with:
- Gaussian KDE broadening
- Regression analysis
- Statistical comparison
- Beautiful scientific visualizations
- HTML report generation
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy.stats import linregress
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("colorblind")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 300


@dataclass
class SpectrumData:
    """Container for spectral data"""
    frequencies: np.ndarray
    intensities: np.ndarray
    labels: List[str]  # e.g., ['fundamental', 'overtone', 'combination']
    

@dataclass
class ComparisonMetrics:
    """Statistical comparison metrics"""
    mae_freq: float  # Mean absolute error in frequencies
    rmse_freq: float  # Root mean square error in frequencies
    r2_freq: float  # R^2 for frequency correlation
    slope_freq: float
    intercept_freq: float
    mae_intensity: float
    r2_intensity: float
    max_error_freq: float
    num_peaks: int


class SpectrumAnalyzer:
    """Main analyzer for IR spectra comparison"""
    
    def __init__(self, 
                 freq_range: Tuple[float, float] = (400, 4000),
                 bandwidth_fwhm: float = 8.0,
                 freq_step: float = 0.5):
        """
        Initialize analyzer
        
        Parameters
        ----------
        freq_range : tuple
            Frequency range in cm^-1 (min, max)
        bandwidth_fwhm : float
            Full width at half maximum for Gaussian broadening in cm^-1
        freq_step : float
            Step size for frequency grid in cm^-1
        """
        self.freq_range = freq_range
        self.bandwidth_fwhm = bandwidth_fwhm
        self.freq_step = freq_step
        
        # Convert FWHM to Gaussian broadening parameter
        # FWHM = 2 * sqrt(2 * ln(2)) * sigma
        # For our exp(-B * (x - x0)^2), B = 1/(2*sigma^2)
        sigma = bandwidth_fwhm / (2 * np.sqrt(2 * np.log(2)))
        self.broad_param = 1.0 / (2 * sigma**2)
        
        # Create frequency grid
        self.freq_grid = np.arange(freq_range[0], freq_range[1], freq_step)
        
        logger.info(f"Initialized analyzer: {freq_range[0]}-{freq_range[1]} cm^-1, "
                   f"FWHM={bandwidth_fwhm} cm^-1, step={freq_step} cm^-1")
    
    def load_results(self, json_path: str) -> Dict:
        """Load results from JSON file"""
        with open(json_path, 'r') as f:
            data = json.load(f)
        return data
    
    def extract_spectrum_data(self, results: Dict, 
                             include_overtones: bool = True,
                             include_combinations: bool = True) -> SpectrumData:
        """
        Extract frequency and intensity data from results
        
        Parameters
        ----------
        results : dict
            Results dictionary from JSON
        include_overtones : bool
            Include overtones in spectrum
        include_combinations : bool
            Include combination bands in spectrum
            
        Returns
        -------
        SpectrumData
            Combined spectrum data
        """
        frequencies = []
        intensities = []
        labels = []
        
        # Get anharmonic fundamentals
        anharmonic = results.get('frequencies', {}).get('anharmonic', [])
        for entry in anharmonic:
            frequencies.append(entry['freq_cm'])
            intensities.append(entry['ir_intensity'])
            labels.append('fundamental')
        
        # Add overtones
        if include_overtones:
            overtones = results.get('frequencies', {}).get('overtones', [])
            for entry in overtones:
                frequencies.append(entry['freq_anharmonic'])
                intensities.append(entry['ir_intensity'])
                labels.append('overtone')
        
        # Add combination bands
        if include_combinations:
            combinations = results.get('frequencies', {}).get('combination_bands', [])
            for entry in combinations:
                frequencies.append(entry['freq_anharmonic'])
                intensities.append(entry['ir_intensity'])
                labels.append('combination')
        
        return SpectrumData(
            frequencies=np.array(frequencies),
            intensities=np.array(intensities),
            labels=labels
        )
    
    def broaden_spectrum(self, spectrum: SpectrumData) -> np.ndarray:
        """
        Apply Gaussian broadening to discrete spectral lines
        
        Parameters
        ----------
        spectrum : SpectrumData
            Discrete spectral lines
            
        Returns
        -------
        np.ndarray
            Broadened spectrum on frequency grid
        """
        broadened = np.zeros_like(self.freq_grid)
        
        for freq, intensity in zip(spectrum.frequencies, spectrum.intensities):
            # Calculate Gaussian contribution from this line
            # I(v) = I_0 * exp(-B * (v - v_0)^2)
            argument = -self.broad_param * (self.freq_grid - freq)**2
            
            # Avoid numerical underflow
            mask = argument > -50  # exp(-50) ~ 2e-22
            broadened[mask] += intensity * np.exp(argument[mask])
        
        return broadened
    
    def match_peaks(self, spectrum1: SpectrumData, 
                   spectrum2: SpectrumData,
                   tolerance: float = 50.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Match peaks between two spectra based on frequency proximity
        
        Parameters
        ----------
        spectrum1, spectrum2 : SpectrumData
            Spectra to match
        tolerance : float
            Maximum frequency difference for matching in cm^-1
            
        Returns
        -------
        freq1, freq2, int1, int2 : np.ndarray
            Matched frequencies and intensities
        """
        freq1_matched = []
        freq2_matched = []
        int1_matched = []
        int2_matched = []
        
        # For each peak in spectrum1, find closest in spectrum2
        used_indices = set()
        
        for i, f1 in enumerate(spectrum1.frequencies):
            # Find closest peak in spectrum2
            differences = np.abs(spectrum2.frequencies - f1)
            min_idx = np.argmin(differences)
            
            if differences[min_idx] <= tolerance and min_idx not in used_indices:
                freq1_matched.append(f1)
                freq2_matched.append(spectrum2.frequencies[min_idx])
                int1_matched.append(spectrum1.intensities[i])
                int2_matched.append(spectrum2.intensities[min_idx])
                used_indices.add(min_idx)
        
        return (np.array(freq1_matched), np.array(freq2_matched),
                np.array(int1_matched), np.array(int2_matched))
    
    def calculate_metrics(self, ml_spectrum: SpectrumData, 
                         dft_spectrum: SpectrumData) -> ComparisonMetrics:
        """
        Calculate comparison metrics between ML and DFT spectra
        
        Parameters
        ----------
        ml_spectrum, dft_spectrum : SpectrumData
            Spectra to compare
            
        Returns
        -------
        ComparisonMetrics
            Statistical comparison metrics
        """
        # Match peaks
        dft_freq, ml_freq, dft_int, ml_int = self.match_peaks(
            dft_spectrum, ml_spectrum
        )
        
        if len(dft_freq) == 0:
            logger.warning("No matching peaks found!")
            return ComparisonMetrics(0, 0, 0, 0, 0, 0, 0, 0, 0)
        
        # Frequency metrics
        freq_errors = ml_freq - dft_freq
        mae_freq = np.mean(np.abs(freq_errors))
        rmse_freq = np.sqrt(np.mean(freq_errors**2))
        max_error_freq = np.max(np.abs(freq_errors))
        
        # Regression for frequencies
        slope_freq, intercept_freq, r_value_freq, _, _ = linregress(dft_freq, ml_freq)
        r2_freq = r_value_freq**2
        
        # Intensity metrics
        int_errors = ml_int - dft_int
        mae_intensity = np.mean(np.abs(int_errors))
        
        # Regression for intensities (if intensities exist)
        if len(dft_int) > 0 and np.any(dft_int > 0):
            _, _, r_value_int, _, _ = linregress(dft_int, ml_int)
            r2_intensity = r_value_int**2
        else:
            r2_intensity = 0.0
        
        return ComparisonMetrics(
            mae_freq=mae_freq,
            rmse_freq=rmse_freq,
            r2_freq=r2_freq,
            slope_freq=slope_freq,
            intercept_freq=intercept_freq,
            mae_intensity=mae_intensity,
            r2_intensity=r2_intensity,
            max_error_freq=max_error_freq,
            num_peaks=len(dft_freq)
        )
    
    def plot_spectra_comparison(self, 
                               ml_spectrum: SpectrumData,
                               dft_spectrum: SpectrumData,
                               ml_name: str,
                               save_path: Optional[str] = None) -> plt.Figure:
        """
        Create comparison plot of broadened spectra
        
        Parameters
        ----------
        ml_spectrum, dft_spectrum : SpectrumData
            Spectra to plot
        ml_name : str
            Name of ML calculator for title
        save_path : str, optional
            Path to save figure
            
        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), 
                                       gridspec_kw={'height_ratios': [2, 1]})
        
        # Broaden spectra
        ml_broadened = self.broaden_spectrum(ml_spectrum)
        dft_broadened = self.broaden_spectrum(dft_spectrum)
        
        # Panel A: Broadened spectra (transmittance style - inverted)
        ax1.plot(self.freq_grid, -dft_broadened, 'k-', linewidth=1.5, 
                label='DFT (anharmonic)', alpha=0.7)
        ax1.plot(self.freq_grid, -ml_broadened, '-', linewidth=1.5, 
                color='#E24A33', label=f'ML ({ml_name})')
        
        ax1.set_ylabel('Intensity (a.u.)', fontsize=11)
        ax1.set_xlim(self.freq_range)
        ax1.legend(loc='upper right', frameon=True)
        ax1.grid(True, alpha=0.3)
        ax1.set_title(f'IR Spectrum Comparison: {ml_name} vs DFT', 
                     fontsize=12, fontweight='bold')
        
        # Panel B: Stick spectra
        # Plot DFT sticks
        for freq, intensity, label in zip(dft_spectrum.frequencies, 
                                          dft_spectrum.intensities,
                                          dft_spectrum.labels):
            color = {'fundamental': '#348ABD', 'overtone': '#988ED5', 
                    'combination': '#8EBA42'}.get(label, 'gray')
            ax2.vlines(freq, 0, -intensity, colors='black', 
                      linewidth=1, alpha=0.4)
        
        # Plot ML sticks
        y_offset = np.max(ml_spectrum.intensities) * 0.1
        for freq, intensity, label in zip(ml_spectrum.frequencies,
                                          ml_spectrum.intensities,
                                          ml_spectrum.labels):
            color = {'fundamental': '#348ABD', 'overtone': '#988ED5',
                    'combination': '#8EBA42'}.get(label, 'gray')
            ax2.vlines(freq, 0, intensity, colors=color, linewidth=1, alpha=0.7)
        
        ax2.axhline(y=0, color='black', linewidth=0.8)
        ax2.set_xlabel('Wavenumber (cm^-1)', fontsize=11)
        ax2.set_ylabel('Intensity', fontsize=11)
        ax2.set_xlim(self.freq_range)
        ax2.grid(True, alpha=0.3, axis='x')
        ax2.text(0.02, 0.95, 'DFT (below) | ML (above)', 
                transform=ax2.transAxes, fontsize=9, 
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved spectrum plot to {save_path}")
        
        return fig
    
    def plot_regression(self,
                       ml_spectrum: SpectrumData,
                       dft_spectrum: SpectrumData,
                       metrics: ComparisonMetrics,
                       ml_name: str,
                       save_path: Optional[str] = None) -> plt.Figure:
        """
        Create regression plot for frequency correlation
        
        Parameters
        ----------
        ml_spectrum, dft_spectrum : SpectrumData
            Spectra to compare
        metrics : ComparisonMetrics
            Calculated metrics
        ml_name : str
            Name of ML calculator
        save_path : str, optional
            Path to save figure
            
        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Match peaks for plotting
        dft_freq, ml_freq, dft_int, ml_int = self.match_peaks(
            dft_spectrum, ml_spectrum
        )
        
        # Panel A: Frequency correlation
        colors = {'fundamental': '#348ABD', 'overtone': '#988ED5',
                 'combination': '#8EBA42'}
        
        # Color by label type (need to track labels through matching)
        scatter_colors = []
        for i, f_ml in enumerate(ml_freq):
            idx = np.argmin(np.abs(ml_spectrum.frequencies - f_ml))
            label = ml_spectrum.labels[idx]
            scatter_colors.append(colors.get(label, 'gray'))
        
        ax1.scatter(dft_freq, ml_freq, c=scatter_colors, s=50, 
                   alpha=0.6, edgecolors='black', linewidth=0.5)
        
        # Perfect agreement line
        lim_min = min(dft_freq.min(), ml_freq.min())
        lim_max = max(dft_freq.max(), ml_freq.max())
        ax1.plot([lim_min, lim_max], [lim_min, lim_max], 
                'k--', linewidth=1, alpha=0.5, label='Perfect agreement')
        
        # Regression line
        regression_line = metrics.slope_freq * dft_freq + metrics.intercept_freq
        ax1.plot(dft_freq, regression_line, 'r-', linewidth=1.5, 
                alpha=0.7, label='Linear fit')
        
        ax1.set_xlabel('DFT Frequency (cm^-1)', fontsize=11)
        ax1.set_ylabel('ML Frequency (cm^-1)', fontsize=11)
        ax1.set_title('Frequency Correlation', fontsize=12, fontweight='bold')
        ax1.legend(loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # Add statistics text box
        textstr = f'R^2 = {metrics.r2_freq:.4f}\n'
        textstr += f'MAE = {metrics.mae_freq:.2f} cm^-1\n'
        textstr += f'RMSE = {metrics.rmse_freq:.2f} cm^-1\n'
        textstr += f'Slope = {metrics.slope_freq:.4f}\n'
        textstr += f'n = {metrics.num_peaks}'
        
        ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Panel B: Intensity correlation
        if len(dft_int) > 0:
            ax2.scatter(dft_int, ml_int, c=scatter_colors, s=50,
                       alpha=0.6, edgecolors='black', linewidth=0.5)
            
            # Perfect agreement line
            lim_min = 0
            lim_max = max(dft_int.max(), ml_int.max()) * 1.1
            ax2.plot([lim_min, lim_max], [lim_min, lim_max],
                    'k--', linewidth=1, alpha=0.5, label='Perfect agreement')
            
            ax2.set_xlabel('DFT Intensity (km/mol)', fontsize=11)
            ax2.set_ylabel('ML Intensity (km/mol)', fontsize=11)
            ax2.set_title('Intensity Correlation', fontsize=12, fontweight='bold')
            ax2.legend(loc='upper left')
            ax2.grid(True, alpha=0.3)
            
            # Add R^2 for intensity
            textstr = f'R^2 = {metrics.r2_intensity:.4f}\n'
            textstr += f'MAE = {metrics.mae_intensity:.2f}'
            ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes,
                    fontsize=9, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.suptitle(f'Regression Analysis: {ml_name} vs DFT',
                    fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved regression plot to {save_path}")
        
        return fig


def main():
    """Example usage"""
    # This will be expanded with full analysis workflow
    analyzer = SpectrumAnalyzer(
        freq_range=(400, 4000),
        bandwidth_fwhm=8.0,
        freq_step=0.5
    )
    
    logger.info("Spectrum analyzer initialized successfully")
    print("\nTo use this analyzer, call the analysis functions with your data paths")
    print("Example: analyzer.load_results('path/to/results.json')")


if __name__ == "__main__":
    main()
