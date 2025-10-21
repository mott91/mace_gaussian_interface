"""
Results Manager for organizing and saving calculation outputs.

Handles directory structure, file naming, and JSON metadata writing
for comparison framework.
"""

import os
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class ResultsManager:
    """Manages output organization and metadata for molecular calculations."""
    
    def __init__(self, base_output_dir: str = "comparison_results"):
        """
        Initialize ResultsManager.
        
        Parameters
        ----------
        base_output_dir : str
            Base directory for all comparison results
        """
        self.base_output_dir = Path(base_output_dir)
        self.base_output_dir.mkdir(exist_ok=True)
        
    def create_molecule_directory(self, molecule_name: str) -> Path:
        """
        Create directory for a specific molecule.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule (e.g., 'acoh', 'water')
            
        Returns
        -------
        Path
            Path to molecule directory
        """
        mol_dir = self.base_output_dir / molecule_name
        mol_dir.mkdir(exist_ok=True)
        return mol_dir
        
    def create_optimization_directory(self, molecule_name: str) -> Path:
        """
        Create directory for geometry optimization results.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
            
        Returns
        -------
        Path
            Path to geometry_opt directory
        """
        mol_dir = self.create_molecule_directory(molecule_name)
        opt_dir = mol_dir / "geometry_opt"
        opt_dir.mkdir(exist_ok=True)
        return opt_dir
        
    def create_frequency_directory(
        self, 
        molecule_name: str,
        energy_calculator: str,
        dipole_calculator: str,
        timestamp: Optional[str] = None
    ) -> Path:
        """
        Create directory for frequency calculation results.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
        energy_calculator : str
            Energy calculator name (e.g., 'mace_mp', 'dft_wB97MV')
        dipole_calculator : str
            Dipole calculator name (e.g., 'espaloma', 'xtb')
        timestamp : str, optional
            Timestamp string, generated if not provided
            
        Returns
        -------
        Path
            Path to frequency calculation directory
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
        mol_dir = self.create_molecule_directory(molecule_name)
        
        # Create directory name
        dir_name = f"{energy_calculator}_{dipole_calculator}_{timestamp}"
        freq_dir = mol_dir / dir_name
        freq_dir.mkdir(exist_ok=True)
        
        return freq_dir
        
    def save_optimization_results(
        self,
        molecule_name: str,
        initial_atoms,
        final_atoms,
        calculator_name: str,
        initial_energy: float,
        final_energy: float,
        converged: bool,
        num_steps: int,
        runtime: float
    ):
        """
        Save geometry optimization results.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
        initial_atoms : ase.Atoms
            Initial molecular structure
        final_atoms : ase.Atoms
            Optimized molecular structure
        calculator_name : str
            Name of calculator used
        initial_energy : float
            Initial energy in eV
        final_energy : float
            Final energy in eV
        converged : bool
            Whether optimization converged
        num_steps : int
            Number of optimization steps
        runtime : float
            Runtime in seconds
        """
        from ase.io import write
        
        opt_dir = self.create_optimization_directory(molecule_name)
        
        # Save geometries
        write(str(opt_dir / "initial.xyz"), initial_atoms)
        write(str(opt_dir / "optimized.xyz"), final_atoms)
        
        # Create metadata
        metadata = {
            "molecule": molecule_name,
            "calculator": calculator_name,
            "initial_energy_eV": float(initial_energy),
            "final_energy_eV": float(final_energy),
            "converged": converged,
            "num_steps": num_steps,
            "runtime_s": float(runtime),
            "files": {
                "initial_geometry": "initial.xyz",
                "optimized_geometry": "optimized.xyz"
            }
        }
        
        # Save JSON
        json_path = opt_dir / "results.json"
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
            
        logger.info(f"✓ Saved optimization results to {opt_dir}")
        
    def save_frequency_results(
        self,
        molecule_name: str,
        energy_calculator: str,
        dipole_calculator: str,
        calculator_type: str,
        frequencies_data: Dict[str, Any],
        energy: float,
        dipole: Optional[Dict[str, Any]],
        runtime: float,
        gaussian_log: Optional[str] = None,
        gaussian_gjf: Optional[str] = None,
        timestamp: Optional[str] = None
    ):
        """
        Save frequency calculation results.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
        energy_calculator : str
            Energy calculator name
        dipole_calculator : str
            Dipole calculator name
        calculator_type : str
            Type of calculator ('ml' or 'dft')
        frequencies_data : dict
            Dictionary with 'harmonic' and 'anharmonic' frequency data
        energy : float
            Final energy in eV
        dipole : dict, optional
            Dipole moment information
        runtime : float
            Runtime in seconds
        gaussian_log : str, optional
            Path to Gaussian log file (should already be in freq_dir)
        gaussian_gjf : str, optional
            Path to Gaussian input file (should already be in freq_dir)
        timestamp : str, optional
            Timestamp for directory name
        """
        freq_dir = self.create_frequency_directory(
            molecule_name, energy_calculator, dipole_calculator, timestamp
        )
        
        # Just record the filenames - files should already be in place
        files = {}
        if gaussian_log:
            # Check if file exists and record relative name
            log_path = Path(gaussian_log)
            if log_path.exists():
                files["gaussian_log"] = log_path.name
                
        if gaussian_gjf:
            # Check if file exists and record relative name
            gjf_path = Path(gaussian_gjf)
            if gjf_path.exists():
                files["gaussian_input"] = gjf_path.name
        
        # Create metadata
        metadata = {
            "molecule": molecule_name,
            "energy_calculator": energy_calculator,
            "dipole_calculator": dipole_calculator,
            "calculator_type": calculator_type,
            "timestamp": timestamp or datetime.now().strftime("%Y%m%d_%H%M%S"),
            "frequencies": frequencies_data,
            "energy_eV": float(energy),
            "dipole": dipole,
            "runtime_s": float(runtime),
            "files": files
        }
        
        # Save JSON
        json_path = freq_dir / "results.json"
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
            
        logger.info(f"✓ Saved frequency results to {freq_dir}")
        
    def load_optimization_results(self, molecule_name: str) -> Dict[str, Any]:
        """
        Load optimization results from JSON.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
            
        Returns
        -------
        dict
            Optimization results metadata
        """
        opt_dir = self.base_output_dir / molecule_name / "geometry_opt"
        json_path = opt_dir / "results.json"
        
        if not json_path.exists():
            raise FileNotFoundError(f"No optimization results found for {molecule_name}")
            
        with open(json_path, 'r') as f:
            return json.load(f)
            
    def get_optimized_geometry_path(self, molecule_name: str) -> Path:
        """
        Get path to optimized geometry file.
        
        Parameters
        ----------
        molecule_name : str
            Name of the molecule
            
        Returns
        -------
        Path
            Path to optimized.xyz file
        """
        opt_dir = self.base_output_dir / molecule_name / "geometry_opt"
        return opt_dir / "optimized.xyz"