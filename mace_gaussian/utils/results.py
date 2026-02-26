"""
Results Manager for organizing and saving calculation outputs.

Handles directory structure, file naming, and JSON metadata writing
for comparison framework. Now includes overtones and combination bands.
"""

import json
import logging
import shutil
from pathlib import Path
from typing import Any, Optional

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
        timestamp: Optional[str] = None,
    ) -> Path:
        """Create directory for frequency calculation results."""
        mol_dir = self.create_molecule_directory(molecule_name)

        # Create directory name - avoid duplication for DFT
        if energy_calculator == dipole_calculator:
            # For DFT where both calculators are the same
            dir_name = energy_calculator
        else:
            # For ML where they're different
            dir_name = f"{energy_calculator}_{dipole_calculator}"

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
        runtime: float,
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

        # Collect version metadata
        try:
            from mace_gaussian.utils.validation import collect_version_metadata

            version_info = collect_version_metadata()
        except ImportError:
            version_info = {}

        # Create metadata
        metadata = {
            "molecule": molecule_name,
            "calculator": calculator_name,
            "initial_energy_eV": float(initial_energy),
            "final_energy_eV": float(final_energy),
            "converged": converged,
            "num_steps": num_steps,
            "runtime_s": float(runtime),
            "files": {"initial_geometry": "initial.xyz", "optimized_geometry": "optimized.xyz"},
            "version_info": version_info,
        }

        # Save JSON
        json_path = opt_dir / "results.json"
        with json_path.open("w") as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"\u2713 Saved optimization results to {opt_dir}")

    def save_frequency_results(
        self,
        molecule_name: str,
        energy_calculator: str,
        dipole_calculator: str,
        calculator_type: str,
        frequencies_data: dict[str, Any],
        energy: float,
        dipole: Optional[dict[str, Any]],
        runtime: float,
        gaussian_log: Optional[str] = None,
        gaussian_gjf: Optional[str] = None,
        timestamp: Optional[str] = None,
        calculation_parameters: Optional[dict[str, Any]] = None,
    ):
        """
        Save frequency calculation results with overtones and combination bands.

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
            Dictionary with 'harmonic', 'anharmonic', 'overtones', and 'combination_bands'
        energy : float
            Final energy in eV
        dipole : dict, optional
            Dipole moment information
        runtime : float
            Runtime in seconds
        gaussian_log : str, optional
            Path to Gaussian log file to copy
        gaussian_gjf : str, optional
            Path to Gaussian input file to copy
        timestamp : str, optional
            Timestamp for directory name
        """
        freq_dir = self.create_frequency_directory(
            molecule_name, energy_calculator, dipole_calculator, timestamp
        )

        # Copy Gaussian files if provided (skip if source and destination are the same)
        files = {}
        if gaussian_log and Path(gaussian_log).exists():
            dest_log = freq_dir / "gaussian_freq.log"
            # Only copy if source and destination are different files
            if Path(gaussian_log).resolve() != dest_log.resolve():
                shutil.copy2(gaussian_log, dest_log)
            files["gaussian_log"] = "gaussian_freq.log"

        if gaussian_gjf and Path(gaussian_gjf).exists():
            dest_gjf = freq_dir / "gaussian_freq.gjf"
            # Only copy if source and destination are different files
            if Path(gaussian_gjf).resolve() != dest_gjf.resolve():
                shutil.copy2(gaussian_gjf, dest_gjf)
            files["gaussian_input"] = "gaussian_freq.gjf"

        # Collect version metadata
        try:
            from mace_gaussian.utils.validation import collect_version_metadata

            version_info = collect_version_metadata()
        except ImportError:
            version_info = {}

        # Create metadata
        metadata = {
            "molecule": molecule_name,
            "energy_calculator": energy_calculator,
            "dipole_calculator": dipole_calculator,
            "calculator_type": calculator_type,
            "timestamp": None,
            "frequencies": frequencies_data,
            "energy_eV": float(energy),
            "dipole": dipole,
            "runtime_s": float(runtime),
            "files": files,
            "version_info": version_info,
            "calculation_parameters": calculation_parameters if calculation_parameters else {},
        }

        # Save JSON
        json_path = freq_dir / "results.json"
        with json_path.open("w") as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"\u2713 Saved frequency results to {freq_dir}")

    def load_optimization_results(self, molecule_name: str) -> dict[str, Any]:
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

        with json_path.open() as f:
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
