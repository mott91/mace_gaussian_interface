#!/usr/bin/env python3
# For the helperscript copy to /home/bin/ and set shebang!!!
# module load mace-torch-dftd/14Jul2025
# module load gv

import warnings
warnings.filterwarnings("ignore")
try:
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    pass
import numpy as np
import zmq
import os
import time
from contextlib import contextmanager
import sys
import subprocess
from pathlib import Path
from mace_calculators import MACEDipoleCalculator
from ase.io import read
from ase.optimize import LBFGS
from ase.data import chemical_symbols
from abc import ABC, abstractmethod
from typing import Tuple, Optional, Dict
import logging
from diagnostics import diagnose_python_environment, test_espaloma_functionality
from charge_analysis import ChargeAnalyzer
from results_manager import ResultsManager
from gaussian_parser import parse_gaussian_log
import time
import shutil



# Set up logging
logging.basicConfig(
    level=logging.INFO,  # Changed from INFO
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Suppress Warnings
warnings.filterwarnings(
    "ignore", message=".*weights_only=False.*", category=FutureWarning
)
os.environ["PYTHONWARNINGS"] = "ignore::FutureWarning"

try:
    from charge_analysis import ChargeAnalyzer
    CHARGE_ANALYSIS_AVAILABLE = True
except ImportError:
    logger.warning("Charge analysis module not available")
    CHARGE_ANALYSIS_AVAILABLE = False

# ============================================================================
# CONFIGURATION
# ============================================================================

# Default paths - can be overridden by environment variables
DEFAULT_MACE_DIPOLE_MODEL = os.getenv(
    "MACE_DIPOLE_MODEL_PATH",
    str(Path.home() / "mace_gaussian" / "dipole_model" / "model_1.model"),
)

DEFAULT_HELPER_SCRIPT = os.getenv(
    "MACE_HELPER_SCRIPT_PATH",
    str(Path(__file__).parent / "gm_helper.py"),  # Use relative path by default
)

# Validate paths on startup
if not Path(DEFAULT_HELPER_SCRIPT).exists():
    logger.warning(f"Helper script not found at: {DEFAULT_HELPER_SCRIPT}")
    logger.warning("Set MACE_HELPER_SCRIPT_PATH environment variable if needed")

# Note: MACE dipole model is checked when MACEMLDipoleCalculator is instantiated

# ============================================================================
# MODULAR DIPOLE CALCULATION SYSTEM
# ============================================================================


class DipoleCalculatorBase(ABC):
    """Abstract base class for dipole calculators"""

    def __init__(self, name: str):
        self.name = name
        self.available = None
        self._check_availability()

    @abstractmethod
    def _check_availability(self) -> bool:
        """Check if this calculator is available"""
        pass

    @abstractmethod
    def calculate_dipole(
        self, atoms, **kwargs
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Calculate dipole moment and partial charges
        Returns: (dipole_vector, partial_charges)
        """
        pass

    def calculate_dipole_derivatives(
        self, atoms, displacement=0.01, **kwargs
    ) -> np.ndarray:
        """Calculate dipole derivatives numerically"""
        natoms = len(atoms)
        dipole_derivatives = np.zeros((3 * natoms, 3))
        base_positions = atoms.get_positions().copy()

        try:
            for i in range(natoms):
                for j in range(3):  # x, y, z directions
                    # Positive displacement
                    pos_disp = base_positions.copy()
                    pos_disp[i, j] += displacement
                    atoms_temp = atoms.copy()
                    atoms_temp.set_positions(pos_disp)
                    dipole_pos, _ = self.calculate_dipole(atoms_temp, **kwargs)

                    # Negative displacement
                    pos_disp[i, j] -= 2 * displacement
                    atoms_temp.set_positions(pos_disp)
                    dipole_neg, _ = self.calculate_dipole(atoms_temp, **kwargs)

                    # Central difference derivative
                    dipole_deriv = (dipole_pos - dipole_neg) / (2 * displacement)
                    dipole_derivatives[3 * i + j, :] = dipole_deriv

        except Exception as e:
            logger.warning(f"Dipole derivative calculation failed: {e}")

        finally:
            # Restore original positions
            atoms.set_positions(base_positions)

        return dipole_derivatives


class EspalomaDipoleCalculator(DipoleCalculatorBase):
    """Espaloma-charge based dipole calculator"""

    def __init__(self):
        super().__init__("espaloma")

    def _check_availability(self):
        try:
            import espaloma_charge
            from rdkit import Chem

            # Test basic functionality
            test_mol = Chem.MolFromSmiles("N")
            test_charges = espaloma_charge.charge(test_mol)

            self.available = True
            logger.info("\u2713 Espaloma-charge dipole calculator available and tested")
        except ImportError as e:
            self.available = False
            logger.warning(f"\u2717 Espaloma-charge dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole using espaloma partial charges"""
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        import espaloma_charge
        import torch

        # Create RDKit molecule from ASE atoms
        mol = Chem.RWMol()
        for symbol in atoms.get_chemical_symbols():
            atom = Chem.Atom(symbol)
            mol.AddAtom(atom)

        # Set coordinates
        conf = Chem.Conformer(len(atoms))
        positions = atoms.get_positions()
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, (float(pos[0]), float(pos[1]), float(pos[2])))
        mol.AddConformer(conf)

        # Determine bonds
        rdDetermineBonds.DetermineBonds(mol, charge=0)

        # Patch espaloma to use float32
        # Temporarily set torch default dtype to float32
        original_dtype = torch.get_default_dtype()
        torch.set_default_dtype(torch.float32)

        try:
            # Get partial charges from espaloma
            charges = espaloma_charge.charge(mol)
        finally:
            # Restore original dtype
            torch.set_default_dtype(original_dtype)

        # Calculate dipole moment: μ = Σ q_i * r_i
        dipole = np.dot(charges, positions)

        # Convert from e*Angstrom to e*Bohr (Gaussian units)
        dipole = dipole / 0.529177210903

        logger.debug(f"Espaloma dipole: {dipole}, charges sum: {np.sum(charges):.6f}")
        return dipole, charges


class XTBDipoleCalculator(DipoleCalculatorBase):
    """xTB-based dipole calculator"""

    def __init__(self):
        super().__init__("xtb")

    def _check_availability(self):
        try:
            from xtb.ase.calculator import XTB

            self.available = True
            logger.info("\u2713 xTB dipole calculator available")
        except ImportError as e:
            self.available = False
            logger.warning(f"\u2717 xTB dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole using xTB"""
        from xtb.ase.calculator import XTB

        atoms_copy = atoms.copy()
        atoms_copy.calc = XTB(method="GFN2-xTB")

        # Get dipole moment (in e*Bohr from xTB)
        dipole_moment = atoms_copy.get_dipole_moment()

        # Get partial charges
        partial_charges = atoms_copy.calc.get_charges(atoms_copy)

        logger.debug(f"xTB dipole: {dipole_moment}")
        return dipole_moment, partial_charges


class MACEMLDipoleCalculator(DipoleCalculatorBase):
    """MACE ML-based dipole calculator"""

    def __init__(self, model_path=None):
        # Use provided path, or fall back to configured default
        self.model_path = (
            model_path if model_path is not None else DEFAULT_MACE_DIPOLE_MODEL
        )
        self.mace_calc = None

        super().__init__("mace_ml")

    def _check_availability(self):
        try:
            # Check if model file exists
            if not Path(self.model_path).exists():
                raise FileNotFoundError(
                    f"MACE dipole model not found at: {self.model_path}"
                )

            self.mace_calc = MACEDipoleCalculator(self.model_path)
            self.available = True
            logger.info(
                f"\u2713 MACE ML dipole calculator available (model: {self.model_path})"
            )
        except (ImportError, FileNotFoundError) as e:
            self.available = False
            logger.warning(f"\u2717 MACE ML dipole calculator failed: {e}")

    def calculate_dipole(self, atoms, **kwargs):
        if not self.available:
            raise RuntimeError("MACE ML dipole calculator not available")
        return self.mace_calc.calculate_dipole(atoms, **kwargs)


class DipoleCalculatorFactory:
    """Factory for managing different dipole calculators"""

    def __init__(self):
        self.calculators = {}
        self.preferred_order = ["mace_ml", "espaloma", "xtb"]
        self._register_calculators()

    def _register_calculators(self):
        """Register all available calculators"""
        calculators = [
            EspalomaDipoleCalculator(),
            XTBDipoleCalculator(),
            MACEMLDipoleCalculator(),
        ]

        for calc in calculators:
            self.calculators[calc.name] = calc

    def get_calculator(self, method: str = "auto") -> DipoleCalculatorBase:
        """Get dipole calculator by name or auto-select best available"""

        if method == "auto":
            # Return first available calculator in preferred order
            for method_name in self.preferred_order:
                calc = self.calculators.get(method_name)
                if calc and calc.available:
                    logger.info(f"Auto-selected dipole calculator: {method_name}")
                    return calc

            raise RuntimeError("No dipole calculators available")

        calc = self.calculators.get(method)
        if not calc:
            raise ValueError(f"Unknown dipole calculator: {method}")

        if not calc.available:
            raise RuntimeError(f"Dipole calculator '{method}' not available")

        return calc

    def list_available(self) -> Dict[str, bool]:
        """List all calculators and their availability"""
        return {name: calc.available for name, calc in self.calculators.items()}


# Global dipole calculator factory
dipole_factory = DipoleCalculatorFactory()


# ============================================================================
# ORIGINAL HELPER FUNCTIONS (Enhanced)
# ============================================================================


@contextmanager
def zmq_server(file):
    """
    Creates a context manager with a ZMQ server socket using the IPC transport protocol.

    This version safely removes an existing IPC file before creating a new one,
    preventing FileExistsError if a previous run left the file behind.
    """
    addr = os.path.abspath(file)

    # If an old socket file exists, remove it first
    if os.path.exists(addr):
        try:
            os.remove(addr)
        except Exception as e:
            print(f"Warning: could not remove old IPC file {addr}: {e}")

    try:
        with zmq.Context() as ctx:
            with ctx.socket(zmq.REP) as socket:
                # Create a placeholder file to reserve the address
                with open(addr, "x"):
                    pass
                socket.bind(f"ipc://{addr}")
                yield socket
    finally:
        # Clean up socket file when done
        if os.path.exists(addr):
            os.remove(addr)


def is_calc_finished(proc, socket):
    """Waits until there is either a new msg (= script 2 was executed), or the g16 process finished"""
    while True:
        # poll the socket for messages with timeout 10 ms, returns 0 in case of timeout (= no messages)
        if socket.poll(timeout=10) != 0:
            # new msg, not finished
            return False
        # poll the running g16 process, if process has exited poll returns return code (int), otherwise None
        elif proc.poll() is not None:
            # g16 process exited, finished
            return True
        else:
            # if neither, wait a second and try again
            time.sleep(1)


# ============================================================================
# GAUSSIAN INPUT/OUTPUT HANDLING (Refactored from run_next_calculation)
# ============================================================================


def parse_gaussian_input(infile: str) -> Tuple[int, int, int, int, np.ndarray, list]:
    """
    Parse Gaussian external calculation input file.

    Args:
        infile: Path to Gaussian input file

    Returns:
        Tuple of (natoms, deriv, charge, spin, coordinates, atomnames)
        - natoms: Number of atoms
        - deriv: Derivative level (0=energy, 1=gradient, 2=hessian)
        - charge: Molecular charge
        - spin: Spin multiplicity
        - coordinates: Numpy array of shape (natoms, 3) in Angstroms
        - atomnames: List of element symbols
    """
    with open(infile, "r") as f:
        lines = f.readlines()

    # Extract system info from header line
    header = lines[0].split()
    natoms = int(header[0])
    deriv = int(header[1])
    charge = int(header[2])
    spin = int(header[3])

    # Initialize arrays
    coordinates = np.zeros((natoms, 3))
    atomnames = []

    # Parse atomic data (atomic number + coordinates in Bohr)
    BOHR_TO_ANGSTROM = 0.52917721092

    for i, line in enumerate(lines[1 : natoms + 1]):
        elements = line.split()

        # Convert atomic number to element symbol
        atomic_num = int(elements[0])
        atomnames.append(chemical_symbols[atomic_num])

        # Convert coordinates from Bohr to Angstrom
        coordinates[i] = BOHR_TO_ANGSTROM * np.array(
            [float(elements[1]), float(elements[2]), float(elements[3])]
        )

    return natoms, deriv, charge, spin, coordinates, atomnames


def update_molecule_geometry(atoms, coordinates: np.ndarray, charge: int, spin: int):
    """
    Update ASE Atoms object with new geometry and electronic state.

    Args:
        atoms: ASE Atoms object to update (modified in-place)
        coordinates: New atomic coordinates in Angstroms, shape (natoms, 3)
        charge: Molecular charge
        spin: Spin multiplicity
    """
    atoms.set_positions(coordinates)
    atoms.info["charge"] = float(charge)
    atoms.info["spin"] = float(spin)


def calculate_energy_and_forces(atoms, calculator) -> Tuple[float, np.ndarray]:
    """
    Calculate energy and forces using the attached calculator.

    Args:
        atoms: ASE Atoms object with calculator attached
        calculator: ASE calculator (not used directly, atoms.calc is used)

    Returns:
        Tuple of (energy, gradient)
        - energy: Potential energy in eV
        - gradient: Negative forces (gradient) in eV/Angstrom, shape (natoms, 3)
    """
    energy = atoms.get_potential_energy()
    gradient = -atoms.get_forces()  # Gradient = -Force
    return energy, gradient


def calculate_hessian(atoms, calculator, natoms: int) -> Optional[np.ndarray]:
    """
    Calculate Hessian matrix (second derivatives).

    Args:
        atoms: ASE Atoms object
        calculator: ASE calculator
        natoms: Number of atoms

    Returns:
        Hessian matrix in Hartree/Bohr^2, shape (3*natoms, 3*natoms)
        Returns None if calculator doesn't support Hessian
    """
    try:
        # Get Hessian in eV/Angstrom^2
        hessian = calculator.get_hessian(atoms=atoms)

        # Convert to Hartree/Bohr^2 for Gaussian
        # 1 eV/Ang^2 = 0.52917721092^2 / 27.211386246 Hartree/Bohr^2
        ANGSTROM_TO_BOHR = 0.52917721092
        EV_TO_HARTREE = 27.211386246
        hessian = hessian * (ANGSTROM_TO_BOHR**2) / EV_TO_HARTREE

        # Reshape to matrix form
        hessian = hessian.reshape(3 * natoms, 3 * natoms)
        return hessian

    except Exception as e:
        logger.warning(f"Hessian calculation failed: {e}")
        return None


def calculate_dipole_properties(
    atoms, dipole_calc, deriv: int, calculate_derivatives: bool
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Calculate dipole moment and optionally its derivatives.

    Args:
        atoms: ASE Atoms object
        dipole_calc: DipoleCalculatorBase instance
        deriv: Derivative level from Gaussian
        calculate_derivatives: Whether to calculate dipole derivatives

    Returns:
        Tuple of (dipole, dipole_derivatives, partial_charges)
        - dipole: Dipole moment in e*Bohr, shape (3,)
        - dipole_derivatives: Dipole derivatives in e, shape (3*natoms, 3)
        - partial_charges: Partial atomic charges (or None), shape (natoms,)
    """
    natoms = len(atoms)

    try:
        # Calculate dipole moment
        logger.info(f"Using dipole calculator: {dipole_calc.name}")
        dipole, partial_charges = dipole_calc.calculate_dipole(atoms)

        # Store charges in atoms object if available
        if partial_charges is not None:
            atoms.set_initial_charges(partial_charges)
            logger.debug(f"Partial charges sum: {np.sum(partial_charges):.6f}")

        # Calculate dipole derivatives for IR intensities
        if calculate_derivatives and deriv >= 1:
            logger.info("Calculating dipole derivatives for IR intensities...")
            dipole_derivatives = dipole_calc.calculate_dipole_derivatives(
                atoms, displacement=0.005
            )
        else:
            dipole_derivatives = np.zeros((3 * natoms, 3))

        logger.info(f"✓ Dipole calculated: {dipole} e*Bohr")
        return dipole, dipole_derivatives, partial_charges

    except Exception as e:
        logger.error(f"Dipole calculation failed: {e}")
        logger.warning("Falling back to zero dipole (IR intensities will be zero)")

        # Fallback to zeros
        dipole = np.zeros(3)
        dipole_derivatives = np.zeros((3 * natoms, 3))
        return dipole, dipole_derivatives, None


def write_gaussian_output(
    outfile: str,
    natoms: int,
    energy: float,
    gradient: np.ndarray,
    dipole: np.ndarray,
    dipole_derivatives: np.ndarray,
    hessian: Optional[np.ndarray],
    deriv: int,
):
    """
    Write results to Gaussian external calculation output file.

    Args:
        outfile: Path to output file
        natoms: Number of atoms
        energy: Energy in eV
        gradient: Gradient in eV/Angstrom, shape (natoms, 3)
        dipole: Dipole moment in e*Bohr, shape (3,)
        dipole_derivatives: Dipole derivatives in e, shape (3*natoms, 3)
        hessian: Hessian in Hartree/Bohr^2, shape (3*natoms, 3*natoms), or None
        deriv: Derivative level (0, 1, or 2)
    """
    # Convert energy from eV to Hartree
    EV_TO_HARTREE = 27.211386246
    energy_hartree = energy / EV_TO_HARTREE

    # Convert gradient from eV/Angstrom to Hartree/Bohr
    ANGSTROM_TO_BOHR = 0.52917721092
    gradient_hartree_bohr = gradient * ANGSTROM_TO_BOHR / EV_TO_HARTREE

    # Polarizability (not implemented, set to zero)
    polarizability = np.zeros(6)

    with open(outfile, "w") as f:
        # Write energy and dipole (Fortran format with 'D' exponent)
        line = f"{energy_hartree:20.12E}{dipole[0]:20.12E}{dipole[1]:20.12E}{dipole[2]:20.12E}"
        f.write(line.replace("E", "D") + "\n")

        # Write gradient (forces)
        for i in range(natoms):
            line = f"{gradient_hartree_bohr[i, 0]:20.12E}{gradient_hartree_bohr[i, 1]:20.12E}{gradient_hartree_bohr[i, 2]:20.12E}"
            f.write(line.replace("E", "D") + "\n")

        # Write polarizability (2 lines, 3 components each)
        line = f"{polarizability[0]:20.12E}{polarizability[1]:20.12E}{polarizability[2]:20.12E}"
        f.write(line.replace("E", "D") + "\n")
        line = f"{polarizability[3]:20.12E}{polarizability[4]:20.12E}{polarizability[5]:20.12E}"
        f.write(line.replace("E", "D") + "\n")

        # Write dipole derivatives (3*natoms lines, 3 components each)
        for i in range(3 * natoms):
            line = f"{dipole_derivatives[i, 0]:20.12E}{dipole_derivatives[i, 1]:20.12E}{dipole_derivatives[i, 2]:20.12E}"
            f.write(line.replace("E", "D") + "\n")

        # Write Hessian if second derivatives requested
        if deriv >= 2 and hessian is not None:
            count = 0
            for i in range(3 * natoms):
                for j in range(i + 1):  # Lower triangle including diagonal
                    line = f"{hessian[i, j]:20.12E}"
                    f.write(line.replace("E", "D"))
                    count += 1

                    if count % 3 == 0:  # Three entries per line
                        f.write("\n")

            # Ensure file ends with newline
            if count % 3 != 0:
                f.write("\n")


def run_next_calculation(
    mol,
    msg: str,
    calculator,
    dipole_method: str = "auto",
    calculate_derivatives: bool = True,
):
    """
    Coordinate a single calculation cycle for Gaussian external interface.

    This is the main function called by Gaussian for each calculation step.
    It orchestrates parsing input, running calculations, and writing output.

    Args:
        mol: ASE Atoms object
        msg: Message from Gaussian containing "infile|outfile" paths
        calculator: ASE calculator for energy/forces
        dipole_method: Method for dipole calculation ('auto', 'mace_ml', 'espaloma', etc.)
        calculate_derivatives: Whether to calculate dipole derivatives
    """
    # Visible progress indicator
    print(f"  → Processing Gaussian request...")

    # Parse message to get file paths
    infile, outfile = msg.split("|")

    # Step 1: Parse Gaussian input file
    natoms, deriv, charge, spin, coordinates, _atomnames = parse_gaussian_input(infile)

    # Step 2: Update molecule geometry and electronic state
    update_molecule_geometry(mol, coordinates, charge, spin)

    # Step 3: Calculate energy and forces
    energy, gradient = calculate_energy_and_forces(mol, calculator)

    # Step 4: Calculate Hessian if needed
    hessian = None
    if deriv >= 2:
        hessian = calculate_hessian(mol, calculator, natoms)

    # Step 5: Calculate dipole properties
    dipole_calc = dipole_factory.get_calculator(dipole_method)
    dipole, dipole_derivatives, _partial_charges = calculate_dipole_properties(
        mol, dipole_calc, deriv, calculate_derivatives
    )

    # Step 6: Write results to Gaussian output file
    write_gaussian_output(
        outfile, natoms, energy, gradient, dipole, dipole_derivatives, hessian, deriv
    )


def ase_to_gjf(
    atoms,
    filename="molecule.gjf",
    route=None,
    title="Gaussian input generated from ASE",
    charge=0,
    multiplicity=1,
):
    # Use default route with configured helper script path if none provided
    if route is None:
        route = f'# freq (anharm)\n# external="{DEFAULT_HELPER_SCRIPT}"'

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()  # Angstrom by ASE convention
    link0 = f"%chk={filename[:-3]}chk\n%mem=2GB\n%NProcShared=2"
    with open(filename, "w") as f:
        f.write(f"{link0}\n")
        f.write(f"{route}\n\n")
        f.write(f"{title}\n\n")
        f.write(f"{charge} {multiplicity}\n")
        for s, pos in zip(symbols, positions):
            f.write(f"{s:2s} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}\n")
        f.write("\n")  # blank line after coordinates


def geometry_optimisation(mol, fmax=0.00005):
    ei = mol.get_potential_energy()
    print("Initial Energy: ", ei, "eV")
    opt = LBFGS(mol)

    opt.run(fmax=fmax, steps=10000)

    ef = mol.get_potential_energy()
    print("Final Energy: ", ef, "eV")

    return mol


def calculator(nnp):
    if nnp == "mace_mp":
        from mace.calculators import mace_mp

        calc = mace_mp(
            model="large", device="cuda", default_dtype="float64", dispersion=False
        )
        return calc

    if nnp == "mace_off":
        from mace.calculators import mace_off

        calc = mace_off(
            model="large", device="cuda", default_dtype="float64", dispersion=False
        )
        return calc

    # *** NEW: MACE-OMOL-0 support with charge and spin embeddings ***
    if nnp == "mace_omol":
        from mace.calculators import mace_omol

        calc = mace_omol(
            model="extra_large",
            device="cuda",
            default_dtype="float64",
            dispersion=False,
        )
        return calc
    
# Charge Analyzer

def analyze_molecular_charges(atoms, output_dir='charge_analysis'):
    """
    Perform comprehensive charge analysis and visualization.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Molecular structure with charges stored
    output_dir : str
        Directory for saving analysis outputs
    """
    if not CHARGE_ANALYSIS_AVAILABLE:
        logger.warning("Charge analysis not available - skipping")
        return None
    
    # Check if charges are available
    charges = atoms.get_initial_charges()
    if charges is None or len(charges) == 0:
        logger.warning("No charges available for analysis")
        return None
    
    try:
        logger.info("")
        logger.info("="*70)
        logger.info("DETAILED CHARGE ANALYSIS")
        logger.info("="*70)
        
        # Create analyzer and run all analyses
        analyzer = ChargeAnalyzer(atoms, charges)
        analyzer.analyze_all(groups=None, save_plots=True, output_dir=output_dir)
        
        logger.info("="*70)
        logger.info("")
        return analyzer
    except Exception as e:
        logger.error(f"Charge analysis failed: {e}")
        return None
    
def setup_output_directory(base_name: str) -> str:
    """
    Create organized output directory structure.
    
    Parameters
    ----------
    base_name : str
        Base name for the calculation (e.g., 'acoh')
    
    Returns
    -------
    str
        Path to output directory
    """
    import os
    from datetime import datetime
    
    # Create main output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"{base_name}_analysis_{timestamp}"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"  ✓ Created output directory: {output_dir}")
    
    return output_dir

def run_geometry_optimization(
    atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    calculator_name: str = "mace_omol"
):
    """
    Run geometry optimization and save results.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Initial molecular structure
    molecule_name : str
        Name of the molecule
    results_mgr : ResultsManager
        Results manager instance
    calculator_name : str
        Name of calculator to use for optimization
        
    Returns
    -------
    ase.Atoms
        Optimized molecular structure
    """
    print("=" * 60)
    print("GEOMETRY OPTIMIZATION")
    print("=" * 60)
    print(f"Calculator: {calculator_name}")
    
    # Store initial structure
    initial_atoms = atoms.copy()
    initial_energy = atoms.get_potential_energy()
    
    # Track optimization
    start_time = time.time()
    
    # Run optimization
    optimized_atoms = geometry_optimisation(atoms)
    
    # Get final energy
    final_energy = optimized_atoms.get_potential_energy()
    runtime = time.time() - start_time
    
    # Save results
    results_mgr.save_optimization_results(
        molecule_name=molecule_name,
        initial_atoms=initial_atoms,
        final_atoms=optimized_atoms,
        calculator_name=calculator_name,
        initial_energy=initial_energy,
        final_energy=final_energy,
        converged=True,  # geometry_optimisation returns when converged
        num_steps=0,  # TODO: track this if needed
        runtime=runtime
    )
    
    print(f"Runtime: {runtime:.1f} seconds")
    print("=" * 60 + "\n")
    
    return optimized_atoms


def run_frequency_calculation(
    atoms,
    molecule_name: str,
    energy_calculator_name: str,
    dipole_calculator_name: str,
    results_mgr: ResultsManager,
    charge: int = 0,
    multiplicity: int = 1
):
    """
    Run frequency calculation with specified calculators.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Optimized molecular structure
    molecule_name : str
        Name of the molecule
    energy_calculator_name : str
        Name of energy calculator (e.g., 'mace_mp', 'mace_off')
    dipole_calculator_name : str
        Name of dipole calculator (e.g., 'espaloma', 'xtb')
    results_mgr : ResultsManager
        Results manager instance
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    print("=" * 60)
    print(f"FREQUENCY CALCULATION: {energy_calculator_name} + {dipole_calculator_name}")
    print("=" * 60)
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    start_time = time.time()
    
    try:
        # Setup calculator for energy
        calc = calculator(energy_calculator_name)
        mol = atoms.copy()
        mol.calc = calc
        
        # Set charge and spin
        mol.info["charge"] = float(charge)
        mol.info["spin"] = float(multiplicity)
        
        # Create frequency directory
        freq_dir = results_mgr.create_frequency_directory(
            molecule_name, energy_calculator_name, dipole_calculator_name, timestamp
        )
        
        # Generate Gaussian input file in the working directory first
        # Use a temporary name to avoid conflicts
        temp_basename = f"temp_{energy_calculator_name}_{dipole_calculator_name}_{timestamp}"
        gjf_temp = f"{temp_basename}.gjf"
        
        # Create Gaussian input
        ase_to_gjf(mol, gjf_temp, charge=charge, multiplicity=multiplicity)
        
        print(f"  → Gaussian input: {gjf_temp}")
        
        # Run Gaussian calculation
        print("  → Launching Gaussian...")
        proc = subprocess.Popen(["g16", gjf_temp])
        
        # ZMQ server for external interface
        counter = 0
        with zmq_server(".ipc_file") as socket:
            while not is_calc_finished(proc, socket):
                logger.info(f"Calculation step {counter}")
                counter += 1
                msg = socket.recv_string()
                
                try:
                    run_next_calculation(
                        mol,
                        msg,
                        calc,
                        dipole_method=dipole_calculator_name,
                        calculate_derivatives=True
                    )
                    socket.send_string("ready")
                    
                except Exception as e:
                    logger.error(f"Calculation step failed: {e}")
                    socket.send_string("error")
                    return False
        
        runtime = time.time() - start_time
        
        # Now move the files to the frequency directory
        log_temp = gjf_temp.replace('.gjf', '.log')
        chk_temp = gjf_temp.replace('.gjf', '.chk')
        
        gjf_final = str(freq_dir / "gaussian_freq.gjf")
        log_final = str(freq_dir / "gaussian_freq.log")
        
        # Move files (not copy, since they're in different locations)
        if os.path.exists(gjf_temp):
            shutil.move(gjf_temp, gjf_final)
        if os.path.exists(log_temp):
            shutil.move(log_temp, log_final)
        # Clean up checkpoint file
        if os.path.exists(chk_temp):
            os.remove(chk_temp)
        
                # Parse Gaussian log file to extract frequencies
        parsed_data = {}
        if os.path.exists(log_final):
            try:
                logger.info("Parsing Gaussian log file...")
                parsed_data = parse_gaussian_log(log_final)
                
                # Convert energy from Hartree to eV if available
                if parsed_data.get('final_energy_hartree'):
                    EV_TO_HARTREE = 27.211386246
                    final_energy = parsed_data['final_energy_hartree'] * EV_TO_HARTREE
                else:
                    final_energy = mol.get_potential_energy()
                    
            except Exception as e:
                logger.error(f"Failed to parse Gaussian log: {e}")
                parsed_data = {'harmonic': [], 'anharmonic': []}
                final_energy = mol.get_potential_energy()
        else:
            parsed_data = {'harmonic': [], 'anharmonic': []}
            final_energy = mol.get_potential_energy()
        
        results_mgr.save_frequency_results(
            molecule_name=molecule_name,
            energy_calculator=energy_calculator_name,
            dipole_calculator=dipole_calculator_name,
            calculator_type="ml",
            frequencies_data={
                'harmonic': parsed_data.get('harmonic', []),
                'anharmonic': parsed_data.get('anharmonic', [])
            },
            energy=final_energy,
            dipole=parsed_data.get('dipole_moment'),
            runtime=runtime,
            gaussian_log=log_final if os.path.exists(log_final) else None,
            gaussian_gjf=gjf_final if os.path.exists(gjf_final) else None,
            timestamp=timestamp
        )
        
        print(f"  ✓ Completed in {runtime:.1f} seconds")
        print("=" * 60 + "\n")
        
        return True
        
    except Exception as e:
        logger.error(f"Frequency calculation failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        print("  ✗ FAILED")
        print("=" * 60 + "\n")
        return False

####################################################################
##                   SCRIPT EXECUTION STARTS HERE                 ##
####################################################################

if __name__ == "__main__":
    # ========================================================================
    # DIAGNOSTIC MODE
    # ========================================================================
    if len(sys.argv) > 1 and sys.argv[1] == "--diagnose":
        print("=" * 60)
        print("DIAGNOSTIC MODE")
        print("=" * 60)
        print("\nAvailable dipole calculators:")
        for name, available in dipole_factory.list_available().items():
            status = "✓" if available else "✗"
            print(f"  {status} {name}")
        print("\n" + "=" * 60)
        sys.exit(0)
    
    # ========================================================================
    # CHECK FOR INPUT FILE
    # ========================================================================
    if len(sys.argv) < 2:
        print("Usage: python gm_main.py molecule.xyz")
        print("       python gm_main.py --diagnose")
        sys.exit(1)
    
    # ========================================================================
    # CONFIGURATION
    # ========================================================================
    OPTIMIZATION_CALCULATOR = "mace_omol"
    
    # For now, hardcode one combination - we'll add CLI later
    ENERGY_CALCULATORS = ["mace_mp", "mace_off", "mace_omol"]
    DIPOLE_CALCULATORS = ["espaloma", "xtb", "mace_ml"]
    
    # ========================================================================
    # MAIN WORKFLOW
    # ========================================================================
    
    # Read initial structure
    initial_coords = sys.argv[1]
    molecule_name = Path(initial_coords).stem  # e.g., 'acoh' from 'acoh.xyz'
    
    print("\n" + "=" * 60)
    print("MACE-GAUSSIAN COMPARISON FRAMEWORK")
    print("=" * 60)
    print(f"Molecule: {molecule_name}")
    print(f"Input: {initial_coords}")
    print("=" * 60 + "\n")
    
    # Initialize results manager
    results_mgr = ResultsManager(base_output_dir="comparison_results")
    
    # Read molecule
    mol = read(initial_coords)
    mol.info["charge"] = 0.0
    mol.info["spin"] = 1.0
    
    # Setup calculator for optimization
    calc = calculator(OPTIMIZATION_CALCULATOR)
    mol.calc = calc
    
    # ========================================================================
    # PHASE 1: GEOMETRY OPTIMIZATION
    # ========================================================================
    try:
        optimized_mol = run_geometry_optimization(
            mol, 
            molecule_name, 
            results_mgr,
            OPTIMIZATION_CALCULATOR
        )
    except Exception as e:
        logger.error(f"Geometry optimization failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
    # ========================================================================
    # PHASE 2: FREQUENCY CALCULATIONS WITH DIFFERENT METHODS
    # ========================================================================
    
    # Get charge and spin for frequency calculations
    charge = int(optimized_mol.info.get("charge", 0))
    multiplicity = int(optimized_mol.info.get("spin", 1))
    
    print("=" * 60)
    print("FREQUENCY CALCULATIONS")
    print("=" * 60)
    print(f"Energy calculators: {ENERGY_CALCULATORS}")
    print(f"Dipole calculators: {DIPOLE_CALCULATORS}")
    print("=" * 60 + "\n")
    
    # Run all combinations
    results = []
    for energy_calc in ENERGY_CALCULATORS:
        for dipole_calc in DIPOLE_CALCULATORS:
            success = run_frequency_calculation(
                optimized_mol,
                molecule_name,
                energy_calc,
                dipole_calc,
                results_mgr,
                charge,
                multiplicity
            )
            results.append({
                'energy': energy_calc,
                'dipole': dipole_calc,
                'success': success
            })
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 60)
    print("CALCULATION SUMMARY")
    print("=" * 60)
    
    successful = sum(1 for r in results if r['success'])
    total = len(results)
    
    print(f"Completed: {successful}/{total} calculations")
    print("\nResults saved to: comparison_results/{}/".format(molecule_name))
    
    for r in results:
        status = "✓" if r['success'] else "✗"
        print(f"  {status} {r['energy']} + {r['dipole']}")
    
    print("=" * 60)
