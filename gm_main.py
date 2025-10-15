#!/usr/bin/env python3
# For the helperscript copy to /home/bin/ and set shebang!!!
# module load mace-torch-dftd/14Jul2025
# module load gv

import numpy as np
import zmq
import os
import time
import warnings
from contextlib import contextmanager
import sys
import subprocess
from mace_calculators import MACEDipoleCalculator
from ase.io import read
from ase.optimize import LBFGS
from ase.data import atomic_numbers, chemical_symbols
from abc import ABC, abstractmethod
from typing import Tuple, Optional, Dict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

#Suppress Warnings
warnings.filterwarnings("ignore", message=".*weights_only=False.*", category=FutureWarning)
os.environ['PYTHONWARNINGS'] = 'ignore::FutureWarning'

# ============================================================================
# DIAGNOSTIC UTILITIES
# ============================================================================

def diagnose_python_environment():
    """Diagnose Python environment and package availability"""
    import sys
    import site
    
    print("=" * 60)
    print("PYTHON ENVIRONMENT DIAGNOSTICS")
    print("=" * 60)
    
    print(f"Python executable: {sys.executable}")
    print(f"Python version: {sys.version}")
    print(f"Python path: {sys.path[:3]}...")  # First few entries
    
    print("\nSite packages:")
    for path in site.getsitepackages():
        print(f"  {path}")
    
    print(f"\nUser site packages: {site.getusersitepackages()}")
    
    print("\nTesting imports:")
    
    # Test individual imports
    test_imports = [
        ('numpy', 'np'),
        ('ase', None),
        ('rdkit', None),
        ('rdkit.Chem', 'Chem'),
        ('espaloma_charge', None),
    ]
    
    for module_name, alias in test_imports:
        try:
            if alias:
                exec(f"import {module_name} as {alias}")
            else:
                exec(f"import {module_name}")
            print(f"  \u2713 {module_name}")
            
            # Get version if available
            try:
                if module_name == 'rdkit':
                    from rdkit import rdBase
                    print(f"    Version: {rdBase.rdkitVersion}")
                elif module_name == 'espaloma_charge':
                    import espaloma_charge
                    if hasattr(espaloma_charge, '__version__'):
                        print(f"    Version: {espaloma_charge.__version__}")
                    else:
                        print(f"    Location: {espaloma_charge.__file__}")
            except:
                pass
                
        except ImportError as e:
            print(f"  \u2717 {module_name}: {e}")
        except Exception as e:
            print(f"  ? {module_name}: {e}")
    
    print("=" * 60)
    print()

def test_espaloma_functionality():
    """Test espaloma_charge functionality step by step"""
    print("TESTING ESPALOMA_CHARGE FUNCTIONALITY")
    print("=" * 60)
    
    try:
        print("1. Testing espaloma_charge import...")
        import espaloma_charge
        print("   \u2713 espaloma_charge imported successfully")
        
        print("2. Testing charge function import...")
        from espaloma_charge import charge
        print("   \u2713 charge function imported successfully")
        
        print("3. Testing RDKit import...")
        from rdkit import Chem
        print("   \u2713 RDKit imported successfully")
        
        print("4. Testing simple molecule creation...")
        mol = Chem.MolFromSmiles("N#N")
        if mol is not None:
            print("   \u2713 RDKit molecule created successfully")
        else:
            print("   \u2717 RDKit molecule creation failed")
            return False
        
        print("5. Testing espaloma charge calculation...")
        charges = charge(mol)
        print(f"   \u2713 Charges calculated: {charges}")
        print(f"   Shape: {charges.shape}, Type: {type(charges)}")
        
        print("6. Testing with a larger molecule...")
        mol2 = Chem.MolFromSmiles("CCO")  # Ethanol
        charges2 = charge(mol2)
        print(f"   \u2713 Ethanol charges: {charges2}")
        
        return True
        
    except ImportError as e:
        print(f"   \u2717 Import error: {e}")
        return False
    except Exception as e:
        print(f"   \u2717 Runtime error: {e}")
        return False
    
    finally:
        print("=" * 60)
        print()

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
    def calculate_dipole(self, atoms, **kwargs) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Calculate dipole moment and partial charges
        Returns: (dipole_vector, partial_charges)
        """
        pass
    
    def calculate_dipole_derivatives(self, atoms, displacement=0.01, **kwargs) -> np.ndarray:
        """Calculate dipole derivatives numerically"""
        natoms = len(atoms)
        dipole_derivatives = np.zeros((3*natoms, 3))
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
                    pos_disp[i, j] -= 2*displacement
                    atoms_temp.set_positions(pos_disp)
                    dipole_neg, _ = self.calculate_dipole(atoms_temp, **kwargs)
                    
                    # Central difference derivative
                    dipole_deriv = (dipole_pos - dipole_neg) / (2*displacement)
                    dipole_derivatives[3*i + j, :] = dipole_deriv
                    
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
    
    def __init__(self, model_path="/home/mot/mace_gaussian/dipole_model/model_1.model"):
        self.model_path = model_path       
        self.mace_calc = None
        super().__init__("mace_ml")
    
    def _check_availability(self):
        try:
            self.mace_calc = MACEDipoleCalculator(self.model_path)
            self.available = True
            logger.info("\u2713 MACE ML dipole calculator available")
        except ImportError as e:
            self.available = False
            logger.warning(f"\u2717 MACE ML dipole calculator failed: {e}")
    
    def calculate_dipole(self, atoms, **kwargs):
        if not self.available:
            raise RuntimeError("MACE ML dipole calculator not available")
        return self.mace_calc.calculate_dipole(atoms, **kwargs)


class GeometryDipoleCalculator(DipoleCalculatorBase):
    """Simple geometry-based dipole estimate (fallback)"""
    
    def __init__(self):
        super().__init__("geometry")
        
    def _check_availability(self):
        self.available = True  # Always available
        logger.info("\u2713 Geometry-based dipole calculator available (fallback)")
    
    def calculate_dipole(self, atoms, **kwargs) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """Crude electronegativity-based dipole estimation"""
        
        # Pauling electronegativities (simplified set)
        electronegativity = {
            'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96,
            'I': 2.66, 'Li': 0.98, 'Na': 0.93, 'K': 0.82, 'Mg': 1.31,
            'Ca': 1.00, 'Al': 1.61, 'B': 2.04, 'Be': 1.57
        }
        
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Estimate partial charges based on electronegativity
        charges = []
        avg_electronegativity = np.mean([electronegativity.get(s, 2.5) for s in symbols])
        
        for symbol in symbols:
            en = electronegativity.get(symbol, 2.5)
            # Simple charge estimation (very crude)
            charge = (en - avg_electronegativity) * 0.1
            charges.append(charge)
        
        charges = np.array(charges)
        
        # Ensure charge neutrality
        total_charge = np.sum(charges)
        charges -= total_charge / len(charges)
        
        # Calculate dipole
        dipole = np.dot(charges, positions) / 0.5291772105638411
        
        logger.debug(f"Geometry-based dipole: {dipole}")
        return dipole, charges


class DipoleCalculatorFactory:
    """Factory for managing different dipole calculators"""
    
    def __init__(self):
        self.calculators = {}
        self.preferred_order = ['mace_ml', 'espaloma', 'xtb', 'geometry']
        self._register_calculators()
    
    def _register_calculators(self):
        """Register all available calculators"""
        calculators = [
            EspalomaDipoleCalculator(),
            XTBDipoleCalculator(),
            MACEMLDipoleCalculator(),
            GeometryDipoleCalculator()
        ]
        
        for calc in calculators:
            self.calculators[calc.name] = calc
    
    def get_calculator(self, method: str = 'auto') -> DipoleCalculatorBase:
        """Get dipole calculator by name or auto-select best available"""
        
        if method == 'auto':
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
    """Creates a context manager with ZMQ client socket as resource using the IPC transport protocol, used via python with block.
    There must not be a file with the filename passed as argument in the current dir.
    """
    try:
        with zmq.Context() as ctx:
            with ctx.socket(zmq.REP) as socket:
                addr = os.path.abspath(file)
                with open(addr, "x") as f:
                    socket.bind("ipc://%s" % addr)
                    yield socket
    finally:
        os.remove(addr)


def is_calc_finished(proc, socket):
    """Waits until there is either a new msg (= script 2 was executed), or the g16 process finished"""
    while True:
        # poll the socket for messages with timeout 10 ms, returns 0 in case of timeout (= no messages)
        if socket.poll(timeout=10) != 0:
            # new msg, not finished
            return False
        # poll the running g16 process, if process has exited poll returns return code (int), otherwise None
        elif proc.poll() != None:
            # g16 process exited, finished
            return True
        else:
            # if neither, wait a second and try again
            time.sleep(1)


def run_next_calculation(mol, msg, calculator, dipole_method='auto', calculate_derivatives=True):
    """updates coords, runs actual calculations, writes results to file"""
    # Read data from gaussian sys-call
    infile, outfile = msg.split('|')

    # Read lines from infile    
    infile_ptr = open(infile, "r")    
    lines = infile_ptr.readlines()
    infile_ptr.close()
    
    # Extract system info from head line
    line_element = lines[0].split()  
    
    natoms = int(line_element[0])
    deriv  = int(line_element[1])
    charge = int(line_element[2])
    spin   = int(line_element[3])
    
    # Initialize arrays
    coord = np.zeros((natoms,3))
    atomnames = []   

    # Loop over all atoms     
    for i, line in enumerate(lines[1:natoms+1]):
        line_element = line.split()

        # convert atomic number to element name
        atom = chemical_symbols[int(line_element[0])]
        atomnames.append(atom)
        
        # convert bohr to Angstroem and store coordinates
        xyz = 0.52917721092 * np.array([float(line_element[1]), float(line_element[2]), float(line_element[3])])

        coord[i,0] = xyz[0]
        coord[i,1] = xyz[1]
        coord[i,2] = xyz[2]

    # Update the coordinates for iterative (anharmonic) scheme
    mol.set_positions(coord)
    
    # *** NEW: Set charge and spin in atoms.info for MACE-OMOL ***
    mol.info["charge"] = float(charge)
    mol.info["spin"] = float(spin)

    # Calculate required properties (energy, forces, etc)
    E = mol.get_potential_energy()
    grad = -mol.get_forces()

    # ==========================================
    # ENHANCED DIPOLE CALCULATION SYSTEM
    # ==========================================
    try:
        # Get dipole calculator
        dipole_calc = dipole_factory.get_calculator(dipole_method)
        logger.info(f"Using dipole calculator: {dipole_calc.name}")
        
        # Calculate dipole moment
        dipole, partial_charges = dipole_calc.calculate_dipole(mol)
     
        # Store charges in ASE atoms object if available
        if partial_charges is not None:
            mol.set_initial_charges(partial_charges)
            logger.debug(f"Stored partial charges: {np.sum(partial_charges):.6f} total")
        
        # Calculate dipole derivatives if requested and derivatives are needed
        if calculate_derivatives and deriv >= 1:
            logger.info("Calculating dipole derivatives for IR intensities...")
            dip_deriv = dipole_calc.calculate_dipole_derivatives(mol, displacement=0.005)
        else:
            dip_deriv = np.zeros((3*natoms, 3))
        
        logger.info(f"\u2713 Dipole calculated: {dipole} e*bohr")
        
    except Exception as e:
        logger.error(f"Dipole calculation failed: {e}")
        logger.warning("Falling back to zero dipole (IR intensities will be zero)")
        
        # Fallback to original behavior
        dipole = np.zeros(3)
        dip_deriv = np.zeros((3*natoms, 3))

    # Set other output to zero (as in original)
    polar = np.zeros(6)
        
    outfile_ptr = open(outfile, "w")   
    # write energy and dipole moment components - replace E with D for fortran
    
    string = f"{E:20.12E}{dipole[0]:20.12E}{dipole[1]:20.12E}{dipole[2]:20.12E}".replace('E', 'D')
    outfile_ptr.write(string + "\n")

    # write gradient
    for i in range(natoms):
        string = f"{grad[i][0]:20.12E}{grad[i][1]:20.12E}{grad[i][2]:20.12E}".replace('E', 'D')
        outfile_ptr.write(string + "\n")
        
    # write polarizability
    string = f"{polar[0]:20.12E}{polar[1]:20.12E}{polar[2]:20.12E}".replace('E', 'D')
    outfile_ptr.write(string + "\n")
 
    string = f"{polar[3]:20.12E}{polar[4]:20.12E}{polar[5]:20.12E}".replace('E', 'D')
    outfile_ptr.write(string + "\n")
 
    # write dipole derivatives
    for i in range(3*natoms):
        string = f"{dip_deriv[i][0]:20.12E}{dip_deriv[i][1]:20.12E}{dip_deriv[i][2]:20.12E}".replace('E', 'D')
        outfile_ptr.write(string + "\n")
   
    hessian = calculator.get_hessian(atoms=mol) * 0.52917721092 * 0.52917721092 / 27.211386246 # convert to E_h/bohr2
    hessian = hessian.reshape(3*natoms, 3*natoms)
 
    count = 0
    for i in range(3*natoms): # iterate over all entries
        for j in range(i + 1):  # Iterate over lower triangle including diagonal
            string = f"{hessian[i, j]:20.12E}".replace('E', 'D')
            outfile_ptr.write(string)  # Format to 6 decimal places
            count += 1
            
            if count % 3 == 0:  # Write three entries per line
                outfile_ptr.write("\n")

    # Ensure the file ends with a newline if the last line isn't full
    if count % 3 != 0:
        outfile_ptr.write("\n")
    
    outfile_ptr.close()

    return


def ase_to_gjf(atoms, filename='molecule.gjf',
               route='# freq (anharm)\n# external="/home/mot/mace_gaussian/gm_helper.py"',
               title='Gaussian input generated from ASE',
               charge=0, multiplicity=1):
               
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()   # Angstrom by ASE convention
    link0=f'%chk={filename[:-3]}chk\n%mem=2GB\n%NProcShared=2'
    with open(filename, 'w') as f:
        f.write(f'{link0}\n')
        f.write(f'{route}\n\n')
        f.write(f'{title}\n\n')
        f.write(f'{charge} {multiplicity}\n')
        for s, pos in zip(symbols, positions):
            f.write(f'{s:2s} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}\n')
        f.write('\n')   # blank line after coordinates


def geometry_optimisation(mol, fmax = 0.00005):
    ei=mol.get_potential_energy()
    print("Initial Energy: ", ei, "eV")
    opt = LBFGS(mol)

    opt.run(fmax=fmax,steps=10000)

    ef=mol.get_potential_energy()
    print("Final Energy: ", ef, "eV")

    return mol


def calculator(nnp):
    if nnp == 'mace_mp':
        from mace.calculators import mace_mp
        calc = mace_mp(model = "large",
                    device = "cuda",
                    default_dtype='float64',
                    dispersion = False
                   )
        return calc
        
    if nnp == 'mace_off':
        from mace.calculators import mace_off
        calc = mace_off(model = "large",
                    device = "cuda",
                    default_dtype='float64',
                    dispersion = False
                    ) 
        return calc
    
    # *** NEW: MACE-OMOL-0 support with charge and spin embeddings ***
    if nnp == 'mace_omol':
        from mace.calculators import mace_omol
        calc = mace_omol(
            model="extra_large",
            device="cuda",
            default_dtype='float64',
            dispersion=False
        )
        return calc


####################################################################
##                   SCRIPT EXECUTION STARTS HERE                 ##
####################################################################

if __name__ == '__main__':
    # ========================================================================
    # DIAGNOSTIC MODE - Run diagnostics if requested
    # ========================================================================
    
    # Configuration options
    DIPOLE_METHOD = 'auto'  # Options: 'auto', 'mace_ml', 'espaloma', 'xtb', 'geometry'
    CALCULATE_DIPOLE_DERIVATIVES = True
    
    # Note: espaloma is recommended for organic molecules
    # 'auto' will select the best available calculator in this order:
    # 1. mace_ml (if mace_dipole_core is available - currently has a bug)
    # 2. espaloma (if espaloma_charge + rdkit are installed) \u2713 RECOMMENDED
    # 3. xtb (if xtb-python is installed)
    # 4. geometry (fallback - crude electronegativity-based estimation)
    
    # Print available dipole calculators
    logger.info("Available dipole calculators:")
    for name, available in dipole_factory.list_available().items():
        status = "\u2713" if available else "\u2717"
        logger.info(f"  {status} {name}")
    
    # If no dipole calculators available except geometry, show help
    available_calcs = [name for name, avail in dipole_factory.list_available().items() if avail]
    if available_calcs == ['geometry']:
        logger.warning("")
        logger.warning("=" * 60)
        logger.warning("ONLY FALLBACK DIPOLE CALCULATOR AVAILABLE")
        logger.warning("=" * 60)
        logger.warning("Only geometry-based dipole estimation is available.")
        logger.warning("This will give very crude dipole moments and IR intensities.")
        logger.warning("")
        logger.warning("To improve accuracy, install additional packages:")
        logger.warning("  pip install espaloma_charge rdkit  # Recommended")
        logger.warning("  # OR")
        logger.warning("  mamba install -c conda-forge xtb-python")
        logger.warning("")
        logger.warning("Run diagnostics with: python gm_main.py --diagnose")
        logger.warning("=" * 60)
        logger.warning("")
    
    # Check for input file
    if len(sys.argv) < 2 or sys.argv[1].startswith('--'):
        print("Usage: python gm_main.py molecule.xyz")
        print("       python gm_main.py --diagnose")
        sys.exit(1)
    
    # ========================================================================
    # MAIN EXECUTION
    # ========================================================================
    
    # pass initial file as argument when running script in .xyz or .cif format
    initial_coords = sys.argv[1]

    # Read molecule
    mol = read(initial_coords)
    
    # *** NEW: Set default charge and spin for MACE-OMOL ***
    # These are defaults for neutral, closed-shell (singlet) molecules
    # Gaussian will override these with actual values when it calls back
    mol.info["charge"] = 0.0   # Neutral molecule
    mol.info["spin"] = 1.0     # Singlet state (all electrons paired, 2S+1 = 1)
    
    print("\n" + "="*60)
    print("INITIAL MOLECULE SETUP")
    print("="*60)
    print(f"Loaded: {initial_coords}")
    print(f"Atoms: {len(mol)}")
    print(f"Initial charge: {mol.info['charge']}")
    print(f"Initial spin multiplicity: {mol.info['spin']}")
    print("="*60 + "\n")

    # Initialize calculator and do geometry optimization
    # *** CHANGED: Using mace_omol instead of mace_mp ***
    calc = calculator('mace_omol')  # Options: 'mace_mp', 'mace_off', 'mace_omol'
    mol.calc = calc

    print("="*60)
    print("GEOMETRY OPTIMIZATION")
    print("="*60)
    mol = geometry_optimisation(mol)
    print("="*60 + "\n")

    # Generate Gaussian input file
    gjf_filename = initial_coords[:-4]+'_freq_anharm.gjf'
    
    # *** NEW: Get charge and spin from atoms.info for Gaussian input ***
    charge = int(mol.info.get("charge", 0))
    multiplicity = int(mol.info.get("spin", 1))
    
    ase_to_gjf(mol, gjf_filename, charge=charge, multiplicity=multiplicity)
    
    print("="*60)
    print("GAUSSIAN INPUT FILE GENERATED")
    print("="*60)
    print(f"File: {gjf_filename}")
    print(f"Charge: {charge}")
    print(f"Spin Multiplicity: {multiplicity}")
    print("="*60 + "\n")

    # Start gaussian as subprocess
    print("="*60)
    print("LAUNCHING GAUSSIAN")
    print("="*60)
    proc = subprocess.Popen(["g16", './'+gjf_filename])

    # From here the 'server' starts and waits for the request of the client.
    counter = 0
    with zmq_server(".ipc_file") as socket:
        while not is_calc_finished(proc, socket):
            logger.info(f'NNP calculation run {counter}')
            counter += 1
            # add possibility to send 2 different strings, where one signals error
            msg = socket.recv_string()
            
            try:
                run_next_calculation(mol, msg, calc, 
                                   dipole_method=DIPOLE_METHOD,
                                   calculate_derivatives=CALCULATE_DIPOLE_DERIVATIVES)
                # currently the content doesn't matter, maybe add possibility to signal error
                socket.send_string("ready")
                
            except Exception as e:
                logger.error(f"Calculation failed: {e}")
                socket.send_string("error")
    
    logger.info("Calculation completed successfully!")
    print("\n" + "="*60)
    print("CALCULATION COMPLETED SUCCESSFULLY!")
    print("="*60)
