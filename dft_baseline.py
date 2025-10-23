"""
DFT Baseline Calculator

Runs pure Gaussian DFT calculations (no external interface) as baselines
for comparison with ML-enhanced calculations.
"""

import subprocess
import time
import logging
import shutil
from pathlib import Path
from typing import Dict, Any, Optional
from ase import Atoms

from gaussian_parser import parse_gaussian_log
from results_manager import ResultsManager

logger = logging.getLogger(__name__)


# DFT methods that the NNPs were trained on
# NOTE: wb97mv requires Gaussian 16 Rev. C.01 or later
# Using wb97xd as alternative - similar range-separated hybrid with dispersion
DFT_BASELINES = {
    'mace_omol': {
        'method': 'wb97xd',  # Alternative to wb97mv for older Gaussian versions
        'basis': 'def2tzvp',
        'description': 'wB97X-D/def2-TZVP (MACE-OMOL alternative)',
        'extra_keywords': ''
    },
    'mace_off': {
        'method': 'wb97xd',
        'basis': 'def2tzvp',
        'description': 'wB97X-D/def2-TZVP (MACE-OFF reference)',
        'extra_keywords': ''
    },
    'mace_mp': {
        'method': 'pbepbe',
        'basis': 'def2svp',
        'description': 'PBE/def2-SVP (MACE-MP reference)',
        'extra_keywords': ''
    }
}

# If you have Gaussian 16 Rev. C.01 or later, you can use the original functional:
# Change mace_omol 'method' to 'wb97mv' for \u03c9B97M-V/def2-TZVP
ALTERNATIVE_BASELINES = {
    'mace_omol_original': {
        'method': 'wb97mv',  # Requires G16 Rev. C.01+
        'basis': 'def2tzvp',
        'description': 'wB97M-V/def2-TZVP (MACE-OMOL original - requires newer Gaussian)',
        'extra_keywords': ''
    },
}


def check_baseline_exists(
    results_mgr: ResultsManager,
    molecule_name: str,
    baseline_name: str
) -> bool:
    """
    Check if a DFT baseline calculation already exists.
    
    Parameters
    ----------
    results_mgr : ResultsManager
        Results manager instance
    molecule_name : str
        Name of the molecule
    baseline_name : str
        Name of the baseline (e.g., 'mace_omol', 'mace_off', 'mace_mp')
        
    Returns
    -------
    bool
        True if results exist and are valid, False otherwise
    """
    config = DFT_BASELINES[baseline_name]
    method = config['method']
    basis = config['basis']
    calculator_name = f"{method}_{basis}"
    
    # Check if directory exists
    mol_dir = results_mgr.create_molecule_directory(molecule_name)
    freq_dir = mol_dir / f"{calculator_name}_{calculator_name}"
    
    if not freq_dir.exists():
        return False
    
    # Check if results.json exists and is valid
    json_file = freq_dir / "results.json"
    if not json_file.exists():
        return False
    
    # Basic validation - check if log file exists
    log_file = freq_dir / "gaussian_freq.log"
    if not log_file.exists():
        return False
    
    logger.info(f"\u2713 DFT baseline {calculator_name} already exists for {molecule_name}")
    return True


def create_gaussian_dft_input(
    atoms: Atoms,
    filename: str,
    method: str,
    basis: str,
    charge: int = 0,
    multiplicity: int = 1,
    title: str = "DFT baseline calculation",
    extra_keywords: str = ""
) -> None:
    """
    Create Gaussian input file for pure DFT calculation.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Molecular structure
    filename : str
        Output filename for .gjf file
    method : str
        DFT method (e.g., 'wb97mv', 'pbepbe')
    basis : str
        Basis set (e.g., 'def2tzvp', 'def2svp')
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity
    title : str
        Title line for Gaussian input
    extra_keywords : str
        Additional keywords for route card
    """
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()  # Angstrom
    
    # Construct route card - pure DFT, no external interface
    route = f"# freq(anharm) {method}/{basis}"
    if extra_keywords:
        route += f" {extra_keywords}"
    
    # Link0 commands
    link0 = f"%chk={filename[:-4]}.chk\n%mem=4GB\n%NProcShared=4"
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"{link0}\n")
        f.write(f"{route}\n\n")
        f.write(f"{title}\n\n")
        f.write(f"{charge} {multiplicity}\n")
        for s, pos in zip(symbols, positions):
            f.write(f"{s:2s} {pos[0]:14.8f} {pos[1]:14.8f} {pos[2]:14.8f}\n")
        f.write("\n")


def run_gaussian_dft(
    gjf_file: str,
    timeout: Optional[int] = None
) -> tuple[bool, str]:
    """
    Run Gaussian calculation and wait for completion.
    
    Parameters
    ----------
    gjf_file : str
        Path to Gaussian input file
    timeout : int, optional
        Timeout in seconds (None = no timeout)
        
    Returns
    -------
    success : bool
        True if calculation completed successfully
    log_file : str
        Path to the log file
    """
    log_file = gjf_file.replace('.gjf', '.log')
    
    try:
        print(f"  \u2192 Launching Gaussian calculation...")
        
        # Run Gaussian (g16 executable)
        proc = subprocess.Popen(
            ['g16', gjf_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # Wait for completion
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            logger.error(f"Gaussian calculation timed out after {timeout}s")
            return False, log_file
        
        # Check return code
        if proc.returncode != 0:
            logger.error(f"Gaussian exited with error code {proc.returncode}")
            return False, log_file
        
        # Check if log file was created
        if not Path(log_file).exists():
            logger.error("Gaussian log file was not created")
            return False, log_file
        
        # Check for normal termination
        with open(log_file, 'r') as f:
            content = f.read()
            if 'Normal termination' not in content:
                logger.error("Gaussian calculation did not terminate normally")
                return False, log_file
        
        print(f"  \u2713 Gaussian calculation completed successfully")
        return True, log_file
        
    except Exception as e:
        logger.error(f"Error running Gaussian: {e}")
        return False, log_file


def run_dft_baseline_calculation(
    atoms: Atoms,
    molecule_name: str,
    baseline_name: str,
    results_mgr: ResultsManager,
    charge: int = 0,
    multiplicity: int = 1,
    skip_if_exists: bool = True
) -> bool:
    """
    Run a DFT baseline calculation.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Optimized molecular structure
    molecule_name : str
        Name of the molecule
    baseline_name : str
        Name of baseline ('mace_omol', 'mace_off', 'mace_mp')
    results_mgr : ResultsManager
        Results manager instance
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity
    skip_if_exists : bool
        If True, skip calculation if results already exist
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if baseline_name not in DFT_BASELINES:
        logger.error(f"Unknown baseline: {baseline_name}")
        return False
    
    config = DFT_BASELINES[baseline_name]
    method = config['method']
    basis = config['basis']
    description = config['description']
    extra_keywords = config.get('extra_keywords', '')
    calculator_name = f"{method}_{basis}"
    
    # Check if already exists
    if skip_if_exists and check_baseline_exists(results_mgr, molecule_name, baseline_name):
        print(f"  \u2192 Skipping {calculator_name} (already exists)")
        return True
    
    print("=" * 60)
    print(f"DFT BASELINE: {description}")
    print("=" * 60)
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    start_time = time.time()
    
    try:
        # Create frequency directory
        freq_dir = results_mgr.create_frequency_directory(
            molecule_name, calculator_name, calculator_name, timestamp
        )
        
        # Generate Gaussian input file
        gjf_file = f"dft_{baseline_name}_{timestamp}.gjf"
        create_gaussian_dft_input(
            atoms,
            gjf_file,
            method,
            basis,
            charge,
            multiplicity,
            title=f"DFT baseline: {description}",
            extra_keywords=extra_keywords
        )
        
        print(f"  \u2192 Created input: {gjf_file}")
        
        # Run Gaussian calculation
        success, log_file = run_gaussian_dft(gjf_file, timeout=7200)  # 2 hour timeout
        
        if not success:
            print(f"  \u2717 Calculation failed")
            print(f"  \u2192 Check log file: {log_file}")
            
            # Check if it's a functional recognition issue
            if Path(log_file).exists():
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    if 'Unrecognized' in log_content or 'Unknown' in log_content:
                        print(f"  \u26a0 Possible issue: Functional '{method}' may not be available in your Gaussian version")
                        print(f"    Try updating DFT_BASELINES in dft_baseline.py to use 'wb97xd' instead")
            
            return False
        
        runtime = time.time() - start_time
        
        # Parse results
        print(f"  \u2192 Parsing results...")
        try:
            parsed_data = parse_gaussian_log(log_file)
        except Exception as e:
            logger.error(f"Failed to parse log file: {e}")
            return False
        
        # Convert energy from Hartree to eV
        HARTREE_TO_EV = 27.211386246
        energy_hartree = parsed_data.get('final_energy_hartree')
        if energy_hartree is not None:
            energy_ev = energy_hartree * HARTREE_TO_EV
        else:
            logger.error("Could not extract energy from log file")
            return False
        
        # Move Gaussian files to output directory
        chk_file = gjf_file.replace('.gjf', '.chk')
        
        # Final paths in output directory
        final_gjf = freq_dir / f"gaussian_dft.gjf"
        final_log = freq_dir / f"gaussian_dft.log"
        final_chk = freq_dir / f"gaussian_dft.chk"
        
        # Move/copy files
        if Path(gjf_file).exists():
            shutil.move(gjf_file, final_gjf)
        if Path(log_file).exists():
            shutil.move(log_file, final_log)
        if Path(chk_file).exists():
            shutil.move(chk_file, final_chk)
        
        # Save results (point to the final file locations)
        results_mgr.save_frequency_results(
            molecule_name=molecule_name,
            energy_calculator=calculator_name,
            dipole_calculator=calculator_name,
            calculator_type="dft",
            frequencies_data={
                'harmonic': parsed_data.get('harmonic', []),
                'anharmonic': parsed_data.get('anharmonic', [])
            },
            energy=energy_ev,
            dipole=parsed_data.get('dipole_moment'),
            runtime=runtime,
            gaussian_log=str(final_log),
            gaussian_gjf=str(final_gjf),
            timestamp=timestamp
        )
        
        print(f"  \u2713 Completed in {runtime:.1f} seconds")
        print("=" * 60 + "\n")
        
        return True
        
    except Exception as e:
        logger.error(f"DFT baseline calculation failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        print(f"  \u2717 FAILED: {e}")
        print("=" * 60 + "\n")
        return False


def run_all_dft_baselines(
    atoms: Atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    charge: int = 0,
    multiplicity: int = 1,
    skip_if_exists: bool = True
) -> Dict[str, bool]:
    """
    Run all DFT baseline calculations.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Optimized molecular structure
    molecule_name : str
        Name of the molecule
    results_mgr : ResultsManager
        Results manager instance
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity
    skip_if_exists : bool
        If True, skip calculations if results already exist
        
    Returns
    -------
    dict
        Dictionary mapping baseline names to success status
    """
    results = {}
    
    print("\n" + "=" * 60)
    print("DFT BASELINE CALCULATIONS")
    print("=" * 60)
    print(f"Methods: {list(DFT_BASELINES.keys())}")
    print("=" * 60 + "\n")
    
    for baseline_name in DFT_BASELINES.keys():
        success = run_dft_baseline_calculation(
            atoms,
            molecule_name,
            baseline_name,
            results_mgr,
            charge,
            multiplicity,
            skip_if_exists
        )
        results[baseline_name] = success
    
    return results
