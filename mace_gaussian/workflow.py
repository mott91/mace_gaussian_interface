"""workflow.py — MACE-Gaussian pipeline orchestrator.

This module is the single authoritative orchestrator for the three-stage pipeline:

  Stage 1 — Geometry optimization (run_geometry_optimization)
  Stage 2 — DFT baseline calculations (run_dft_baselines)
  Stage 3 — ML frequency calculations (run_ml_calculations)

The top-level coordinator run_pipeline() calls all three stages in order.
Stages 2 and 3 each depend only on the geometry output from stage 1 (the
cluster seam): run DFT on cluster → copy results → resume with
include_dft_baselines=False.

Low-level helpers (update_molecule_geometry, calculate_energy_and_forces,
calculate_hessian, calculate_dipole_properties, run_next_calculation,
geometry_optimisation, calculator) are also provided here for use by the
Gaussian ZMQ callback.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import time
import warnings
from pathlib import Path

import numpy as np
from ase.io import read
from ase.optimize import LBFGS

from .calculators import dipole_factory
from .gaussian.fchk import convert_chk_to_fchk
from .gaussian.io import ase_to_gjf, parse_gaussian_input, write_gaussian_output
from .gaussian.parser import parse_gaussian_log
from .gaussian.runner import DEFAULT_TIMEOUT_SECONDS, run_gaussian_with_zmq
from .utils.results import ResultsManager
from .utils.units import BOHR_TO_ANGSTROM, HARTREE_TO_EV

# Conversion factor for polarizability: Angstrom^3 -> Bohr^3
ANGSTROM3_TO_BOHR3 = (1.0 / BOHR_TO_ANGSTROM) ** 3

warnings.filterwarnings("ignore", message=".*weights_only=False.*", category=FutureWarning)
os.environ["PYTHONWARNINGS"] = "ignore::FutureWarning"

# Elements supported by mace_anicc (HCNO molecules only)
_MACE_ANICC_SUPPORTED_ELEMENTS = frozenset({"H", "C", "N", "O"})

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# ============================================================================
# LOW-LEVEL HELPERS
# ============================================================================


def update_molecule_geometry(atoms, coordinates: np.ndarray, charge: int, spin: int):
    """Update ASE Atoms object with new geometry and electronic state.

    Args:
        atoms: ASE Atoms object to update (modified in-place)
        coordinates: New atomic coordinates in Angstroms, shape (natoms, 3)
        charge: Molecular charge
        spin: Spin multiplicity
    """
    atoms.set_positions(coordinates)
    atoms.info["charge"] = float(charge)
    atoms.info["spin"] = float(spin)


def calculate_energy_and_forces(atoms, calculator) -> tuple[float, np.ndarray]:
    """Calculate energy and forces using the attached calculator."""
    try:
        energy = atoms.get_potential_energy()
        gradient = -atoms.get_forces()
        return energy, gradient
    except Exception as e:
        logger.error(f"Energy/forces calculation failed: {e}")
        import traceback

        logger.error(f"Full traceback:\n{traceback.format_exc()}")
        raise


def calculate_hessian(atoms, calculator, natoms: int) -> np.ndarray | None:
    """Calculate Hessian matrix (second derivatives)."""
    try:
        hessian = calculator.get_hessian(atoms=atoms)

        # Convert to Hartree/Bohr^2 for Gaussian
        hessian = hessian * (BOHR_TO_ANGSTROM**2) / HARTREE_TO_EV
        hessian = hessian.reshape(3 * natoms, 3 * natoms)

        return hessian
    except Exception as e:
        logger.error(f"Hessian calculation failed: {e}")
        import traceback

        logger.error(f"Full traceback:\n{traceback.format_exc()}")
        return None


def calculate_dipole_properties(
    atoms, dipole_calc, deriv: int, calculate_derivatives: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]:
    """Calculate dipole moment and optionally its derivatives.

    Args:
        atoms: ASE Atoms object
        dipole_calc: DipoleCalculatorBase instance
        deriv: Derivative level from Gaussian
        calculate_derivatives: Whether to calculate dipole derivatives

    Returns:
        Tuple of (dipole, dipole_derivatives, partial_charges, polarizability_voigt6)
        - dipole: Dipole moment in e*Bohr, shape (3,)
        - dipole_derivatives: Dipole derivatives in e, shape (3*natoms, 3)
        - partial_charges: Partial atomic charges (or None), shape (natoms,)
        - polarizability_voigt6: shape (6,) in Bohr^3 [axx,axy,ayy,axz,ayz,azz], or zeros
    """
    natoms = len(atoms)

    try:
        # Calculate dipole moment
        logger.info(f"Using dipole calculator: {dipole_calc.name}")
        td1 = time.perf_counter()
        dipole, partial_charges = dipole_calc.calculate_dipole(atoms)
        td_dipole = time.perf_counter() - td1

        # Store charges in atoms object if available
        if partial_charges is not None:
            atoms.set_initial_charges(partial_charges)
            logger.debug(f"Partial charges sum: {np.sum(partial_charges):.6f}")

        # Calculate dipole derivatives for IR intensities
        td_derivs = 0.0
        if calculate_derivatives and deriv >= 1:
            logger.info("Calculating dipole derivatives for IR intensities...")
            td1 = time.perf_counter()
            dipole_derivatives = dipole_calc.calculate_dipole_derivatives(atoms, displacement=0.005)
            td_derivs = time.perf_counter() - td1
        else:
            dipole_derivatives = np.zeros((3 * natoms, 3))

        # Extract polarizability if calculator supports it (MACEDipoleCalculator)
        td_polar = 0.0
        if hasattr(dipole_calc, "calculate_polarizability"):
            td1 = time.perf_counter()
            polar_3x3 = dipole_calc.calculate_polarizability(atoms)  # (3,3) Angstrom^3
            td_polar = time.perf_counter() - td1
            p = polar_3x3 * ANGSTROM3_TO_BOHR3
            polarizability_voigt6 = np.array(
                [p[0, 0], p[0, 1], p[1, 1], p[0, 2], p[1, 2], p[2, 2]]
            )
        else:
            polarizability_voigt6 = np.zeros(6)

        logger.info(f"Dipole calculated: {dipole} e*Bohr")
        logger.debug(f"Dipole timing: dipole={td_dipole:.2f}s derivs={td_derivs:.2f}s polar={td_polar:.2f}s")
        return dipole, dipole_derivatives, partial_charges, polarizability_voigt6

    except Exception as e:
        logger.error(f"Dipole calculation failed: {e}")
        logger.warning("Falling back to zero dipole (IR intensities will be zero)")

        # Fallback to zeros
        dipole = np.zeros(3)
        dipole_derivatives = np.zeros((3 * natoms, 3))
        return dipole, dipole_derivatives, None, np.zeros(6)


def run_next_calculation(
    mol,
    msg: str,
    calculator,
    dipole_method: str = "auto",
    calculate_derivatives: bool = True,
):
    """Coordinate a single calculation cycle for Gaussian external interface.

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
    print("  -> Processing Gaussian request...")
    t0 = time.perf_counter()

    # Parse message to get file paths
    infile, outfile = msg.split("|")

    # Step 1: Parse Gaussian input file
    natoms, deriv, charge, spin, coordinates, _atomnames = parse_gaussian_input(infile)

    # Step 2: Update molecule geometry and electronic state
    update_molecule_geometry(mol, coordinates, charge, spin)

    # Step 3: Calculate energy and forces
    t1 = time.perf_counter()
    energy, gradient = calculate_energy_and_forces(mol, calculator)
    t_energy = time.perf_counter() - t1

    # Step 4: Calculate Hessian if needed
    hessian = None
    t_hessian = 0.0
    if deriv >= 2:
        t1 = time.perf_counter()
        hessian = calculate_hessian(mol, calculator, natoms)
        t_hessian = time.perf_counter() - t1

    # Step 5: Calculate dipole properties
    t1 = time.perf_counter()
    dipole_calc = dipole_factory.get_calculator(dipole_method)
    dipole, dipole_derivatives, _partial_charges, polarizability_voigt6 = (
        calculate_dipole_properties(mol, dipole_calc, deriv, calculate_derivatives)
    )
    t_dipole = time.perf_counter() - t1

    # Thread dalpha_dr through pipeline (Python-only; NOT written to Gaussian file)
    dalpha_dr = getattr(dipole_calc, "_last_dalpha_dr", None)
    if dalpha_dr is not None:
        mol.info["dalpha_dr"] = dalpha_dr
        logger.debug("dalpha_dr shape: %s", dalpha_dr.shape)

    # Step 6: Write results to Gaussian output file
    write_gaussian_output(
        outfile,
        natoms,
        energy,
        gradient,
        dipole,
        dipole_derivatives,
        hessian,
        deriv,
        polarizability=polarizability_voigt6,
    )

    elapsed = time.perf_counter() - t0
    print(f"  -> Gaussian call completed in {elapsed:.2f}s (deriv={deriv})")
    logger.debug(
        f"Call timing: energy={t_energy:.2f}s hessian={t_hessian:.2f}s dipole={t_dipole:.2f}s"
    )


def geometry_optimisation(mol, fmax=0.000001):
    """Run LBFGS geometry optimisation on mol in-place and return (mol, steps)."""
    ei = mol.get_potential_energy()
    print("Initial Energy: ", ei, "eV")
    opt = LBFGS(mol)

    opt.run(fmax=fmax, steps=10000)

    num_steps = opt.get_number_of_steps()
    print(f"Optimization steps: {num_steps}")

    ef = mol.get_potential_energy()
    print("Final Energy: ", ef, "eV")

    return mol, num_steps


def _check_mace_anicc_elements(atoms) -> None:
    """Raise ValueError if atoms contains elements outside H/C/N/O (mace_anicc restriction)."""
    symbols = set(atoms.get_chemical_symbols())
    unsupported = sorted(symbols - _MACE_ANICC_SUPPORTED_ELEMENTS)
    if unsupported:
        raise ValueError(
            f"mace_anicc only supports H/C/N/O — "
            f"molecule contains: {unsupported}"
        )


def calculator(nnp):
    """Return an ASE calculator for the named neural network potential."""
    if nnp == "mace_mp":
        from mace.calculators import mace_mp

        calc = mace_mp(model="large", device="cuda", default_dtype="float64", dispersion=False)
        return calc

    if nnp == "mace_off":
        from mace.calculators import mace_off

        calc = mace_off(model="large", device="cuda", default_dtype="float64", dispersion=False)
        return calc

    if nnp == "mace_omol":
        from mace.calculators import mace_omol

        calc = mace_omol(
            model="extra_large",
            device="cuda",
            default_dtype="float64",
            dispersion=False,
        )
        return calc

    if nnp == "mace_anicc":
        from mace.calculators import mace_anicc

        calc = mace_anicc(device="cuda")
        return calc


# ============================================================================
# STAGE FUNCTIONS
# ============================================================================


def run_geometry_optimization(
    atoms, molecule_name: str, results_mgr: ResultsManager, calculator_name: str = "mace_omol"
):
    """Run geometry optimization and save results.

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
    optimized_atoms, num_steps = geometry_optimisation(atoms)

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
        num_steps=num_steps,
        runtime=runtime,
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
    multiplicity: int = 1,
):
    """Run frequency calculation with specified calculators.

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

    # Element guard for mace_anicc (must run before model load)
    if energy_calculator_name == "mace_anicc":
        _check_mace_anicc_elements(atoms)

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

        print(f"  -> Gaussian input: {gjf_temp}")

        # Run Gaussian calculation
        print("  -> Launching Gaussian...")

        def _on_request(msg: str) -> str:
            run_next_calculation(
                mol,
                msg,
                calc,
                dipole_method=dipole_calculator_name,
                calculate_derivatives=True,
            )
            return "ready"

        run_gaussian_with_zmq(
            gjf_temp,
            on_request=_on_request,
            timeout_seconds=DEFAULT_TIMEOUT_SECONDS,
        )

        runtime = time.time() - start_time

        # Now move the files to the frequency directory
        log_temp = gjf_temp.replace(".gjf", ".log")
        chk_temp = gjf_temp.replace(".gjf", ".chk")

        gjf_final = str(freq_dir / "gaussian_freq.gjf")
        log_final = str(freq_dir / "gaussian_freq.log")
        chk_final = str(freq_dir / "gaussian_freq.chk")

        # Move files (not copy, since they're in different locations)
        if Path(gjf_temp).exists():
            shutil.move(gjf_temp, gjf_final)
        if Path(log_temp).exists():
            shutil.move(log_temp, log_final)
        # Save checkpoint file for mode matching
        if Path(chk_temp).exists():
            shutil.move(chk_temp, chk_final)
            # Automatically convert to .fchk for mode matching
            try:
                fchk_final = chk_final.replace(".chk", ".fchk")
                logger.info("Converting .chk to .fchk for mode matching...")
                convert_chk_to_fchk(chk_final, fchk_final)
                logger.info(f"Created {fchk_final}")
            except Exception as e:
                logger.warning(f"Could not convert .chk to .fchk: {e}")
                logger.warning("Mode matching will not be available for this calculation")

        # Parse Gaussian log file to extract frequencies
        parsed_data = {}
        if Path(log_final).exists():
            try:
                logger.info("Parsing Gaussian log file...")
                parsed_data = parse_gaussian_log(log_final)

                # Convert energy from Hartree to eV if available
                if parsed_data.get("final_energy_hartree"):
                    final_energy = parsed_data["final_energy_hartree"] * HARTREE_TO_EV
                else:
                    final_energy = mol.get_potential_energy()

            except Exception as e:
                logger.error(f"Failed to parse Gaussian log: {e}")
                parsed_data = {"harmonic": [], "anharmonic": []}
                final_energy = mol.get_potential_energy()
        else:
            parsed_data = {"harmonic": [], "anharmonic": []}
            final_energy = mol.get_potential_energy()

        calculation_parameters = {
            "energy_calculator": energy_calculator_name,
            "dipole_calculator": dipole_calculator_name,
        }
        results_mgr.save_frequency_results(
            molecule_name=molecule_name,
            energy_calculator=energy_calculator_name,
            dipole_calculator=dipole_calculator_name,
            calculator_type="ml",
            frequencies_data={
                "harmonic": parsed_data.get("harmonic", []),
                "anharmonic": parsed_data.get("anharmonic", []),
                "overtones": parsed_data.get("overtones", []),
                "combination_bands": parsed_data.get("combination_bands", []),
            },
            energy=final_energy,
            dipole=parsed_data.get("dipole_moment"),
            runtime=runtime,
            gaussian_log=log_final if Path(log_final).exists() else None,
            gaussian_gjf=gjf_final if Path(gjf_final).exists() else None,
            timestamp=timestamp,
            calculation_parameters=calculation_parameters,
        )

        print(f"  Completed in {runtime:.1f} seconds")
        print("=" * 60 + "\n")

        return True

    except Exception as e:
        logger.error(f"Frequency calculation failed: {e}")
        import traceback

        logger.error(traceback.format_exc())
        print("  FAILED")
        print("=" * 60 + "\n")
        return False


def run_dft_baselines(
    optimized_atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    charge: int = 0,
    multiplicity: int = 1,
) -> dict[str, bool]:
    """Run DFT baseline calculations. Wraps dft_baseline.run_all_dft_baselines().

    The dft_baseline import is lazy (inside this function body) to avoid
    heavy DGL/espaloma imports at module load time.
    """
    from .dft_baseline import run_all_dft_baselines

    return run_all_dft_baselines(
        optimized_atoms, molecule_name, results_mgr, charge, multiplicity, skip_if_exists=True
    )


def run_ml_calculations(
    optimized_atoms,
    molecule_name: str,
    results_mgr: ResultsManager,
    energy_calculators: list[str],
    dipole_calculators: list[str],
    charge: int = 0,
    multiplicity: int = 1,
) -> list[dict]:
    """Run all energy+dipole calculator combinations via Gaussian ZMQ interface."""
    print("=" * 60)
    print("ML FREQUENCY CALCULATIONS")
    print("=" * 60)
    print(f"Energy calculators: {energy_calculators}")
    print(f"Dipole calculators: {dipole_calculators}")
    print("=" * 60 + "\n")

    results = []
    for energy_calc in energy_calculators:
        for dipole_calc in dipole_calculators:
            success = run_frequency_calculation(
                optimized_atoms,
                molecule_name,
                energy_calc,
                dipole_calc,
                results_mgr,
                charge,
                multiplicity,
            )
            results.append({"energy": energy_calc, "dipole": dipole_calc, "success": success})
    return results


# ============================================================================
# TOP-LEVEL COORDINATOR
# ============================================================================


def run_pipeline(
    input_file: str,
    optimization_calculator: str = "mace_omol",
    energy_calculators: list | None = None,
    dipole_calculators: list | None = None,
    force_optimization: bool = False,
    include_dft_baselines: bool = True,
    base_output_dir: str = "comparison_results",
) -> dict:
    """Run the complete MACE-Gaussian comparison pipeline.

    Three stages in sequence:
      1. Geometry optimization (run_geometry_optimization)
      2. DFT baseline calculations (run_dft_baselines)  -- cluster seam
      3. ML frequency calculations (run_ml_calculations) -- cluster seam

    Stages 2 and 3 each depend only on the geometry output from stage 1.
    This is the cluster seam: run DFT on cluster, copy results back, then
    resume with include_dft_baselines=False to run stage 3 only.

    Parameters
    ----------
    input_file : str
        Path to input XYZ file
    optimization_calculator : str
        Calculator to use for geometry optimization (default: 'mace_omol')
    energy_calculators : list, optional
        List of energy calculators to use (default: ['mace_mp', 'mace_omol'])
    dipole_calculators : list, optional
        List of dipole calculators to use (default: ['espaloma', 'mace_ml'])
    force_optimization : bool
        If True, force re-optimization even if optimized geometry exists
    include_dft_baselines : bool
        If True, run DFT baseline calculations (default: True)
    base_output_dir : str
        Base directory for results (default: 'comparison_results')

    Returns
    -------
    dict
        Summary of results with success status for each calculation
    """
    from .utils.validation import detect_device

    # Detect and log compute device
    detect_device()

    # Set defaults
    if energy_calculators is None:
        energy_calculators = ["mace_mp", "mace_omol"]
    if dipole_calculators is None:
        dipole_calculators = ["espaloma", "mace_ml"]

    # Extract molecule name
    molecule_name = Path(input_file).stem

    print("\n" + "=" * 60)
    print("MACE-GAUSSIAN COMPARISON FRAMEWORK")
    print("=" * 60)
    print(f"Molecule: {molecule_name}")
    print(f"Input: {input_file}")
    print("=" * 60 + "\n")

    # Initialize results manager
    results_mgr = ResultsManager(base_output_dir=base_output_dir)

    # ========================================================================
    # STAGE 1: GEOMETRY OPTIMIZATION
    # ========================================================================

    # Check if optimized geometry already exists
    opt_dir = results_mgr.create_optimization_directory(molecule_name)
    opt_file = opt_dir / "optimized.xyz"

    if opt_file.exists() and not force_optimization:
        print("=" * 60)
        print("GEOMETRY OPTIMIZATION")
        print("=" * 60)
        print("Found existing optimized geometry")
        print(f"  Loading from: {opt_file}")
        print("  (use --force-optimization to re-run)")

        optimized_mol = read(str(opt_file))

        # Load metadata to get charge and spin
        json_file = opt_dir / "results.json"
        if json_file.exists():
            with json_file.open() as f:
                json.load(f)
            # Note: charge/spin might not be in old metadata, use defaults
            optimized_mol.info["charge"] = 0.0
            optimized_mol.info["spin"] = 1.0
        else:
            optimized_mol.info["charge"] = 0.0
            optimized_mol.info["spin"] = 1.0

        print("=" * 60 + "\n")
    else:
        if force_optimization and opt_file.exists():
            print("=" * 60)
            print("GEOMETRY OPTIMIZATION")
            print("=" * 60)
            print("Forcing re-optimization (existing geometry will be overwritten)")
            print("=" * 60 + "\n")

        # Read initial structure
        mol = read(input_file)
        mol.info["charge"] = 0.0
        mol.info["spin"] = 1.0

        # Element guard for mace_anicc (must run before model load)
        if optimization_calculator == "mace_anicc":
            _check_mace_anicc_elements(mol)
        # Setup calculator for optimization
        calc = calculator(optimization_calculator)
        mol.calc = calc

        try:
            optimized_mol = run_geometry_optimization(
                mol, molecule_name, results_mgr, optimization_calculator
            )
        except Exception as e:
            logger.error(f"Geometry optimization failed: {e}")
            import traceback

            logger.error(traceback.format_exc())
            raise

    # Get charge and spin for subsequent calculations
    charge = int(optimized_mol.info.get("charge", 0))
    multiplicity = int(optimized_mol.info.get("spin", 1))

    # ========================================================================
    # STAGE 2: DFT BASELINE CALCULATIONS
    # ========================================================================

    dft_results = {}
    if include_dft_baselines:
        dft_results = run_dft_baselines(
            optimized_mol, molecule_name, results_mgr, charge, multiplicity
        )

    # ========================================================================
    # STAGE 3: ML FREQUENCY CALCULATIONS
    # ========================================================================

    ml_results = run_ml_calculations(
        optimized_mol,
        molecule_name,
        results_mgr,
        energy_calculators,
        dipole_calculators,
        charge,
        multiplicity,
    )

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 60)
    print("CALCULATION SUMMARY")
    print("=" * 60)

    # DFT baselines summary
    if dft_results:
        print("\nDFT Baselines:")
        for baseline, success in dft_results.items():
            status = "+" if success else "-"
            print(f"  {status} {baseline}")

    # ML calculations summary
    print("\nML Calculations:")
    successful_ml = sum(1 for r in ml_results if r["success"])
    total_ml = len(ml_results)

    print(f"Completed: {successful_ml}/{total_ml}")

    for r in ml_results:
        status = "+" if r["success"] else "-"
        print(f"  {status} {r['energy']} + {r['dipole']}")

    print(f"\nResults saved to: {base_output_dir}/{molecule_name}/")
    print("=" * 60)

    return {
        "dft_baselines": dft_results,
        "ml_calculations": ml_results,
        "molecule_name": molecule_name,
    }
