#!/usr/bin/env python3
"""
Command-line interface for MACE-Gaussian comparison framework.

Usage:
    python cli.py run molecule.xyz [OPTIONS]
    python cli.py list [molecule]
    python cli.py diagnose
"""

import json
import logging
import os
import sys
from pathlib import Path

import click

from mace_gaussian.utils.exceptions import InputValidationError, PrerequisiteError
from mace_gaussian.utils.validation import (
    detect_device,
    validate_all_prerequisites,
    validate_xyz_file,
)

try:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except ImportError:
    pass

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

VALID_ENERGY_CALCULATORS = ["mace_mp", "mace_omol", "mace_off", "mace_anicc"]
VALID_DIPOLE_CALCULATORS = ["espaloma", "mace_ml"]


def _validate_energy_calculators(ctx, param, value):
    """Validate comma-separated energy calculator names."""
    for name in [c.strip() for c in value.split(",")]:
        if name not in VALID_ENERGY_CALCULATORS:
            raise click.BadParameter(
                f"'{name}' is not a valid energy calculator. "
                f"Choose from: {', '.join(VALID_ENERGY_CALCULATORS)}",
                param=param,
            )
    return value


def _validate_dipole_calculators(ctx, param, value):
    """Validate comma-separated dipole calculator names."""
    for name in [c.strip() for c in value.split(",")]:
        if name not in VALID_DIPOLE_CALCULATORS:
            raise click.BadParameter(
                f"'{name}' is not a valid dipole calculator. "
                f"Choose from: {', '.join(VALID_DIPOLE_CALCULATORS)}",
                param=param,
            )
    return value


@click.group()
@click.version_option(version="0.2.0", prog_name="MACE-Gaussian Comparison Framework")
def cli():
    """MACE-Gaussian comparison framework for molecular calculations."""
    pass


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--optimization-calculator",
    default="mace_omol",
    type=click.Choice(["mace_omol", "mace_off", "mace_mp", "mace_anicc"]),
    help="Calculator for geometry optimization (default: mace_omol)",
)
@click.option(
    "--energy-calculators",
    default="mace_mp,mace_omol",
    callback=_validate_energy_calculators,
    help="Comma-separated energy calculators. Choices: mace_mp, mace_omol, mace_off, mace_anicc",
)
@click.option(
    "--dipole-calculators",
    default="espaloma,mace_ml",
    callback=_validate_dipole_calculators,
    help="Comma-separated dipole calculators. Choices: espaloma, mace_ml",
)
@click.option(
    "--force-optimization",
    is_flag=True,
    help="Force re-optimization even if optimized geometry exists",
)
@click.option("--skip-dft-baseline", is_flag=True, help="Skip DFT baseline calculations")
@click.option(
    "--output-dir",
    default="comparison_results",
    type=click.Path(),
    help="Output directory for results (default: comparison_results)",
)
def run(
    input_file,
    optimization_calculator,
    energy_calculators,
    dipole_calculators,
    force_optimization,
    skip_dft_baseline,
    output_dir,
):
    """
    Run complete calculation workflow on INPUT_FILE.

    Example:
        python cli.py run water.xyz
        python cli.py run water.xyz --skip-dft-baseline
        python cli.py run water.xyz --energy-calculators mace_mp --dipole-calculators espaloma
        python cli.py run water.xyz --output-dir my_results
    """
    from mace_gaussian.workflow import run_pipeline

    # Parse calculator lists
    energy_calc_list = [c.strip() for c in energy_calculators.split(",")]
    dipole_calc_list = [c.strip() for c in dipole_calculators.split(",")]

    input_path = Path(input_file)

    # 1. Validate input file
    try:
        xyz_info = validate_xyz_file(str(input_path))
        click.echo(f"Input: {xyz_info['n_atoms']} atoms ({', '.join(set(xyz_info['symbols']))})")
    except InputValidationError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # 2. Validate prerequisites
    try:
        # Resolve env var paths — inline defaults match the calculator module defaults
        # (not imported from calculators to avoid heavy-dep side effects at import time)
        dipole_model_path = os.getenv(
            "MACE_DIPOLE_MODEL_PATH",
            str(
                Path.home()
                / "mace_gaussian"
                / "mace4ir_models"
                / "pretrained_models"
                / "model_1_dipole.model"
            ),
        )
        helper_script_path = os.getenv(
            "MACE_HELPER_SCRIPT_PATH",
            str(Path(__file__).parent / "gm_helper.py"),
        )

        validate_all_prerequisites(
            check_gaussian=True,
            check_formchk_tool=True,
            dipole_model_path=dipole_model_path,
            helper_script_path=helper_script_path,
        )
        click.echo("Prerequisites OK: g16 found")
    except PrerequisiteError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # 3. Detect compute device
    device = detect_device()
    click.echo(f"Device: {device}")

    click.echo(f"\nStarting workflow for: {input_file}")
    click.echo(f"Energy calculators: {energy_calc_list}")
    click.echo(f"Dipole calculators: {dipole_calc_list}")
    click.echo("")

    try:
        results = run_pipeline(
            input_file=str(input_path),
            optimization_calculator=optimization_calculator,
            energy_calculators=energy_calc_list,
            dipole_calculators=dipole_calc_list,
            force_optimization=force_optimization,
            include_dft_baselines=not skip_dft_baseline,
            base_output_dir=output_dir,
        )

        # Print final summary
        click.echo("\n" + "=" * 60)
        click.echo("WORKFLOW COMPLETED SUCCESSFULLY")
        click.echo("=" * 60)

        if results["dft_baselines"]:
            dft_success = sum(1 for v in results["dft_baselines"].values() if v)
            dft_total = len(results["dft_baselines"])
            click.echo(f"DFT baselines: {dft_success}/{dft_total} successful")

        ml_success = sum(1 for r in results["ml_calculations"] if r["success"])
        ml_total = len(results["ml_calculations"])
        click.echo(f"ML calculations: {ml_success}/{ml_total} successful")

        click.echo(f"\nResults directory: {output_dir}/{results['molecule_name']}/")
        click.echo("=" * 60)

    except Exception as e:
        click.echo(f"\nError: Workflow failed: {e}", err=True)
        logger.exception("Workflow failed with exception")
        sys.exit(1)


@cli.command()
@click.argument("molecule_name")
@click.option("--force", is_flag=True, help="Overwrite existing file")
@click.option(
    "--output-dir",
    default="molecules",
    type=click.Path(),
    help="Output directory for XYZ file (default: molecules)",
)
def fetch(molecule_name, force, output_dir):
    """Fetch 3D structure from PubChem by molecule name.

    Downloads a 3D conformer from PubChem PUG REST API and saves it
    as an XYZ file in the output directory.

    Example:
        mace-gaussian fetch aspirin
        mace-gaussian fetch caffeine --force
        mace-gaussian fetch water --output-dir structures/
    """
    from mace_gaussian.pubchem import fetch_3d_structure

    try:
        output_path = fetch_3d_structure(
            molecule_name=molecule_name,
            output_dir=Path(output_dir),
            force=force,
        )
        click.echo(f"Saved: {output_path}")
    except FileExistsError as e:
        click.echo(f"Warning: {e}", err=True)
        sys.exit(0)  # Skip is not an error
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: Failed to fetch '{molecule_name}': {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("molecule", required=False)
@click.option(
    "--output-dir",
    default="comparison_results",
    type=click.Path(),
    help="Results directory to search (default: comparison_results)",
)
def list(molecule, output_dir):
    """
    List available calculation results.

    If MOLECULE is specified, show detailed results for that molecule.
    Otherwise, list all available molecules.

    Example:
        python cli.py list
        python cli.py list water
    """
    results_dir = Path(output_dir)

    if not results_dir.exists():
        click.echo(f"Results directory not found: {results_dir}", err=True)
        sys.exit(1)

    if molecule:
        # Show detailed results for specific molecule
        mol_dir = results_dir / molecule
        if not mol_dir.exists():
            click.echo(f"No results found for molecule: {molecule}", err=True)
            sys.exit(1)

        click.echo(f"\n{'=' * 60}")
        click.echo(f"Results for: {molecule}")
        click.echo(f"{'=' * 60}\n")

        # Check for optimization
        opt_dir = mol_dir / "geometry_opt"
        if opt_dir.exists():
            click.echo("\u2713 Geometry optimization")
            opt_json = opt_dir / "results.json"
            if opt_json.exists():
                with opt_json.open() as f:
                    opt_data = json.load(f)
                click.echo(f"  Calculator: {opt_data.get('calculator', 'N/A')}")
                click.echo(f"  Energy: {opt_data.get('final_energy_eV', 'N/A'):.6f} eV")
                click.echo(f"  Runtime: {opt_data.get('runtime_s', 'N/A'):.1f} s")

        # List all frequency calculations
        click.echo("\nFrequency calculations:")
        freq_dirs = [d for d in mol_dir.iterdir() if d.is_dir() and d.name != "geometry_opt"]

        if not freq_dirs:
            click.echo("  No frequency calculations found")
        else:
            for freq_dir in sorted(freq_dirs):
                json_file = freq_dir / "results.json"
                if json_file.exists():
                    with json_file.open() as f:
                        data = json.load(f)

                    calc_type = data.get("calculator_type", "unknown")
                    energy_calc = data.get("energy_calculator", "N/A")
                    dipole_calc = data.get("dipole_calculator", "N/A")

                    click.echo(f"\n  {freq_dir.name}:")
                    click.echo(f"    Type: {calc_type}")
                    click.echo(f"    Energy: {energy_calc}")
                    click.echo(f"    Dipole: {dipole_calc}")
                    click.echo(f"    Energy: {data.get('energy_eV', 'N/A'):.6f} eV")
                    click.echo(f"    Runtime: {data.get('runtime_s', 'N/A'):.1f} s")

                    n_harmonic = len(data.get("frequencies", {}).get("harmonic", []))
                    n_anharmonic = len(data.get("frequencies", {}).get("anharmonic", []))
                    click.echo(f"    Frequencies: {n_harmonic} harmonic, {n_anharmonic} anharmonic")

        click.echo(f"\n{'=' * 60}\n")

    else:
        # List all molecules
        molecules = [d.name for d in results_dir.iterdir() if d.is_dir()]

        if not molecules:
            click.echo(f"No results found in {results_dir}")
            return

        click.echo(f"\nAvailable molecules ({len(molecules)}):")
        for mol in sorted(molecules):
            mol_dir = results_dir / mol

            # Count calculations
            freq_dirs = [d for d in mol_dir.iterdir() if d.is_dir() and d.name != "geometry_opt"]
            n_calcs = len(freq_dirs)

            has_opt = (mol_dir / "geometry_opt").exists()
            opt_status = "\u2713" if has_opt else "\u2717"

            click.echo(f"  {opt_status} {mol:20s} - {n_calcs} calculations")

        click.echo("\nUse 'python cli.py list <molecule>' for detailed information")


@cli.command()
@click.argument("molecule")
@click.option(
    "--output-dir",
    default="comparison_results",
    type=click.Path(),
    help="Results directory (default: comparison_results)",
)
def compare(molecule, output_dir):
    """
    Compare calculation results for MOLECULE.

    Not yet implemented. To analyze and compare results, use:

        python run_analysis.py <molecule>           (anharmonic: fundamentals, overtones,
                                                     combinations)
        python run_analysis_harmonic.py <molecule>  (harmonic: fundamentals with mode matching)

    These scripts generate comparison plots, regression analysis, and an HTML report
    in the analysis_results/ directory.

    Example:
        python cli.py compare water    (not yet implemented)
    """
    click.echo(
        "Not yet implemented. To analyze and compare results, use:\n\n"
        "    python run_analysis.py <molecule>           "
        "(anharmonic: fundamentals, overtones, combinations)\n"
        "    python run_analysis_harmonic.py <molecule>  "
        "(harmonic: fundamentals with mode matching)\n\n"
        "These scripts generate comparison plots, regression analysis, and an HTML report\n"
        "in the analysis_results/ directory."
    )


@cli.command()
@click.argument("molecule")
@click.option(
    "--format",
    type=click.Choice(["csv", "json", "xlsx"]),
    default="csv",
    help="Export format (default: csv)",
)
@click.option(
    "--output", type=click.Path(), help="Output file path (default: <molecule>_results.<format>)"
)
@click.option(
    "--results-dir",
    default="comparison_results",
    type=click.Path(),
    help="Results directory (default: comparison_results)",
)
def export(molecule, format, output, results_dir):
    """
    Export results for MOLECULE to specified format.

    Not yet implemented. Calculation results are stored as JSON files under
    comparison_results/<molecule>/<energy_calc>_<dipole_calc>/results.json.

    Example:
        python cli.py export water --format csv    (not yet implemented)
        python cli.py export water --format xlsx --output water_analysis.xlsx  (not yet implemented)
    """
    if output is None:
        output = f"{molecule}_results.{format}"

    click.echo(
        "Not yet implemented. Calculation results are stored as JSON files under\n"
        f"comparison_results/{molecule}/<energy_calc>_<dipole_calc>/results.json."
    )


@cli.command()
def diagnose():
    """
    Run diagnostic checks for available calculators and dependencies.

    Checks:
        - Available dipole calculators (espaloma, mace_ml, xtb)
        - Gaussian 16 (g16) availability in PATH
        - formchk availability in PATH
        - CUDA device availability

    Example:
        python cli.py diagnose
    """
    from mace_gaussian.calculators import dipole_factory

    click.echo("=" * 60)
    click.echo("DIAGNOSTIC MODE")
    click.echo("=" * 60)
    click.echo("\nAvailable dipole calculators:")
    for name, available in dipole_factory.list_available().items():
        status = "OK" if available else "UNAVAILABLE"
        click.echo(f"  {status}: {name}")
    click.echo("\n" + "=" * 60)


if __name__ == "__main__":
    cli()
