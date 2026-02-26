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
import sys
from pathlib import Path

import click

try:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except ImportError:
    pass

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


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
    type=click.Choice(["mace_omol", "mace_off", "mace_mp"]),
    help="Calculator for geometry optimization (default: mace_omol)",
)
@click.option(
    "--energy-calculators",
    default="mace_mp,mace_omol",
    help="Comma-separated list of energy calculators (default: mace_mp,mace_omol)",
)
@click.option(
    "--dipole-calculators",
    default="espaloma,mace_ml",
    help="Comma-separated list of dipole calculators (default: espaloma,mace_ml)",
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
    """
    from mace_gaussian.workflow import run_pipeline

    # Parse calculator lists
    energy_calc_list = [c.strip() for c in energy_calculators.split(",")]
    dipole_calc_list = [c.strip() for c in dipole_calculators.split(",")]

    input_path = Path(input_file)

    # Validate prerequisites before expensive imports
    from mace_gaussian.utils.exceptions import InputValidationError, PrerequisiteError
    from mace_gaussian.utils.validation import (
        detect_device,
        validate_all_prerequisites,
        validate_xyz_file,
    )

    # 1. Validate input file
    try:
        xyz_info = validate_xyz_file(str(input_path))
        click.echo(f"Input: {xyz_info['n_atoms']} atoms ({', '.join(set(xyz_info['symbols']))})")
    except InputValidationError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # 2. Validate prerequisites
    try:
        validate_all_prerequisites(
            check_gaussian=True,
            check_formchk_tool=True,
            dipole_model_path=None,  # Checked later by calculators
            helper_script_path=None,  # Checked later by workflow
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

    [PLACEHOLDER - Coming soon]

    This will generate comparison tables and plots for:
    - Energy comparisons
    - Frequency comparisons
    - Runtime analysis
    - Accuracy metrics

    Example:
        python cli.py compare water
    """
    click.echo(f"\n{'=' * 60}")
    click.echo("COMPARE COMMAND - COMING SOON")
    click.echo(f"{'=' * 60}\n")
    click.echo(f"This command will compare results for: {molecule}")
    click.echo(f"Results directory: {output_dir}")
    click.echo("\nPlanned features:")
    click.echo("  - Energy comparison table")
    click.echo("  - Frequency comparison (harmonic vs anharmonic)")
    click.echo("  - Runtime analysis")
    click.echo("  - Visualization plots")
    click.echo("\nTo implement this feature, add analysis code to a new module.")
    click.echo(f"{'=' * 60}\n")


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

    [PLACEHOLDER - Coming soon]

    Exports all calculation results to a single file for further analysis.

    Example:
        python cli.py export water --format csv
        python cli.py export water --format xlsx --output water_analysis.xlsx
    """
    if output is None:
        output = f"{molecule}_results.{format}"

    click.echo(f"\n{'=' * 60}")
    click.echo("EXPORT COMMAND - COMING SOON")
    click.echo(f"{'=' * 60}\n")
    click.echo(f"Molecule: {molecule}")
    click.echo(f"Format: {format}")
    click.echo(f"Output: {output}")
    click.echo(f"Results directory: {results_dir}")
    click.echo("\nPlanned features:")
    click.echo("  - Export all frequencies to CSV/Excel")
    click.echo("  - Export energies and dipole moments")
    click.echo("  - Export metadata (runtimes, parameters)")
    click.echo("  - JSON export for programmatic access")
    click.echo("\nTo implement this feature, add export code to a new module.")
    click.echo(f"{'=' * 60}\n")


@cli.command()
def diagnose():
    """
    Run diagnostic checks for available calculators and dependencies.

    Checks:
        - Available dipole calculators
        - Python environment
        - Required packages

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
