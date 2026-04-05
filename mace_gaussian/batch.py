"""Batch runner with manifest-driven per-calculator restart.

Runs the full MACE-Gaussian pipeline over a list of molecules with:
- Per-calculator failure isolation (one bad combo doesn't halt the batch)
- Manifest-based restart (skip already-complete combinations)
- Atomic manifest writes (crash-safe via tempfile + os.replace)
- Summary table at end with molecule status and timing
"""

from __future__ import annotations

import contextlib
import json
import os
import tempfile
import time
from datetime import datetime
from pathlib import Path

import click
from ase.io import read

from .utils.results import ResultsManager

STATUS_PENDING = "pending"
STATUS_COMPLETE = "complete"
STATUS_FAILED = "failed"


def load_manifest(path: Path) -> dict:
    """Load JSON manifest or return empty skeleton.

    Parameters
    ----------
    path : Path
        Path to batch_manifest.json

    Returns
    -------
    dict
        Manifest data with at least {"molecules": {}}
    """
    if path.exists():
        with path.open() as f:
            return json.load(f)
    return {"molecules": {}}


def save_manifest(manifest: dict, path: Path) -> None:
    """Write manifest atomically to prevent corruption on interrupt.

    Uses tempfile + os.replace for POSIX atomic write safety.

    Parameters
    ----------
    manifest : dict
        Manifest data to write
    path : Path
        Target path for the manifest file
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(manifest, f, indent=2)
        os.replace(tmp, str(path))  # noqa: PTH105
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)  # noqa: PTH108
        raise


def parse_batch_file(batch_file: Path) -> list[Path]:
    """Read molecule paths from a batch file.

    Skips blank lines and lines starting with #. Resolves relative
    paths against os.getcwd().

    Parameters
    ----------
    batch_file : Path
        Path to text file with one .xyz path per line

    Returns
    -------
    list[Path]
        List of resolved, absolute paths to .xyz files

    Raises
    ------
    FileNotFoundError
        If any resolved path does not exist
    """
    paths = []
    cwd = Path.cwd()
    with batch_file.open() as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            p = Path(stripped)
            if not p.is_absolute():
                p = cwd / p
            p = p.resolve()
            if not p.exists():
                raise FileNotFoundError(f"Molecule file not found: {p}")
            paths.append(p)
    return paths


def _combination_key(energy_calc: str, dipole_calc: str) -> str:
    """Return manifest key for an energy+dipole calculator combination."""
    return f"{energy_calc}_{dipole_calc}"


def run_batch(
    batch_file: Path,
    optimization_calculator: str,
    energy_calculators: list[str],
    dipole_calculators: list[str],
    skip_dft_baseline: bool,
    output_dir: str = "comparison_results",
    keep_scratch: bool = False,
    dft_on_cluster: str | None = None,
    slurm_template: str | None = None,
) -> dict:
    """Run the full pipeline for multiple molecules with manifest-based restart.

    For each molecule, runs geometry optimization, optional DFT baselines,
    and all energy x dipole calculator combinations. Each combination is
    tracked independently in the manifest for fine-grained restart.

    Parameters
    ----------
    batch_file : Path
        Path to text file listing .xyz files (one per line)
    optimization_calculator : str
        Calculator for geometry optimization
    energy_calculators : list[str]
        Energy calculators to run
    dipole_calculators : list[str]
        Dipole calculators to run
    skip_dft_baseline : bool
        If True, skip DFT baseline calculations
    output_dir : str
        Output directory for results
    keep_scratch : bool
        If True, preserve scratch directories
    dft_on_cluster : str or None
        SSH target for SLURM DFT offloading (e.g. ``user@hostname``).
        When set, DFT baselines are submitted as SLURM jobs instead of
        running locally.
    slurm_template : str or None
        Path to custom SLURM job template. Defaults to
        ``templates/slurm_dft.sh`` relative to the package root.

    Returns
    -------
    dict
        Summary with keys: complete, failed, skipped, molecules
    """
    molecules = parse_batch_file(batch_file)
    manifest_path = Path(output_dir) / "batch_manifest.json"
    manifest = load_manifest(manifest_path)

    # Store run options in manifest for drift detection
    manifest["options"] = {
        "energy_calculators": energy_calculators,
        "dipole_calculators": dipole_calculators,
        "skip_dft_baseline": skip_dft_baseline,
        "optimization_calculator": optimization_calculator,
    }
    manifest["started"] = manifest.get("started", datetime.now().isoformat())
    manifest["updated"] = datetime.now().isoformat()

    results_mgr = ResultsManager(base_output_dir=output_dir)
    total = len(molecules)
    summary = {"complete": 0, "failed": 0, "skipped": 0}

    for i, xyz_path in enumerate(molecules, 1):
        molecule_name = xyz_path.stem
        mol_start = time.time()
        click.echo(f"[{i}/{total}] Running {molecule_name}...")

        # Initialize molecule entry in manifest if not present
        if molecule_name not in manifest["molecules"]:
            manifest["molecules"][molecule_name] = {
                "xyz_path": str(xyz_path),
                "geometry_opt": STATUS_PENDING,
                "dft_baseline": STATUS_PENDING if not skip_dft_baseline else "skipped",
                "combinations": {},
            }
        mol_manifest = manifest["molecules"][molecule_name]

        try:
            # Stage 1: Geometry optimization (per-molecule, run once)
            atoms = read(str(xyz_path))
            atoms.info["charge"] = 0.0
            atoms.info["spin"] = 1.0
            opt_geom_path = results_mgr.get_optimized_geometry_path(molecule_name)

            if opt_geom_path.exists() and mol_manifest.get("geometry_opt") == STATUS_COMPLETE:
                optimized_atoms = read(str(opt_geom_path))
            else:
                from .workflow import calculator, run_geometry_optimization

                calc = calculator(optimization_calculator)
                atoms.calc = calc
                optimized_atoms = run_geometry_optimization(
                    atoms,
                    molecule_name,
                    results_mgr,
                    calculator_name=optimization_calculator,
                )
                mol_manifest["geometry_opt"] = STATUS_COMPLETE
                save_manifest(manifest, manifest_path)

            # Stage 2: DFT baselines (per-molecule, run once)
            if not skip_dft_baseline and mol_manifest.get("dft_baseline") != STATUS_COMPLETE:
                if dft_on_cluster:
                    # Submit DFT to SLURM immediately after geom opt
                    slurm_info = mol_manifest.get("slurm", {})
                    if slurm_info.get("status") not in ("SUBMITTED", "RUNNING", "COMPLETED"):
                        from .dft_baseline import create_gaussian_dft_input
                        from .slurm import submit_dft_jobs

                        template = (
                            Path(slurm_template)
                            if slurm_template
                            else Path(__file__).parent.parent / "templates" / "slurm_dft.sh"
                        )

                        gjf_dir = Path(output_dir) / molecule_name / "b3lyp_6-31Gdp"
                        gjf_dir.mkdir(parents=True, exist_ok=True)
                        gjf_path = gjf_dir / f"{molecule_name}_freq_anharm.gjf"
                        create_gaussian_dft_input(
                            optimized_atoms,
                            str(gjf_path),
                            method="b3lyp",
                            basis="6-31G(d,p)",
                            title=molecule_name,
                            output_dir=str(gjf_dir),
                            nproc=4,
                            mem="4GB",
                        )
                        job_ids = submit_dft_jobs(
                            [{"name": molecule_name, "gjf_path": str(gjf_path)}],
                            dft_on_cluster,
                            template,
                            output_dir,
                        )
                        if molecule_name in job_ids:
                            mol_manifest["slurm"] = {
                                "job_id": job_ids[molecule_name],
                                "host": dft_on_cluster,
                                "remote_dir": f"~/mace_gaussian_dft/{molecule_name}",
                                "status": "SUBMITTED",
                            }
                            click.echo(
                                f"  Submitted DFT for {molecule_name} "
                                f"(job {job_ids[molecule_name]})"
                            )
                        save_manifest(manifest, manifest_path)
                else:
                    from .workflow import run_dft_baselines

                    run_dft_baselines(
                        optimized_atoms,
                        molecule_name,
                        results_mgr,
                    )
                    mol_manifest["dft_baseline"] = STATUS_COMPLETE
                    save_manifest(manifest, manifest_path)

            # Stage 3: ML combinations (per-calculator granularity)
            from .workflow import run_frequency_calculation

            mol_failed = False
            for energy_calc in energy_calculators:
                for dipole_calc in dipole_calculators:
                    combo_key = _combination_key(energy_calc, dipole_calc)
                    existing = mol_manifest["combinations"].get(combo_key, {})

                    if existing.get("status") == STATUS_COMPLETE:
                        summary["skipped"] += 1
                        continue

                    combo_start = time.time()
                    try:
                        success = run_frequency_calculation(
                            optimized_atoms,
                            molecule_name,
                            energy_calc,
                            dipole_calc,
                            results_mgr,
                        )
                        combo_runtime = time.time() - combo_start
                        if success:
                            mol_manifest["combinations"][combo_key] = {
                                "status": STATUS_COMPLETE,
                                "runtime_s": round(combo_runtime, 1),
                            }
                            summary["complete"] += 1
                        else:
                            mol_manifest["combinations"][combo_key] = {
                                "status": STATUS_FAILED,
                                "error": "run_frequency_calculation returned False",
                                "runtime_s": round(combo_runtime, 1),
                            }
                            summary["failed"] += 1
                            mol_failed = True
                    except Exception as e:
                        combo_runtime = time.time() - combo_start
                        mol_manifest["combinations"][combo_key] = {
                            "status": STATUS_FAILED,
                            "error": str(e),
                            "runtime_s": round(combo_runtime, 1),
                        }
                        summary["failed"] += 1
                        mol_failed = True
                    save_manifest(manifest, manifest_path)

            # Stage 4: Harmonic analysis (if any combos succeeded)
            complete_count = sum(
                1
                for c in mol_manifest["combinations"].values()
                if c.get("status") == STATUS_COMPLETE
            )
            if complete_count > 0:
                try:
                    from .analysis import analyze_molecule_harmonic

                    click.echo(f"  Running harmonic analysis for {molecule_name}...")
                    analyze_molecule_harmonic(molecule_name, base_results_dir=output_dir)
                    mol_manifest["analysis_harmonic"] = STATUS_COMPLETE
                except Exception as e:
                    click.echo(f"  Warning: Harmonic analysis failed: {e}", err=True)
                    mol_manifest["analysis_harmonic"] = STATUS_FAILED
                save_manifest(manifest, manifest_path)

            mol_runtime = time.time() - mol_start
            status_str = "done" if not mol_failed else "done (with failures)"
            click.echo(
                f"[{i}/{total}] {molecule_name} {status_str} "
                f"({mol_runtime / 60:.0f}m {mol_runtime % 60:.0f}s)"
            )

        except Exception as e:
            mol_runtime = time.time() - mol_start
            click.echo(f"[{i}/{total}] {molecule_name} FAILED: {e}")
            # Mark all pending combinations as failed
            for energy_calc in energy_calculators:
                for dipole_calc in dipole_calculators:
                    combo_key = _combination_key(energy_calc, dipole_calc)
                    if combo_key not in mol_manifest.get("combinations", {}):
                        mol_manifest.setdefault("combinations", {})[combo_key] = {
                            "status": STATUS_FAILED,
                            "error": f"Molecule-level failure: {e}",
                        }
                        summary["failed"] += 1
            save_manifest(manifest, manifest_path)

    # SLURM DFT offload: poll submitted jobs, retrieve results
    if dft_on_cluster and not skip_dft_baseline:
        from .slurm import TERMINAL_STATES as _SLURM_TERMINAL
        from .slurm import poll_jobs, retrieve_results

        # Collect submitted jobs that haven't reached a terminal state
        pending_jobs: dict[str, str] = {}
        for mol_name, mol_data in manifest["molecules"].items():
            slurm_info = mol_data.get("slurm", {})
            job_id = slurm_info.get("job_id")
            if job_id and slurm_info.get("status") not in _SLURM_TERMINAL:
                pending_jobs[mol_name] = job_id

        if pending_jobs:
            click.echo(
                f"\nWaiting for {len(pending_jobs)} DFT job(s) on {dft_on_cluster}..."
            )
            final_states = poll_jobs(dft_on_cluster, pending_jobs)

            # Update manifest with final states
            completed_mols = []
            for mol_name, state in final_states.items():
                manifest["molecules"][mol_name]["slurm"]["status"] = state
                if state == "COMPLETED":
                    completed_mols.append(mol_name)
                else:
                    manifest["molecules"][mol_name]["dft_baseline"] = "dft_failed"
                    click.echo(f"  SLURM job for {mol_name} ended with state: {state}")
            save_manifest(manifest, manifest_path)

            # Retrieve results for completed molecules
            if completed_mols:
                click.echo(f"Retrieving results for {len(completed_mols)} molecule(s)...")
                retrieve_results(dft_on_cluster, completed_mols, output_dir)
                for mol_name in completed_mols:
                    manifest["molecules"][mol_name]["dft_baseline"] = STATUS_COMPLETE
                save_manifest(manifest, manifest_path)

                # Re-run harmonic analysis now that DFT baseline is available
                from .analysis import analyze_molecule_harmonic

                for mol_name in completed_mols:
                    try:
                        click.echo(f"  Re-running harmonic analysis for {mol_name}...")
                        analyze_molecule_harmonic(mol_name, base_results_dir=output_dir)
                        manifest["molecules"][mol_name]["analysis_harmonic"] = STATUS_COMPLETE
                    except Exception as e:
                        click.echo(f"  Warning: Harmonic analysis failed for {mol_name}: {e}", err=True)
                        manifest["molecules"][mol_name]["analysis_harmonic"] = STATUS_FAILED
                save_manifest(manifest, manifest_path)

    # Print summary table
    click.echo("\n" + "=" * 60)
    click.echo("BATCH SUMMARY")
    click.echo("=" * 60)
    click.echo(f"{'Molecule':<25} {'Status':<15} {'Time':>10}")
    click.echo("-" * 50)
    for mol_name, mol_data in manifest["molecules"].items():
        combos = mol_data.get("combinations", {})
        n_complete = sum(1 for c in combos.values() if c.get("status") == STATUS_COMPLETE)
        n_failed = sum(1 for c in combos.values() if c.get("status") == STATUS_FAILED)
        total_time = sum(c.get("runtime_s", 0) for c in combos.values())
        status = f"{n_complete} ok, {n_failed} fail" if n_failed > 0 else "complete"
        click.echo(f"{mol_name:<25} {status:<15} {total_time:>8.0f}s")
    click.echo("-" * 50)
    click.echo(
        f"Complete: {summary['complete']}, Failed: {summary['failed']}, "
        f"Skipped: {summary['skipped']}"
    )
    click.echo(f"Manifest: {manifest_path}")
    click.echo("=" * 60)

    manifest["updated"] = datetime.now().isoformat()
    save_manifest(manifest, manifest_path)
    return summary
