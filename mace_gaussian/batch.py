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
            opt_geom_path = results_mgr.get_optimized_geometry_path(molecule_name)

            if opt_geom_path.exists() and mol_manifest.get("geometry_opt") == STATUS_COMPLETE:
                optimized_atoms = read(str(opt_geom_path))
            else:
                from .workflow import run_geometry_optimization

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
                    # SLURM DFT: handled after all molecules' geometry opt
                    pass
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

    # SLURM DFT offload: submit all, poll, retrieve (per D-02, D-07)
    if dft_on_cluster and not skip_dft_baseline:
        from .gaussian.io import ase_to_gjf
        from .slurm import TERMINAL_STATES as _SLURM_TERMINAL
        from .slurm import poll_jobs, retrieve_results, submit_dft_jobs

        # Collect molecules needing DFT
        need_dft = []
        for mol_name, mol_data in manifest["molecules"].items():
            if mol_data.get("dft_baseline") not in (STATUS_COMPLETE, "dft_failed"):
                slurm_info = mol_data.get("slurm", {})
                if slurm_info.get("status") not in _SLURM_TERMINAL:
                    need_dft.append(mol_name)

        if need_dft:
            # Resolve template path
            template = (
                Path(slurm_template)
                if slurm_template
                else Path(__file__).parent.parent / "templates" / "slurm_dft.sh"
            )

            # Generate .gjf files for molecules needing DFT
            molecules_for_slurm = []
            for mol_name in need_dft:
                mol_data = manifest["molecules"][mol_name]
                opt_geom_path = results_mgr.get_optimized_geometry_path(mol_name)
                if opt_geom_path.exists():
                    atoms = read(str(opt_geom_path))
                    gjf_dir = Path(output_dir) / mol_name / "b3lyp_6-31Gdp"
                    gjf_dir.mkdir(parents=True, exist_ok=True)
                    gjf_path = gjf_dir / f"{mol_name}_freq_anharm.gjf"
                    ase_to_gjf(atoms, str(gjf_path), molecule_name=mol_name)
                    molecules_for_slurm.append(
                        {"name": mol_name, "gjf_path": str(gjf_path)}
                    )

            if molecules_for_slurm:
                click.echo(
                    f"\nSubmitting {len(molecules_for_slurm)} DFT job(s) "
                    f"to {dft_on_cluster}..."
                )
                job_ids = submit_dft_jobs(
                    molecules_for_slurm, dft_on_cluster, template, output_dir
                )

                # Record job IDs in manifest
                for mol_name, job_id in job_ids.items():
                    manifest["molecules"][mol_name]["slurm"] = {
                        "job_id": job_id,
                        "host": dft_on_cluster,
                        "remote_dir": f"~/mace_gaussian_dft/{mol_name}",
                        "status": "SUBMITTED",
                    }
                save_manifest(manifest, manifest_path)

                # Poll until all complete
                click.echo("Polling SLURM jobs (hourly updates)...")
                final_states = poll_jobs(dft_on_cluster, job_ids)

                # Update manifest with final states
                completed_mols = []
                for mol_name, state in final_states.items():
                    manifest["molecules"][mol_name]["slurm"]["status"] = state
                    if state == "COMPLETED":
                        completed_mols.append(mol_name)
                    else:
                        # Per D-17: mark dft_failed, continue
                        manifest["molecules"][mol_name]["dft_baseline"] = "dft_failed"
                        click.echo(
                            f"  SLURM job for {mol_name} ended with state: {state}"
                        )
                save_manifest(manifest, manifest_path)

                # Retrieve results for completed molecules
                if completed_mols:
                    click.echo(
                        f"Retrieving results for {len(completed_mols)} molecule(s)..."
                    )
                    retrieve_results(dft_on_cluster, completed_mols, output_dir)
                    for mol_name in completed_mols:
                        manifest["molecules"][mol_name]["dft_baseline"] = STATUS_COMPLETE
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
