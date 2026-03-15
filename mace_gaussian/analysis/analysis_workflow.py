"""
Comprehensive Analysis Workflow

Changes:
- PNG instead of PDF for HTML display
- Combined plots with all ML methods vs B3LYP DFT baseline
- Filters B3LYP DFT as the reference baseline
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

from .analyze_spectra import SpectrumAnalyzer, SpectrumData
from .mode_matching import (
    create_alignment_matrix,
    extract_mode_data_from_checkpoint,
    match_modes,
    plot_mode_overlap_heatmap,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ComparisonWorkflow:
    """Orchestrates full comparison analysis workflow"""

    def __init__(
        self,
        molecule_name: str,
        base_results_dir: str = "comparison_results",
        output_dir: str = "analysis_results",
        use_harmonic: bool = False,
    ):
        """
        Initialize workflow

        Parameters
        ----------
        molecule_name : str
            Name of molecule to analyze
        base_results_dir : str
            Base directory containing comparison_results
        output_dir : str
            Directory for analysis outputs
        use_harmonic : bool
            If True, analyze only harmonic fundamental frequencies.
            Creates separate output directory with '_harmonic' suffix.
        """
        self.molecule_name = molecule_name
        self.base_results_dir = Path(base_results_dir)
        self.molecule_dir = self.base_results_dir / molecule_name
        self.use_harmonic = use_harmonic

        # Create output directories (separate folder for harmonic analysis)
        if use_harmonic:
            self.output_dir = Path(output_dir + "_harmonic") / molecule_name
        else:
            self.output_dir = Path(output_dir) / molecule_name

        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"

        for dir_path in [self.output_dir, self.plots_dir, self.data_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        # Initialize analyzer
        self.analyzer = SpectrumAnalyzer(freq_range=(400, 4000), bandwidth_fwhm=8.0, freq_step=0.5)

        mode_str = "harmonic" if use_harmonic else "anharmonic"
        logger.info(f"Initialized workflow for {molecule_name} ({mode_str} mode)")

    def find_dft_baseline(self, prefer_b3lyp: bool = True) -> Optional[Path]:
        """
        Find DFT anharmonic baseline results

        Scans all subdirectories for results.json files with calculator_type='dft'
        Prefers B3LYP if multiple DFT results exist

        Parameters
        ----------
        prefer_b3lyp : bool
            If True, prefer B3LYP DFT method when multiple DFT results exist

        Returns
        -------
        Path or None
            Path to DFT baseline results.json
        """
        dft_results = []

        # Scan all subdirectories for DFT results
        for item in self.molecule_dir.iterdir():
            if item.is_dir() and item.name != "geometry_opt":
                json_path = item / "results.json"
                if json_path.exists():
                    try:
                        with json_path.open() as f:
                            data = json.load(f)
                            if data.get("calculator_type") == "dft":
                                dft_results.append((item.name, json_path))
                    except (json.JSONDecodeError, KeyError) as e:
                        logger.warning(f"Error reading {json_path}: {e}")
                        continue

        if not dft_results:
            logger.warning("No DFT anharmonic baseline found!")
            return None

        # Prefer B3LYP if requested and available
        if prefer_b3lyp:
            b3lyp_results = [(name, path) for name, path in dft_results if "b3lyp" in name.lower()]
            if b3lyp_results:
                # If multiple B3LYP results, prefer one with .fchk files
                results_with_fchk = []
                results_without_fchk = []

                for name, path in b3lyp_results:
                    dft_dir = path.parent
                    fchk_files = list(dft_dir.glob("*.fchk"))
                    if fchk_files:
                        results_with_fchk.append((name, path))
                    else:
                        results_without_fchk.append((name, path))

                # Prefer results with .fchk files
                if results_with_fchk:
                    # If multiple with .fchk, use newest
                    best_path = max(results_with_fchk, key=lambda x: x[1].stat().st_mtime)[1]
                    logger.info(f"Found DFT baseline (B3LYP with .fchk): {best_path}")
                    return best_path
                else:
                    # No .fchk files, use newest
                    best_path = max(b3lyp_results, key=lambda x: x[1].stat().st_mtime)[1]
                    logger.info(f"Found DFT baseline (B3LYP, no .fchk): {best_path}")
                    return best_path

        # Otherwise return newest DFT result
        best_path = max(dft_results, key=lambda x: x[1].stat().st_mtime)[1]
        logger.info(f"Found DFT baseline: {best_path}")
        return best_path

    def find_ml_results(self) -> list[tuple[str, Path]]:
        """
        Find all ML calculation results

        Returns
        -------
        list of (name, path) tuples
            List of ML calculator combinations and their result paths
        """
        ml_results = []

        # Scan molecule directory for ML results
        for item in self.molecule_dir.iterdir():
            if item.is_dir() and item.name != "geometry_opt":
                json_path = item / "results.json"
                if json_path.exists():
                    try:
                        with json_path.open() as f:
                            data = json.load(f)
                            # Check if it's an ML calculation
                            if data.get("calculator_type") == "ml":
                                name = f"{data['energy_calculator']}_{data['dipole_calculator']}"
                                ml_results.append((name, json_path))
                                logger.info(f"Found ML result: {name}")
                    except (json.JSONDecodeError, KeyError) as e:
                        logger.warning(f"Error reading {json_path}: {e}")
                        continue

        return ml_results

    def _select_best_fchk(self, candidates: list[Path]) -> Optional[Path]:
        """
        Select best .fchk file from candidates.

        Prefers standard names (gaussian_freq.fchk, gaussian_dft.fchk),
        excludes temporary files, and falls back to newest file.

        Parameters
        ----------
        candidates : list of Path
            List of candidate .fchk files

        Returns
        -------
        Path or None
            Best .fchk file, or None if no candidates
        """
        if not candidates:
            return None

        # Filter out temporary files
        non_temp = [f for f in candidates if not f.name.startswith("temp_")]
        if non_temp:
            candidates = non_temp

        # Prefer standard names
        for preferred_name in ["gaussian_freq.fchk", "gaussian_dft.fchk"]:
            for f in candidates:
                if f.name == preferred_name:
                    return f

        # Otherwise take the newest file
        return max(candidates, key=lambda f: f.stat().st_mtime)

    def create_comparison_table(
        self,
        dft_spectrum: SpectrumData,
        ml_spectrum: SpectrumData,
        ml_name: str,
        mode_mapping: Optional[dict[int, int]] = None,
        mode_overlaps: Optional[dict[int, float]] = None,
    ) -> pd.DataFrame:
        """
        Create detailed comparison table

        Parameters
        ----------
        dft_spectrum, ml_spectrum : SpectrumData
            Spectra to compare
        ml_name : str
            Name of ML calculator
        mode_mapping : dict, optional
            Mapping from ML mode index to DFT mode index
        mode_overlaps : dict, optional
            Mapping from ML mode index to eigenvector overlap value

        Returns
        -------
        pd.DataFrame
            Comparison table
        """
        # Match modes with eigenvector mapping
        dft_freq, ml_freq, dft_int, ml_int, _ = self.analyzer.match_by_mode(
            dft_spectrum, ml_spectrum, mode_mapping=mode_mapping
        )

        # Build overlap column by determining which ML mode indices were used
        overlaps = []
        if mode_mapping is not None and mode_overlaps is not None:
            # Remap ML mode IDs using mode_mapping (same logic as match_by_mode)
            ml_mode_ids_remapped = []
            original_ml_indices = []
            for _, mode_id in enumerate(ml_spectrum.mode_ids):
                if mode_id.startswith("F"):
                    ml_mode_num = int(mode_id[1:])
                    ml_idx = ml_mode_num - 1
                    if ml_idx in mode_mapping:
                        dft_idx = mode_mapping[ml_idx]
                        dft_mode_num = dft_idx + 1
                        remapped_id = f"F{dft_mode_num}"
                        ml_mode_ids_remapped.append(remapped_id)
                        original_ml_indices.append(ml_idx)
                    else:
                        ml_mode_ids_remapped.append(mode_id)
                        original_ml_indices.append(None)
                else:
                    ml_mode_ids_remapped.append(mode_id)
                    original_ml_indices.append(None)

            # Find matches (same logic as match_by_mode)
            dft_modes_set = set(dft_spectrum.mode_ids)
            ml_modes_set = set(ml_mode_ids_remapped)
            matched_modes = sorted(dft_modes_set & ml_modes_set)

            # For each matched mode, find the overlap
            ml_mode_dict = {mode_id: i for i, mode_id in enumerate(ml_mode_ids_remapped)}
            for mode_id in matched_modes:
                ml_list_idx = ml_mode_dict[mode_id]
                ml_original_idx = original_ml_indices[ml_list_idx]
                if ml_original_idx is not None and ml_original_idx in mode_overlaps:
                    overlaps.append(mode_overlaps[ml_original_idx])
                else:
                    overlaps.append(None)
        else:
            overlaps = [None] * len(dft_freq)

        # Create DataFrame
        df_data = {
            "DFT_Frequency_cm": dft_freq,
            "ML_Frequency_cm": ml_freq,
            "Freq_Difference_cm": ml_freq - dft_freq,
            "Percent_Error": 100 * (ml_freq - dft_freq) / dft_freq,
            "Mode_Overlap": overlaps,
            "DFT_Intensity_km_mol": dft_int,
            "ML_Intensity_km_mol": ml_int,
            "Intensity_Difference": ml_int - dft_int,
        }
        df = pd.DataFrame(df_data)

        # Sort by DFT frequency
        df = df.sort_values("DFT_Frequency_cm").reset_index(drop=True)

        return df

    def extract_mode_mapping(
        self, ml_path: Path, dft_path: Path
    ) -> Optional[tuple[dict[int, int], dict[int, float]]]:
        """
        Extract mode mapping between ML and DFT calculations using eigenvector matching.

        Parameters
        ----------
        ml_path : Path
            Path to ML results.json
        dft_path : Path
            Path to DFT results.json

        Returns
        -------
        tuple of (dict, dict) or None
            (mode_mapping, mode_overlaps) where:
            - mode_mapping: ML mode index -> DFT mode index
            - mode_overlaps: ML mode index -> overlap value
            Returns None if .fchk files not found
        """
        try:
            # Find .fchk files
            ml_dir = ml_path.parent
            dft_dir = dft_path.parent

            ml_fchk_candidates = list(ml_dir.glob("*.fchk"))
            dft_fchk_candidates = list(dft_dir.glob("*.fchk"))

            ml_fchk = self._select_best_fchk(ml_fchk_candidates)
            dft_fchk = self._select_best_fchk(dft_fchk_candidates)

            if not ml_fchk or not dft_fchk:
                logger.warning("  No .fchk files found for mode matching")
                return None

            logger.debug(f"  Using ML .fchk: {ml_fchk.name}")
            logger.debug(f"  Using DFT .fchk: {dft_fchk.name}")

            # Extract modes from checkpoints - USE HARMONIC MODES
            # For mode matching, we use harmonic eigenvectors since they are the
            # fundamental normal modes. Anharmonic modes include coupling/perturbations.
            modes_ml, _, _, _, n_atoms_ml = extract_mode_data_from_checkpoint(
                str(ml_fchk), force_harmonic=True
            )
            modes_dft, _, _, _, n_atoms_dft = extract_mode_data_from_checkpoint(
                str(dft_fchk), force_harmonic=True
            )

            # # COMMENTED OUT: Old code that used anharmonic modes (kept for reference)
            # # This was using whatever modes were available (anharmonic if present)
            # modes_ml, _, _, _, n_atoms_ml = extract_mode_data_from_checkpoint(str(ml_fchk))
            # modes_dft, _, _, _, n_atoms_dft = extract_mode_data_from_checkpoint(str(dft_fchk))

            if n_atoms_ml != n_atoms_dft:
                logger.warning(f"  Different number of atoms ({n_atoms_ml} vs {n_atoms_dft})")
                return None

            # Match modes using eigenvector overlap
            matches = match_modes(modes_ml, modes_dft, threshold=0.5)

            # Convert to mapping dicts (ml_idx -> dft_idx) and (ml_idx -> overlap)
            mode_mapping = {
                ml_idx: dft_idx
                for ml_idx, (dft_idx, overlap) in matches.items()
                if dft_idx is not None
            }
            mode_overlaps = {
                ml_idx: overlap
                for ml_idx, (dft_idx, overlap) in matches.items()
                if dft_idx is not None
            }

            logger.info(f"  Extracted mode mapping for {len(mode_mapping)} modes")
            return (mode_mapping, mode_overlaps)

        except Exception as e:
            logger.warning(f"  Failed to extract mode mapping: {e}")
            return None

    def run_single_comparison(self, ml_name: str, ml_path: Path, dft_path: Path) -> dict:
        """
        Run complete comparison for one ML calculator

        Parameters
        ----------
        ml_name : str
            Name of ML calculator
        ml_path : Path
            Path to ML results.json
        dft_path : Path
            Path to DFT results.json

        Returns
        -------
        dict
            Comparison results including metrics and file paths
        """
        logger.info(f"Comparing {ml_name} vs DFT...")

        # Extract mode mapping from eigenvector matching
        mapping_result = self.extract_mode_mapping(ml_path, dft_path)
        if mapping_result is not None:
            mode_mapping, mode_overlaps = mapping_result
        else:
            mode_mapping, mode_overlaps = None, None

        # Extract spectra
        # HARMONIC MODE: Extract directly from .fchk files to get ALL degenerate modes
        # ANHARMONIC MODE: Extract from results.json (collapsed degenerates)
        if self.use_harmonic:
            # Find .fchk files
            ml_dir = ml_path.parent
            dft_dir = dft_path.parent

            ml_fchk = self._select_best_fchk(list(ml_dir.glob("*.fchk")))
            dft_fchk = self._select_best_fchk(list(dft_dir.glob("*.fchk")))

            if not ml_fchk or not dft_fchk:
                logger.error("  Missing .fchk files for harmonic analysis")
                logger.error(f"  ML .fchk: {ml_fchk}")
                logger.error(f"  DFT .fchk: {dft_fchk}")
                raise FileNotFoundError("Harmonic analysis requires .fchk files")

            logger.info("  Using .fchk-based extraction for harmonic analysis")
            logger.info(f"    ML .fchk: {ml_fchk.name}")
            logger.info(f"    DFT .fchk: {dft_fchk.name}")

            # Extract spectrum directly from .fchk (includes all degenerate modes)
            ml_spectrum = self.analyzer.extract_spectrum_from_fchk(ml_fchk)
            dft_spectrum = self.analyzer.extract_spectrum_from_fchk(dft_fchk)

            # Still load JSON results for runtime info
            ml_results = self.analyzer.load_results(ml_path)
            dft_results = self.analyzer.load_results(dft_path)

        else:
            # Anharmonic mode: use JSON extraction
            ml_results = self.analyzer.load_results(ml_path)
            dft_results = self.analyzer.load_results(dft_path)

            ml_spectrum = self.analyzer.extract_spectrum_data(
                ml_results, include_overtones=True, include_combinations=True, use_harmonic=False
            )
            dft_spectrum = self.analyzer.extract_spectrum_data(
                dft_results, include_overtones=True, include_combinations=True, use_harmonic=False
            )

        logger.info(
            f"  ML peaks: {len(ml_spectrum.frequencies)}, "
            f"DFT peaks: {len(dft_spectrum.frequencies)}"
        )

        # Calculate metrics with mode mapping
        metrics = self.analyzer.calculate_metrics(
            ml_spectrum, dft_spectrum, mode_mapping=mode_mapping
        )

        # Create plots - SAVE AS PNG NOT PDF!
        spectrum_plot_path = self.plots_dir / f"spectrum_{ml_name}.png"
        regression_plot_path = self.plots_dir / f"regression_{ml_name}.png"

        self.analyzer.plot_spectra_comparison(
            ml_spectrum,
            dft_spectrum,
            ml_name,
            molecule_name=self.molecule_name,
            save_path=str(spectrum_plot_path),
        )

        self.analyzer.plot_regression(
            ml_spectrum,
            dft_spectrum,
            metrics,
            ml_name,
            molecule_name=self.molecule_name,
            save_path=str(regression_plot_path),
            mode_mapping=mode_mapping,
        )

        # Create comparison table (also needs mode mapping and overlaps)
        comparison_df = self.create_comparison_table(
            dft_spectrum,
            ml_spectrum,
            ml_name,
            mode_mapping=mode_mapping,
            mode_overlaps=mode_overlaps,
        )

        table_path = self.data_dir / f"comparison_{ml_name}.csv"
        comparison_df.to_csv(table_path, index=False, float_format="%.4f")

        # Get runtime info
        ml_runtime = ml_results.get("runtime_s", 0)
        dft_runtime = dft_results.get("runtime_s", 0)
        speedup = dft_runtime / ml_runtime if ml_runtime > 0 else 0

        return {
            "name": ml_name,
            "metrics": metrics,
            "ml_peaks": len(ml_spectrum.frequencies),
            "dft_peaks": len(dft_spectrum.frequencies),
            "ml_runtime": ml_runtime,
            "dft_runtime": dft_runtime,
            "speedup": speedup,
            "spectrum_plot": spectrum_plot_path.name,
            "regression_plot": regression_plot_path.name,
            "table_file": table_path.name,
            "comparison_df": comparison_df,
            "ml_spectrum": ml_spectrum,  # Store for combined plots
            "dft_spectrum": dft_spectrum,
            "mode_mapping": mode_mapping,  # Store mode mapping for combined plots
        }

    def create_combined_plots(self, comparisons: list[dict], dft_spectrum: SpectrumData):
        """
        Create combined plots with all ML methods vs DFT

        Parameters
        ----------
        comparisons : list
            List of comparison results from run_single_comparison
        dft_spectrum : SpectrumData
            DFT reference spectrum
        """
        if not comparisons:
            return

        # Create combined spectrum plot (400-4000 cm⁻¹)
        combined_spectrum_path = self.plots_dir / "spectrum_combined.png"
        self.analyzer.plot_combined_spectra(
            ml_spectra=[c["ml_spectrum"] for c in comparisons],
            ml_names=[c["name"] for c in comparisons],
            dft_spectrum=dft_spectrum,
            molecule_name=self.molecule_name,
            save_path=str(combined_spectrum_path),
        )

        # Create extended combined spectrum plot (400-8000 cm⁻¹, includes overtones)
        combined_spectrum_extended_path = self.plots_dir / "spectrum_combined_extended.png"
        self.analyzer.plot_combined_spectra_extended(
            ml_spectra=[c["ml_spectrum"] for c in comparisons],
            ml_names=[c["name"] for c in comparisons],
            dft_spectrum=dft_spectrum,
            molecule_name=self.molecule_name,
            save_path=str(combined_spectrum_extended_path),
        )

        # Create combined regression plot with mode mappings
        combined_regression_path = self.plots_dir / "regression_combined.png"
        self.analyzer.plot_combined_regression(
            ml_spectra=[c["ml_spectrum"] for c in comparisons],
            ml_names=[c["name"] for c in comparisons],
            dft_spectrum=dft_spectrum,
            metrics_list=[c["metrics"] for c in comparisons],
            molecule_name=self.molecule_name,
            save_path=str(combined_regression_path),
            mode_mappings=[c.get("mode_mapping") for c in comparisons],  # Pass mode mappings
        )

        logger.info("Created combined plots")

    def generate_mode_overlap_heatmaps(self, ml_results: list[tuple[str, Path]], dft_path: Path):
        """
        Generate mode overlap heatmaps for all ML calculations vs DFT baseline.

        Parameters
        ----------
        ml_results : list
            List of (ml_name, ml_results_path) tuples
        dft_path : Path
            Path to DFT baseline results.json
        """
        logger.info("Generating mode overlap heatmaps...")

        # Find DFT .fchk file from the DFT directory
        dft_dir = dft_path.parent
        dft_fchk_candidates = list(dft_dir.glob("*.fchk"))
        dft_fchk = self._select_best_fchk(dft_fchk_candidates)

        if not dft_fchk:
            logger.warning(f"No .fchk file found for DFT baseline in {dft_dir}")
            return

        logger.info(f"Using DFT .fchk: {dft_fchk.name}")

        # Generate heatmap for each ML calculation
        for ml_name, ml_path in ml_results:
            try:
                # Find ML .fchk file
                ml_dir = ml_path.parent
                ml_fchk_candidates = list(ml_dir.glob("*.fchk"))
                ml_fchk = self._select_best_fchk(ml_fchk_candidates)

                if not ml_fchk:
                    logger.warning(f"No .fchk file found for {ml_name} in {ml_dir}")
                    continue

                logger.info(f"  Using ML .fchk: {ml_fchk.name}")

                # Extract modes - USE HARMONIC MODES ONLY FOR HEATMAPS
                # Harmonic modes are the true eigenvectors of the Hessian
                # Anharmonic modes are perturbed versions and not suitable for mode overlap
                logger.info(f"  Generating heatmap for {ml_name}...")
                modes_ml, freqs_ml, *_ = extract_mode_data_from_checkpoint(
                    str(ml_fchk), force_harmonic=True
                )
                modes_dft, freqs_dft, *_ = extract_mode_data_from_checkpoint(
                    str(dft_fchk), force_harmonic=True
                )

                # Check if same number of atoms
                if modes_ml.shape[1] != modes_dft.shape[1]:
                    logger.warning(f"  Skipping {ml_name}: different number of atoms")
                    continue

                # Create alignment matrix and match modes
                alignment_matrix = create_alignment_matrix(modes_ml, modes_dft)
                matches = match_modes(modes_ml, modes_dft)

                # Generate heatmap
                dft_baseline_name = dft_dir.name
                output_file = self.plots_dir / f"mode_overlap_{ml_name}_vs_{dft_baseline_name}.png"
                plot_mode_overlap_heatmap(
                    alignment_matrix,
                    output_file=str(output_file),
                    calc_label="ML Calculation",
                    ref_label="DFT Reference",
                    freqs_calc=freqs_ml,
                    freqs_ref=freqs_dft,
                    matches=matches,
                )

                logger.info(f"  ✓ Saved heatmap: {output_file.name}")

            except Exception as e:
                logger.warning(f"  Failed to generate heatmap for {ml_name}: {e}")
                import traceback

                traceback.print_exc()

        logger.info("Mode overlap heatmaps complete")

    def run_full_analysis(self) -> dict:
        """
        Run complete analysis workflow

        Returns
        -------
        dict
            Complete analysis results for all comparisons
        """
        logger.info("=" * 60)
        logger.info(f"STARTING ANALYSIS FOR {self.molecule_name.upper()}")
        logger.info("=" * 60)

        # Find DFT baseline (prefer B3LYP)
        dft_path = self.find_dft_baseline(prefer_b3lyp=True)
        if dft_path is None:
            logger.error("=" * 60)
            logger.error(f"NO DFT BASELINE FOUND FOR {self.molecule_name.upper()}")
            logger.error("=" * 60)
            logger.error(f"Searched in: {self.molecule_dir}")
            logger.error("Could not find any directory with calculator_type='dft' in results.json")
            logger.error("=" * 60)

            # Return empty result instead of raising exception
            return {
                "molecule": self.molecule_name,
                "comparisons": [],
                "output_dir": self.output_dir,
                "error": "No DFT baseline found",
            }

        # Find all ML results
        ml_results = self.find_ml_results()
        if not ml_results:
            logger.warning("=" * 60)
            logger.warning(f"NO ML RESULTS FOUND FOR {self.molecule_name.upper()}")
            logger.warning("=" * 60)
            logger.warning(f"Searched in: {self.molecule_dir}")
            logger.warning("No directories with calculator_type='ml' found.")
            logger.warning("=" * 60)

            return {
                "molecule": self.molecule_name,
                "comparisons": [],
                "output_dir": self.output_dir,
                "error": "No ML results found",
            }

        logger.info(f"Found DFT baseline: {dft_path.parent.name}")
        logger.info(f"Found {len(ml_results)} ML calculations to compare")

        # Load DFT spectrum once for combined plots
        if self.use_harmonic:
            # Harmonic: extract from .fchk to get all degenerate modes
            dft_dir = dft_path.parent
            dft_fchk = self._select_best_fchk(list(dft_dir.glob("*.fchk")))
            if not dft_fchk:
                raise FileNotFoundError(f"No .fchk file found for DFT baseline in {dft_dir}")
            logger.info(f"Using DFT .fchk for harmonic analysis: {dft_fchk.name}")
            dft_spectrum = self.analyzer.extract_spectrum_from_fchk(dft_fchk)
        else:
            # Anharmonic: extract from JSON
            dft_results = self.analyzer.load_results(dft_path)
            dft_spectrum = self.analyzer.extract_spectrum_data(
                dft_results, include_overtones=True, include_combinations=True, use_harmonic=False
            )

        # Run comparisons
        comparisons = []
        for ml_name, ml_path in ml_results:
            try:
                result = self.run_single_comparison(ml_name, ml_path, dft_path)
                comparisons.append(result)
                logger.info(
                    f"  [OK] {ml_name}: MAE={result['metrics'].mae_freq:.2f} cm^-1, "
                    f"R^2={result['metrics'].r2_freq:.4f}"
                )
            except Exception as e:
                logger.error(f"  [FAIL] {ml_name} failed: {e}")
                import traceback

                traceback.print_exc()

        # Create combined plots
        self.create_combined_plots(comparisons, dft_spectrum)

        # Generate mode overlap heatmaps using the same DFT path
        self.generate_mode_overlap_heatmaps(ml_results, dft_path)

        # Only save summary if we have comparisons
        if comparisons:
            # Save summary metrics
            metrics_summary = {
                "molecule": self.molecule_name,
                "analysis_date": datetime.now().isoformat(),
                "num_ml_calculators": len(comparisons),
                "dft_method": dft_path.parent.name,
                "comparisons": [
                    {
                        "name": c["name"],
                        "mae_freq": c["metrics"].mae_freq,
                        "rmse_freq": c["metrics"].rmse_freq,
                        "r2_freq": c["metrics"].r2_freq,
                        "r2_intensity": c["metrics"].r2_intensity,
                        "speedup": c["speedup"],
                    }
                    for c in comparisons
                ],
            }

            with (self.data_dir / "metrics_summary.json").open("w") as f:
                json.dump(metrics_summary, f, indent=2)

        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info("=" * 60)

        return {
            "molecule": self.molecule_name,
            "comparisons": comparisons,
            "output_dir": self.output_dir,
            "has_combined_plots": True,
        }

    def generate_html_report(self, analysis_results: dict):
        """
        Generate comprehensive HTML report

        Parameters
        ----------
        analysis_results : dict
            Results from run_full_analysis()
        """
        from .html_report_generator import HTMLReportGenerator

        generator = HTMLReportGenerator(
            molecule_name=self.molecule_name, output_dir=self.output_dir
        )

        generator.generate_report(analysis_results)

        # Only show success message if we have comparisons
        if analysis_results.get("comparisons"):
            logger.info(f"HTML report generated: {self.output_dir}/report.html")
        else:
            logger.warning(f"Error report generated: {self.output_dir}/report.html")


def analyze_molecule(
    molecule_name: str,
    base_results_dir: str = "comparison_results",
    output_dir: str = "analysis_results",
    use_harmonic: bool = False,
) -> dict:
    """
    Convenience function to run complete analysis

    Parameters
    ----------
    molecule_name : str
        Name of molecule to analyze
    base_results_dir : str
        Base directory containing comparison_results
    output_dir : str
        Directory for analysis outputs
    use_harmonic : bool
        If True, analyze only harmonic fundamental frequencies

    Returns
    -------
    dict
        Analysis results
    """
    workflow = ComparisonWorkflow(
        molecule_name=molecule_name,
        base_results_dir=base_results_dir,
        output_dir=output_dir,
        use_harmonic=use_harmonic,
    )

    results = workflow.run_full_analysis()
    workflow.generate_html_report(results)

    return results


def analyze_molecule_harmonic(
    molecule_name: str,
    base_results_dir: str = "comparison_results",
    output_dir: str = "analysis_results",
) -> dict:
    """
    Convenience function to run harmonic-only analysis.

    This analyzes only harmonic fundamental frequencies with proper
    eigenvector-based mode matching. Overtones and combinations are
    excluded since they don't have proper mode matching.

    Parameters
    ----------
    molecule_name : str
        Name of molecule to analyze
    base_results_dir : str
        Base directory containing comparison_results
    output_dir : str
        Directory for analysis outputs (will create output_dir_harmonic)

    Returns
    -------
    dict
        Analysis results
    """
    return analyze_molecule(
        molecule_name=molecule_name,
        base_results_dir=base_results_dir,
        output_dir=output_dir,
        use_harmonic=True,
    )


def run_analysis_main() -> None:
    """Entry point for anharmonic IR spectral analysis (run_analysis.py shim target)."""
    import sys

    if len(sys.argv) < 2:
        print("=" * 60)
        print("IR SPECTRAL ANALYSIS FRAMEWORK")
        print("=" * 60)
        print("\nUsage:")
        print("  python run_analysis.py <molecule_name>")
        print("\nExample:")
        print("  python run_analysis.py acoh")
        print("\nThis will:")
        print("  1. Find DFT anharmonic baseline")
        print("  2. Find all ML calculation results")
        print("  3. Generate comparison plots")
        print("  4. Calculate statistical metrics")
        print("  5. Create HTML report")
        print("\n" + "=" * 60)
        sys.exit(1)

    molecule_name = sys.argv[1]

    print("\n" + "=" * 60)
    print(f"ANALYZING: {molecule_name.upper()}")
    print("=" * 60 + "\n")

    try:
        results = analyze_molecule(molecule_name)

        print("\n" + "=" * 60)
        print("SUCCESS!")
        print("=" * 60)
        print(f"\nResults directory: {results['output_dir']}")
        print(f"HTML Report: {results['output_dir']}/report.html")
        print("\nOpen the HTML report in your browser to view the analysis!")
        print("=" * 60 + "\n")

    except FileNotFoundError as e:
        print("\n" + "=" * 60)
        print("ERROR: Files not found")
        print("=" * 60)
        print(f"\n{e}")
        print("\nMake sure you have:")
        print(f"  - comparison_results/{molecule_name}/freq_anharm/results.json (DFT baseline)")
        print(f"  - comparison_results/{molecule_name}/<calculator>/results.json (ML results)")
        print("=" * 60 + "\n")
        sys.exit(1)

    except Exception as e:
        print("\n" + "=" * 60)
        print("ERROR: Analysis failed")
        print("=" * 60)
        print(f"\n{e}")
        print("\nCheck the logs above for details.")
        print("=" * 60 + "\n")
        sys.exit(1)


def run_analysis_harmonic_main() -> None:
    """Entry point for harmonic IR spectral analysis (run_analysis_harmonic.py shim target)."""
    import sys

    if len(sys.argv) < 2:
        print("=" * 60)
        print("HARMONIC-ONLY IR SPECTRAL ANALYSIS")
        print("=" * 60)
        print("\nUsage:")
        print("  python run_analysis_harmonic.py <molecule_name>")
        print("\nExample:")
        print("  python run_analysis_harmonic.py water")
        print("\nThis will:")
        print("  1. Analyze ONLY harmonic fundamental frequencies")
        print("  2. Use eigenvector dot product mode matching")
        print("  3. Generate comparison plots (fundamentals only)")
        print("  4. Create HTML report")
        print("  5. Save to analysis_results_harmonic/")
        print("\nNote: Overtones and combinations are excluded")
        print("=" * 60)
        sys.exit(1)

    molecule_name = sys.argv[1]

    print("\n" + "=" * 60)
    print(f"HARMONIC ANALYSIS: {molecule_name.upper()}")
    print("=" * 60 + "\n")

    try:
        results = analyze_molecule_harmonic(molecule_name)

        print("\n" + "=" * 60)
        print("SUCCESS!")
        print("=" * 60)
        print(f"\nResults directory: {results['output_dir']}")
        print(f"HTML Report: {results['output_dir']}/report.html")
        print("\nOpen the HTML report in your browser to view the analysis!")
        print("=" * 60 + "\n")

    except FileNotFoundError as e:
        print("\n" + "=" * 60)
        print("ERROR: Files not found")
        print("=" * 60)
        print(f"\n{e}")
        print("\nMake sure you have:")
        print(f"  - comparison_results/{molecule_name}/<dft>/results.json (DFT baseline)")
        print(f"  - comparison_results/{molecule_name}/<calculator>/results.json (ML results)")
        print("=" * 60 + "\n")
        sys.exit(1)

    except Exception as e:
        print("\n" + "=" * 60)
        print("ERROR: Analysis failed")
        print("=" * 60)
        print(f"\n{e}")
        print("\nCheck the logs above for details.")
        print("=" * 60 + "\n")
        sys.exit(1)


if __name__ == "__main__":
    run_analysis_main()
