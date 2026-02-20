#!/usr/bin/env python3
"""
Batch convert all .chk files to .fchk in comparison_results directory.

Useful for converting existing calculations that don't have .fchk files yet.
"""

import logging
from pathlib import Path

from gaussian.fchk import convert_chk_to_fchk

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def find_all_chk_files(base_dir: str = "comparison_results") -> list:
    """
    Find all .chk files in the comparison_results directory.

    Parameters
    ----------
    base_dir : str
        Base directory to search (default: comparison_results)

    Returns
    -------
    chk_files : list
        List of Path objects for .chk files
    """
    base_path = Path(base_dir)
    if not base_path.exists():
        logger.error(f"Directory not found: {base_dir}")
        return []

    chk_files = list(base_path.rglob("*.chk"))
    return chk_files


def convert_all(base_dir: str = "comparison_results", force: bool = False):
    """
    Convert all .chk files to .fchk.

    Parameters
    ----------
    base_dir : str
        Base directory to search
    force : bool
        If True, reconvert even if .fchk already exists
    """
    chk_files = find_all_chk_files(base_dir)

    if not chk_files:
        logger.info(f"No .chk files found in {base_dir}")
        return

    logger.info(f"Found {len(chk_files)} .chk files")
    logger.info("="*60)

    success_count = 0
    skip_count = 0
    fail_count = 0

    for chk_file in chk_files:
        fchk_file = chk_file.with_suffix('.fchk')

        # Skip if .fchk exists and is newer (unless force=True)
        if not force and fchk_file.exists():
            if fchk_file.stat().st_mtime > chk_file.stat().st_mtime:
                logger.info(f"⏭  Skip (already exists): {chk_file.name}")
                skip_count += 1
                continue

        try:
            logger.info(f"🔄 Converting: {chk_file.relative_to(base_dir)}")
            convert_chk_to_fchk(str(chk_file), str(fchk_file))
            logger.info(f"✓  Success: {fchk_file.name}")
            success_count += 1
        except Exception as e:
            logger.error(f"✗  Failed: {chk_file.name}")
            logger.error(f"   Error: {e}")
            fail_count += 1

    logger.info("="*60)
    logger.info("SUMMARY:")
    logger.info(f"  Converted: {success_count}")
    logger.info(f"  Skipped:   {skip_count}")
    logger.info(f"  Failed:    {fail_count}")
    logger.info(f"  Total:     {len(chk_files)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Batch convert all .chk files to .fchk"
    )
    parser.add_argument(
        "--dir",
        default="comparison_results",
        help="Base directory to search (default: comparison_results)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force reconversion even if .fchk exists"
    )

    args = parser.parse_args()

    convert_all(args.dir, args.force)
