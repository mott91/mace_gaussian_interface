"""Create trimmed test fixture files from full comparison result logs.

This script carefully extracts the sections needed for parser testing while
preserving exact line formatting (column alignment, spaces) so the parser
produces identical results.
"""

import shutil
import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
COMP = ROOT / "comparison_results"
FIX = ROOT / "tests" / "fixtures"


def extract_sections(lines, sections_spec):
    """Extract specified line ranges from file, preserving exact formatting.

    sections_spec: list of (start_line, end_line) tuples (1-indexed, inclusive)
    """
    result = []
    for start, end in sorted(sections_spec):
        # Convert to 0-indexed
        for i in range(start - 1, min(end, len(lines))):
            result.append(lines[i])
    return result


def find_line(lines, pattern, start=0):
    """Find first line containing pattern, starting from given index."""
    for i in range(start, len(lines)):
        if pattern in lines[i]:
            return i
    return -1


def find_line_from_end(lines, pattern):
    """Find last line containing pattern."""
    for i in range(len(lines) - 1, -1, -1):
        if pattern in lines[i]:
            return i
    return -1


def trim_gaussian_log(input_path, output_path, include_energy_archive=True):
    """Trim a Gaussian log file keeping only parser-relevant sections.

    Keeps:
    - First 5 header lines
    - 3 lines context before each "Frequencies --" block through displacement lines
    - The "Anharmonic Infrared Spectroscopy" section (both I and DS blocks)
    - Energy= archive lines (one representative)
    - Dipole= archive line
    - Normal termination line
    """
    with open(input_path) as f:
        lines = f.readlines()

    sections = []

    # 1. Header (first 5 lines)
    sections.append((1, 5))

    # 2. Find all "Frequencies --" blocks (with context)
    i = 0
    while i < len(lines):
        if "Frequencies --" in lines[i]:
            # Include 3 lines before (symmetry labels, mode numbers, etc.)
            start = max(0, i - 3)
            # Find end of this block: look for Atom AN + displacement lines + blank
            block_end = i
            found_atom_an = False
            for j in range(i, min(i + 50, len(lines))):
                if "Atom  AN" in lines[j]:
                    found_atom_an = True
                if found_atom_an and lines[j].strip() == "":
                    block_end = j
                    break
                block_end = j
            sections.append((start + 1, block_end + 1))  # 1-indexed
            i = block_end + 1
        else:
            i += 1

    # 3. Thermochemistry section after frequency blocks (for Energy= parsing)
    for i_line in range(len(lines)):
        if "Thermal correction to Energy=" in lines[i_line]:
            # Include from 2 lines before to 2 lines after Gibbs
            start = max(0, i_line - 2)
            for j in range(i_line, min(i_line + 10, len(lines))):
                if "Gibbs Free Energy=" in lines[j]:
                    sections.append((start + 1, j + 2))  # 1-indexed
                    break
            break  # Only need first occurrence

    # 4. Anharmonic Infrared Spectroscopy section (the full thing)
    anharm_start = find_line(lines, "Anharmonic Infrared Spectroscopy")
    if anharm_start >= 0:
        # Go back to find "====" separator
        sep_start = anharm_start
        for j in range(anharm_start, max(anharm_start - 5, -1), -1):
            if "====" in lines[j]:
                sep_start = j
                break

        # Find end: after DS(anharm) Combination Bands section
        # Look for the blank line after the last data block
        anharm_end = anharm_start
        found_ds_combo = False
        for j in range(anharm_start, len(lines)):
            if "DS(anharm)" in lines[j]:
                # This is in the DS section
                pass
            if found_ds_combo and lines[j].strip() == "":
                anharm_end = j
                break
            if "Combination Bands" in lines[j] and find_line(lines, "DS(anharm)", anharm_start) < j:
                found_ds_combo = True
            anharm_end = j

        # More robust: find the second "Combination Bands" after Anharmonic Infrared
        combo_count = 0
        for j in range(anharm_start, len(lines)):
            if "Combination Bands" in lines[j]:
                combo_count += 1
                if combo_count == 2:
                    # Find end of this combo section (next blank or === line)
                    for k in range(j + 1, len(lines)):
                        line = lines[k].strip()
                        if line == "" or line.startswith("===") or "GradGrad" in lines[k]:
                            anharm_end = k - 1
                            break
                        anharm_end = k
                    break

        sections.append((sep_start + 1, anharm_end + 1))

    # 5. Vibrational Energies at Anharmonic Level section (for acoh and others)
    vib_start = find_line(lines, "Vibrational Energies at Anharmonic Level")
    if vib_start >= 0 and anharm_start < 0:
        # Only if no Anharmonic Infrared section (like acoh)
        # Go back to find "====" separator
        sep_start = vib_start
        for j in range(vib_start, max(vib_start - 5, -1), -1):
            if "====" in lines[j]:
                sep_start = j
                break

        # Find end of this section: after Combination Bands data
        vib_end = vib_start
        for j in range(vib_start, len(lines)):
            if "====" in lines[j] and j > vib_start + 5:
                vib_end = j - 1
                break
            vib_end = j

        sections.append((sep_start + 1, vib_end + 1))

    # 6. Energy= lines (one representative from archive)
    if include_energy_archive:
        for i_line in range(len(lines)):
            if "Energy=" in lines[i_line] and "NIter=" in lines[i_line]:
                sections.append((i_line + 1, i_line + 1))
                break

    # 7. Dipole= archive line
    dipole_line = find_line_from_end(lines, "Dipole=")
    if dipole_line >= 0:
        sections.append((dipole_line + 1, dipole_line + 1))

    # 8. Normal termination (last few lines)
    term_line = find_line_from_end(lines, "Normal termination")
    if term_line >= 0:
        sections.append((term_line + 1, term_line + 1))

    # Now merge overlapping sections and extract
    sections.sort()
    merged = []
    for start, end in sections:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    result = extract_sections(lines, merged)

    with open(output_path, 'w') as f:
        f.writelines(result)

    return len(result)


def trim_acoh_log(input_path, output_path):
    """Trim acoh log file, preserving the bug-demonstrating sections.

    Key bug features:
    1. No "Anharmonic Infrared Spectroscopy" section (no I(anharm)/DS(anharm) columns)
    2. H/L prefixed lines in "Vibrational Energies at Anharmonic Level" section
    """
    with open(input_path) as f:
        lines = f.readlines()

    sections = []

    # 1. Header
    sections.append((1, 5))

    # 2. First occurrence of Frequencies -- blocks (all of them)
    i = 0
    first_freq_done = False
    while i < len(lines):
        if "Frequencies --" in lines[i]:
            start = max(0, i - 3)
            block_end = i
            found_atom_an = False
            for j in range(i, min(i + 50, len(lines))):
                if "Atom  AN" in lines[j]:
                    found_atom_an = True
                if found_atom_an and lines[j].strip() == "":
                    block_end = j
                    break
                block_end = j
            sections.append((start + 1, block_end + 1))
            i = block_end + 1
            # Only get first set of frequency blocks
            if lines[i:i+5] and any("Thermochemistry" in l or "---" in l for l in lines[i:i+5]):
                first_freq_done = True
            if first_freq_done:
                break
        else:
            i += 1

    # 3. Vibrational Energies at Anharmonic Level section (with H/L prefixed lines)
    vib_start = find_line(lines, "Vibrational Energies at Anharmonic Level")
    if vib_start >= 0:
        # Go back to find "====" separator
        sep_start = vib_start
        for j in range(vib_start, max(vib_start - 5, -1), -1):
            if "====" in lines[j]:
                sep_start = j
                break

        # Find end: look for next "====" separator or "Predicted change"
        vib_end = vib_start
        in_combo = False
        for j in range(vib_start + 1, len(lines)):
            if "Combination Bands" in lines[j]:
                in_combo = True
            if in_combo and (lines[j].strip() == "" or "Predicted change" in lines[j] or ("====" in lines[j])):
                # Check if it's just a blank between data lines
                if lines[j].strip() == "":
                    # Check next non-blank line
                    for k in range(j + 1, min(j + 3, len(lines))):
                        if lines[k].strip():
                            if any(c.isdigit() for c in lines[k][:20]):
                                continue  # Still data
                            vib_end = j - 1
                            break
                    else:
                        vib_end = j - 1
                    break
                else:
                    vib_end = j - 1
                    break
            vib_end = j

        sections.append((sep_start + 1, vib_end + 1))

    # Also include the "Vibrational Energies (cm^-1)" section if present
    vib_cm_start = find_line(lines, "Vibrational Energies (cm^-1)")
    if vib_cm_start >= 0 and vib_cm_start < vib_start:
        # Include up to the Anharmonic Zero Point Energy section
        vib_cm_end = vib_cm_start
        for j in range(vib_cm_start + 1, len(lines)):
            if "====" in lines[j]:
                vib_cm_end = j - 1
                break
            vib_cm_end = j
        sections.append((vib_cm_start + 1, vib_cm_end + 1))

    # 4. One Energy= line
    for i_line in range(len(lines)):
        if "Energy=" in lines[i_line] and "NIter=" in lines[i_line]:
            sections.append((i_line + 1, i_line + 1))
            break

    # 5. Normal termination
    term_line = find_line_from_end(lines, "Normal termination")
    if term_line >= 0:
        sections.append((term_line + 1, term_line + 1))

    # Merge and extract
    sections.sort()
    merged = []
    for start, end in sections:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    result = extract_sections(lines, merged)

    with open(output_path, 'w') as f:
        f.writelines(result)

    return len(result)


def main():
    # Create directories
    for d in ["water", "CH4_ase", "acoh"]:
        (FIX / d).mkdir(parents=True, exist_ok=True)

    # ===== Water fixtures =====
    print("=== Water fixtures ===")

    # Water DFT log (trimmed)
    n = trim_gaussian_log(
        COMP / "water" / "b3lyp_6-31Gdp" / "gaussian_dft.log",
        FIX / "water" / "dft_b3lyp.log"
    )
    print(f"  dft_b3lyp.log: {n} lines")

    # Water DFT fchk (full copy)
    shutil.copy2(
        COMP / "water" / "b3lyp_6-31Gdp" / "gaussian_dft.fchk",
        FIX / "water" / "dft_b3lyp.fchk"
    )
    print(f"  dft_b3lyp.fchk: copied")

    # Water ML log (trimmed)
    n = trim_gaussian_log(
        COMP / "water" / "mace_mp_espaloma" / "gaussian_freq.log",
        FIX / "water" / "ml_mace_mp_esp.log"
    )
    print(f"  ml_mace_mp_esp.log: {n} lines")

    # Water ML fchk (full copy)
    shutil.copy2(
        COMP / "water" / "mace_mp_espaloma" / "gaussian_freq.fchk",
        FIX / "water" / "ml_mace_mp_esp.fchk"
    )
    print(f"  ml_mace_mp_esp.fchk: copied")

    # Water results.json (copy)
    shutil.copy2(
        COMP / "water" / "mace_mp_espaloma" / "results.json",
        FIX / "water" / "results.json"
    )
    print(f"  results.json: copied")

    # ===== CH4 fixtures =====
    print("\n=== CH4 fixtures ===")

    # CH4 DFT log (trimmed)
    n = trim_gaussian_log(
        COMP / "CH4_ase" / "b3lyp_6-31Gdp" / "gaussian_dft.log",
        FIX / "CH4_ase" / "dft_b3lyp.log"
    )
    print(f"  dft_b3lyp.log: {n} lines")

    # CH4 DFT fchk (full copy)
    shutil.copy2(
        COMP / "CH4_ase" / "b3lyp_6-31Gdp" / "gaussian_dft.fchk",
        FIX / "CH4_ase" / "dft_b3lyp.fchk"
    )
    print(f"  dft_b3lyp.fchk: copied")

    # CH4 ML log (trimmed)
    n = trim_gaussian_log(
        COMP / "CH4_ase" / "mace_mp_espaloma" / "gaussian_freq.log",
        FIX / "CH4_ase" / "ml_mace_mp_esp.log"
    )
    print(f"  ml_mace_mp_esp.log: {n} lines")

    # CH4 ML fchk (full copy)
    shutil.copy2(
        COMP / "CH4_ase" / "mace_mp_espaloma" / "gaussian_freq.fchk",
        FIX / "CH4_ase" / "ml_mace_mp_esp.fchk"
    )
    print(f"  ml_mace_mp_esp.fchk: copied")

    # CH4 results.json (copy)
    shutil.copy2(
        COMP / "CH4_ase" / "mace_mp_espaloma" / "results.json",
        FIX / "CH4_ase" / "results.json"
    )
    print(f"  results.json: copied")

    # ===== Acoh fixtures =====
    print("\n=== Acoh fixtures ===")

    # Acoh ML log (trimmed)
    n = trim_acoh_log(
        COMP / "acoh" / "mace_mp_espaloma" / "gaussian_freq.log",
        FIX / "acoh" / "ml_mace_mp_esp.log"
    )
    print(f"  ml_mace_mp_esp.log: {n} lines")

    print("\nDone!")


if __name__ == "__main__":
    main()
