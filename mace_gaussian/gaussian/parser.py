"""
Gaussian log file parser for extracting frequencies and IR intensities.

This module lives at gaussian/parser.py (relocated from top-level gaussian_parser.py).
The top-level gaussian_parser.py is kept intact for now; callers are updated in plan 05.
"""

import logging
import re
from pathlib import Path
from typing import Optional

from ..utils.exceptions import GaussianParseError

logger = logging.getLogger(__name__)

# A VPT2 state label: ``3(1)`` for asymmetric tops, ``3(1,+1)`` for linear molecules and
# symmetric tops, where the second number is the vibrational angular momentum l of a
# degenerate mode. Groups: mode number, quanta, optional l.
_STATE = r"(\d+)\((\d+)(?:,([+-]?\d+))?\)"


class GaussianLogParser:
    """Parser for Gaussian 16 log files."""

    def __init__(self, log_file: str):
        """
        Initialize parser with log file path.

        Parameters
        ----------
        log_file : str
            Path to Gaussian log file
        """
        self.log_file = Path(log_file)
        if not self.log_file.exists():
            raise FileNotFoundError(f"Log file not found: {log_file}")

        with self.log_file.open() as f:
            self.content = f.read()

        self._components: dict[tuple[int, int], int] | None = None

    def _degenerate_components(self) -> dict[tuple[int, int], int]:
        """Number the fundamental components of a linear or symmetric-top VPT2 run.

        Gaussian labels a degenerate pair once, ``3(1,-1)`` and ``3(1,+1)``, while the
        harmonic table and the .fchk list it as two modes. Each component gets its own
        mode number (order of appearance) so that downstream code, which keys everything
        by ``mode``, sees 3N-5/3N-6 distinct fundamentals. Empty for asymmetric tops,
        whose labels carry no l.
        """
        if self._components is None:
            self._components = {}
            start = self.content.find("Fundamental Bands")
            if start >= 0:
                for line in self.content[start:].split("\n")[3:]:
                    m = re.match(r"^\s*(?:[HL]\s+)?" + _STATE, line)
                    if not m:
                        break
                    if m.group(3) is not None and m.group(2) == "1":
                        key = (int(m.group(1)), int(m.group(3)))
                        self._components.setdefault(key, len(self._components) + 1)
        return self._components

    def _fundamental_number(self, n: int, ell: str | None) -> int:
        """Mode number of fundamental component ``n(1,l)``; ``n`` itself without l."""
        if ell is None:
            return n
        return self._degenerate_components().get((n, int(ell)), n)

    def _fundamental_label(self, match: re.Match) -> dict:
        """``mode`` for a fundamental row; symmetric-top rows also keep Gaussian's n and l."""
        n, ell = int(match.group(1)), match.group(3)
        if ell is None:
            return {"mode": n}
        return {"mode": self._fundamental_number(n, ell), "mode_gaussian": n, "l": int(ell)}

    def _is_degenerate(self, n: int) -> bool:
        return sum(1 for mode, _ in self._degenerate_components() if mode == n) > 1

    def _base_number(self, n: int) -> int:
        """Lowest component number of mode ``n`` (the overtone's label)."""
        numbers = [k for (mode, _), k in self._degenerate_components().items() if mode == n]
        return min(numbers) if numbers else n

    def parse_harmonic_frequencies(self) -> list[dict[str, float]]:
        """
        Parse harmonic frequencies and IR intensities.

        Deduplicates repeated frequency blocks (Gaussian prints frequencies twice
        for anharmonic calculations) while preserving degenerate modes within a
        single block that share identical (freq, intensity) values.

        Returns
        -------
        list of dict
            List of dictionaries with 'freq_cm' and 'ir_intensity' keys
        """
        # Frequencies may be negative (imaginary modes printed as -xxx.x)
        freq_pattern = r"Frequencies\s+--\s+([-\d\.\s]+)"
        ir_pattern = r"IR Inten\s+--\s+([\d\.\s]+)"

        # Collect all frequency blocks (each "Frequencies --" line is one block)
        all_blocks: list[list[tuple[float, float]]] = []
        lines = self.content.split("\n")

        i = 0
        while i < len(lines):
            line = lines[i]

            if "Frequencies --" in line:
                freq_match = re.search(freq_pattern, line)
                if freq_match:
                    freqs = [float(x) for x in freq_match.group(1).split()]

                    for j in range(i + 1, min(i + 10, len(lines))):
                        if "IR Inten" in lines[j]:
                            ir_match = re.search(ir_pattern, lines[j])
                            if ir_match:
                                intensities = [float(x) for x in ir_match.group(1).split()]
                                block = list(zip(freqs, intensities))
                                all_blocks.append(block)
                            break

            i += 1

        # Deduplicate repeated blocks while preserving degenerate modes.
        # A block is identified by its rounded (freq, intensity) tuples.
        # Identical blocks from repeated Gaussian output sections are dropped.
        seen_block_keys: set[tuple[tuple[float, float], ...]] = set()
        frequencies: list[dict[str, float]] = []

        for block in all_blocks:
            block_key = tuple((round(f, 4), round(it, 4)) for f, it in block)
            if block_key not in seen_block_keys:
                seen_block_keys.add(block_key)
                for freq, intensity in block:
                    frequencies.append({"freq_cm": freq, "ir_intensity": intensity})

        if not frequencies:
            raise GaussianParseError(
                f"No harmonic frequencies found in {self.log_file}. "
                "Check that the calculation completed normally "
                "with a frequency calculation."
            )

        logger.info(f"Parsed {len(frequencies)} harmonic frequencies")
        return frequencies

    def parse_anharmonic_frequencies(self, strict: bool = False) -> list[dict[str, float]]:
        """
        Parse anharmonic frequencies and IR intensities (Fundamental Bands only).

        Returns
        -------
        list of dict
            List of dictionaries with 'mode', 'freq_cm', 'ir_intensity', and 'freq_harmonic' keys
        """
        frequencies = []

        # Look for the section with IR intensities
        # Format A (DFT logs — "Anharmonic Infrared Spectroscopy"):
        #   Mode(n)  E(harm)  E(anharm)  I(harm)  I(anharm)
        #      1(1)  3764.146  3579.741           625.83031627
        #
        # Format B (ML external calc logs — "Vibrational Energies at Anharmonic Level"):
        #   Mode(n)  Status  E(harm)  E(anharm)  Aa(x)  Ba(y)  Ca(z)
        #      1(1)  active  3796.914  3544.674  ...
        #   H  4(1)  active  1828.929  1812.980  ...   <- H/L prefix = high/low overlap

        lines = self.content.split("\n")
        in_fundamental_section = False
        in_format_b = False

        for i, line in enumerate(lines):
            # Check if we're in the Fundamental Bands section.
            # Only the "Anharmonic Infrared Spectroscopy" table (I(anharm), km/mol)
            # is accepted; the later "Dipole strengths" table (DS(anharm),
            # 10^-40 esu^2.cm^2) is a different physical quantity and must NOT
            # overwrite the intensities.
            if "Fundamental Bands" in line:
                # Lookahead: determine which format this section is
                for j in range(i, min(i + 10, len(lines))):
                    if "I(anharm)" in lines[j]:
                        in_fundamental_section = True
                        in_format_b = False  # Format A: has intensity column
                        break
                    if "E(anharm)" in lines[j] and "Status" in lines[j]:
                        in_fundamental_section = True
                        in_format_b = True  # Format B: no intensity column
                        break
                continue

            # Check if we're leaving the fundamental section
            if in_fundamental_section and ("Overtones" in line or "Combination Bands" in line):
                in_fundamental_section = False
                in_format_b = False
                continue

            if in_fundamental_section:
                if in_format_b:
                    # Format B: optional H/L overlap prefix, mode(1), status word, two floats
                    # Matches: "   1(1)  active  3796.914  3544.674 ..."
                    # Matches: "H  4(1)  active  1828.929  1812.980 ..."
                    # Matches: "  18(1)  active  -255.103  -250.532 ..."  (imaginary/negative modes)
                    match = re.match(
                        r"^\s*(?:[HL]\s+)?" + _STATE + r"\s+\w+\s+(-?[\d\.]+)\s+(-?[\d\.]+)",
                        line,
                    )
                    if match and match.group(2) == "1":
                        frequencies.append(
                            {
                                **self._fundamental_label(match),
                                "freq_cm": float(match.group(5)),
                                "ir_intensity": 0.0,  # Format B logs have no IR intensity column
                                "freq_harmonic": float(match.group(4)),
                            }
                        )
                else:
                    # Format A: mode(1), harm freq, anharm freq, optional harm intensity,
                    # then anharm intensity
                    # Match lines like:
                    #    1(1)                  3764.146   3579.741                    625.83031627
                    # or with I(harm) value:
                    #    1(1)                  3764.146   3579.741    653.06339135    625.83031627
                    match = re.match(
                        r"^\s*"
                        + _STATE
                        + r"\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(?:([\d\.]+)\s+)?([\d\.]+)\s*$",
                        line,
                    )
                    if match and match.group(2) == "1":
                        # Group 6 is the optional harmonic intensity
                        entry = {
                            **self._fundamental_label(match),
                            "freq_cm": float(match.group(5)),
                            "ir_intensity": float(match.group(7)),
                            "freq_harmonic": float(match.group(4)),
                        }
                        mode = entry["mode"]
                        # Update existing entry from Format B, or append new
                        existing = next((f for f in frequencies if f["mode"] == mode), None)
                        if existing:
                            existing.update(entry)
                        else:
                            frequencies.append(entry)

        if strict and not frequencies:
            raise GaussianParseError(
                f"No Fundamental Bands section found in {self.log_file}. "
                "The log file may only contain harmonic frequency data."
            )

        logger.info(f"Parsed {len(frequencies)} anharmonic frequencies")
        return frequencies

    def parse_overtones(self, strict: bool = False) -> list[dict[str, float]]:
        """
        Parse overtones frequencies and IR intensities.

        Returns
        -------
        list of dict
            List of dictionaries with 'mode', 'overtone_level', 'freq_harmonic',
            'freq_anharmonic', and 'ir_intensity' keys
        """
        overtones = []

        lines = self.content.split("\n")
        in_overtones_section = False

        for i, line in enumerate(lines):
            # Check if we're entering the Overtones section.
            # Only accept the km/mol intensity table (I(anharm)); the later
            # dipole-strength table (DS(anharm)) is a different quantity.
            if "Overtones" in line and "---" in lines[i + 1]:
                # Look ahead to see if this section has intensities
                for j in range(i, min(i + 10, len(lines))):
                    if "I(anharm)" in lines[j]:
                        in_overtones_section = True
                        break
                continue

            # Check if we're leaving the overtones section
            if in_overtones_section and "Combination Bands" in line:
                break

            if in_overtones_section:
                # Match lines like:
                #    1(2)                  7528.291   6994.185                     11.18668104
                # Pattern: mode(overtone_level), harmonic freq, anharmonic freq, intensity
                match = re.match(
                    r"^\s*" + _STATE + r"\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+([\d\.]+)\s*$", line
                )

                if match:
                    n = int(match.group(1))
                    entry = {
                        "mode": self._base_number(n),
                        "overtone_level": int(match.group(2)),
                        "freq_harmonic": float(match.group(4)),
                        "freq_anharmonic": float(match.group(5)),
                        "ir_intensity": float(match.group(6)),
                    }
                    # A degenerate mode's overtone splits by l (HCN 2v3: l = 0, +-2);
                    # l keeps those states apart.
                    if match.group(3) is not None and self._is_degenerate(n):
                        entry["l"] = int(match.group(3))
                    overtones.append(entry)

        if strict and not overtones:
            raise GaussianParseError(f"No Overtones section found in {self.log_file}.")

        logger.info(f"Parsed {len(overtones)} overtones")
        return overtones

    def parse_combination_bands(self, strict: bool = False) -> list[dict[str, float]]:
        """
        Parse combination bands frequencies and IR intensities.

        Returns
        -------
        list of dict
            List of dictionaries with 'mode1', 'mode2', 'freq_harmonic',
            'freq_anharmonic', and 'ir_intensity' keys
        """
        combination_bands = []

        lines = self.content.split("\n")
        in_combination_section = False

        for i, line in enumerate(lines):
            # Check if we're entering the Combination Bands section.
            # Only accept the km/mol intensity table (I(anharm)); the later
            # dipole-strength table (DS(anharm)) is a different quantity.
            if "Combination Bands" in line and "---" in lines[i + 1]:
                # Look ahead to see if this section has intensities
                for j in range(i, min(i + 10, len(lines))):
                    if "I(anharm)" in lines[j]:
                        in_combination_section = True
                        break
                continue

            # Check if we're leaving the combination bands section
            # (next major section, or the start of the dipole-strengths tables)
            if in_combination_section and (
                "Electric dipole :" in line
                or "Rotational Constants" in line
                or "Dipole strengths" in line
                or line.strip().startswith("==")
            ):
                break

            if in_combination_section:
                # Match lines like:
                #    2(1)        1(1)      6953.940   6650.547                      0.03741575
                # Pattern: mode1(1), mode2(1), harmonic freq, anharmonic freq, intensity
                match = re.match(
                    r"^\s*"
                    + _STATE
                    + r"\s+"
                    + _STATE
                    + r"\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+([\d\.]+)\s*$",
                    line,
                )

                if match and match.group(2) == "1" and match.group(5) == "1":
                    mode1 = self._fundamental_number(int(match.group(1)), match.group(3))
                    mode2 = self._fundamental_number(int(match.group(4)), match.group(6))
                    freq_harm = float(match.group(7))
                    freq_anharm = float(match.group(8))
                    ir_intensity = float(match.group(9))

                    combination_bands.append(
                        {
                            "mode1": mode1,
                            "mode2": mode2,
                            "freq_harmonic": freq_harm,
                            "freq_anharmonic": freq_anharm,
                            "ir_intensity": ir_intensity,
                        }
                    )

        if strict and not combination_bands:
            raise GaussianParseError(f"No Combination Bands section found in {self.log_file}.")

        logger.info(f"Parsed {len(combination_bands)} combination bands")
        return combination_bands

    def parse_final_energy(self) -> Optional[float]:
        """
        Parse final energy from log file.

        Returns
        -------
        float or None
            Final energy in Hartrees, or None if not found
        """
        # Two sources, both anchored (review finding H4). An unanchored "Energy="
        # pattern also matched "Thermal correction to Gibbs Free Energy=" in the
        # thermochemistry block and returned that instead of the energy.
        #   External runs:  " Energy=   -76.4312   NIter=  1"
        #   DFT runs:       " SCF Done:  E(RB3LYP) =  -76.4196339639     A.U. after ..."
        patterns = (
            r"^\s*Energy=\s+(-?\d+\.\d+)\s+NIter",
            r"SCF Done:\s+E\([^)]*\)\s+=\s+(-?\d+\.\d+)",
        )
        last = None
        for pattern in patterns:
            for m in re.finditer(pattern, self.content, re.MULTILINE):
                if last is None or m.start() > last.start():
                    last = m

        if last is not None:
            energy_hartree = float(last.group(1))
            logger.info(f"Parsed final energy: {energy_hartree} Hartree")
            return energy_hartree

        logger.warning("Could not find final energy in log file")
        return None

    def _unwrapped_archive(self) -> Optional[str]:
        """Return the Gaussian archive section with its 70-column line wrapping removed.

        The archive hard-wraps mid-token (continuation lines start with a single
        space), so regexes must run on the unwrapped text or numbers get truncated
        at line breaks.
        """
        start = self.content.find("1\\1\\GINC")
        if start == -1:
            return None
        end = self.content.find("@", start)
        if end == -1:
            end = len(self.content)
        lines = self.content[start:end].splitlines()
        return "".join(line[1:] if line.startswith(" ") else line for line in lines)

    def parse_dipole_moment(self) -> Optional[dict[str, float]]:
        """
        Parse dipole moment from log file.

        Returns
        -------
        dict or None
            Dictionary with dipole moment components and magnitude in Debye
        """
        # Look for "Dipole=" line in archive section
        # Format: Dipole=x,y,z (may have '-' for undefined components in linear molecules)
        archive_pattern = r"Dipole=([-\d\.]+),([-\d\.]+),([-\d\.]+)"

        # The archive wraps at 70 columns and can split the entry mid-number,
        # so unwrap it first; fall back to raw content for truncated logs.
        search_text = self._unwrapped_archive() or self.content
        match = re.search(archive_pattern, search_text)
        if match:
            try:
                x = float(match.group(1))
                y = float(match.group(2))
                z = float(match.group(3))
                magnitude = (x**2 + y**2 + z**2) ** 0.5

                # Convert from a.u. to Debye (1 a.u. = 2.54174 Debye)
                AU_TO_DEBYE = 2.54174623

                dipole = {
                    "x": x * AU_TO_DEBYE,
                    "y": y * AU_TO_DEBYE,
                    "z": z * AU_TO_DEBYE,
                    "magnitude": magnitude * AU_TO_DEBYE,
                }
                logger.info(f"Parsed dipole moment: {dipole['magnitude']:.4f} Debye")
                return dipole
            except ValueError:
                # Handle linear molecules where some components may be '-' or undefined
                logger.warning("Could not parse dipole moment values (may be linear molecule)")
                return None

        logger.warning("Could not find dipole moment in log file")
        return None

    def parse_timing(self) -> list[dict[str, float]]:
        """Parse wall-clock and CPU timing from Gaussian log file.

        Gaussian prints timing per job step. A DFT opt+freq log has two
        sections; an ML freq-only log has one.

        Returns
        -------
        list[dict]
            One dict per job step with keys ``cpu_s`` and ``elapsed_s``.
        """
        content = self.log_file.read_text()
        time_re = re.compile(r"(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+([\d.]+) seconds")
        sections: list[dict[str, float]] = []
        current: dict[str, float] = {}

        for line in content.splitlines():
            m = time_re.search(line)
            if not m:
                continue
            total_s = (
                int(m.group(1)) * 86400
                + int(m.group(2)) * 3600
                + int(m.group(3)) * 60
                + float(m.group(4))
            )
            if "Job cpu time" in line:
                current["cpu_s"] = total_s
            elif "Elapsed time" in line:
                current["elapsed_s"] = total_s
                sections.append(current)
                current = {}

        return sections

    def parse_timing_summary(self) -> dict[str, float]:
        """Parse timing and return a summary with total and per-stage times.

        Returns
        -------
        dict
            Keys: ``total_elapsed_s``, ``total_cpu_s``, and ``stages`` list.
        """
        stages = self.parse_timing()
        return {
            "total_elapsed_s": sum(s.get("elapsed_s", 0) for s in stages),
            "total_cpu_s": sum(s.get("cpu_s", 0) for s in stages),
            "stages": stages,
        }

    def parse_vpt2_diagnostics(self) -> dict:
        """Health indicators of the VPT2 run (report review item 4, 2026-09-18).

        Returns
        -------
        dict
            ``unreliable_cubic`` and ``unreliable_quartic``: number of
            ``WARNING: Unreliable ... force constant`` lines. Gaussian derives each
            cubic constant from several displacement pairs and warns when they
            disagree; ML surfaces trigger this more often than DFT (observation O1).
            ``fermi_resonances`` / ``darling_dennison_resonances``: True when Gaussian
            applied a deperturbation, i.e. the log does NOT contain the corresponding
            "No ... resonance found" line. ``has_anharmonic``: whether a VPT2 section
            exists at all.
        """
        c = self.content
        has_anharm = "Second-order Perturbative Anharmonic Analysis" in c
        return {
            "has_anharmonic": has_anharm,
            "unreliable_cubic": len(re.findall(r"Unreliable CUBIC force constant", c)),
            "unreliable_quartic": len(re.findall(r"Unreliable QUARTIC force constant", c)),
            "fermi_resonances": has_anharm and "No Fermi resonance found" not in c,
            "darling_dennison_resonances": has_anharm
            and not (
                "No 2-2 Darling-Dennison resonance found" in c
                and "No 1-1 Darling-Dennison resonance found" in c
            ),
        }

    def parse_x_matrix(self) -> dict[tuple[int, int], float]:
        """Total anharmonic constants x_ij (cm^-1) from the "Total Anharmonic X Matrix" block.

        Gaussian prints the lower triangle in blocks of five columns, Fortran ``D``
        exponents, indexed by its own ``Mode(n)`` numbers (the anharmonic-table numbering,
        not the ascending-frequency .fchk order). Returns ``{(i, j): x}`` with i >= j, or
        an empty dict when the log has no VPT2 section. This is the deperturbed X matrix:
        where Gaussian treats a resonance variationally, the band positions in the
        spectroscopy tables no longer follow from it alone.
        """
        start = self.content.find("Total Anharmonic X Matrix")
        if start < 0:
            return {}
        x: dict[tuple[int, int], float] = {}
        cols: list[int] = []
        for line in self.content[start:].split("\n")[2:]:
            tokens = line.split()
            if not tokens:
                break
            if all(t.isdigit() for t in tokens):
                cols = [int(t) for t in tokens]
                continue
            if not tokens[0].isdigit():
                break
            row = int(tokens[0])
            for col, value in zip(cols, tokens[1:]):
                x[(row, col)] = float(value.replace("D", "E"))
        return x

    def parse_all(self) -> dict:
        """
        Parse all available data from log file.

        Returns
        -------
        dict
            Dictionary with all parsed data including overtones and combination bands
        """
        return {
            "harmonic": self.parse_harmonic_frequencies(),
            "anharmonic": self.parse_anharmonic_frequencies(),
            "overtones": self.parse_overtones(),
            "combination_bands": self.parse_combination_bands(),
            "final_energy_hartree": self.parse_final_energy(),
            "dipole_moment": self.parse_dipole_moment(),
            "timing": self.parse_timing_summary(),
            "vpt2_diagnostics": self.parse_vpt2_diagnostics(),
        }


def parse_gaussian_log(log_file: str) -> dict:
    """
    Convenience function to parse Gaussian log file.

    Parameters
    ----------
    log_file : str
        Path to Gaussian log file

    Returns
    -------
    dict
        Dictionary with parsed data
    """
    parser = GaussianLogParser(log_file)
    return parser.parse_all()
