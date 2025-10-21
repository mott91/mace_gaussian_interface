"""
Gaussian log file parser for extracting frequencies and IR intensities.
"""

import re
import logging
from typing import Dict, List, Tuple, Optional
from pathlib import Path

logger = logging.getLogger(__name__)


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
            
        with open(self.log_file, 'r') as f:
            self.content = f.read()
            
    def parse_harmonic_frequencies(self) -> List[Dict[str, float]]:
        """
        Parse harmonic frequencies and IR intensities.
        
        Returns
        -------
        list of dict
            List of dictionaries with 'freq_cm' and 'ir_intensity' keys
        """
        frequencies = []
        seen_freqs = set()  # To avoid duplicates
        
        # Look for the frequency section
        # Pattern: "Frequencies --" followed by frequencies
        # Then "IR Inten    --" followed by intensities
        freq_pattern = r'Frequencies\s+--\s+([\d\.\s]+)'
        ir_pattern = r'IR Inten\s+--\s+([\d\.\s]+)'
        
        # Find all frequency blocks
        lines = self.content.split('\n')
        
        i = 0
        while i < len(lines):
            line = lines[i]
            
            if 'Frequencies --' in line:
                # Extract frequencies from this line
                freq_match = re.search(freq_pattern, line)
                if freq_match:
                    freqs = [float(x) for x in freq_match.group(1).split()]
                    
                    # Look for corresponding IR intensities (usually a few lines below)
                    for j in range(i+1, min(i+10, len(lines))):
                        if 'IR Inten' in lines[j]:
                            ir_match = re.search(ir_pattern, lines[j])
                            if ir_match:
                                intensities = [float(x) for x in ir_match.group(1).split()]
                                
                                # Pair them up, avoiding duplicates
                                for freq, intensity in zip(freqs, intensities):
                                    freq_key = (round(freq, 4), round(intensity, 4))
                                    if freq_key not in seen_freqs:
                                        seen_freqs.add(freq_key)
                                        frequencies.append({
                                            'freq_cm': freq,
                                            'ir_intensity': intensity
                                        })
                                break
                    
            i += 1
        
        logger.info(f"Parsed {len(frequencies)} harmonic frequencies")
        return frequencies
    
    def parse_anharmonic_frequencies(self) -> List[Dict[str, float]]:
        """
        Parse anharmonic frequencies and IR intensities.

        Returns
        -------
        list of dict
            List of dictionaries with 'freq_cm' and 'ir_intensity' keys
        """
        frequencies = []
        
        # Look for the section with IR intensities
        # Format:
        # Mode(n)                  E(harm)   E(anharm)       I(harm)        I(anharm)
        #    1(1)                  3764.146   3579.741                    625.83031627
        
        lines = self.content.split('\n')
        in_fundamental_section = False
        
        for i, line in enumerate(lines):
            # Check if we're in the Fundamental Bands section with intensities
            if 'Fundamental Bands' in line:
                # Look ahead to see if this section has intensities
                for j in range(i, min(i+10, len(lines))):
                    if 'I(anharm)' in lines[j] or 'DS(anharm)' in lines[j]:
                        in_fundamental_section = True
                        break
                continue
            
            # Check if we're leaving the fundamental section
            if in_fundamental_section and ('Overtones' in line or 'Combination Bands' in line):
                break
            
            if in_fundamental_section:
                # Match lines like:
                #    1(1)                  3764.146   3579.741                    625.83031627
                # or with I(harm) value:
                #    1(1)                  3764.146   3579.741    653.06339135    625.83031627
                
                # Pattern: mode number, harmonic freq, anharmonic freq, optional harm intensity, anharm intensity
                match = re.match(r'^\s*(\d+)\(1\)\s+([\d\.]+)\s+([\d\.]+)\s+(?:([\d\.]+)\s+)?([\d\.]+)\s*$', line)
                
                if match:
                    mode = int(match.group(1))
                    freq_harm = float(match.group(2))
                    freq_anharm = float(match.group(3))
                    # Group 4 is optional harmonic intensity
                    ir_anharm = float(match.group(5))
                    
                    frequencies.append({
                        'mode': mode,
                        'freq_cm': freq_anharm,
                        'ir_intensity': ir_anharm,
                        'freq_harmonic': freq_harm
                    })
        
        logger.info(f"Parsed {len(frequencies)} anharmonic frequencies")
        return frequencies
    
    def parse_final_energy(self) -> Optional[float]:
        """
        Parse final energy from log file.
        
        Returns
        -------
        float or None
            Final energy in Hartrees, or None if not found
        """
        # Look for the last "Energy=" line from external calculation
        pattern = r'Energy=\s+([-\d\.]+)'
        
        matches = re.findall(pattern, self.content)
        if matches:
            # Return the last one (final energy)
            energy_hartree = float(matches[-1])
            logger.info(f"Parsed final energy: {energy_hartree} Hartree")
            return energy_hartree
            
        logger.warning("Could not find final energy in log file")
        return None
    
    def parse_dipole_moment(self) -> Optional[Dict[str, float]]:
        """
        Parse dipole moment from log file.
        
        Returns
        -------
        dict or None
            Dictionary with dipole moment components and magnitude in Debye
        """
        # Look for "Dipole=" line in archive section
        # Format: Dipole=x,y,z
        archive_pattern = r'Dipole=([-\d\.]+),([-\d\.]+),([-\d\.]+)'
        
        match = re.search(archive_pattern, self.content)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))
            magnitude = (x**2 + y**2 + z**2)**0.5
            
            # Convert from a.u. to Debye (1 a.u. = 2.54174 Debye)
            AU_TO_DEBYE = 2.54174623
            
            dipole = {
                'x': x * AU_TO_DEBYE,
                'y': y * AU_TO_DEBYE,
                'z': z * AU_TO_DEBYE,
                'magnitude': magnitude * AU_TO_DEBYE
            }
            logger.info(f"Parsed dipole moment: {dipole['magnitude']:.4f} Debye")
            return dipole
            
        logger.warning("Could not find dipole moment in log file")
        return None
    
    def parse_all(self) -> Dict:
        """
        Parse all available data from log file.
        
        Returns
        -------
        dict
            Dictionary with all parsed data
        """
        return {
            'harmonic': self.parse_harmonic_frequencies(),
            'anharmonic': self.parse_anharmonic_frequencies(),
            'final_energy_hartree': self.parse_final_energy(),
            'dipole_moment': self.parse_dipole_moment()
        }


def parse_gaussian_log(log_file: str) -> Dict:
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