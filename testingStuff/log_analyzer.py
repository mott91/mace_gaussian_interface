import re
import numpy as np

def extract_frequencies_from_gaussian(log_file, output_file="Vib_freq.txt", verbose=True):
    """
    Extract vibrational frequencies and IR intensities from Gaussian output files.
    
    Parameters:
    -----------
    log_file : str
        Path to Gaussian output file (.log or .out)
    output_file : str
        Path to output frequency file (default: "Vib_freq.txt")
    verbose : bool
        Print extraction progress
    
    Returns:
    --------
    frequencies : list
        List of frequencies in cm^-1
    intensities : list
        List of IR intensities in KM/Mole
    """
    
    frequencies = []
    intensities = []
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find Gaussian output file: {log_file}")
    
    # Look for the harmonic frequency section
    # The section starts with "Harmonic frequencies (cm**-1), IR intensities"
    freq_section_pattern = r"Harmonic frequencies \(cm\*\*-1\), IR intensities.*?(?=\n\s*\n|\Z)"
    freq_section_match = re.search(freq_section_pattern, content, re.DOTALL)
    
    if not freq_section_match:
        raise ValueError("Could not find harmonic frequency section in Gaussian output")
    
    freq_section = freq_section_match.group(0)
    
    # Extract frequencies and intensities
    # Pattern matches lines like: "Frequencies --     68.0098               398.2991               474.9680"
    freq_lines = re.findall(r'Frequencies\s+--\s+([\d\.\s\-]+)', freq_section)
    inten_lines = re.findall(r'IR Inten\s+--\s+([\d\.\s\-]+)', freq_section)
    
    if verbose:
        print(f"Found {len(freq_lines)} frequency blocks and {len(inten_lines)} intensity blocks")
    
    # Parse each block (each line can contain multiple frequencies/intensities)
    for freq_line, inten_line in zip(freq_lines, inten_lines):
        # Split and convert to float
        freqs_in_line = [float(f) for f in freq_line.split()]
        intens_in_line = [float(i) for i in inten_line.split()]
        
        frequencies.extend(freqs_in_line)
        intensities.extend(intens_in_line)
    
    if len(frequencies) != len(intensities):
        raise ValueError(f"Mismatch: {len(frequencies)} frequencies vs {len(intensities)} intensities")
    
    if verbose:
        print(f"Successfully extracted {len(frequencies)} vibrational modes")
        print(f"Frequency range: {min(frequencies):.1f} - {max(frequencies):.1f} cm^-1")
        print(f"Intensity range: {min(intensities):.4f} - {max(intensities):.4f} KM/Mole")
    
    # Convert to format expected by plotting function
    # Mode  meV  cm^-1  Intensity
    write_vib_freq_file(frequencies, intensities, output_file, verbose)
    
    return frequencies, intensities


def write_vib_freq_file(frequencies, intensities, output_file="Vib_freq.txt", verbose=True):
    """
    Write frequencies and intensities to file in the format expected by the plotting function.
    
    Parameters:
    -----------
    frequencies : list
        Frequencies in cm^-1
    intensities : list
        IR intensities in KM/Mole
    output_file : str
        Output filename
    verbose : bool
        Print writing progress
    """
    
    # Conversion factor: 1 cm^-1 = 0.124 meV (approximately)
    cm_to_mev = 0.124
    
    with open(output_file, 'w') as f:
        # Write header (similar to what ASE vibrations produces)
        f.write("-------------------------------------\n")
        f.write(" Mode    Frequency        Intensity\n") 
        f.write("  #    meV     cm^-1   (KM/Mole)\n")
        f.write("-------------------------------------\n")
        
        for i, (freq, inten) in enumerate(zip(frequencies, intensities)):
            # Convert cm^-1 to meV
            freq_mev = freq * cm_to_mev
            
            # Handle imaginary frequencies (negative frequencies)
            if freq < 0:
                freq_str = f"{abs(freq):8.1f}i"
                mev_str = f"{abs(freq_mev):6.1f}i"
            else:
                freq_str = f"{freq:8.1f}"
                mev_str = f"{freq_mev:6.1f}"
                
            f.write(f"{i:3d}  {mev_str:>8s}  {freq_str:>8s}    {inten:8.4f}\n")
    
    if verbose:
        print(f"Vibrational frequencies written to: {output_file}")


def extract_frequencies_advanced(log_file, output_file="Vib_freq.txt", 
                                anharmonic=False, verbose=True):
    """
    Advanced extractor that can handle both harmonic and anharmonic frequencies.
    
    Parameters:
    -----------
    log_file : str
        Path to Gaussian output file
    output_file : str
        Output filename
    anharmonic : bool
        If True, extract anharmonic frequencies (if available)
    verbose : bool
        Print extraction progress
    
    Returns:
    --------
    frequencies : list
        Extracted frequencies in cm^-1
    intensities : list
        Extracted IR intensities
    """
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find file: {log_file}")
    
    frequencies = []
    intensities = []
    
    if anharmonic:
        # Try to extract anharmonic frequencies first
        if verbose:
            print("Searching for anharmonic frequency data...")
            
        # Look for anharmonic fundamental bands section
        anharm_pattern = r"Fundamental Bands\s*\-+\s*Mode\(n\)\s+E\(harm\)\s+E\(anharm\)\s+I\(harm\)\s+I\(anharm\).*?(?=\n\s*\n|\Z)"
        anharm_match = re.search(anharm_pattern, content, re.DOTALL)
        
        if anharm_match:
            anharm_section = anharm_match.group(0)
            # Extract anharmonic data
            lines = anharm_section.split('\n')
            for line in lines:
                # Look for lines with frequency data
                if re.match(r'\s*\d+\(1\)', line):
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            freq_anharm = float(parts[2])  # E(anharm)
                            inten_anharm = float(parts[4])  # I(anharm)
                            frequencies.append(freq_anharm)
                            intensities.append(inten_anharm)
                        except (ValueError, IndexError):
                            continue
            
            if frequencies and verbose:
                print(f"Found {len(frequencies)} anharmonic frequencies")
        
        if not frequencies and verbose:
            print("No anharmonic data found, falling back to harmonic frequencies")
    
    # If no anharmonic data found or not requested, extract harmonic
    if not frequencies:
        frequencies, intensities = extract_frequencies_from_gaussian(
            log_file, output_file, verbose=False)
    
    # Write output file
    write_vib_freq_file(frequencies, intensities, output_file, verbose)
    
    return frequencies, intensities


def compare_frequencies(gaussian_file, ml_file="Vib_freq.txt", tolerance=50):
    """
    Compare frequencies from Gaussian and ML calculations.
    
    Parameters:
    -----------
    gaussian_file : str
        Gaussian output file
    ml_file : str
        ML frequency file (from ASE vibrations)
    tolerance : float
        Matching tolerance in cm^-1
    """
    
    # Extract Gaussian frequencies
    print("Extracting Gaussian frequencies...")
    gauss_freq, gauss_inten = extract_frequencies_from_gaussian(gaussian_file, "gauss_freq.txt")
    
    # Read ML frequencies
    print("Reading ML frequencies...")
    ml_freq = []
    ml_inten = []
    
    try:
        with open(ml_file, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            parts = line.split()
            if len(parts) >= 4 and not line.startswith('-') and 'Mode' not in line:
                try:
                    freq = float(parts[2].replace('i', ''))
                    inten = float(parts[3])
                    if freq > 0:  # Only positive frequencies
                        ml_freq.append(freq)
                        ml_inten.append(inten)
                except (ValueError, IndexError):
                    continue
                    
    except FileNotFoundError:
        print(f"ML frequency file {ml_file} not found")
        return
    
    # Compare frequencies
    print(f"\nComparison (tolerance: {tolerance} cm^-1):")
    print("Gaussian (cm^-1)  |  ML (cm^-1)  |  Difference")
    print("-" * 45)
    
    matched = 0
    for gf, gi in zip(gauss_freq, gauss_inten):
        if gf > 0:  # Only positive frequencies
            # Find closest ML frequency
            diffs = [abs(gf - mf) for mf in ml_freq]
            if diffs:
                min_diff_idx = np.argmin(diffs)
                min_diff = diffs[min_diff_idx]
                closest_ml = ml_freq[min_diff_idx]
                
                if min_diff <= tolerance:
                    print(f"{gf:>10.1f}      |  {closest_ml:>8.1f}   |  {gf-closest_ml:>8.1f}")
                    matched += 1
                else:
                    print(f"{gf:>10.1f}      |    ---     |    ---")
    
    print(f"\nMatched: {matched}/{len([f for f in gauss_freq if f > 0])} frequencies")


# Example usage and testing
if __name__ == "__main__":
    
    # Example 1: Basic extraction
    try:
        frequencies, intensities = extract_frequencies_from_gaussian(
            "acoh_freq_anharm.log", 
            "gaussian_vib_freq.txt"
        )
        print(f"Extracted {len(frequencies)} modes")
    except Exception as e:
        print(f"Error: {e}")
    
    # Example 2: Advanced extraction with anharmonic data
    try:
        frequencies, intensities = extract_frequencies_advanced(
            "acoh_freq_anharm.log",
            "vib_freq.txt", 
            anharmonic=True
        )
    except Exception as e:
        print(f"Error: {e}")
    
    # Example 3: Compare with ML results
    # compare_frequencies("acoh_freq_anharm.log", "Vib_freq.txt")
