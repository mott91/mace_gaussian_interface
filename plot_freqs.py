import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_freq_file(filename):
    frequencies = []
    intensities = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    freq = float(parts[2].replace('i', ''))
                    intensity = float(parts[3])
                    if freq > 0:
                        frequencies.append(freq)
                        intensities.append(intensity)
                except:
                    continue
    return frequencies, intensities

def make_spectrum(frequencies, intensities, freq_range, sigma=40):
    spectrum = np.zeros_like(freq_range)
    for f, I in zip(frequencies, intensities):
        spectrum += I * np.exp(-((freq_range - f) ** 2) / (2 * sigma ** 2))
    if spectrum.max() > spectrum.min():
        spectrum = (spectrum - spectrum.min()) / (spectrum.max() - spectrum.min())
    return spectrum

def write_comparison_file(all_data, output_file="frequency_comparison.txt"):
    """Write a comparison file with all frequencies and intensities"""
    with open(output_file, 'w') as f:
        f.write("IR Frequency Comparison\n")
        f.write("=" * 60 + "\n\n")
        
        # Write summary
        f.write("SUMMARY:\n")
        f.write("-" * 30 + "\n")
        for label, frequencies, intensities in all_data:
            f.write(f"{label}:\n")
            f.write(f"  Number of modes: {len(frequencies)}\n")
            f.write(f"  Frequency range: {min(frequencies):.1f} - {max(frequencies):.1f} cm-1\n")
            f.write(f"  Max intensity: {max(intensities):.4f}\n\n")
        
        # Write detailed frequencies
        f.write("\nDETAILED FREQUENCIES (cm-1):\n")
        f.write("-" * 50 + "\n\n")
        
        for label, frequencies, intensities in all_data:
            f.write(f"{label.upper()}:\n")
            f.write("Freq (cm-1)    Intensity\n")
            f.write("-" * 25 + "\n")
            
            # Sort by frequency
            sorted_data = sorted(zip(frequencies, intensities))
            for freq, intensity in sorted_data:
                f.write(f"{freq:8.1f}      {intensity:8.4f}\n")
            f.write("\n")
        
        # Find and write top 10 peaks for each dataset
        f.write("\nTOP 10 STRONGEST PEAKS:\n")
        f.write("-" * 40 + "\n\n")
        
        for label, frequencies, intensities in all_data:
            f.write(f"{label.upper()}:\n")
            f.write("Rank  Freq (cm-1)  Intensity\n")
            f.write("-" * 30 + "\n")
            
            # Sort by intensity (highest first)
            sorted_by_intensity = sorted(zip(frequencies, intensities), 
                                       key=lambda x: x[1], reverse=True)
            
            for i, (freq, intensity) in enumerate(sorted_by_intensity[:10], 1):
                f.write(f"{i:2d}    {freq:8.1f}    {intensity:8.4f}\n")
            f.write("\n")
        
        # Simple frequency matching (within 50 cm-1)
        if len(all_data) > 1:
            f.write("\nFREQUENCY MATCHING (tolerance: 50 cm-1):\n")
            f.write("-" * 50 + "\n\n")
            
            ref_label, ref_freqs, ref_intensities = all_data[0]
            f.write(f"Reference: {ref_label}\n\n")
            
            f.write(f"{'Ref Freq':>10} ")
            for label, _, _ in all_data[1:]:
                f.write(f"{label:>12} ")
            f.write("Difference\n")
            f.write("-" * (15 + 15 * len(all_data)) + "\n")
            
            for ref_freq in sorted(ref_freqs):
                matches = [f"{ref_freq:8.1f}  "]
                differences = []
                
                for label, frequencies, intensities in all_data[1:]:
                    # Find closest frequency
                    diffs = [abs(ref_freq - f) for f in frequencies]
                    min_diff_idx = np.argmin(diffs)
                    min_diff = diffs[min_diff_idx]
                    closest_freq = frequencies[min_diff_idx]
                    
                    if min_diff <= 50:  # within tolerance
                        matches.append(f"{closest_freq:8.1f}    ")
                        differences.append(f"{ref_freq - closest_freq:+6.1f}")
                    else:
                        matches.append("   ---      ")
                        differences.append("  ---  ")
                
                # Write the line
                f.write("".join(matches))
                f.write("  " + "  ".join(differences) + "\n")

    print(f"Comparison file written to: {output_file}")

# Get files from command line
files = []
labels = []
reading_labels = False

for i, arg in enumerate(sys.argv[1:]):
    if arg == '-l':
        reading_labels = True
    elif reading_labels:
        labels.append(arg)
    else:
        reading_labels = False
        files.append(arg)

if not files:
    print("Usage: python plot_freqs.py file1.txt file2.txt file3.txt -l label1 label2 label3")
    sys.exit()

# Plot setup
freq_range = np.linspace(500, 4000, 1500)
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
plt.figure(figsize=(10, 8))

# Store all data for comparison file
all_data = []

offset = 1.2
for i, filename in enumerate(files):
    print(f"Processing {filename}")
    
    frequencies, intensities = read_freq_file(filename)
    if not frequencies:
        print(f"No data in {filename}")
        continue
        
    spectrum = make_spectrum(frequencies, intensities, freq_range)
    offset_spectrum = spectrum + i * offset
    
    label = labels[i] if i < len(labels) else f"File {i+1}"
    color = colors[i % len(colors)]
    
    plt.fill_between(freq_range, i * offset, offset_spectrum, alpha=0.3, color=color)
    plt.plot(freq_range, offset_spectrum, color=color, label=label, linewidth=1.5)
    
    # Store data for comparison file
    all_data.append((label, frequencies, intensities))

plt.xlabel('Wavenumber (cm-1)', fontsize=16)
plt.ylabel('Absorbance', fontsize=16)
plt.legend()
plt.yticks([])
plt.show()

# Write comparison file
if all_data:
    write_comparison_file(all_data)
