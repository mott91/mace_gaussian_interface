#!/usr/bin/env python3
"""
Generate flowcharts using matplotlib only (no Graphviz dependency).
Creates professional-looking diagrams for the presentation.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Polygon
from matplotlib.lines import Line2D
import numpy as np


def add_box(ax, x, y, width, height, text, color='#E3F2FD', edge_color='#424242',
            text_color='black', fontsize=10, fontweight='normal', style='round'):
    """Add a styled box to the axes."""
    if style == 'round':
        box = FancyBboxPatch(
            (x - width/2, y - height/2), width, height,
            boxstyle="round,pad=0.05",
            edgecolor=edge_color,
            facecolor=color,
            linewidth=2
        )
    elif style == 'hexagon':
        # Create hexagon points
        angles = np.linspace(0, 2*np.pi, 7)
        hw = width / 2
        hh = height / 2
        points = [(x + hw * np.cos(a) * 0.9, y + hh * np.sin(a)) for a in angles]
        box = Polygon(points, edgecolor=edge_color, facecolor=color, linewidth=2)
    else:  # square
        box = mpatches.Rectangle(
            (x - width/2, y - height/2), width, height,
            edgecolor=edge_color,
            facecolor=color,
            linewidth=2
        )

    ax.add_patch(box)

    # Add text (handle multiline)
    lines = text.split('\n')
    if len(lines) == 1:
        ax.text(x, y, text, ha='center', va='center',
               fontsize=fontsize, color=text_color, weight=fontweight,
               zorder=10)
    else:
        y_offset = (len(lines) - 1) * 0.08
        for i, line in enumerate(lines):
            ax.text(x, y + y_offset - i*0.16, line, ha='center', va='center',
                   fontsize=fontsize-1, color=text_color, weight=fontweight,
                   zorder=10)


def add_arrow(ax, x1, y1, x2, y2, label='', style='->', color='#424242',
              curvature=0, label_offset=(0, 0)):
    """Add an arrow between two points."""
    if curvature == 0:
        connectionstyle = "arc3,rad=0"
    else:
        connectionstyle = f"arc3,rad={curvature}"

    # Handle different arrow styles
    if style == '-' or style == '--':
        # Simple line without arrowhead
        line = Line2D([x1, x2], [y1, y2], color=color, linewidth=2,
                     linestyle=style, zorder=5)
        ax.add_line(line)
    else:
        # Arrow with head
        arrowstyle = style + ',head_width=0.3,head_length=0.4'
        arrow = FancyArrowPatch(
            (x1, y1), (x2, y2),
            arrowstyle=arrowstyle,
            color=color,
            linewidth=2,
            connectionstyle=connectionstyle,
            zorder=5
        )
        ax.add_patch(arrow)

    if label:
        mid_x = (x1 + x2) / 2 + label_offset[0]
        mid_y = (y1 + y2) / 2 + label_offset[1]
        ax.text(mid_x, mid_y, label, fontsize=8,
               style='italic', color='#666666',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                        edgecolor='none', alpha=0.8),
               ha='center', va='center', zorder=15)


def create_main_workflow(output_file, dpi=300):
    """Create the main MACE-Gaussian workflow diagram."""
    fig, ax = plt.subplots(figsize=(10, 12), dpi=dpi)
    ax.set_xlim(-1, 11)
    ax.set_ylim(0, 14)
    ax.axis('off')

    # Title
    ax.text(5, 13.5, 'MACE-Gaussian Interface Pipeline',
           ha='center', fontsize=16, weight='bold', color='#1565C0')

    # Workflow nodes
    y = 12.5

    # Input
    add_box(ax, 5, y, 2.5, 0.6, 'XYZ Input File', '#C5E1A5', fontweight='bold')
    add_arrow(ax, 5, y - 0.3, 5, y - 0.9)

    # Geometry optimization
    y -= 1.5
    add_box(ax, 5, y, 3.5, 0.6, 'Geometry Optimization\n(MACE-MP/OMOL)', '#FFE082')
    add_arrow(ax, 5, y - 0.3, 5, y - 0.9, label='optimized\natoms')

    # Frequency calculation start
    y -= 1.5
    add_box(ax, 5, y, 3.5, 0.6, 'Frequency Calculation\n(Gaussian + ML)', '#FFE082')

    # Branch to ZMQ system
    add_arrow(ax, 3.5, y - 0.3, 2, y - 1.4, curvature=-0.3, label='geometry\nrequest')

    # ZMQ Bridge subsystem
    y_zmq = y - 2.0
    add_box(ax, 2, y_zmq, 2.8, 0.6, 'ZMQ Bridge\n(IPC)', '#CE93D8', style='hexagon')
    add_arrow(ax, 2, y_zmq - 0.35, 2, y_zmq - 0.95)

    y_zmq -= 1.2
    add_box(ax, 2, y_zmq, 3.2, 0.6, 'ML Dipole Calculator\n(Espaloma/MACE-ML)', '#90CAF9')
    add_arrow(ax, 2, y_zmq + 0.35, 2, y_zmq + 0.95)

    # Return arrow
    add_arrow(ax, 3.2, y_zmq + 1.3, 6.5, y - 0.3, curvature=0.4,
             label='dipoles +\nderivatives', color='#1976D2')

    # Continue main flow
    y -= 2.5
    add_arrow(ax, 5, y + 0.5, 5, y - 0.2)

    add_box(ax, 5, y, 3.5, 0.6, 'Anharmonic Analysis\n(Gaussian VPT2)', '#FFE082')
    add_arrow(ax, 5, y - 0.3, 5, y - 0.9, label='overtones +\ncombinations')

    # Parsing
    y -= 1.5
    add_box(ax, 5, y, 3.0, 0.6, 'Parse Results\n(JSON Export)', '#FFAB91')
    add_arrow(ax, 5, y - 0.3, 5, y - 0.9)

    # Analysis
    y -= 1.5
    add_box(ax, 5, y, 3.5, 0.6, 'Statistical Analysis\n(Regression, KDE)', '#FFAB91')
    add_arrow(ax, 5, y - 0.3, 5, y - 0.9)

    # Output
    y -= 1.5
    add_box(ax, 5, y, 3.0, 0.6, 'HTML Report\n(Plots + Tables)', '#C5E1A5', fontweight='bold')

    # Add legend
    legend_y = 1.0
    ax.text(1, legend_y, 'Legend:', fontsize=9, weight='bold')
    add_box(ax, 0.8, legend_y - 0.4, 0.8, 0.3, 'Input/Output', '#C5E1A5')
    add_box(ax, 1.8, legend_y - 0.4, 0.8, 0.3, 'Gaussian', '#FFE082')
    add_box(ax, 2.8, legend_y - 0.4, 0.8, 0.3, 'ML Engine', '#90CAF9')
    add_box(ax, 3.8, legend_y - 0.4, 0.8, 0.3, 'Analysis', '#FFAB91')

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white', dpi=dpi)
    plt.close()
    print(f"Main workflow saved: {output_file}")


def create_zmq_architecture(output_file, dpi=300):
    """Create detailed ZMQ architecture diagram."""
    fig, ax = plt.subplots(figsize=(12, 8), dpi=dpi)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Title
    ax.text(6, 7.5, 'ZMQ Inter-Process Communication Architecture',
           ha='center', fontsize=14, weight='bold', color='#1565C0')

    # Main Python process
    ax.add_patch(mpatches.FancyBboxPatch((0.3, 4.5), 3.5, 2.5,
                 boxstyle="round,pad=0.1", edgecolor='#1976D2',
                 facecolor='#E3F2FD', linewidth=2, linestyle='--'))
    ax.text(2.05, 6.8, 'gm_main.py (Python)', ha='center', fontsize=10,
           weight='bold', color='#1565C0')

    add_box(ax, 2, 6.2, 2.5, 0.5, 'Launch Gaussian\n(subprocess)', '#BBDEFB')
    add_box(ax, 2, 5.4, 2.5, 0.5, 'ZMQ Server\n(REP socket)', '#90CAF9')
    add_box(ax, 2, 4.8, 2.5, 0.5, 'Dipole Calculator', '#C5E1A5')

    # Helper script
    ax.add_patch(mpatches.FancyBboxPatch((4.5, 5.5), 2.8, 1.5,
                 boxstyle="round,pad=0.1", edgecolor='#388E3C',
                 facecolor='#E8F5E9', linewidth=2, linestyle='--'))
    ax.text(5.9, 6.9, 'gm_helper.py (Bridge)', ha='center', fontsize=10,
           weight='bold', color='#388E3C')

    add_box(ax, 5.9, 6.3, 2.2, 0.5, 'ZMQ Client\n(REQ socket)', '#A5D6A7')

    # Gaussian process
    ax.add_patch(mpatches.FancyBboxPatch((8.2, 4.5), 3.5, 2.5,
                 boxstyle="round,pad=0.1", edgecolor='#D32F2F',
                 facecolor='#FFEBEE', linewidth=2, linestyle='--'))
    ax.text(9.95, 6.8, 'Gaussian 16 (Fortran)', ha='center', fontsize=10,
           weight='bold', color='#D32F2F')

    add_box(ax, 10, 6.2, 2.5, 0.5, 'Gaussian Engine', '#FFCDD2')
    add_box(ax, 10, 5.4, 2.5, 0.5, 'External Interface', '#EF9A9A')

    # IPC file (center)
    add_box(ax, 6, 3.5, 2.0, 0.8, '.ipc_file\n(Unix Socket)', '#FFE0B2',
           style='hexagon', fontweight='bold')

    # Arrows showing communication flow
    # 1. Gaussian startup
    add_arrow(ax, 3.2, 6.2, 8.8, 6.2, label='1. spawn', curvature=0.2, color='#F44336')

    # 2. External call
    add_arrow(ax, 8.8, 5.4, 7.0, 6.1, label='2. call\nExternal=', curvature=-0.15)

    # 3. Connect to socket
    add_arrow(ax, 5.9, 5.8, 6, 4.3, label='3. connect', curvature=0)

    # 4. ZMQ communication
    add_arrow(ax, 5.2, 3.5, 3.2, 5.0, label='4. REQ-REP\nprotocol', curvature=0.2, color='#2196F3')

    # 5. Calculate
    add_arrow(ax, 2, 5.15, 2, 5.0, label='', style='-')

    # 6. Response
    add_arrow(ax, 2.8, 5.4, 5.2, 3.9, label='6. response', curvature=-0.2, color='#4CAF50')

    # 7. Return
    add_arrow(ax, 6.7, 6.3, 8.8, 5.7, label='7. return', curvature=0.15, color='#4CAF50')

    # 8. Continue
    add_arrow(ax, 10, 5.1, 10, 4.5, label='8. continue', style='-')

    # Add sequence box
    seq_text = "Sequence:\n1. Main spawns Gaussian\n2. Gaussian calls helper via External=\n3. Helper connects to ZMQ\n4. Request sent (geometry)\n5. Main calculates dipoles\n6. Response sent\n7. Helper returns to Gaussian\n8. Gaussian continues calculation"
    ax.text(1, 2.5, seq_text, fontsize=8, verticalalignment='top',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF9C4',
                    edgecolor='#F57F17', linewidth=1.5))

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white', dpi=dpi)
    plt.close()
    print(f"ZMQ architecture saved: {output_file}")


def create_module_swapping(output_file, dpi=300):
    """Create module swapping mechanism diagram."""
    fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Title
    ax.text(5, 7.5, 'Dynamic Module Swapping Architecture',
           ha='center', fontsize=14, weight='bold', color='#1565C0')

    # Two packages at top
    add_box(ax, 2.5, 6.5, 2.5, 0.8, 'mace_ML_pkg/\n(Standard MACE)', '#E1BEE7',
           fontweight='bold')
    add_box(ax, 7.5, 6.5, 2.5, 0.8, 'mace_dipole_pkg/\n(Dipole MACE)', '#B2DFDB',
           fontweight='bold')

    # sys.modules in middle
    add_box(ax, 5, 4.5, 3.5, 0.8, 'sys.modules\n(Python Import Cache)', '#FFECB3',
           style='hexagon', fontweight='bold')

    # Swapping functions
    add_box(ax, 2.5, 3.0, 2.8, 0.6, 'load_standard_mace()', '#B3E5FC')
    add_box(ax, 7.5, 3.0, 2.8, 0.6, 'load_dipole_mace()', '#B3E5FC')

    # Use cases at bottom
    add_box(ax, 2.5, 1.2, 2.5, 0.6, 'Geometry\nOptimization', '#FFCCBC')
    add_box(ax, 7.5, 1.2, 2.5, 0.6, 'Dipole\nCalculation', '#C5E1A5')

    # Arrows
    # Standard MACE flow
    add_arrow(ax, 2.5, 5.9, 2.5, 3.6, color='#9C27B0')
    add_arrow(ax, 2.5, 3.6, 3.8, 4.4, curvature=-0.2, color='#9C27B0')
    add_arrow(ax, 2.5, 2.4, 2.5, 1.8, color='#9C27B0')
    add_arrow(ax, 2.5, 1.8, 2.5, 3.0, curvature=0.3,
             label='need\nenergy/forces', color='#9C27B0', style='<-')

    # Dipole MACE flow
    add_arrow(ax, 7.5, 5.9, 7.5, 3.6, color='#00796B')
    add_arrow(ax, 7.5, 3.6, 6.2, 4.4, curvature=0.2, color='#00796B')
    add_arrow(ax, 7.5, 2.4, 7.5, 1.8, color='#00796B')
    add_arrow(ax, 7.5, 1.8, 7.5, 3.0, curvature=-0.3,
             label='need\ndipoles', color='#00796B', style='<-')

    # Reload arrows
    add_arrow(ax, 3.5, 4.5, 2.5, 3.6, label='reload', style='--', color='#F57C00')
    add_arrow(ax, 6.5, 4.5, 7.5, 3.6, label='reload', style='--', color='#F57C00')

    # Add explanation
    explanation = ("Key Mechanism:\n"
                  "• Two incompatible MACE versions coexist\n"
                  "• fake_module_from_real() copies module\n"
                  "• Inject into sys.modules dynamically\n"
                  "• importlib.reload() activates new version\n"
                  "• cleanup_mace_modules() prevents contamination")
    ax.text(5, 0.3, explanation, fontsize=8, ha='center',
           bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF9C4',
                    edgecolor='#F57F17', linewidth=1.5))

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white', dpi=dpi)
    plt.close()
    print(f"Module swapping diagram saved: {output_file}")


def main():
    """Generate all matplotlib-based flowcharts."""
    output_dir = "../assets/flowcharts"

    print("Generating matplotlib flowcharts...")
    print("=" * 60)

    print("Creating main workflow...")
    create_main_workflow(f"{output_dir}/workflow_matplotlib.png", dpi=300)

    print("Creating ZMQ architecture...")
    create_zmq_architecture(f"{output_dir}/zmq_architecture_matplotlib.png", dpi=300)

    print("Creating module swapping diagram...")
    create_module_swapping(f"{output_dir}/module_swapping_matplotlib.png", dpi=300)

    print("\n" + "=" * 60)
    print("Flowchart generation complete!")
    print(f"Check the '{output_dir}' directory for results.")


if __name__ == "__main__":
    main()
