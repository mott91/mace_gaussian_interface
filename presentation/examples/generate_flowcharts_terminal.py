#!/usr/bin/env python3
"""
Generate clean, terminal-style flowcharts with minimalist aesthetic.
Uses monospace fonts and terminal color scheme.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.font_manager import FontProperties
import numpy as np


# Terminal color schemes
TERMINAL_DARK = {
    'bg': '#0D1117',           # GitHub dark background
    'fg': '#C9D1D9',           # Main text
    'accent1': '#58A6FF',      # Bright blue
    'accent2': '#3FB950',      # Bright green
    'accent3': '#F778BA',      # Bright pink
    'accent4': '#D2A8FF',      # Bright purple
    'dim': '#8B949E',          # Dim text
    'border': '#30363D',       # Subtle borders
}

TERMINAL_MINIMAL = {
    'bg': '#FAFAFA',           # Off-white
    'fg': '#2E3440',           # Dark gray text
    'accent1': '#5E81AC',      # Muted blue
    'accent2': '#88C0D0',      # Light blue
    'accent3': '#81A1C1',      # Soft blue
    'accent4': '#B48EAD',      # Muted purple
    'dim': '#D8DEE9',          # Light gray
    'border': '#4C566A',       # Medium gray
}


def get_mono_font(size=10):
    """Get monospace font."""
    for font_name in ['DejaVu Sans Mono', 'Courier New', 'Consolas', 'Monaco']:
        try:
            return FontProperties(family=font_name, size=size)
        except:
            continue
    return FontProperties(family='monospace', size=size)


def add_terminal_box(ax, x, y, width, height, text, theme, style='normal'):
    """Add a minimal terminal-style box."""
    colors = TERMINAL_DARK if theme == 'dark' else TERMINAL_MINIMAL

    # Choose color based on style
    if style == 'input':
        fill_color = colors['accent2']
        text_color = colors['bg'] if theme == 'dark' else '#FFFFFF'
    elif style == 'process':
        fill_color = colors['accent1']
        text_color = colors['bg'] if theme == 'dark' else '#FFFFFF'
    elif style == 'ml':
        fill_color = colors['accent3']
        text_color = colors['bg'] if theme == 'dark' else '#FFFFFF'
    elif style == 'special':
        fill_color = colors['accent4']
        text_color = colors['bg'] if theme == 'dark' else '#FFFFFF'
    else:
        fill_color = colors['border']
        text_color = colors['fg']

    # Draw box with subtle corners
    box = Rectangle(
        (x - width/2, y - height/2), width, height,
        facecolor=fill_color,
        edgecolor=colors['border'],
        linewidth=1.5,
        alpha=0.9
    )
    ax.add_patch(box)

    # Add text
    font = get_mono_font(9)
    ax.text(x, y, text, ha='center', va='center',
           fontproperties=font, color=text_color, weight='normal')


def add_simple_arrow(ax, x1, y1, x2, y2, theme, label=''):
    """Add a clean, simple arrow."""
    colors = TERMINAL_DARK if theme == 'dark' else TERMINAL_MINIMAL

    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle='->',
        color=colors['dim'],
        linewidth=1.5,
        mutation_scale=15,
        alpha=0.7
    )
    ax.add_patch(arrow)

    if label:
        font = get_mono_font(7)
        mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mid_x + 0.3, mid_y, label, fontproperties=font,
               color=colors['dim'], fontsize=7, style='italic')


def create_workflow_terminal(output_file, theme='dark', dpi=300):
    """Create minimal terminal-style workflow."""
    colors = TERMINAL_DARK if theme == 'dark' else TERMINAL_MINIMAL

    fig, ax = plt.subplots(figsize=(8, 10), dpi=dpi, facecolor=colors['bg'])
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 12)
    ax.axis('off')
    ax.set_facecolor(colors['bg'])

    # Title
    font_title = get_mono_font(12)
    ax.text(4, 11.5, '// MACE-Gaussian Pipeline', ha='center',
           fontproperties=font_title, color=colors['fg'], weight='bold')

    # Fixed grid layout (no overlaps)
    y = 10.5
    dy = 1.3

    # 1. Input
    add_terminal_box(ax, 4, y, 2.5, 0.6, 'XYZ Input', theme, 'input')
    add_simple_arrow(ax, 4, y-0.3, 4, y-dy+0.3, theme)

    # 2. Optimization
    y -= dy
    add_terminal_box(ax, 4, y, 3.0, 0.6, 'Geometry Opt [MACE]', theme, 'ml')
    add_simple_arrow(ax, 4, y-0.3, 4, y-dy+0.3, theme)

    # 3. Frequency calc
    y -= dy
    add_terminal_box(ax, 4, y, 3.2, 0.6, 'Frequency [Gaussian]', theme, 'process')

    # Branch to ZMQ
    add_simple_arrow(ax, 2.8, y-0.15, 1.5, y-0.8, theme)

    # 4. ZMQ system (side branch)
    y_zmq = y - 1.1
    add_terminal_box(ax, 1.5, y_zmq, 2.0, 0.5, 'ZMQ Bridge', theme, 'special')
    add_simple_arrow(ax, 1.5, y_zmq-0.25, 1.5, y_zmq-0.7, theme)

    y_zmq -= 0.9
    add_terminal_box(ax, 1.5, y_zmq, 2.2, 0.5, 'ML Dipoles', theme, 'ml')

    # Return arrow
    add_simple_arrow(ax, 2.6, y_zmq+0.25, 4, y-0.5, theme, label='inject')

    # Continue main flow
    y -= 2.6
    add_simple_arrow(ax, 4, y+1.3, 4, y+0.3, theme)

    # 5. Anharmonic
    add_terminal_box(ax, 4, y, 3.2, 0.6, 'Anharmonic [VPT2]', theme, 'process')
    add_simple_arrow(ax, 4, y-0.3, 4, y-dy+0.3, theme)

    # 6. Parse
    y -= dy
    add_terminal_box(ax, 4, y, 2.5, 0.6, 'Parse → JSON', theme, 'normal')
    add_simple_arrow(ax, 4, y-0.3, 4, y-dy+0.3, theme)

    # 7. Analysis
    y -= dy
    add_terminal_box(ax, 4, y, 3.0, 0.6, 'Analysis [Stats]', theme, 'normal')
    add_simple_arrow(ax, 4, y-0.3, 4, y-dy+0.3, theme)

    # 8. Output
    y -= dy
    add_terminal_box(ax, 4, y, 2.5, 0.6, 'HTML Report', theme, 'input')

    # Add minimal footer
    font_small = get_mono_font(7)
    ax.text(0.3, 0.3, f'// theme: {theme}', fontproperties=font_small,
           color=colors['dim'], alpha=0.5)

    plt.tight_layout(pad=0.5)
    plt.savefig(output_file, bbox_inches='tight', facecolor=colors['bg'], dpi=dpi)
    plt.close()
    print(f"Terminal workflow ({theme}) saved: {output_file}")


def create_zmq_terminal(output_file, theme='dark', dpi=300):
    """Create minimal ZMQ architecture diagram."""
    colors = TERMINAL_DARK if theme == 'dark' else TERMINAL_MINIMAL

    fig, ax = plt.subplots(figsize=(10, 6), dpi=dpi, facecolor=colors['bg'])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')
    ax.set_facecolor(colors['bg'])

    # Title
    font_title = get_mono_font(12)
    ax.text(5, 5.5, '// ZMQ Architecture', ha='center',
           fontproperties=font_title, color=colors['fg'], weight='bold')

    # Three columns: Main | Helper | Gaussian
    x1, x2, x3 = 1.5, 5.0, 8.5

    # Main process
    font = get_mono_font(8)
    ax.text(x1, 4.8, 'gm_main.py', ha='center', fontproperties=font,
           color=colors['accent1'], weight='bold')
    add_terminal_box(ax, x1, 4.2, 2.0, 0.4, 'ZMQ Server', theme, 'process')
    add_terminal_box(ax, x1, 3.6, 2.0, 0.4, 'Calculate', theme, 'ml')
    add_terminal_box(ax, x1, 3.0, 2.0, 0.4, 'Spawn G16', theme, 'normal')

    # Helper process
    ax.text(x2, 4.8, 'gm_helper.py', ha='center', fontproperties=font,
           color=colors['accent2'], weight='bold')
    add_terminal_box(ax, x2, 4.2, 2.0, 0.4, 'ZMQ Client', theme, 'process')
    add_terminal_box(ax, x2, 3.6, 2.0, 0.4, 'Bridge', theme, 'special')

    # Gaussian process
    ax.text(x3, 4.8, 'Gaussian 16', ha='center', fontproperties=font,
           color=colors['accent4'], weight='bold')
    add_terminal_box(ax, x3, 4.2, 2.0, 0.4, 'Engine', theme, 'process')
    add_terminal_box(ax, x3, 3.6, 2.0, 0.4, 'External=', theme, 'normal')

    # IPC socket (center bottom)
    add_terminal_box(ax, 5, 1.8, 2.5, 0.5, '.ipc_file', theme, 'special')

    # Arrows showing flow
    # 1. Spawn
    add_simple_arrow(ax, 2.5, 3.0, 7.5, 3.6, theme, label='1.spawn')

    # 2. Call helper
    add_simple_arrow(ax, 7.5, 3.6, 6.0, 3.6, theme, label='2.call')

    # 3. Connect to socket
    add_simple_arrow(ax, x2, 3.4, x2, 2.3, theme, label='3.connect')
    add_simple_arrow(ax, x1, 3.4, x1, 2.3, theme)

    # 4. Request/Response
    add_simple_arrow(ax, 4.0, 1.8, 2.5, 1.8, theme, label='4.REQ/REP')

    # 5. Return
    add_simple_arrow(ax, 6.0, 4.0, 7.5, 4.0, theme, label='5.return')

    # Add sequence
    seq_font = get_mono_font(7)
    seq_text = "1→spawn G16\n2→call helper\n3→connect IPC\n4→request/respond\n5→return to G16"
    ax.text(0.3, 1.5, seq_text, fontproperties=seq_font,
           color=colors['dim'], verticalalignment='top', linespacing=1.5)

    plt.tight_layout(pad=0.5)
    plt.savefig(output_file, bbox_inches='tight', facecolor=colors['bg'], dpi=dpi)
    plt.close()
    print(f"ZMQ architecture ({theme}) saved: {output_file}")


def create_module_swap_terminal(output_file, theme='dark', dpi=300):
    """Create minimal module swapping diagram."""
    colors = TERMINAL_DARK if theme == 'dark' else TERMINAL_MINIMAL

    fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi, facecolor=colors['bg'])
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 6)
    ax.axis('off')
    ax.set_facecolor(colors['bg'])

    # Title
    font_title = get_mono_font(12)
    ax.text(4, 5.5, '// Module Swapping', ha='center',
           fontproperties=font_title, color=colors['fg'], weight='bold')

    # Two packages at top
    add_terminal_box(ax, 2, 4.5, 2.5, 0.5, 'mace_ML_pkg', theme, 'ml')
    add_terminal_box(ax, 6, 4.5, 2.5, 0.5, 'mace_dipole_pkg', theme, 'ml')

    # sys.modules in middle
    add_terminal_box(ax, 4, 3.2, 3.0, 0.6, 'sys.modules', theme, 'special')

    # Swapper functions
    add_terminal_box(ax, 2, 2.0, 2.2, 0.4, 'load_standard()', theme, 'process')
    add_terminal_box(ax, 6, 2.0, 2.2, 0.4, 'load_dipole()', theme, 'process')

    # Use cases
    add_terminal_box(ax, 2, 0.8, 2.0, 0.4, 'Optimize', theme, 'input')
    add_terminal_box(ax, 6, 0.8, 2.0, 0.4, 'Dipoles', theme, 'input')

    # Arrows
    # Standard path
    add_simple_arrow(ax, 2, 4.2, 2, 3.8, theme)
    add_simple_arrow(ax, 2.7, 3.2, 2, 2.4, theme)
    add_simple_arrow(ax, 2, 1.6, 2, 1.2, theme)

    # Dipole path
    add_simple_arrow(ax, 6, 4.2, 6, 3.8, theme)
    add_simple_arrow(ax, 5.3, 3.2, 6, 2.4, theme)
    add_simple_arrow(ax, 6, 1.6, 6, 1.2, theme)

    # Explanation
    font_small = get_mono_font(7)
    explanation = "// Swap MACE implementations\n// via sys.modules injection\n// cleanup_mace_modules() after use"
    ax.text(4, 0.2, explanation, ha='center', fontproperties=font_small,
           color=colors['dim'], linespacing=1.4)

    plt.tight_layout(pad=0.5)
    plt.savefig(output_file, bbox_inches='tight', facecolor=colors['bg'], dpi=dpi)
    plt.close()
    print(f"Module swapping ({theme}) saved: {output_file}")


def main():
    """Generate all terminal-style flowcharts."""
    output_dir = "../assets/flowcharts"

    print("Generating terminal-style flowcharts...")
    print("=" * 60)

    # Dark theme versions
    print("\nDark theme:")
    create_workflow_terminal(f"{output_dir}/workflow_terminal_dark.png", 'dark', dpi=300)
    create_zmq_terminal(f"{output_dir}/zmq_terminal_dark.png", 'dark', dpi=300)
    create_module_swap_terminal(f"{output_dir}/module_swap_terminal_dark.png", 'dark', dpi=300)

    # Light/minimal theme versions
    print("\nMinimal theme:")
    create_workflow_terminal(f"{output_dir}/workflow_terminal_minimal.png", 'minimal', dpi=300)
    create_zmq_terminal(f"{output_dir}/zmq_terminal_minimal.png", 'minimal', dpi=300)
    create_module_swap_terminal(f"{output_dir}/module_swap_terminal_minimal.png", 'minimal', dpi=300)

    print("\n" + "=" * 60)
    print("Terminal-style flowcharts complete!")


if __name__ == "__main__":
    main()
