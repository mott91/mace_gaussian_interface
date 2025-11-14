#!/usr/bin/env python3
"""
Test different terminal-style code snippet generation methods.

This script generates code snippets in various terminal themes:
- VS Code Dark
- Monokai
- Solarized Dark
- GitHub Dark
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from pygments import highlight
from pygments.lexers import PythonLexer, BashLexer
from pygments.formatters import ImageFormatter
from pygments.styles import get_style_by_name
import io
from PIL import Image
import numpy as np


# Color schemes for different terminal themes
THEMES = {
    "vscode_dark": {
        "background": "#1E1E1E",
        "text": "#D4D4D4",
        "comment": "#6A9955",
        "keyword": "#569CD6",
        "string": "#CE9178",
        "function": "#DCDCAA",
        "number": "#B5CEA8",
        "prompt": "#4EC9B0",
    },
    "monokai": {
        "background": "#272822",
        "text": "#F8F8F2",
        "comment": "#75715E",
        "keyword": "#F92672",
        "string": "#E6DB74",
        "function": "#A6E22E",
        "number": "#AE81FF",
        "prompt": "#66D9EF",
    },
    "solarized_dark": {
        "background": "#002B36",
        "text": "#839496",
        "comment": "#586E75",
        "keyword": "#268BD2",
        "string": "#2AA198",
        "function": "#B58900",
        "number": "#D33682",
        "prompt": "#859900",
    },
    "github_dark": {
        "background": "#0D1117",
        "text": "#C9D1D9",
        "comment": "#8B949E",
        "keyword": "#FF7B72",
        "string": "#A5D6FF",
        "function": "#D2A8FF",
        "number": "#79C0FF",
        "prompt": "#7EE787",
    },
}


def create_code_snippet_matplotlib(code, theme_name="vscode_dark", language="python",
                                   title="", show_prompt=True, dpi=300):
    """
    Create a terminal-style code snippet using matplotlib.

    Args:
        code: The code text to display
        theme_name: Name of the theme to use
        language: 'python' or 'bash'
        title: Optional title/comment at the top
        show_prompt: Whether to show a shell prompt
        dpi: DPI for the output image
    """
    theme = THEMES[theme_name]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6), dpi=dpi)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Background
    bg_rect = mpatches.Rectangle((0, 0), 1, 1,
                                  facecolor=theme["background"],
                                  edgecolor='none')
    ax.add_patch(bg_rect)

    # Add border for polish
    border = mpatches.Rectangle((0.01, 0.01), 0.98, 0.98,
                                facecolor='none',
                                edgecolor=theme["text"],
                                linewidth=1,
                                alpha=0.3)
    ax.add_patch(border)

    # Font properties
    try:
        # Try to use common monospace fonts
        font = FontProperties(family='monospace', size=11)
        # Test if we can use specific fonts
        for font_name in ['Consolas', 'Monaco', 'DejaVu Sans Mono', 'Courier New']:
            try:
                font = FontProperties(family=font_name, size=11)
                break
            except:
                continue
    except:
        font = FontProperties(family='monospace', size=11)

    # Starting y position
    y_pos = 0.92
    line_height = 0.04

    # Title/comment if provided
    if title:
        ax.text(0.04, y_pos, f"# {title}",
               fontproperties=font,
               color=theme["comment"],
               verticalalignment='top',
               fontsize=10)
        y_pos -= line_height

    # Process code lines
    lines = code.strip().split('\n')

    for line in lines:
        if y_pos < 0.08:  # Don't go below bottom margin
            break

        # Handle different line types
        if language == "bash" and show_prompt and not line.startswith('#'):
            # Add prompt for bash commands
            ax.text(0.04, y_pos, "$ ",
                   fontproperties=font,
                   color=theme["prompt"],
                   verticalalignment='top',
                   fontweight='bold',
                   fontsize=11)
            x_offset = 0.07
        else:
            x_offset = 0.04

        # Color the line based on content
        if line.strip().startswith('#'):
            color = theme["comment"]
        elif any(kw in line for kw in ['def ', 'class ', 'import ', 'from ']):
            color = theme["keyword"]
        elif '"' in line or "'" in line:
            color = theme["string"]
        else:
            color = theme["text"]

        ax.text(x_offset, y_pos, line,
               fontproperties=font,
               color=color,
               verticalalignment='top',
               fontsize=11)

        y_pos -= line_height

    # Add theme label in corner
    ax.text(0.96, 0.04, theme_name.replace('_', ' ').title(),
           fontproperties=font,
           color=theme["text"],
           alpha=0.3,
           fontsize=8,
           horizontalalignment='right')

    plt.tight_layout(pad=0)
    return fig


def create_code_snippet_pygments(code, style="monokai", language="python"):
    """
    Create a code snippet using Pygments directly.

    Args:
        code: The code text
        style: Pygments style name
        language: 'python' or 'bash'
    """
    lexer = PythonLexer() if language == "python" else BashLexer()

    # Generate image using Pygments
    formatter = ImageFormatter(
        style=style,
        font_name='DejaVu Sans Mono',
        font_size=14,
        line_numbers=False,
        line_pad=10
    )

    # Render to bytes
    img_bytes = highlight(code, lexer, formatter)

    # Load as PIL image
    img = Image.open(io.BytesIO(img_bytes))

    return img


# Example code snippets
PYTHON_EXAMPLE = """
# Import MACE calculator
from mace_calculators import load_dipole_mace_calculator

# Load custom dipole-enabled MACE model
calc = load_dipole_mace_calculator(
    model_paths="model_1.model",
    device="cuda",
    model_type="DipolePolarizabilityMACE"
)

# Calculate dipoles for geometry
dipoles = calc.calculate_dipoles(atoms)
"""

BASH_EXAMPLE = """
python cli.py run water.xyz --energy-calculators mace_mp --dipole-calculators espaloma
python cli.py list water
python run_analysis.py water
"""

ZMQ_EXAMPLE = """
# ZMQ Bridge in gm_helper.py
context = zmq.Context()
socket = context.socket(zmq.REQ)
socket.connect(f"ipc://{ipc_file}")

# Send geometry files to main process
message = f"{infile}|{outfile}"
socket.send_string(message)

# Wait for calculation to complete
response = socket.recv_string()
"""


def main():
    """Generate all test styles."""
    output_dir = "../assets/code_snippets"

    print("Generating terminal-style code snippets...")
    print("=" * 60)

    # Test 1: Python code in different themes (matplotlib)
    for theme_name in THEMES.keys():
        print(f"Creating {theme_name} theme (Python)...")
        fig = create_code_snippet_matplotlib(
            PYTHON_EXAMPLE,
            theme_name=theme_name,
            language="python",
            title="MACE Dipole Calculator Setup",
            show_prompt=False
        )
        filename = f"{output_dir}/{theme_name}_python_matplotlib.png"
        fig.savefig(filename, bbox_inches='tight', facecolor=THEMES[theme_name]["background"])
        plt.close(fig)
        print(f"  Saved: {filename}")

    # Test 2: Bash commands in different themes
    for theme_name in ["vscode_dark", "monokai"]:
        print(f"Creating {theme_name} theme (Bash)...")
        fig = create_code_snippet_matplotlib(
            BASH_EXAMPLE,
            theme_name=theme_name,
            language="bash",
            title="CLI Workflow Commands",
            show_prompt=True
        )
        filename = f"{output_dir}/{theme_name}_bash_matplotlib.png"
        fig.savefig(filename, bbox_inches='tight', facecolor=THEMES[theme_name]["background"])
        plt.close(fig)
        print(f"  Saved: {filename}")

    # Test 3: ZMQ bridge code
    theme_name = "vscode_dark"
    print(f"Creating ZMQ example ({theme_name})...")
    fig = create_code_snippet_matplotlib(
        ZMQ_EXAMPLE,
        theme_name=theme_name,
        language="python",
        title="ZMQ Inter-Process Communication",
        show_prompt=False
    )
    filename = f"{output_dir}/{theme_name}_zmq.png"
    fig.savefig(filename, bbox_inches='tight', facecolor=THEMES[theme_name]["background"])
    plt.close(fig)
    print(f"  Saved: {filename}")

    # Test 4: Pygments-based rendering
    print("\nGenerating Pygments-based snippets...")
    for style in ["monokai", "github-dark", "material"]:
        try:
            print(f"Creating {style} theme (Pygments)...")
            img = create_code_snippet_pygments(PYTHON_EXAMPLE, style=style)
            filename = f"{output_dir}/pygments_{style}.png"
            img.save(filename)
            print(f"  Saved: {filename}")
        except Exception as e:
            print(f"  Skipped {style}: {e}")

    print("\n" + "=" * 60)
    print(f"Code snippet generation complete!")
    print(f"Check the '{output_dir}' directory for results.")


if __name__ == "__main__":
    main()
