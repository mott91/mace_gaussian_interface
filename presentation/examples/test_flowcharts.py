#!/usr/bin/env python3
"""
Test different flowchart generation methods.

Methods tested:
1. Graphviz - Programmatic graph layout
2. Mermaid - Text-based diagrams (exported via CLI if available)
3. Matplotlib - Custom manual drawing
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import graphviz
from pathlib import Path


def create_flowchart_graphviz(output_file, format='png', dpi='300'):
    """
    Create the main MACE-Gaussian workflow using Graphviz.

    Args:
        output_file: Output filename (without extension)
        format: Output format (png, svg, pdf)
        dpi: Resolution for raster formats
    """
    dot = graphviz.Digraph(
        name='MACE_Gaussian_Pipeline',
        comment='MACE-Gaussian Interface Workflow',
        format=format,
        graph_attr={
            'rankdir': 'TB',
            'splines': 'ortho',
            'nodesep': '0.5',
            'ranksep': '0.8',
            'bgcolor': 'white',
            'dpi': dpi,
            'fontname': 'Helvetica',
            'fontsize': '12'
        },
        node_attr={
            'shape': 'box',
            'style': 'rounded,filled',
            'fillcolor': '#E3F2FD',
            'fontname': 'Helvetica',
            'fontsize': '11',
            'margin': '0.3,0.2'
        },
        edge_attr={
            'fontname': 'Helvetica',
            'fontsize': '10',
            'color': '#424242'
        }
    )

    # Add nodes
    dot.node('input', 'XYZ Input File', fillcolor='#C5E1A5')
    dot.node('opt', 'Geometry Optimization\n(MACE-MP/OMOL)', fillcolor='#FFE082')
    dot.node('freq', 'Frequency Calculation\n(Gaussian + ML Dipoles)', fillcolor='#FFE082')
    dot.node('zmq', 'ZMQ Bridge\n(Real-time Communication)', fillcolor='#CE93D8', shape='hexagon')
    dot.node('dipole', 'ML Dipole Calculator\n(Espaloma/MACE-ML)', fillcolor='#90CAF9')
    dot.node('anharm', 'Anharmonic Analysis\n(Gaussian)', fillcolor='#FFE082')
    dot.node('parse', 'Parse Results\n(JSON Export)', fillcolor='#FFAB91')
    dot.node('analyze', 'Statistical Analysis\n(Regression, KDE)', fillcolor='#FFAB91')
    dot.node('report', 'HTML Report\n(Plots + Tables)', fillcolor='#C5E1A5')

    # Add edges
    dot.edge('input', 'opt', label='atoms')
    dot.edge('opt', 'freq', label='optimized\ngeometry')
    dot.edge('freq', 'zmq', label='geometry\nrequest')
    dot.edge('zmq', 'dipole', label='IPC')
    dot.edge('dipole', 'zmq', label='dipoles +\nderivatives')
    dot.edge('zmq', 'freq', label='inject\ndata')
    dot.edge('freq', 'anharm', label='harmonic\nfrequencies')
    dot.edge('anharm', 'parse', label='log file')
    dot.edge('parse', 'analyze', label='frequencies +\nintensities')
    dot.edge('analyze', 'report', label='metrics +\nplots')

    # Render
    dot.render(output_file, cleanup=True)
    print(f"Graphviz flowchart saved: {output_file}.{format}")


def create_zmq_architecture_graphviz(output_file, format='png', dpi='300'):
    """
    Create detailed ZMQ architecture diagram.
    """
    dot = graphviz.Digraph(
        name='ZMQ_Architecture',
        format=format,
        graph_attr={
            'rankdir': 'LR',
            'splines': 'true',
            'bgcolor': 'white',
            'dpi': dpi,
            'fontname': 'Helvetica'
        },
        node_attr={
            'shape': 'box',
            'style': 'rounded,filled',
            'fontname': 'Helvetica',
            'fontsize': '11'
        },
        edge_attr={
            'fontname': 'Helvetica',
            'fontsize': '9'
        }
    )

    # Main process
    with dot.subgraph(name='cluster_main') as c:
        c.attr(label='gm_main.py (Python)', style='dashed', color='blue')
        c.node('server', 'ZMQ Server\n(REP socket)', fillcolor='#BBDEFB')
        c.node('calc', 'Dipole Calculator\n(MACE/Espaloma)', fillcolor='#C5E1A5')
        c.node('gaussian_start', 'Launch Gaussian\n(subprocess)', fillcolor='#FFE082')

    # Helper process
    with dot.subgraph(name='cluster_helper') as c:
        c.attr(label='gm_helper.py (Bridge)', style='dashed', color='green')
        c.node('client', 'ZMQ Client\n(REQ socket)', fillcolor='#C8E6C9')

    # Gaussian process
    with dot.subgraph(name='cluster_gaussian') as c:
        c.attr(label='Gaussian 16 (Fortran)', style='dashed', color='red')
        c.node('g16', 'Gaussian Engine', fillcolor='#FFCDD2')
        c.node('external', 'External Interface', fillcolor='#F8BBD0')

    # IPC file
    dot.node('ipc', '.ipc_file\n(Unix Socket)', shape='cylinder', fillcolor='#FFE0B2')

    # Connections
    dot.edge('gaussian_start', 'g16', label='spawn')
    dot.edge('g16', 'external', label='needs\ndipoles')
    dot.edge('external', 'client', label='call\nExternal=')
    dot.edge('client', 'ipc', label='connect')
    dot.edge('ipc', 'server', label='REQ-REP')
    dot.edge('server', 'calc', label='request')
    dot.edge('calc', 'server', label='dipoles')
    dot.edge('server', 'ipc', label='response')
    dot.edge('ipc', 'client', label='ready')
    dot.edge('client', 'external', label='return')
    dot.edge('external', 'g16', label='continue')

    # Render
    dot.render(output_file, cleanup=True)
    print(f"ZMQ architecture diagram saved: {output_file}.{format}")


def create_module_swapping_graphviz(output_file, format='png', dpi='300'):
    """
    Create module swapping mechanism diagram.
    """
    dot = graphviz.Digraph(
        name='Module_Swapping',
        format=format,
        graph_attr={
            'rankdir': 'TB',
            'bgcolor': 'white',
            'dpi': dpi,
            'fontname': 'Helvetica'
        },
        node_attr={
            'fontname': 'Helvetica',
            'fontsize': '11'
        }
    )

    # Packages
    dot.node('pkg1', 'mace_ML_pkg/\n(Standard MACE)', shape='folder', fillcolor='#E1BEE7', style='filled')
    dot.node('pkg2', 'mace_dipole_pkg/\n(Dipole MACE)', shape='folder', fillcolor='#B2DFDB', style='filled')

    # sys.modules
    dot.node('sysmod', 'sys.modules\n(Python Import Cache)', shape='box3d', fillcolor='#FFECB3', style='filled')

    # States
    dot.node('state1', 'State 1:\nGeometry Optimization', shape='ellipse', fillcolor='#FFCCBC', style='filled')
    dot.node('state2', 'State 2:\nDipole Calculation', shape='ellipse', fillcolor='#C5E1A5', style='filled')

    # Swap operations
    dot.node('swap1', 'load_standard_mace_calculator()', shape='hexagon', fillcolor='#B3E5FC', style='filled')
    dot.node('swap2', 'load_dipole_mace_calculator()', shape='hexagon', fillcolor='#B3E5FC', style='filled')

    # Connections
    dot.edge('state1', 'swap1', label='need\nenergy/forces')
    dot.edge('swap1', 'pkg1', label='import')
    dot.edge('pkg1', 'sysmod', label='inject into\nsys.modules')

    dot.edge('state2', 'swap2', label='need\ndipoles')
    dot.edge('swap2', 'pkg2', label='import')
    dot.edge('pkg2', 'sysmod', label='inject into\nsys.modules')

    dot.edge('sysmod', 'swap1', label='reload', style='dashed')
    dot.edge('sysmod', 'swap2', label='reload', style='dashed')

    # Render
    dot.render(output_file, cleanup=True)
    print(f"Module swapping diagram saved: {output_file}.{format}")


def create_flowchart_matplotlib(output_file, dpi=300):
    """
    Create a custom flowchart using matplotlib (manual drawing).
    This gives complete control over styling.
    """
    fig, ax = plt.subplots(figsize=(10, 12), dpi=dpi)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')

    # Helper function to create boxes
    def add_box(x, y, width, height, text, color='#E3F2FD', text_color='black'):
        box = FancyBboxPatch(
            (x - width/2, y - height/2), width, height,
            boxstyle="round,pad=0.1",
            edgecolor='#424242',
            facecolor=color,
            linewidth=2
        )
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center',
               fontsize=10, color=text_color, weight='bold')

    # Helper function to create arrows
    def add_arrow(x1, y1, x2, y2, label=''):
        arrow = FancyArrowPatch(
            (x1, y1), (x2, y2),
            arrowstyle='->,head_width=0.4,head_length=0.8',
            color='#424242',
            linewidth=2,
            connectionstyle="arc3,rad=0"
        )
        ax.add_patch(arrow)
        if label:
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(mid_x + 0.3, mid_y, label, fontsize=8,
                   style='italic', color='#666666')

    # Title
    ax.text(5, 13.5, 'MACE-Gaussian Interface Pipeline',
           ha='center', fontsize=14, weight='bold')

    # Build flowchart
    y = 12
    add_box(5, y, 3, 0.8, 'XYZ Input', '#C5E1A5')
    add_arrow(5, y - 0.4, 5, y - 1.1)

    y -= 1.5
    add_box(5, y, 4, 0.8, 'Geometry Optimization\n(MACE-MP/OMOL)', '#FFE082')
    add_arrow(5, y - 0.4, 5, y - 1.1)

    y -= 1.5
    add_box(5, y, 4, 0.8, 'Frequency Calculation\n(Gaussian)', '#FFE082')

    # Branch to ZMQ
    add_arrow(5, y - 0.4, 2, y - 1.9)
    add_arrow(8, y - 1.9, 5, y - 2.6)

    y -= 2.5
    add_box(2, y, 3.5, 0.8, 'ZMQ Bridge\n(IPC)', '#CE93D8')
    add_arrow(2, y - 0.4, 2, y - 1.1)

    y -= 1.5
    add_box(2, y, 3.5, 0.8, 'ML Dipole Calculator\n(Espaloma/MACE)', '#90CAF9')

    # Continue main flow
    y = 6.5
    add_box(5, y, 4, 0.8, 'Anharmonic Analysis\n(Gaussian)', '#FFE082')
    add_arrow(5, y - 0.4, 5, y - 1.1)

    y -= 1.5
    add_box(5, y, 3.5, 0.8, 'Parse Results\n(JSON)', '#FFAB91')
    add_arrow(5, y - 0.4, 5, y - 1.1)

    y -= 1.5
    add_box(5, y, 4, 0.8, 'Statistical Analysis\n(Regression, KDE)', '#FFAB91')
    add_arrow(5, y - 0.4, 5, y - 1.1)

    y -= 1.5
    add_box(5, y, 3.5, 0.8, 'HTML Report', '#C5E1A5')

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white', dpi=dpi)
    plt.close()
    print(f"Matplotlib flowchart saved: {output_file}")


def create_mermaid_diagram(output_file):
    """
    Create a Mermaid diagram (text format).
    Note: This creates the .mmd file. To render to PNG, you need mermaid-cli:
    npm install -g @mermaid-js/mermaid-cli
    mmdc -i diagram.mmd -o diagram.png
    """
    mermaid_code = """
graph TB
    Start[XYZ Input File] --> Opt[Geometry Optimization<br/>MACE-MP/OMOL]
    Opt --> Freq[Frequency Calculation<br/>Gaussian + ML Dipoles]
    Freq --> ZMQ{ZMQ Bridge<br/>IPC Communication}
    ZMQ --> Dipole[ML Dipole Calculator<br/>Espaloma/MACE-ML]
    Dipole --> ZMQ
    ZMQ --> Anharm[Anharmonic Analysis<br/>Gaussian]
    Anharm --> Parse[Parse Results<br/>JSON Export]
    Parse --> Analyze[Statistical Analysis<br/>Regression, KDE]
    Analyze --> Report[HTML Report<br/>Plots + Tables]

    style Start fill:#C5E1A5
    style Opt fill:#FFE082
    style Freq fill:#FFE082
    style ZMQ fill:#CE93D8
    style Dipole fill:#90CAF9
    style Anharm fill:#FFE082
    style Parse fill:#FFAB91
    style Analyze fill:#FFAB91
    style Report fill:#C5E1A5
    """

    with open(output_file, 'w') as f:
        f.write(mermaid_code.strip())

    print(f"Mermaid diagram saved: {output_file}")
    print("To render: mmdc -i {output_file} -o {output_file.replace('.mmd', '.png')}")


def main():
    """Generate all test flowcharts."""
    output_dir = Path("../assets/flowcharts")

    print("Generating flowcharts...")
    print("=" * 60)

    # Test 1: Main workflow (Graphviz)
    print("Creating main workflow (Graphviz)...")
    create_flowchart_graphviz(
        str(output_dir / "workflow_graphviz"),
        format='png',
        dpi='300'
    )

    # Test 2: ZMQ architecture (Graphviz)
    print("Creating ZMQ architecture (Graphviz)...")
    create_zmq_architecture_graphviz(
        str(output_dir / "zmq_architecture_graphviz"),
        format='png',
        dpi='300'
    )

    # Test 3: Module swapping (Graphviz)
    print("Creating module swapping diagram (Graphviz)...")
    create_module_swapping_graphviz(
        str(output_dir / "module_swapping_graphviz"),
        format='png',
        dpi='300'
    )

    # Test 4: Custom matplotlib
    print("Creating custom flowchart (Matplotlib)...")
    create_flowchart_matplotlib(
        str(output_dir / "workflow_matplotlib.png"),
        dpi=300
    )

    # Test 5: Mermaid (text output)
    print("Creating Mermaid diagram...")
    create_mermaid_diagram(str(output_dir / "workflow.mmd"))

    print("\n" + "=" * 60)
    print("Flowchart generation complete!")
    print(f"Check the '{output_dir}' directory for results.")


if __name__ == "__main__":
    main()
