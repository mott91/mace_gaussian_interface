#!/usr/bin/env python3
"""
Test detailed workflow flowcharts with inputs, outputs, and side branches.
Multiple styles to choose from.
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor


# Terminal Dark color scheme
BG = RGBColor(13, 17, 23)          # #0D1117 GitHub dark
TEXT = RGBColor(201, 209, 217)     # #C9D1D9
ACCENT = RGBColor(88, 166, 255)    # #58A6FF bright blue
DIM = RGBColor(139, 148, 158)      # #8B949E
GREEN = RGBColor(87, 171, 90)      # #57AB5A
ORANGE = RGBColor(219, 109, 40)    # #DB6D28
FONT = 'Consolas'


def create_detailed_vertical_workflow(prs):
    """Detailed vertical workflow with inputs/outputs and branches."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat workflow_detailed_v1.txt"
    p = tf.paragraphs[0]
    p.font.size = Pt(24)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Flowchart
    flow_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.0), Inches(9), Inches(6.2))
    tf = flow_box.text_frame

    lines = [
        ("INPUT: molecule.xyz (initial geometry)", GREEN),
        ("    │", DIM),
        ("    ▼", DIM),
        ("┌───────────────────────────────────────────────────────────────┐", ACCENT),
        ("│  Phase 1: Geometry Optimization                              │", ACCENT),
        ("│  Calculator: MACE-OMOL (default) or MACE-OFF/MACE-MP         │", TEXT),
        ("│  Output: optimized.xyz, results.json                         │", TEXT),
        ("└───────────────────────────────────────────────────────────────┘", ACCENT),
        ("    │", DIM),
        ("    ├──────────────────┐", DIM),
        ("    │                  │ (optional, default=True)", DIM),
        ("    ▼                  ▼", DIM),
        ("┌────────────────┐   ┌──────────────────────────────────────────┐", ACCENT),
        ("│  Phase 2:      │   │  Phase 3: ML Frequency Calculations      │", ACCENT),
        ("│  DFT Baseline  │   │  Energy: MACE-MP, MACE-OMOL              │", GREEN),
        ("│  B3LYP/        │   │  Dipole: Espaloma, MACE-ML               │", TEXT),
        ("│  6-31G(d,p)    │   │  ┌────────────────────────────────────┐  │", ORANGE),
        ("│  Output:       │   │  │ ZMQ Bridge                         │  │", ORANGE),
        ("│  reference     │   │  │ • Main: ZMQ server + ML dipoles    │  │", ORANGE),
        ("│  frequencies   │   │  │ • Helper: connects processes       │  │", ORANGE),
        ("└────────────────┘   │  │ • Gaussian: External calculation   │  │", ORANGE),
        ("    │                │  └────────────────────────────────────┘  │", ORANGE),
        ("    │                │  Output: {energy}_{dipole}/results.json  │", TEXT),
        ("    │                └──────────────────────────────────────────┘", ACCENT),
        ("    │                  │", DIM),
        ("    └──────────────────┘", DIM),
        ("    │", DIM),
        ("    ▼", DIM),
        ("┌───────────────────────────────────────────────────────────────┐", ACCENT),
        ("│  Phase 4: Analysis & Comparison                               │", ACCENT),
        ("│  • Compare ML vs DFT baseline                                 │", TEXT),
        ("│  • Metrics: MAE, RMSE, R², slope/intercept                    │", TEXT),
        ("│  • Generate regression plots                                  │", TEXT),
        ("│  • KDE broadening (FWHM=8 cm⁻¹)                               │", TEXT),
        ("│  Output: {molecule}_spectral_analysis.html                    │", TEXT),
        ("└───────────────────────────────────────────────────────────────┘", ACCENT),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(9)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(1)
        p.line_spacing = 1.0


def create_swimlane_workflow(prs):
    """Swim lane style showing different components."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat workflow_swimlane.txt"
    p = tf.paragraphs[0]
    p.font.size = Pt(24)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Flowchart
    flow_box = slide.shapes.add_textbox(Inches(0.4), Inches(1.0), Inches(9.2), Inches(6.2))
    tf = flow_box.text_frame

    lines = [
        ("═══════════════════════════════════════════════════════════════════", ACCENT),
        ("LANE 1: GEOMETRY    molecule.xyz → [MACE-OMOL] → optimized.xyz", ACCENT),
        ("═══════════════════════════════════════════════════════════════════", ACCENT),
        ("                                      │", DIM),
        ("═══════════════════════════════════════════════════════════════════", GREEN),
        ("LANE 2: DFT         optimized.xyz → [Gaussian B3LYP] → ref_freq", GREEN),
        ("═══════════════════════════════════════════════════════════════════", GREEN),
        ("                                      │", DIM),
        ("═══════════════════════════════════════════════════════════════════", ORANGE),
        ("LANE 3: ML FREQ     optimized.xyz ──┬→ [MACE-MP + Espaloma]", ORANGE),
        ("                                    ├→ [MACE-OMOL + Espaloma]", ORANGE),
        ("                                    ├→ [MACE-MP + MACE-ML]", ORANGE),
        ("                                    └→ [MACE-OMOL + MACE-ML]", ORANGE),
        ("", TEXT),
        ("                    ┌─────────────────────────────────────┐", ORANGE),
        ("                    │  ZMQ Bridge (Real-time)             │", ORANGE),
        ("                    │  Main ↔ Helper ↔ Gaussian           │", ORANGE),
        ("                    └─────────────────────────────────────┘", ORANGE),
        ("", TEXT),
        ("                    Each → {energy}_{dipole}/results.json", TEXT),
        ("═══════════════════════════════════════════════════════════════════", ORANGE),
        ("                                      │", DIM),
        ("═══════════════════════════════════════════════════════════════════", ACCENT),
        ("LANE 4: ANALYSIS    All results → [Compare] → HTML report", ACCENT),
        ("                    • Statistical metrics (MAE, RMSE, R²)", TEXT),
        ("                    • Regression plots", TEXT),
        ("                    • KDE spectra", TEXT),
        ("═══════════════════════════════════════════════════════════════════", ACCENT),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(9)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(1)
        p.line_spacing = 1.0


def create_tree_workflow(prs):
    """Tree-style workflow showing parallel branches."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat workflow_tree.txt"
    p = tf.paragraphs[0]
    p.font.size = Pt(24)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Flowchart
    flow_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.0), Inches(9), Inches(6.2))
    tf = flow_box.text_frame

    lines = [
        ("molecule.xyz", GREEN),
        ("    │", DIM),
        ("    ├─▶ [Phase 1: Optimize with MACE-OMOL]", ACCENT),
        ("    │       → optimized.xyz", TEXT),
        ("    │", DIM),
        ("    ├─────┬─▶ [Phase 2: DFT Baseline]", GREEN),
        ("    │     │       B3LYP/6-31G(d,p)", TEXT),
        ("    │     │       → reference_frequencies.json", TEXT),
        ("    │     │", DIM),
        ("    │     └─▶ [Phase 3: ML Frequencies] ──┬─▶ MACE-MP + Espaloma", ORANGE),
        ("    │             Energy × Dipole         ├─▶ MACE-MP + MACE-ML", ORANGE),
        ("    │             combinations            ├─▶ MACE-OMOL + Espaloma", ORANGE),
        ("    │                                     └─▶ MACE-OMOL + MACE-ML", ORANGE),
        ("    │             │", DIM),
        ("    │             ├─▶ [ZMQ Bridge]", ORANGE),
        ("    │             │       Main Process (Python)", TEXT),
        ("    │             │       ├─▶ ZMQ Server", TEXT),
        ("    │             │       └─▶ ML Dipole Calc", TEXT),
        ("    │             │       │", DIM),
        ("    │             │       Helper Script", TEXT),
        ("    │             │       ├─▶ Connect Main ↔ Gaussian", TEXT),
        ("    │             │       └─▶ File path exchange", TEXT),
        ("    │             │       │", DIM),
        ("    │             │       Gaussian 16", TEXT),
        ("    │             │       └─▶ Frequency calculation", TEXT),
        ("    │             │", DIM),
        ("    │             └─▶ → {energy}_{dipole}/results.json", TEXT),
        ("    │", DIM),
        ("    └─────────────────▶ [Phase 4: Analysis]", ACCENT),
        ("                            Compare all methods", TEXT),
        ("                            ├─▶ Statistics (MAE, RMSE, R²)", TEXT),
        ("                            ├─▶ Regression plots", TEXT),
        ("                            ├─▶ KDE spectra (FWHM=8)", TEXT),
        ("                            └─▶ HTML report", TEXT),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(10)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(1)
        p.line_spacing = 1.0


def create_compact_detailed_workflow(prs):
    """Compact but detailed workflow with annotations."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat workflow_compact.txt"
    p = tf.paragraphs[0]
    p.font.size = Pt(24)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Flowchart
    flow_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.0), Inches(8.4), Inches(6.2))
    tf = flow_box.text_frame

    lines = [
        ("┌─────────────────────────────────────────────────────┐", GREEN),
        ("│  INPUT: molecule.xyz                                │", GREEN),
        ("└─────────────────────────────────────────────────────┘", GREEN),
        ("    ↓", DIM),
        ("┌─────────────────────────────────────────────────────┐", ACCENT),
        ("│  [1] Optimize → MACE-OMOL → optimized.xyz          │", ACCENT),
        ("└─────────────────────────────────────────────────────┘", ACCENT),
        ("    ↓", DIM),
        ("    ├──────────────────────┬─────────────────────────┐", DIM),
        ("    ▼                      ▼                         ▼", DIM),
        ("┌─────────────┐  ┌──────────────────────────────────────┐", GREEN),
        ("│  [2] DFT    │  │  [3] ML Frequencies (4 combos)       │", GREEN),
        ("│  Baseline   │  │  ┌────────────────────────────────┐  │", ORANGE),
        ("│  B3LYP      │  │  │  Energy: MACE-MP, MACE-OMOL    │  │", ORANGE),
        ("│  6-31G**    │  │  │  Dipole: Espaloma, MACE-ML     │  │", ORANGE),
        ("│             │  │  │  ────────────────────────────  │  │", TEXT),
        ("│  Output:    │  │  │  ZMQ: Main↔Helper↔Gaussian    │  │", ORANGE),
        ("│  reference  │  │  └────────────────────────────────┘  │", ORANGE),
        ("└─────────────┘  │  Output: 4 × results.json            │", TEXT),
        ("    │            └──────────────────────────────────────┘", GREEN),
        ("    │                      │", DIM),
        ("    └──────────────────────┘", DIM),
        ("    ↓", DIM),
        ("┌─────────────────────────────────────────────────────┐", ACCENT),
        ("│  [4] Compare & Analyze                              │", ACCENT),
        ("│  • Metrics: MAE, RMSE, R²                           │", TEXT),
        ("│  • Plots: regression, KDE spectra (FWHM=8)          │", TEXT),
        ("│  • Output: {molecule}_spectral_analysis.html        │", TEXT),
        ("└─────────────────────────────────────────────────────┘", ACCENT),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(10)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(1)
        p.line_spacing = 1.0


def main():
    """Generate detailed workflow test slides."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    print("Generating detailed workflow diagrams...")

    create_detailed_vertical_workflow(prs)
    print("  ✓ Detailed vertical with branches")

    create_swimlane_workflow(prs)
    print("  ✓ Swim lane style")

    create_tree_workflow(prs)
    print("  ✓ Tree-style with parallel branches")

    create_compact_detailed_workflow(prs)
    print("  ✓ Compact detailed workflow")

    # Save presentation
    output_path = "../test_detailed_workflow.pptx"
    prs.save(output_path)
    print(f"\n✓ Test presentation created: {output_path}")
    print(f"  → 4 slides with detailed workflow variations")
    print(f"  → Shows inputs, outputs, branches, and ZMQ details")
    print(f"  → Terminal Dark styling")


if __name__ == "__main__":
    main()
