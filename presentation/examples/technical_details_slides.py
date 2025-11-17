#!/usr/bin/env python3
"""
Generate technical detail slides for MACE-Gaussian Interface presentation.

Covers:
1. Dipole calculation methods (MACE-ML and Espaloma)
2. Eigenvector mode matching algorithm
3. Output file parsing requirements
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor

# Terminal Dark Color Scheme (GitHub Dark)
BG = RGBColor(13, 17, 23)        # #0D1117 - background
TEXT = RGBColor(201, 209, 217)   # #C9D1D9 - main text
ACCENT = RGBColor(88, 166, 255)  # #58A6FF - blue accent
DIM = RGBColor(139, 148, 158)    # #8B949E - dimmed text
GREEN = RGBColor(87, 171, 90)    # #57AB5A - success/input
ORANGE = RGBColor(219, 109, 40)  # #DB6D28 - warning/highlight
PURPLE = RGBColor(188, 140, 255) # #BC8CFF - special

FONT = "Consolas"


def create_dipole_calculation_slide(prs):
    """Slide explaining dipole calculation methods."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.4), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat dipole_methods.md"
    p = tf.paragraphs[0]
    p.font.size = Pt(28)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Content
    content_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(4.8))
    tf = content_box.text_frame

    lines = [
        ("# MACE-ML (Direct Prediction)", ACCENT),
        ("", TEXT),
        ("  Input → E(3)-equivariant NN → μ = (μx, μy, μz)", DIM),
        ("", TEXT),
        ("  → Direct dipole prediction", TEXT),
        ("  → Fast (~0.1s per geometry)", TEXT),
        ("  → Requires custom model", TEXT),
        ("", TEXT),
        ("", TEXT),
        ("# Espaloma (Charge-Based)", ACCENT),
        ("", TEXT),
        ("  Predict charges → μ = Σ qᵢ·rᵢ", DIM),
        ("", TEXT),
        ("  → ML partial charges (q₁, q₂, ..., qₙ)", TEXT),
        ("  → Sum over atomic positions", TEXT),
        ("  → Widely available, interpretable", TEXT),
        ("", TEXT),
        ("", TEXT),
        ("# Fallback Chain", ACCENT),
        ("  MACE-ML → Espaloma → xTB → Geometry", DIM),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(6)
        p.line_spacing = 1.0


def create_eigenvector_matching_slide(prs):
    """Slide explaining eigenvector mode matching."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.4), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ cat mode_matching.md"
    p = tf.paragraphs[0]
    p.font.size = Pt(28)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Content
    content_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(4.8))
    tf = content_box.text_frame

    lines = [
        ("# Problem: Frequency Order Mismatch", ACCENT),
        ("", TEXT),
        ("  DFT:  [500, 1200, 1500, 2900, 3100] cm⁻¹", DIM),
        ("  ML:   [498, 3095, 1205, 2895, 1498] cm⁻¹", DIM),
        ("", TEXT),
        ("  → Cannot compare by index", TEXT),
        ("", TEXT),
        ("", TEXT),
        ("# Solution: Eigenvector Dot Product", ACCENT),
        ("", TEXT),
        ("  similarity = |v_ML · v_DFT|", DIM),
        ("", TEXT),
        ("  → Compare atomic displacement patterns", TEXT),
        ("  → Find best match (threshold > 0.8)", TEXT),
        ("  → Robust to frequency ordering", TEXT),
        ("", TEXT),
        ("", TEXT),
        ("# Example", ACCENT),
        ("  ML 498 cm⁻¹ → DFT 500 cm⁻¹  [0.996]", DIM),
        ("  ML 3095 cm⁻¹ → DFT 3100 cm⁻¹ [0.989]", DIM),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(6)
        p.line_spacing = 1.0


def create_parsing_slide(prs):
    """Slide explaining output file parsing requirements."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.4), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = "$ ls -la parsing/"
    p = tf.paragraphs[0]
    p.font.size = Pt(28)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Content
    content_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(4.8))
    tf = content_box.text_frame

    lines = [
        ("# Gaussian .log Files (2-10 MB)", ACCENT),
        ("", TEXT),
        ("  → Harmonic frequencies", TEXT),
        ("  → Anharmonic frequencies", TEXT),
        ("  → Overtones & combination bands", TEXT),
        ("  → IR intensities", TEXT),
        ("  → Eigenvectors (3N-6 modes)", TEXT),
        ("", TEXT),
        ("  Challenges:", DIM),
        ("  • Multi-column format (3 modes/block)", TEXT),
        ("  • 10,000-50,000 lines to search", TEXT),
        ("  • Regex pattern matching required", TEXT),
        ("", TEXT),
        ("", TEXT),
        ("# JSON Results (ML)", ACCENT),
        ("", TEXT),
        ("  {", DIM),
        ("    \"frequencies\": [...],", DIM),
        ("    \"ir_intensities\": [...],", DIM),
        ("    \"eigenvectors\": [[...]]", DIM),
        ("  }", DIM),
        ("", TEXT),
        ("  → Structured, easy to parse", TEXT),
    ]

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(6)
        p.line_spacing = 1.0


def main():
    """Generate technical details presentation."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    print("Generating technical detail slides...")

    create_dipole_calculation_slide(prs)
    print("  ✓ Dipole calculation methods (MACE-ML & Espaloma)")

    create_eigenvector_matching_slide(prs)
    print("  ✓ Eigenvector mode matching algorithm")

    create_parsing_slide(prs)
    print("  ✓ Output file parsing requirements")

    # Save presentation
    output_path = "../technical_details.pptx"
    prs.save(output_path)
    print(f"\n✓ Technical details presentation created: {output_path}")
    print(f"  → 3 slides covering implementation details")
    print(f"  → Terminal Dark styling")


if __name__ == "__main__":
    main()
