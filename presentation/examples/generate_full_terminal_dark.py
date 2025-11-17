#!/usr/bin/env python3
"""
Generate complete 10-slide presentation using Terminal Dark template styling.
Based on the original template_terminal_dark.pptx design.
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
FONT = 'Consolas'


def add_footer(slide):
    """Add footer to slide (consistent across all slides)."""
    info_box = slide.shapes.add_textbox(Inches(1), Inches(6.5), Inches(8), Inches(0.6))
    tf = info_box.text_frame
    tf.text = "Your Name · Research Group · 2025-11-15"
    p = tf.paragraphs[0]
    p.font.size = Pt(12)
    p.font.name = FONT
    p.font.color.rgb = DIM
    p.alignment = PP_ALIGN.CENTER


def create_title_slide(prs):
    """Slide 0: Title slide."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Prompt-style title
    title_box = slide.shapes.add_textbox(Inches(1), Inches(2.5), Inches(8), Inches(1.5))
    tf = title_box.text_frame
    tf.text = "$ ./presentation.sh\n> MACE-Gaussian Interface"
    p = tf.paragraphs[0]
    p.font.size = Pt(32)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True
    p.alignment = PP_ALIGN.LEFT

    # Subtitle
    subtitle_box = slide.shapes.add_textbox(Inches(1), Inches(4.2), Inches(8), Inches(0.6))
    tf = subtitle_box.text_frame
    tf.text = "// ML-accelerated anharmonic spectroscopy"
    p = tf.paragraphs[0]
    p.font.size = Pt(16)
    p.font.name = FONT
    p.font.color.rgb = DIM

    # Footer
    add_footer(slide)


def create_content_slide(prs, command, lines):
    """Create a content slide with command header and bullet points."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Title with prompt (28pt as requested)
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.4), Inches(9), Inches(0.5))
    tf = title_box.text_frame
    tf.text = command
    p = tf.paragraphs[0]
    p.font.size = Pt(28)  # Increased from 20pt to 28pt
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    # Content box
    content_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(4.8))
    tf = content_box.text_frame

    for i, (line, color) in enumerate(lines):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(6)

    # Footer
    add_footer(slide)


def main():
    """Generate complete presentation."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    # Slide 0: Title
    create_title_slide(prs)

    # Slide 1: Motivation
    create_content_slide(prs, "$ cat motivation.md", [
        ("# Problem", ACCENT),
        ("  → DFT anharmonic: computationally expensive", TEXT),
        ("  → Hours to days per molecule", TEXT),
        ("  → Limited throughput for screening", TEXT),
        ("", TEXT),
        ("# Solution", ACCENT),
        ("  → ML potentials: 10-100× faster", TEXT),
        ("  → Maintain accuracy via hybrid approach", TEXT),
        ("  → Enable high-throughput predictions", TEXT),
        ("", TEXT),
        ("# Impact", ACCENT),
        ("  → Rapid spectral database generation", TEXT),
        ("  → Practical molecular screening", TEXT),
    ])

    # Slide 2: Architecture Overview
    create_content_slide(prs, "$ cat architecture.md", [
        ("# System Architecture", ACCENT),
        ("", TEXT),
        ("┌─────────────┐      ┌──────────────┐", DIM),
        ("│ MACE        │──────│ Gaussian 16  │", DIM),
        ("│ (Python)    │ ZMQ  │ (Fortran)    │", DIM),
        ("└─────────────┘      └──────────────┘", DIM),
        ("", TEXT),
        ("# Components", ACCENT),
        ("  → Geometry optimization: MACE-OMOL", TEXT),
        ("  → Energy calculation: Gaussian/MACE-MP", TEXT),
        ("  → Dipole derivatives: ML calculators", TEXT),
        ("  → Anharmonic analysis: Gaussian VPT2", TEXT),
    ])

    # Slide 3: ZMQ Bridge
    create_content_slide(prs, "$ cat zmq_bridge.md", [
        ("# Inter-Process Communication", ACCENT),
        ("", TEXT),
        ("  1. Gaussian calls external program", TEXT),
        ("  2. Helper script connects via ZMQ socket", TEXT),
        ("  3. Main process receives geometry", TEXT),
        ("  4. ML calculator computes dipoles", TEXT),
        ("  5. Results written to Gaussian format", TEXT),
        ("  6. Gaussian continues calculation", TEXT),
        ("", TEXT),
        ("# Key Innovation", ACCENT),
        ("  → Real-time injection of ML predictions", TEXT),
        ("  → No modification to Gaussian source", TEXT),
        ("  → File-based IPC for reliability", TEXT),
    ])

    # Slide 4: Dipole Calculators
    create_content_slide(prs, "$ ls -la calculators/", [
        ("# Modular Dipole System", ACCENT),
        ("", TEXT),
        ("drwxr-xr-x  mace_ml/          # Primary", DIM),
        ("drwxr-xr-x  espaloma/         # Fallback", DIM),
        ("drwxr-xr-x  xtb/              # Fast", DIM),
        ("drwxr-xr-x  geometry_based/   # Last resort", DIM),
        ("", TEXT),
        ("# Features", ACCENT),
        ("  → Automatic fallback chain", TEXT),
        ("  → Custom MACE dipole model", TEXT),
        ("  → Charge-based dipole (Espaloma)", TEXT),
        ("  → Semi-empirical xTB backup", TEXT),
    ])

    # Slide 5: Module Swapping
    create_content_slide(prs, "$ cat module_swap.py", [
        ("# Dual MACE Architecture", ACCENT),
        ("", TEXT),
        ("mace_ML_pkg/          # Standard MACE", DIM),
        ("mace_dipole_pkg/      # Custom fork", DIM),
        ("", TEXT),
        ("# Implementation", ACCENT),
        ("  → Dynamic module replacement", TEXT),
        ("  → Cache invalidation required", TEXT),
        ("  → Cleanup after dipole calculations", TEXT),
        ("", TEXT),
        ("# Challenge", ACCENT),
        ("  → Prevent module contamination", TEXT),
        ("  → Thread-safe switching mechanism", TEXT),
    ])

    # Slide 6: Validation Method
    create_content_slide(prs, "$ cat validation.md", [
        ("# Eigenvector-Based Matching", ACCENT),
        ("", TEXT),
        ("  → Compare DFT vs ML normal modes", TEXT),
        ("  → Dot product similarity threshold", TEXT),
        ("  → Handle mode ordering differences", TEXT),
        ("", TEXT),
        ("# Metrics", ACCENT),
        ("  → MAE (Mean Absolute Error)", TEXT),
        ("  → RMSE (Root Mean Square Error)", TEXT),
        ("  → R² (Coefficient of Determination)", TEXT),
        ("  → Slope/intercept analysis", TEXT),
        ("", TEXT),
        ("# Spectral Analysis", ACCENT),
        ("  → Gaussian KDE broadening (8 cm⁻¹)", TEXT),
    ])

    # Slide 7: Results - Fundamentals
    create_content_slide(prs, "$ python analyze_spectra.py", [
        ("# Frequency Predictions", ACCENT),
        ("", TEXT),
        ("Test molecules:", TEXT),
        ("  → Water (H₂O)", TEXT),
        ("  → Acetic acid (CH₃COOH)", TEXT),
        ("", TEXT),
        ("# Accuracy vs DFT Baseline", ACCENT),
        ("  → Fundamentals:    MAE < 15 cm⁻¹", TEXT),
        ("  → Overtones:       MAE < 30 cm⁻¹", TEXT),
        ("  → Combinations:    MAE < 25 cm⁻¹", TEXT),
        ("", TEXT),
        ("# Performance", ACCENT),
        ("  → 50-100× speedup over pure DFT", TEXT),
    ])

    # Slide 8: Results - Overtones
    create_content_slide(prs, "$ grep 'overtone' results.json", [
        ("# Anharmonic Analysis", ACCENT),
        ("", TEXT),
        ("  → Gaussian VPT2 framework", TEXT),
        ("  → Overtones: 2ν₁, 2ν₂, ...", TEXT),
        ("  → Combination bands: ν₁+ν₂, ν₁+ν₃, ...", TEXT),
        ("", TEXT),
        ("# Parser Implementation", ACCENT),
        ("  → Regex matching for 2(X) patterns", TEXT),
        ("  → Handle X(1) + Y(1) combinations", TEXT),
        ("  → Separate harmonic/anharmonic sections", TEXT),
        ("", TEXT),
        ("# Challenges", ACCENT),
        ("  → Missing anharmonic data in some calcs", TEXT),
        ("  → Mode numbering consistency", TEXT),
    ])

    # Slide 9: Workflow Pipeline
    create_content_slide(prs, "$ ./run_analysis.py water", [
        ("# End-to-End Workflow", ACCENT),
        ("", TEXT),
        ("  1. Geometry optimization (MACE)", TEXT),
        ("  2. Frequency calculation (Gaussian + ML)", TEXT),
        ("  3. Anharmonic VPT2 (Gaussian)", TEXT),
        ("  4. Parse results (JSON)", TEXT),
        ("  5. Statistical analysis", TEXT),
        ("  6. HTML report generation", TEXT),
        ("", TEXT),
        ("# Automation", ACCENT),
        ("  → Hierarchical directory structure", TEXT),
        ("  → Automatic DFT baseline detection", TEXT),
        ("  → Publication-ready visualizations", TEXT),
    ])

    # Slide 10: Future Directions
    create_content_slide(prs, "$ git log --oneline --graph", [
        ("# Current Status", ACCENT),
        ("  → Proof of concept validated", TEXT),
        ("  → Water & acetic acid benchmarks", TEXT),
        ("  → Automated analysis pipeline", TEXT),
        ("", TEXT),
        ("# Next Steps", ACCENT),
        ("  → Expand molecular test set", TEXT),
        ("  → Support ORCA quantum chemistry", TEXT),
        ("  → Additional ML potentials", TEXT),
        ("  → Optimize dipole model training", TEXT),
        ("", TEXT),
        ("# Vision", ACCENT),
        ("  → Universal ML-QM spectroscopy platform", TEXT),
    ])

    # Save presentation
    output_path = "../mace_gaussian_presentation.pptx"
    prs.save(output_path)
    print(f"✓ Full presentation created: {output_path}")
    print(f"  → 1 title slide + 10 content slides")
    print(f"  → Headers: 28pt")
    print(f"  → Footer on all slides")
    print(f"  → Terminal Dark styling")


if __name__ == "__main__":
    main()
