#!/usr/bin/env python3
"""
Simplified Terminal Dark presentation - clean and minimal.
- No macOS title bar
- No borders
- Meslo Nerd Font
- Cleaner layout
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE


# Terminal color scheme
BG = RGBColor(13, 17, 23)          # #0D1117 GitHub dark
TEXT = RGBColor(201, 209, 217)     # #C9D1D9 light gray/white
ACCENT = RGBColor(88, 166, 255)    # #58A6FF bright blue
GREEN = RGBColor(63, 185, 80)      # #3FB950 bright green
PURPLE = RGBColor(210, 168, 255)   # #D2A8FF bright purple
YELLOW = RGBColor(229, 192, 123)   # #E5C07B bright yellow
CYAN = RGBColor(86, 182, 194)      # #56B6C2 bright cyan
RED = RGBColor(248, 81, 73)        # #F85149 bright red
DIM = RGBColor(139, 148, 158)      # #8B949E dim gray

FONT_NAME = 'Meslo Nerd Font'


def add_p10k_prompt(slide, y_pos, command="", show_command=True):
    """Add Powerlevel10k-style prompt."""
    # Top line of prompt with segments
    prompt_top = slide.shapes.add_textbox(Inches(0.5), Inches(y_pos), Inches(9), Inches(0.25))
    tf = prompt_top.text_frame
    tf.clear()
    p = tf.add_paragraph()

    # Corner
    run = p.add_run()
    run.text = "╭─"
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = DIM

    # User@host segment
    run = p.add_run()
    run.text = "  mot@rs1-c724 "
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = CYAN
    run.font.bold = True

    # Separator
    run = p.add_run()
    run.text = " "
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = DIM

    # Directory segment
    run = p.add_run()
    run.text = " ~/presentation "
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = PURPLE

    # Separator
    run = p.add_run()
    run.text = " "
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = DIM

    # Git branch
    run = p.add_run()
    run.text = " master "
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = GREEN

    # Status icon
    run = p.add_run()
    run.text = "⚡"
    run.font.name = FONT_NAME
    run.font.size = Pt(11)
    run.font.color.rgb = YELLOW

    # Bottom line with arrow prompt
    if show_command:
        prompt_bottom = slide.shapes.add_textbox(Inches(0.5), Inches(y_pos + 0.25), Inches(9), Inches(0.25))
        tf = prompt_bottom.text_frame
        tf.clear()
        p = tf.add_paragraph()

        # Corner and arrow
        run = p.add_run()
        run.text = "╰─"
        run.font.name = FONT_NAME
        run.font.size = Pt(11)
        run.font.color.rgb = DIM

        run = p.add_run()
        run.text = "❯ "
        run.font.name = FONT_NAME
        run.font.size = Pt(12)
        run.font.color.rgb = GREEN
        run.font.bold = True

        # Command
        if command:
            run = p.add_run()
            run.text = command
            run.font.name = FONT_NAME
            run.font.size = Pt(11)
            run.font.color.rgb = TEXT


def generate_ascii_art(text, style='banner'):
    """Generate ASCII art for titles."""
    if style == 'banner':
        art_map = {
            'M': ['███╗   ███╗', '████╗ ████║', '██╔████╔██║', '██║╚██╔╝██║', '██║ ╚═╝ ██║'],
            'A': [' █████╗ ', '██╔══██╗', '███████║', '██╔══██║', '██║  ██║'],
            'C': ['██████╗ ', '██╔════╝', '██║     ', '██║     ', '╚██████╗'],
            'E': ['███████╗', '██╔════╝', '█████╗  ', '██╔══╝  ', '███████╗'],
            'G': [' ██████╗ ', '██╔════╝ ', '██║  ███╗', '██║   ██║', '╚██████╔╝'],
            'U': ['██╗   ██╗', '██║   ██║', '██║   ██║', '██║   ██║', '╚██████╔╝'],
            'S': ['███████╗', '██╔════╝', '███████╗', '╚════██║', '███████║'],
            'I': ['██╗', '██║', '██║', '██║', '██║'],
            'N': ['███╗   ██╗', '████╗  ██║', '██╔██╗ ██║', '██║╚██╗██║', '██║ ╚████║'],
            ' ': ['     ', '     ', '     ', '     ', '     '],
        }

        lines = ['', '', '', '', '']
        for char in text.upper():
            if char in art_map:
                for i, part in enumerate(art_map[char]):
                    lines[i] += part + ' '

        return '\n'.join(lines)

    elif style == 'simple':
        return f"""╔══════════════════════════════════════╗
║  {text.upper().center(36)}  ║
╚══════════════════════════════════════╝"""


def create_clean_terminal_dark():
    """Create clean, simplified terminal dark presentation."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    # ========================================================================
    # Slide 1: Title - Simple and clean
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Just the prompt, no command
    add_p10k_prompt(slide, Inches(0.5), show_command=False)

    # ASCII art title - ALL LINES COLORED
    ascii_art = generate_ascii_art("MACE", 'banner')
    ascii_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.3), Inches(8.4), Inches(1.8))
    tf = ascii_box.text_frame
    tf.text = ascii_art
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = FONT_NAME
    p.font.color.rgb = ACCENT
    p.line_spacing = 0.9

    # Subtitle
    subtitle_box = slide.shapes.add_textbox(Inches(1.5), Inches(3.2), Inches(7), Inches(0.5))
    tf = subtitle_box.text_frame
    tf.text = ">>> Gaussian Interface <<<"
    p = tf.paragraphs[0]
    p.font.size = Pt(20)
    p.font.name = FONT_NAME
    p.font.color.rgb = PURPLE
    p.font.bold = True
    p.alignment = PP_ALIGN.CENTER

    # Description with icons
    desc_box = slide.shapes.add_textbox(Inches(1.8), Inches(4.0), Inches(6.4), Inches(1.2))
    tf = desc_box.text_frame
    tf.clear()

    lines_data = [
        ("⚡ ML-Accelerated Anharmonic IR Spectroscopy", YELLOW),
        ("🔬 Hybrid ML-QM Approach", CYAN),
        ("🚀 10-100× Speedup", GREEN),
    ]

    for i, (text, color) in enumerate(lines_data):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = text
        p.font.size = Pt(13)
        p.font.name = FONT_NAME
        p.font.color.rgb = color
        p.line_spacing = 1.4

    # Footer
    footer_box = slide.shapes.add_textbox(Inches(1.5), Inches(6.0), Inches(7), Inches(0.5))
    tf = footer_box.text_frame
    tf.text = "Your Name  │  Research Group  │  2025-11-15"
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = FONT_NAME
    p.font.color.rgb = DIM
    p.alignment = PP_ALIGN.CENTER

    # Blinking cursor
    cursor_box = slide.shapes.add_textbox(Inches(4.9), Inches(6.8), Inches(0.2), Inches(0.3))
    tf = cursor_box.text_frame
    tf.text = "█"
    p = tf.paragraphs[0]
    p.font.size = Pt(14)
    p.font.name = FONT_NAME
    p.font.color.rgb = GREEN

    # ========================================================================
    # Slide 2: cat motivation.md
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # P10k prompt with command
    add_p10k_prompt(slide, Inches(0.5), "cat motivation.md", show_command=True)

    # ASCII title
    ascii_title = generate_ascii_art("MOTIVATION", 'simple')
    ascii_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(0.6))
    tf = ascii_box.text_frame
    tf.text = ascii_title
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = FONT_NAME
    p.font.color.rgb = ACCENT
    p.line_spacing = 0.9

    # Content
    content_box = slide.shapes.add_textbox(Inches(0.8), Inches(2.0), Inches(8.4), Inches(5.0))
    tf = content_box.text_frame
    tf.clear()

    content_data = [
        ("# Problem Statement", PURPLE),
        ("  → DFT anharmonic calculations: expensive", TEXT),
        ("  → Hours to days per molecule", TEXT),
        ("  → Limits high-throughput applications", TEXT),
        ("", TEXT),
        ("# Our Solution", PURPLE),
        ("  → ML potentials for fast components", GREEN),
        ("  → Gaussian for anharmonic corrections", GREEN),
        ("  → ZMQ bridge for real-time communication", CYAN),
        ("", TEXT),
        ("# Impact", PURPLE),
        ("  ⚡ 10-100× computational speedup", YELLOW),
        ("  ✓ Comparable accuracy to pure DFT", GREEN),
    ]

    for i, (text, color) in enumerate(content_data):
        if i > 0:
            tf.add_paragraph()
        p = tf.paragraphs[i]
        p.text = text
        p.font.size = Pt(12)
        p.font.name = FONT_NAME
        p.font.color.rgb = color

    # Status bar (keep on content slides)
    status_bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
        Inches(0), Inches(7.25), Inches(10), Inches(0.25))
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.5), Inches(7.28), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = "[2/10]  NORMAL  motivation.md"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = FONT_NAME
    p.font.color.rgb = TEXT

    # ========================================================================
    # Slide 3: ls -la architecture/
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # P10k prompt
    add_p10k_prompt(slide, Inches(0.5), "ls -la architecture/", show_command=True)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(0.6))
    tf = title_box.text_frame
    tf.text = "╔═══════════════════════════════╗\n║     SYSTEM ARCHITECTURE       ║\n╚═══════════════════════════════╝"
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = FONT_NAME
    p.font.color.rgb = CYAN
    p.line_spacing = 0.9

    # Listing
    listing_box = slide.shapes.add_textbox(Inches(0.8), Inches(2.0), Inches(8.4), Inches(4.8))
    tf = listing_box.text_frame
    listing_text = """total 5
drwxr-xr-x  5 mot  group   160  Nov 15 12:00  .
drwxr-xr-x 12 mot  group   384  Nov 15 11:00  ..
-rw-r--r--  1 mot  group  2.1K  Nov 15 10:30  workflow.png
-rw-r--r--  1 mot  group  3.4K  Nov 15 10:31  zmq_bridge.png
-rw-r--r--  1 mot  group  1.8K  Nov 15 10:32  module_swap.png

📁 Input → 🔧 MACE Opt → ⚛️  Gaussian
                            ↓
                      🔌 ZMQ Bridge
                            ↓
                        🤖 ML Dipole
                            ↓
          📊 Analysis → 📄 Report"""
    tf.text = listing_text
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = FONT_NAME
    p.font.color.rgb = TEXT
    p.line_spacing = 1.3

    # Status bar
    status_bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
        Inches(0), Inches(7.25), Inches(10), Inches(0.25))
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.5), Inches(7.28), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = "[3/10] 📁 NORMAL  architecture/"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = FONT_NAME
    p.font.color.rgb = TEXT

    # ========================================================================
    # Slide 4: git diff --stat
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # P10k prompt
    add_p10k_prompt(slide, Inches(0.5), "git diff --stat dft ml", show_command=True)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.8), Inches(1.2), Inches(8.4), Inches(0.6))
    tf = title_box.text_frame
    tf.text = "═══════════════════════════════\n     RESULTS & COMPARISON\n═══════════════════════════════"
    p = tf.paragraphs[0]
    p.font.size = Pt(12)
    p.font.name = FONT_NAME
    p.font.color.rgb = PURPLE
    p.line_spacing = 0.9
    p.alignment = PP_ALIGN.CENTER

    # Results
    results_box = slide.shapes.add_textbox(Inches(0.9), Inches(2.0), Inches(8.2), Inches(4.8))
    tf = results_box.text_frame
    results_text = """Comparing: DFT (B3LYP) vs ML (MACE-MP + Espaloma)

 computation_time    | ████████████████████ -95%
 frequency_accuracy  | ██ +2% MAE
 intensity_accuracy  | ███ +5% RMSE
 mode_matching       | ████████████████████ 98%

 Summary:
 ──────────────────────────────────────
  Metric            DFT        ML
 ──────────────────────────────────────
  Runtime         8.5 hrs    6.2 min
  Freq MAE           --      12 cm⁻¹
  R²                 --      0.996
  Speed Factor       1×       ~100×
 ──────────────────────────────────────

 ✓ Validated against DFT baseline
 ✓ All fundamental modes matched"""
    tf.text = results_text
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = FONT_NAME
    p.font.color.rgb = TEXT
    p.line_spacing = 1.15

    # Status bar
    status_bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
        Inches(0), Inches(7.25), Inches(10), Inches(0.25))
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.5), Inches(7.28), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = "[4/10]  master  results.md"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = FONT_NAME
    p.font.color.rgb = TEXT

    # Save
    prs.save("../template_terminal_enhanced.pptx")
    print("✓ Clean Terminal Dark template created (4 slides)")
    print("\nChanges:")
    print("  • Removed macOS title bar (no traffic lights)")
    print("  • Removed slide borders")
    print("  • Changed font to 'Meslo Nerd Font'")
    print("  • Title slide: no command, no status bar")
    print("  • Content slides: keep status bar for context")
    print("\nBlinking cursor: Select █ → Animations → Pulse → Repeat")


def main():
    """Generate clean terminal presentation."""
    print("Generating clean terminal dark presentation...")
    print("=" * 60)
    create_clean_terminal_dark()
    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
