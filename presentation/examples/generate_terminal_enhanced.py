#!/usr/bin/env python3
"""
Enhanced Terminal Dark presentation with ASCII art, Nerd Font, and animations.
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE


# Terminal color scheme
BG = RGBColor(13, 17, 23)          # #0D1117 GitHub dark
TEXT = RGBColor(201, 209, 217)     # #C9D1D9
ACCENT = RGBColor(88, 166, 255)    # #58A6FF bright blue
GREEN = RGBColor(63, 185, 80)      # #3FB950 bright green
PURPLE = RGBColor(210, 168, 255)   # #D2A8FF bright purple
RED = RGBColor(248, 81, 73)        # #F85149 bright red
DIM = RGBColor(139, 148, 158)      # #8B949E


def generate_ascii_art(text, style='banner'):
    """Generate ASCII art for titles."""

    if style == 'banner':
        # Simple block letters
        art_map = {
            'M': ['███╗   ███╗', '████╗ ████║', '██╔████╔██║', '██║╚██╔╝██║', '██║ ╚═╝ ██║'],
            'A': [' █████╗ ', '██╔══██╗', '███████║', '██╔══██║', '██║  ██║'],
            'C': [ '██████╗', '██╔════╝', '██║     ', '██║     ', '╚██████╗'],
            'E': ['███████╗', '██╔════╝', '█████╗  ', '██╔══╝  ', '███████╗'],
            'G': [' ██████╗ ', '██╔════╝ ', '██║  ███╗', '██║   ██║', '╚██████╔╝'],
            'U': ['██╗   ██╗', '██║   ██║', '██║   ██║', '██║   ██║', '╚██████╔╝'],
            'S': ['███████╗', '██╔════╝', '███████╗', '╚════██║', '███████║'],
            'I': ['██╗', '██║', '██║', '██║', '██║'],
            'N': ['███╗   ██╗', '████╗  ██║', '██╔██╗ ██║', '██║╚██╗██║', '██║ ╚████║'],
            'R': ['██████╗ ', '██╔══██╗', '██████╔╝', '██╔══██╗', '██║  ██║'],
            'F': ['███████╗', '██╔════╝', '█████╗  ', '██╔══╝  ', '██║     '],
            'O': [' ██████╗ ', '██╔═══██╗', '██║   ██║', '██║   ██║', '╚██████╔╝'],
            'T': ['████████╗', '╚══██╔══╝', '   ██║   ', '   ██║   ', '   ██║   '],
            'W': ['██╗    ██╗', '██║    ██║', '██║ █╗ ██║', '██║███╗██║', '╚███╔███╔╝'],
            'L': ['██╗     ', '██║     ', '██║     ', '██║     ', '███████╗'],
            ' ': ['   ', '   ', '   ', '   ', '   '],
        }

        lines = ['', '', '', '', '']
        for char in text.upper():
            if char in art_map:
                for i, part in enumerate(art_map[char]):
                    lines[i] += part + ' '

        return '\n'.join(lines)

    elif style == 'simple':
        # Simpler ASCII style
        return f"""
╔══════════════════════════════════════╗
║  {text.upper().center(36)}  ║
╚══════════════════════════════════════╝
"""

    elif style == 'double':
        # Double line box
        return f"""
╔═══════════════════════════════════════╗
║                                       ║
║  {text.upper().center(35)}  ║
║                                       ║
╚═══════════════════════════════════════╝
"""


def add_blinking_cursor(slide, left, top):
    """Add a blinking cursor with animation."""
    cursor_box = slide.shapes.add_textbox(left, top, Inches(0.2), Inches(0.3))
    tf = cursor_box.text_frame
    tf.text = "█"
    p = tf.paragraphs[0]
    p.font.size = Pt(14)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN

    # Add blink animation
    # PowerPoint animations: 0=appear, 1=fly_in, etc.
    # We'll use emphasis > blink
    try:
        from pptx.enum.action import PP_ACTION
        # Add fade in/out animation to simulate blinking
        # Note: python-pptx has limited animation support, but we can add the shape
        # and the user can manually add blink animation in PowerPoint
    except:
        pass

    return cursor_box


def create_terminal_slide(prs, title_cmd, ascii_title, content_lines=None, slide_num=1):
    """Create a terminal-style slide."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Top bar (like terminal title bar)
    top_bar = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(0), Inches(0), Inches(10), Inches(0.3)
    )
    top_bar.fill.solid()
    top_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    top_bar.line.fill.background()

    # Terminal title bar text
    bar_text = slide.shapes.add_textbox(Inches(0.3), Inches(0.05), Inches(9), Inches(0.2))
    tf = bar_text.text_frame
    tf.text = f"  presentation.sh — user@research-group"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = 'Menlo'
    p.font.color.rgb = DIM

    # Command prompt
    prompt_box = slide.shapes.add_textbox(Inches(0.4), Inches(0.5), Inches(9), Inches(0.3))
    tf = prompt_box.text_frame
    tf.text = f"user@research-group:~/presentation$ {title_cmd}"
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN
    p.font.bold = True

    # ASCII art title
    if ascii_title:
        ascii_box = slide.shapes.add_textbox(Inches(0.6), Inches(1.0), Inches(8.8), Inches(1.5))
        tf = ascii_box.text_frame
        tf.text = ascii_title
        p = tf.paragraphs[0]
        p.font.size = Pt(11)
        p.font.name = 'Menlo'
        p.font.color.rgb = ACCENT
        p.line_spacing = 0.8

    # Content area
    if content_lines:
        y_pos = 2.7
        for line_type, text, color_name in content_lines:
            color = {
                'green': GREEN,
                'blue': ACCENT,
                'purple': PURPLE,
                'text': TEXT,
                'dim': DIM,
                'red': RED,
            }.get(color_name, TEXT)

            text_box = slide.shapes.add_textbox(Inches(0.6), Inches(y_pos), Inches(8.8), Inches(0.35))
            tf = text_box.text_frame
            tf.text = text
            p = tf.paragraphs[0]
            p.font.size = Pt(13)
            p.font.name = 'Menlo'
            p.font.color.rgb = color

            y_pos += 0.4

    # Bottom status line (like vim/tmux)
    status_bar = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(0), Inches(7.2), Inches(10), Inches(0.3)
    )
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.3), Inches(7.25), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = f"[{slide_num}/{10}] -- INSERT -- MACE-Gaussian Interface"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT

    return slide


def create_enhanced_terminal_dark():
    """Create enhanced terminal dark presentation."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    # ========================================================================
    # Slide 1: Title with big ASCII art
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Terminal header
    header_box = slide.shapes.add_textbox(Inches(0.3), Inches(0.3), Inches(9), Inches(0.4))
    tf = header_box.text_frame
    tf.text = "user@research-group:~/presentation$ ./run_presentation.sh"
    p = tf.paragraphs[0]
    p.font.size = Pt(13)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN
    p.font.bold = True

    # ASCII art title
    ascii_art = generate_ascii_art("MACE", 'banner')
    ascii_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.2), Inches(9), Inches(2))
    tf = ascii_box.text_frame
    tf.text = ascii_art
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = 'Menlo'
    p.font.color.rgb = ACCENT
    p.line_spacing = 0.9

    # Subtitle
    subtitle_box = slide.shapes.add_textbox(Inches(0.8), Inches(3.3), Inches(8), Inches(0.5))
    tf = subtitle_box.text_frame
    tf.text = ">>> Gaussian Interface <<<"
    p = tf.paragraphs[0]
    p.font.size = Pt(20)
    p.font.name = 'Menlo'
    p.font.color.rgb = PURPLE
    p.font.bold = True
    p.alignment = PP_ALIGN.CENTER

    # Description with icons (using unicode)
    desc_box = slide.shapes.add_textbox(Inches(1.5), Inches(4.2), Inches(7), Inches(1.5))
    tf = desc_box.text_frame
    tf.text = """⚡ ML-Accelerated Anharmonic IR Spectroscopy
🔬 Hybrid ML-QM Approach
🚀 10-100× Speedup"""
    p = tf.paragraphs[0]
    p.font.size = Pt(14)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT
    p.line_spacing = 1.4

    # Footer
    footer_box = slide.shapes.add_textbox(Inches(1), Inches(6.5), Inches(8), Inches(0.6))
    tf = footer_box.text_frame
    tf.text = "Your Name  │  Research Group  │  2025-11-15"
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = 'Menlo'
    p.font.color.rgb = DIM
    p.alignment = PP_ALIGN.CENTER

    # Blinking cursor
    cursor_box = slide.shapes.add_textbox(Inches(4.8), Inches(6.9), Inches(0.2), Inches(0.3))
    tf = cursor_box.text_frame
    tf.text = "█"
    p = tf.paragraphs[0]
    p.font.size = Pt(14)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN

    # ========================================================================
    # Slide 2: cat motivation.md
    # ========================================================================
    ascii_title = generate_ascii_art("MOTIVATION", 'simple')
    content = [
        ('header', '# Problem Statement', 'purple'),
        ('text', '  → DFT anharmonic calculations: expensive', 'text'),
        ('text', '  → Hours to days per molecule', 'text'),
        ('text', '  → Limits high-throughput applications', 'text'),
        ('', '', 'text'),
        ('header', '# Our Solution', 'purple'),
        ('text', '  → ML potentials for fast components', 'green'),
        ('text', '  → Gaussian for anharmonic corrections', 'green'),
        ('text', '  → ZMQ bridge for real-time communication', 'blue'),
        ('', '', 'text'),
        ('header', '# Impact', 'purple'),
        ('text', '  ⚡ 10-100× computational speedup', 'green'),
        ('text', '  ✓ Comparable accuracy to pure DFT', 'green'),
    ]
    create_terminal_slide(prs, 'cat motivation.md', ascii_title, content, slide_num=2)

    # ========================================================================
    # Slide 3: ls -la architecture/
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Command
    cmd_box = slide.shapes.add_textbox(Inches(0.4), Inches(0.4), Inches(9), Inches(0.3))
    tf = cmd_box.text_frame
    tf.text = "user@research-group:~/presentation$ ls -la architecture/"
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN
    p.font.bold = True

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.6), Inches(0.9), Inches(8.8), Inches(0.4))
    tf = title_box.text_frame
    tf.text = "╔═══════════════════════════════╗\n║     SYSTEM ARCHITECTURE       ║\n╚═══════════════════════════════╝"
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = 'Menlo'
    p.font.color.rgb = ACCENT
    p.line_spacing = 0.9

    # File listing
    listing_box = slide.shapes.add_textbox(Inches(0.6), Inches(1.8), Inches(8.8), Inches(4.5))
    tf = listing_box.text_frame
    listing_text = """total 5
drwxr-xr-x  5 user  group   160  Nov 15 12:00  .
drwxr-xr-x 12 user  group   384  Nov 15 11:00  ..
-rw-r--r--  1 user  group  2.1K  Nov 15 10:30  workflow.png
-rw-r--r--  1 user  group  3.4K  Nov 15 10:31  zmq_bridge.png
-rw-r--r--  1 user  group  1.8K  Nov 15 10:32  module_swap.png

📁 XYZ Input → 🔧 MACE Optimization → ⚛️  Gaussian Freq
                                            ↓
                                      🔌 ZMQ Bridge
                                            ↓
                                      🤖 ML Dipoles
                                            ↓
                        📊 Analysis → 📄 HTML Report"""
    tf.text = listing_text
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT
    p.line_spacing = 1.3

    # Status bar
    status_bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(7.2), Inches(10), Inches(0.3))
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.3), Inches(7.25), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = "[3/10] -- NORMAL -- architecture/"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT

    # ========================================================================
    # Slide 4: git diff --stat results
    # ========================================================================
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG

    # Command
    cmd_box = slide.shapes.add_textbox(Inches(0.4), Inches(0.4), Inches(9), Inches(0.3))
    tf = cmd_box.text_frame
    tf.text = "user@research-group:~/presentation$ git diff --stat dft ml"
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = 'Menlo'
    p.font.color.rgb = GREEN
    p.font.bold = True

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.6), Inches(1.0), Inches(8.8), Inches(0.6))
    tf = title_box.text_frame
    tf.text = "═══════════════════════════════\n     RESULTS & COMPARISON\n═══════════════════════════════"
    p = tf.paragraphs[0]
    p.font.size = Pt(12)
    p.font.name = 'Menlo'
    p.font.color.rgb = PURPLE
    p.line_spacing = 0.9
    p.alignment = PP_ALIGN.CENTER

    # Git diff style results
    results_box = slide.shapes.add_textbox(Inches(0.8), Inches(2.0), Inches(8.4), Inches(4.5))
    tf = results_box.text_frame
    results_text = """Comparing: DFT (B3LYP) vs ML (MACE-MP + Espaloma)

 computation_time    | ████████████████████ -95%
 frequency_accuracy  | ██ +2% MAE
 intensity_accuracy  | ███ +5% RMSE
 mode_matching       | ████████████████████ 98%

 Summary:
 ───────────────────────────────────────────
  Metric               DFT          ML
 ───────────────────────────────────────────
  Runtime            8.5 hrs      6.2 min
  Freq MAE              --        12 cm⁻¹
  R²                    --        0.996
  Speed Factor          1×         ~100×
 ───────────────────────────────────────────

 ✓ Successfully validated against DFT baseline
 ✓ All fundamental modes matched
 ✓ Overtones and combinations included"""
    tf.text = results_text
    p = tf.paragraphs[0]
    p.font.size = Pt(10)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT
    p.line_spacing = 1.2

    # Status bar
    status_bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(7.2), Inches(10), Inches(0.3))
    status_bar.fill.solid()
    status_bar.fill.fore_color.rgb = RGBColor(22, 27, 34)
    status_bar.line.fill.background()

    status_text = slide.shapes.add_textbox(Inches(0.3), Inches(7.25), Inches(9), Inches(0.2))
    tf = status_text.text_frame
    tf.text = "[4/10] -- NORMAL -- results.md"
    p = tf.paragraphs[0]
    p.font.size = Pt(9)
    p.font.name = 'Menlo'
    p.font.color.rgb = TEXT

    # Save
    prs.save("../template_terminal_enhanced.pptx")
    print("✓ Enhanced Terminal Dark template created (4 slides)")
    print("\nFeatures:")
    print("  • ASCII art titles")
    print("  • Menlo Nerd Font")
    print("  • Terminal commands (cat, ls -la, git diff)")
    print("  • Unicode icons (⚡🔬🚀📁🔧)")
    print("  • Status bars (vim-style)")
    print("  • Blinking cursor placeholders")
    print("\nNote: Add blink animation manually in PowerPoint:")
    print("  1. Select the cursor (█)")
    print("  2. Animations → Add Animation → Emphasis → Pulse")
    print("  3. Set to 'Repeat Until End of Slide'")


def main():
    """Generate enhanced terminal presentation."""
    print("Generating enhanced terminal dark presentation...")
    print("=" * 60)
    create_enhanced_terminal_dark()
    print("\n" + "=" * 60)
    print("Done! Check template_terminal_enhanced.pptx")


if __name__ == "__main__":
    main()
