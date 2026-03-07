#!/usr/bin/env python3
"""
Title slide banner workshop — 4 style variations.
Run: python3.8 presentation/banners/generate_banners.py
Output: presentation/banners/title_banners.pptx
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor

BG     = RGBColor(13,  17,  23)
TEXT   = RGBColor(201, 209, 217)
ACCENT = RGBColor(88,  166, 255)
GREEN  = RGBColor(63,  185, 80)
DIM    = RGBColor(139, 148, 158)
FONT   = "Consolas"
DATE   = "Manuel Ott · Hofer Lab · 2026-03-26"


def blank(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG
    return slide


def footer(slide):
    box = slide.shapes.add_textbox(Inches(1), Inches(6.85), Inches(8), Inches(0.4))
    tf = box.text_frame
    tf.text = DATE
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = FONT
    p.font.color.rgb = DIM
    p.alignment = PP_ALIGN.CENTER


def label(slide, text, color=DIM):
    """Small top-left style label."""
    box = slide.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(4), Inches(0.3))
    tf = box.text_frame
    tf.text = text
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = FONT
    p.font.color.rgb = color


def text_block(slide, lines, x, y, w, h, size=15, align=PP_ALIGN.CENTER):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = False
    for i, (text, color) in enumerate(lines):
        p = tf.paragraphs[i] if i == 0 else tf.add_paragraph()
        p.text = text
        p.font.size = Pt(size)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(0)
        p.alignment = align
    return box


# ── Style 1: Minimal — v1 spirit, centered, breathing room ───────────────────

def style_1(prs):
    slide = blank(prs)
    label(slide, "style 1 — minimal")

    text_block(slide, [
        ("$ ./presentation.sh",         ACCENT),
        ("",                             DIM),
        ("// ML-accelerated IR spectroscopy", DIM),
    ], x=1, y=2.8, w=8, h=1.5, size=28)

    footer(slide)


# ── Style 2: Big ASCII block title ───────────────────────────────────────────

def style_2(prs):
    slide = blank(prs)
    label(slide, "style 2 — block letters")

    # MACE in figlet-style block letters (hand-crafted)
    banner = [
        (r" __  __   _   ___ ___    __ _   _   _   _ ___ ___ _   _   _  _ ",  ACCENT),
        (r"|  \/  | /_\ / __| __|  / _` | /_\ | | | / __/ __| | /_\ | \| |", ACCENT),
        (r"| |\/| |/ _ \ (__| _|  | (_| |/ _ \| |_| \__ \__ \ |/ _ \| .` |", ACCENT),
        (r"|_|  |_/_/ \_\___|___|  \__, /_/ \_\\___/|___/___/_/_/ \_\_|\_|", ACCENT),
        (r"                        |___/                                    ", DIM),
    ]

    text_block(slide, banner, x=0.2, y=1.6, w=9.6, h=1.5, size=10, align=PP_ALIGN.CENTER)

    text_block(slide, [
        ("",                                        DIM),
        ("// ML-accelerated IR spectroscopy",       DIM),
        ("// MACE × Gaussian · ZMQ · VPT2",         DIM),
    ], x=1, y=3.5, w=8, h=1.2, size=16, align=PP_ALIGN.CENTER)

    footer(slide)


# ── Style 3: Terminal box frame ───────────────────────────────────────────────

def style_3(prs):
    slide = blank(prs)
    label(slide, "style 3 — box frame")

    W = 52
    inner = "mace-gaussian  //  ML-accelerated IR spectroscopy"
    sub   = "MACE × Gaussian 16 · ZMQ bridge · VPT2 pipeline"
    pad   = lambda s: s.center(W)

    box_lines = [
        ("┌" + "─" * W + "┐",       DIM),
        ("│" + " " * W + "│",       DIM),
        ("│" + pad("$ ./mace_gaussian.sh") + "│",  ACCENT),
        ("│" + " " * W + "│",       DIM),
        ("│" + pad(inner) + "│",    TEXT),
        ("│" + pad(sub)   + "│",    DIM),
        ("│" + " " * W + "│",       DIM),
        ("└" + "─" * W + "┘",       DIM),
    ]

    text_block(slide, box_lines, x=0.5, y=2.4, w=9, h=2.2, size=14, align=PP_ALIGN.CENTER)

    footer(slide)


# ── Style 4: Molecule-first, command below ────────────────────────────────────

def style_4(prs):
    slide = blank(prs)
    label(slide, "style 4 — molecule-first")

    molecules = [
        ("      O          O   OH   ",       GREEN),
        ("     / \\         ‖   |    ",       GREEN),
        ("    H   H    CH₃-C-O-C₆H₄-COOH",  GREEN),
        ("    H₂O           aspirin  ",      DIM),
    ]

    text_block(slide, molecules, x=1, y=1.6, w=8, h=1.4, size=17, align=PP_ALIGN.CENTER)

    text_block(slide, [
        ("",                                     DIM),
        ("$ python3 mace_gaussian.py water.xyz", ACCENT),
        ("",                                     DIM),
        ("// ML potentials + Gaussian VPT2",     DIM),
    ], x=1, y=3.3, w=8, h=1.5, size=20, align=PP_ALIGN.CENTER)

    footer(slide)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    prs = Presentation()
    prs.slide_width  = Inches(10)
    prs.slide_height = Inches(7.5)

    style_1(prs)
    style_2(prs)
    style_3(prs)
    style_4(prs)

    out = "presentation/banners/title_banners.pptx"
    prs.save(out)
    print(f"✓ Saved: {out}  (4 slides)")
    print("  1 — minimal (v1 spirit)")
    print("  2 — ASCII block letters")
    print("  3 — terminal box frame")
    print("  4 — molecule-first")


if __name__ == "__main__":
    main()
