#!/usr/bin/env python3
"""
MACE-Gaussian Interface — Research Presentation v2
March 2026 · ~15 minutes · 13 slides
Run: python generate_pptx_v2.py
Output: presentation_v2.pptx
"""

import json
import os

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor

# ── Color scheme (GitHub Dark) ────────────────────────────────────────────────
BG     = RGBColor(13,  17,  23)   # #0D1117
TEXT   = RGBColor(201, 209, 217)  # #C9D1D9
ACCENT = RGBColor(88,  166, 255)  # #58A6FF  bright blue
GREEN  = RGBColor(63,  185, 80)   # #3FB950
YELLOW = RGBColor(210, 153, 34)   # #D29922
DIM    = RGBColor(139, 148, 158)  # #8B949E
RED    = RGBColor(248, 81,  73)   # #F85149
FONT   = "Consolas"

DATE   = "Manuel Ott · Hofer Lab · 2026-03-11"


# ── Helpers ───────────────────────────────────────────────────────────────────

def blank_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = BG
    return slide


def add_footer(slide, text=DATE):
    box = slide.shapes.add_textbox(Inches(1), Inches(6.85), Inches(8), Inches(0.4))
    tf = box.text_frame
    tf.text = text
    p = tf.paragraphs[0]
    p.font.size = Pt(11)
    p.font.name = FONT
    p.font.color.rgb = DIM
    p.alignment = PP_ALIGN.CENTER


def add_prompt(slide, command, y=0.35):
    box = slide.shapes.add_textbox(Inches(0.5), Inches(y), Inches(9), Inches(0.6))
    tf = box.text_frame
    tf.text = command
    p = tf.paragraphs[0]
    p.font.size = Pt(28)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True


def content_slide(prs, command, lines):
    """Standard slide: command header + list of (text, color) lines."""
    slide = blank_slide(prs)
    add_prompt(slide, command)

    box = slide.shapes.add_textbox(Inches(0.75), Inches(1.15), Inches(8.5), Inches(5.4))
    tf = box.text_frame
    tf.word_wrap = True

    for i, (line, color) in enumerate(lines):
        p = tf.paragraphs[i] if i == 0 else tf.add_paragraph()
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(4)

    add_footer(slide)
    return slide


def two_col_slide(prs, command, left_lines, right_lines, split=4.3):
    """Slide with two column panels."""
    slide = blank_slide(prs)
    add_prompt(slide, command)

    lbox = slide.shapes.add_textbox(Inches(0.5), Inches(1.15), Inches(split - 0.3), Inches(5.4))
    ltf = lbox.text_frame
    ltf.word_wrap = True
    for i, (line, color) in enumerate(left_lines):
        p = ltf.paragraphs[i] if i == 0 else ltf.add_paragraph()
        p.text = line
        p.font.size = Pt(13)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(3)

    rbox = slide.shapes.add_textbox(Inches(split), Inches(1.15), Inches(9.5 - split), Inches(5.4))
    rtf = rbox.text_frame
    rtf.word_wrap = True
    for i, (line, color) in enumerate(right_lines):
        p = rtf.paragraphs[i] if i == 0 else rtf.add_paragraph()
        p.text = line
        p.font.size = Pt(13)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(3)

    add_footer(slide)
    return slide


# ── Data helpers ──────────────────────────────────────────────────────────────

def load_combo_metrics(molecule):
    """Load combo results from analysis_results_harmonic/{molecule}/data/metrics_summary.json.

    Data source: pre-aggregated metrics_summary.json (not comparison_results/ directory walk).
    Sorted by mae_freq ascending so the best-performing combo (lowest error) appears first.
    Per CONTEXT.md line 28: overrides original comparison_results/ approach from discuss-phase.
    """
    summary_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "analysis_results_harmonic",
        molecule,
        "data",
        "metrics_summary.json",
    )
    if not os.path.isfile(summary_path):
        return []
    with open(summary_path) as f:
        data = json.load(f)
    combos = []
    for entry in data.get("comparisons", []):
        combos.append({
            "name": entry.get("name", "unknown"),
            "mae_freq": entry.get("mae_freq", float("inf")),
            "r2_freq": entry.get("r2_freq", 0.0),
            "r2_intensity": entry.get("r2_intensity", 0.0),
            "speedup": entry.get("speedup"),
        })
    # Sort by mae_freq ascending (lowest MAE = best = shown first)
    # Per CONTEXT.md line 26: overrides original "R² descending" from discuss-phase.
    combos.sort(key=lambda c: c["mae_freq"], reverse=False)
    return combos


def parse_combo_name(name):
    """Split combo name into (energy_model, dipole_model)."""
    for suffix in ("_mace_ml", "_espaloma"):
        if name.endswith(suffix):
            return name[: -len(suffix)], suffix[1:]
    return name, "?"


# ── Slides ────────────────────────────────────────────────────────────────────

def slide_title(prs):
    """Slide 0: neofetch-style title."""
    slide = blank_slide(prs)

    # ASCII art — left panel (molecule + project logo)
    ascii_art = [
        "  __  __   _   ___ ___    ",
        " |  \\/  | /_\\ / __| __|   ",
        " | |\\/| |/ _ \\ (__| _|    ",
        " |_|  |_/_/ \\_\\___|___|   ",
        "       GAUSSIAN            ",
        "                           ",
        "  // ML-accelerated        ",
        "  // IR spectroscopy       ",
        "                           ",
        "     O                     ",
        "    / \\                    ",
        "   H   H   H₂O             ",
        "                           ",
        "     O   OH                ",
        "     ‖   |                 ",
        " CH₃-C-O-C₆H₄-COOH        ",
        "   aspirin C₉H₈O₄          ",
    ]

    abox = slide.shapes.add_textbox(Inches(0.4), Inches(1.0), Inches(4.4), Inches(5.6))
    atf = abox.text_frame
    for i, line in enumerate(ascii_art):
        p = atf.paragraphs[i] if i == 0 else atf.add_paragraph()
        p.text = line
        p.font.size = Pt(15)
        p.font.name = FONT
        p.font.color.rgb = GREEN if i < 5 else (ACCENT if i >= 9 else DIM)
        p.space_after = Pt(0)

    # Divider
    div = slide.shapes.add_textbox(Inches(4.9), Inches(1.0), Inches(0.1), Inches(5.6))
    dtf = div.text_frame
    for i in range(17):
        p = dtf.paragraphs[i] if i == 0 else dtf.add_paragraph()
        p.text = "│"
        p.font.size = Pt(15)
        p.font.name = FONT
        p.font.color.rgb = DIM
        p.space_after = Pt(0)

    # Info panel — right
    info = [
        ("mot@hofer-hpc", ACCENT),
        ("─" * 28, DIM),
        ("OS       Hofer Lab · TU München", TEXT),
        ("Shell    zsh + oh-my-zsh", TEXT),
        ("─" * 28, DIM),
        ("Project  mace-gaussian v1.1", GREEN),
        ("Method   VPT2 + ML Dipoles", TEXT),
        ("Engine   Gaussian 16", TEXT),
        ("Bridge   ZMQ (IPC socket)", TEXT),
        ("─" * 28, DIM),
        ("ML       MACE-OMOL-0, MACE-MP", TEXT),
        ("Dipoles  MACE4IR / Espaloma", TEXT),
        ("Combos   4 energy × 2 dipole", TEXT),
        ("─" * 28, DIM),
        ("Tests    131 passing ✓", GREEN),
        ("CI       GitHub Actions ✓", GREEN),
        ("─" * 28, DIM),
        ("▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓", ACCENT),
    ]

    ibox = slide.shapes.add_textbox(Inches(5.1), Inches(1.0), Inches(4.6), Inches(5.6))
    itf = ibox.text_frame
    for i, (line, color) in enumerate(info):
        p = itf.paragraphs[i] if i == 0 else itf.add_paragraph()
        p.text = line
        p.font.size = Pt(14)
        p.font.name = FONT
        p.font.color.rgb = color
        p.space_after = Pt(1)

    add_footer(slide)


def slide_motivation(prs):
    """Slide 1: Why we built this."""
    content_slide(prs, "$ cat motivation.md", [
        ("# Problem", ACCENT),
        ("  \u2192 DFT frequency calculations: slow", TEXT),
        ("  \u2192 Minutes to hours per molecule", TEXT),
        ("  \u2192 Bottleneck for systematic benchmark campaigns", TEXT),
        ("", DIM),
        ("# Solution", ACCENT),
        ("  \u2192 Hybrid ML-QM approach", TEXT),
        ("  \u2192 ML dipoles and energies", TEXT),
        ("  \u2192 Possibly 10\u2013100\u00d7 faster than pure DFT", TEXT),
        ("", DIM),
        ("# Impact", ACCENT),
        ("  \u2192 Enable rapid IR spectral predictions", TEXT),
        ("  \u2192 Maintain DFT-level accuracy", TEXT),
        ("  \u2192 Systematic comparison: energy surface vs dipole model quality", TEXT),
    ])


def slide_ir_theory(prs):
    """Slide 2: What the software computes — computational framing."""
    two_col_slide(prs, "$ man ir_spectroscopy",
        left_lines=[
            ("# What the software computes", ACCENT),
            ("", DIM),
            ("  Harmonic approximation:", TEXT),
            ("  V(x) \u2248 V\u2080 + \u00bdk\u00b7x\u00b2", TEXT),
            ("", DIM),
            ("  \u2192 Diagonalise Hessian \u2192 normal modes", TEXT),
            ("  \u2192 3N\u22126 frequencies \u03c9\u1d62", TEXT),
            ("  \u2192 Intensities \u221d |\u0394\u03bc/\u0394Q|\u00b2", TEXT),
            ("", DIM),
            ("  For each displaced geometry:", TEXT),
            ("  \u2192 energy (forces \u2192 Hessian)", TEXT),
            ("  \u2192 dipole \u03bc\u20d7 (intensity)", TEXT),
            ("  \u2192 Gaussian displaces ~hundreds of times", TEXT),
            ("", DIM),
            ("# Anharmonic (VPT2): ongoing work", ACCENT),
            ("  \u2192 Already running. Overtones + combination bands.", TEXT),
            ("  \u2192 Results pending full validation \u2014 next direction.", DIM),
        ],
        right_lines=[
            ("# What ML provides", ACCENT),
            ("", DIM),
            ("  Each displaced geometry \u2192 ML query:", TEXT),
            ("", DIM),
            ("  [geometry]", GREEN),
            ("      \u2193  MACE energy model", TEXT),
            ("  E, F, Hessian contribution", DIM),
            ("", DIM),
            ("  [geometry]", GREEN),
            ("      \u2193  MACE dipole model", TEXT),
            ("  \u03bc\u20d7 = (\u03bc_x, \u03bc_y, \u03bc_z)", DIM),
            ("", DIM),
            ("  Both injected into Gaussian", TEXT),
            ("  via ZMQ bridge \u2014 next slide.", TEXT),
            ("", DIM),
            ("  IR range: 400\u20134000 cm\u207b\u00b9", DIM),
            ("  1 cm\u207b\u00b9 \u2248 0.124 meV", DIM),
        ],
        split=5.0,
    )


def slide_architecture(prs):
    content_slide(prs, "$ python3 main.py  # system architecture", [
        ("  [molecule.xyz]", GREEN),
        ("       ↓", DIM),
        ("  ┌─────────────────────────────────────────────────────┐", DIM),
        ("  │  GEOMETRY OPTIMIZER  (MACE-OMOL-0 + ASE)           │", DIM),
        ("  └──────────────┬──────────────────────────────────────┘", DIM),
        ("                 ↓", DIM),
        ("  ┌──────────────┴──────────────┐", DIM),
        ("  │ DFT BASELINE                │   ┌─────────────────────────────┐", DIM),
        ("  │ B3LYP/6-31G(d,p)           │   │ ML FREQ CALC                │", DIM),
        ("  │ Gaussian 16                 │   │ MACE energy + forces        │", DIM),
        ("  │ freq + dipoles              │   │ MACE dipole model           │", DIM),
        ("  └──────────────┬──────────────┘   │ ZMQ bridge (zmq_server.py) │", DIM),
        ("                 │                  │                             │", DIM),
        ("                 └──────────────────┤ → fort.7 → Gaussian        │", DIM),
        ("                                    └──────────────┬──────────────┘", DIM),
        ("                                                   ↓", DIM),
        ("  ┌──────────────────────────────────────────────────────────────┐", DIM),
        ("  │ ANALYSIS  eigenvector matching → KDE broadening → HTML      │", DIM),
        ("  └──────────────────────────────────────────────────────────────┘", DIM),
    ])


def slide_dipole_methods(prs):
    two_col_slide(prs, "$ cat dipole_methods.md",
        left_lines=[
            ("# MACE4IR (direct dipole)", ACCENT),
            ("", DIM),
            ("  Input: atomic positions + species", TEXT),
            ("       ↓  E(3)-equivariant GNN", TEXT),
            ("  Output: μ = (μx, μy, μz)", TEXT),
            ("", DIM),
            ("  → Equivariant: respects rotations", TEXT),
            ("  → Trained specifically on dipoles", TEXT),
            ("  → 78-element coverage", TEXT),
            ("  → Custom MACE4IR model (78 el.)", TEXT),
            ("  → Fast — single forward pass", TEXT),
            ("", DIM),
            ("# Why E(3)-equivariance matters", ACCENT),
            ("  Dipole is a vector — it must rotate", TEXT),
            ("  with the molecule. Standard NNs", TEXT),
            ("  don't guarantee this.", TEXT),
        ],
        right_lines=[
            ("# Espaloma (charge-based dipole)", ACCENT),
            ("", DIM),
            ("  Input: molecular graph", TEXT),
            ("       ↓  graph neural network", TEXT),
            ("  Output: partial charges q₁…qₙ", TEXT),
            ("       ↓  μ = Σ qᵢ · rᵢ", TEXT),
            ("  Output: dipole vector", TEXT),
            ("", DIM),
            ("  → Physics-based construction", TEXT),
            ("  → Interpretable (charges visible)", TEXT),
            ("  → Available as pip package", TEXT),
            ("  → No custom model needed", TEXT),
            ("", DIM),
            ("# Energy calculators", ACCENT),
            ("  MACE-OMOL-0  molecules, default", TEXT),
            ("  MACE-MP-0    universal (130 el.)", TEXT),
            ("  MACE-OFF     organic molecules", TEXT),
            ("  MACE-ANI     organic H,C,N,O,S,F", TEXT),
        ],
        split=5.1,
    )


def slide_zmq(prs):
    """Slide 5: ZMQ bridge — single-column pedagogical flow."""
    content_slide(prs, "$ cat zmq_bridge.md", [
        ("# Problem", ACCENT),
        ("  Gaussian needs dipole derivatives for every displaced geometry", TEXT),
        ("  → DFT dipoles: expensive. This is the bottleneck.", TEXT),
        ("", DIM),
        ("# Solution: ZMQ bridge", ACCENT),
        ("  Gaussian 'External' keyword calls gm_helper.py per geometry", TEXT),
        ("  → gm_helper.py sends geometry over ZMQ IPC socket", TEXT),
        ("  → zmq_server.py receives, queries ML, returns results", TEXT),
        ("", DIM),
        ("# How it works", ACCENT),
        ("  [Gaussian]  → writes geometry → calls gm_helper.py", DIM),
        ("                                        ↓  ZMQ (IPC)", DIM),
        ("                                  zmq_server.py", DIM),
        ("                                        ↓  MACE energy + forces", DIM),
        ("                                        ↓  MACE dipole model", DIM),
        ("  [Gaussian]  ← reads fort.7  ←         ↑  returns results", DIM),
        ("", DIM),
        ("# Why it's non-trivial", ACCENT),
        ("  Socket cleanup: LINGER=0 + SIGKILL timeout (no orphan processes)", TEXT),
        ("  Absolute paths: Gaussian needs full path resolved at CLI startup", TEXT),
        ("  Model loading: dipole class remapping via torch.load()", TEXT),
    ])


def slide_mode_matching(prs):
    content_slide(prs, "$ cat mode_matching.md", [
        ("# Problem: ML and DFT don't output modes in the same order", ACCENT),
        ("", DIM),
        ("  DFT:  [500, 1200, 1500, 2900, 3100] cm⁻¹", TEXT),
        ("  ML:   [498, 3095, 1205, 2895, 1498] cm⁻¹   ← different order!", YELLOW),
        ("", DIM),
        ("  → Cannot compare by index. Need to identify which mode is which.", TEXT),
        ("", DIM),
        ("# Solution: Eigenvector dot product", ACCENT),
        ("", DIM),
        ("  Each normal mode has an eigenvector: 3N-dimensional displacement pattern.", TEXT),
        ("  Modes with the same physical motion will have parallel eigenvectors.", TEXT),
        ("", DIM),
        ("  similarity(i, j) = |v_ML[i] · v_DFT[j]|        (0 = orthogonal, 1 = same)", TEXT),
        ("", DIM),
        ("  → Assign each ML mode to its best-matching DFT mode", TEXT),
        ("  → Robust to reordering, robust to small frequency errors", TEXT),
        ("", DIM),
        ("# Example (water)", ACCENT),
        ("  ML 3095 cm⁻¹ → DFT 1200 cm⁻¹  [overlap: 0.027]   ✗ wrong match", RED),
        ("  ML 3095 cm⁻¹ → DFT 3100 cm⁻¹  [overlap: 0.989]   ✓ correct", GREEN),
    ])


def slide_results_overview(prs):
    """Slide 7: HTML report listing with ASCII molecule art."""
    content_slide(prs, "$ ls analysis_results_harmonic/", [
        ("# Spectrum analysis reports", ACCENT),
        ("", DIM),
        ("       O                       O   OH", GREEN),
        ("      / \\                      ‖   |", GREEN),
        ("     H   H   H₂O          CH₃-C-O-C₆H₄-COOH", GREEN),
        ("                               aspirin C₉H₈O₄", DIM),
        ("", DIM),
        ("  water/report.html       3 atoms  · 3 modes  · 8 combos", GREEN),
        ("  aspirin/report.html    19 atoms  · 51 modes · 4 combos", GREEN),
        ("", DIM),
        ("  Each report: regression plots · KDE spectra · mode matching", TEXT),
        ("  Open: $ open analysis_results_harmonic/water/report.html", DIM),
        ("         $ open analysis_results_harmonic/aspirin/report.html", DIM),
    ])


def slide_results_table(prs):
    """Slides 8–9: Harmonic benchmark — water (8 combos) and aspirin (4 combos).

    Rows sorted by MAE(freq) ascending (lowest error = best = shown first).
    Uniform TEXT color — no per-row color logic.
    Data loaded from analysis_results_harmonic/{molecule}/data/metrics_summary.json.
    """
    water = load_combo_metrics("water")
    aspirin = load_combo_metrics("aspirin")

    def table_lines(combos, mol_label):
        lines = [
            (f"# Harmonic benchmark — {mol_label}", ACCENT),
            ("  Ref: B3LYP/6-31G(d,p)   |   sorted by MAE(freq) ascending", DIM),
            ("", DIM),
            ("  Energy model    Dipole      MAE(freq)   R²(freq)    R²(intens)", DIM),
            ("  " + "─" * 60, DIM),
        ]
        for c in combos:
            energy, dipole = parse_combo_name(c["name"])
            row = (
                f"  {energy:<14}  {dipole:<10}"
                f"  {c['mae_freq']:>6.1f} cm⁻¹"
                f"   {c['r2_freq']:.7f}"
                f"   {c['r2_intensity']:.2f}"
            )
            lines.append((row, TEXT))
        return lines

    water_lines = table_lines(water, "water (H₂O, 3 modes)")
    aspirin_lines = table_lines(aspirin, "aspirin (C₉H₈O₄, 51 modes)")

    # Water slide
    content_slide(prs, "$ cat analysis_results_harmonic/water/data/metrics_summary.json | sort-by-mae", water_lines + [
        ("", DIM),
        ("  → Frequencies: all combos excellent (R² > 0.9999)", TEXT),
        ("  → Intensities: mace_ml consistently beats espaloma", TEXT),
    ])

    # Aspirin slide
    content_slide(prs, "$ cat analysis_results_harmonic/aspirin/data/metrics_summary.json | sort-by-mae", aspirin_lines + [
        ("", DIM),
        ("  → mace_omol: best accuracy (MAE < 9 cm⁻¹, R² 0.9998)", TEXT),
        ("  → mace_mp: most speedup (~29×) but worse accuracy", TEXT),
    ])


def slide_scaling(prs):
    """Slide 10: Molecule size vs speedup — ASCII bar chart."""
    MAX_BAR = 20
    molecules = [
        ("water",   3,  1,  "~1×"),
        ("glycine", 10, 6,  "~7–10×"),
        ("aspirin", 21, 20, "~18–29×"),
    ]
    bar_lines = []
    for name, atoms, bar_val, label in molecules:
        bar = "=" * bar_val + " " * (MAX_BAR - bar_val)
        bar_lines.append((
            f"  {name:<10}  ({atoms:>2} atoms)  |{bar}|  {label}",
            TEXT,
        ))

    content_slide(prs, "$ python benchmark.py --scaling", [
        ("# ML speedup scales with molecule size", ACCENT),
        ("  Dipole step: ML vs DFT (B3LYP/6-31G(d,p))", DIM),
        ("", DIM),
        ("  molecule      atoms         speedup", DIM),
        ("  " + "─" * 50, DIM),
        *bar_lines,
        ("", DIM),
        ("  → Water (3 atoms): ML overhead visible — no speedup at this scale", TEXT),
        ("  → Glycine (10 atoms): ~7–10× — target use case", TEXT),
        ("  → Aspirin (21 atoms): ~18–29× — strong benefit", TEXT),
        ("", DIM),
        ("  Speedup source: DFT dipoles O(N³); ML dipoles O(N)", DIM),
        ("  (estimated ranges — not benchmarked to the second)", DIM),
    ])


def slide_status(prs):
    """Slide 9: Research journey — what we built, ran, found, and what's next."""
    two_col_slide(prs, "$ git log --oneline --graph",
        left_lines=[
            ("# What we built", ACCENT),
            ("  ZMQ bridge: ML dipoles \u2192 Gaussian VPT2", TEXT),
            ("  CLI: mace-gaussian run molecule.xyz", DIM),
            ("  4 energy \u00d7 2 dipole models = 8 combinations", TEXT),
            ("", DIM),
            ("# What we ran", ACCENT),
            ("  Harmonic benchmark across 8 molecules", TEXT),
            ("  Water, aspirin, glycine, methane, ethane, NH\u2083, CO", DIM),
            ("", DIM),
            ("# What we found", ACCENT),
            ("  Frequencies: R\u00b2 > 0.999 across all combos \u2713", GREEN),
            ("  Intensities: mace_ml >> espaloma (0.72 vs 0.52)", GREEN),
            ("  Speedup scales with molecule size (aspirin: ~18\u00d7)", GREEN),
        ],
        right_lines=[
            ("# Still open", ACCENT),
            ("", DIM),
            ("  Anharmonic validation:", TEXT),
            ("  \u2192 VPT2 pipeline runs; accuracy TBD", DIM),
            ("  \u2192 Overtones + combination bands", DIM),
            ("", DIM),
            ("  25-molecule benchmark campaign:", TEXT),
            ("  \u2192 Systematic across functional groups", DIM),
            ("  \u2192 Batch runner in progress", DIM),
            ("", DIM),
            ("# Thesis question", ACCENT),
            ("  Energy surface quality", TEXT),
            ("  vs dipole model quality \u2014", TEXT),
            ("  which dominates IR accuracy?", TEXT),
            ("", DIM),
            ("  $ mace-gaussian run molecule.xyz", GREEN),
        ],
        split=5.1,
    )


def slide_questions(prs):
    slide = blank_slide(prs)

    box = slide.shapes.add_textbox(Inches(1), Inches(2.5), Inches(8), Inches(2.0))
    tf = box.text_frame

    p = tf.paragraphs[0]
    p.text = "$ ./questions.sh"
    p.font.size = Pt(36)
    p.font.name = FONT
    p.font.color.rgb = ACCENT
    p.font.bold = True

    p2 = tf.add_paragraph()
    p2.text = ""

    p3 = tf.add_paragraph()
    p3.text = "// Thank you"
    p3.font.size = Pt(18)
    p3.font.name = FONT
    p3.font.color.rgb = DIM

    add_footer(slide)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    prs = Presentation()
    prs.slide_width  = Inches(10)
    prs.slide_height = Inches(7.5)

    slide_title(prs)            # 0  — neofetch header
    slide_motivation(prs)       # 1  — problem / solution / impact
    slide_ir_theory(prs)        # 2  — what is IR, harmonic limit
    slide_architecture(prs)     # 3  — system architecture diagram
    slide_dipole_methods(prs)   # 4  — MACE4IR + Espaloma + energy models
    slide_zmq(prs)              # 5  — ZMQ bridge — pedagogical flow
    slide_mode_matching(prs)    # 6  — eigenvector dot product
    slide_results_overview(prs) # 7  — HTML report listing + ASCII art
    slide_results_table(prs)    # 8–9 — harmonic table (water + aspirin, 2 slides)
    slide_scaling(prs)          # 10 — molecule size vs speedup bar chart
    slide_status(prs)           # 11 — research journey
    slide_questions(prs)        # 12 — questions

    out = "presentation/presentation_v2.pptx"
    prs.save(out)

    print(f"✓ Saved: {out}")
    print(f"  13 slides  (~15 min)")
    print(f"  Slides: title · motivation · IR theory ·")
    print(f"          architecture · dipole methods · ZMQ bridge ·")
    print(f"          mode matching · results overview ·")
    print(f"          results table (water) · results table (aspirin) ·")
    print(f"          scaling · status · questions")


if __name__ == "__main__":
    main()
