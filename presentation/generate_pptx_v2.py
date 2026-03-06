#!/usr/bin/env python3
"""
MACE-Gaussian Interface — Research Presentation v2
March 2026 · ~15 minutes · 12 slides
Run: python generate_pptx_v2.py
Output: presentation_v2.pptx
"""

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
    content_slide(prs, "$ cat motivation.md", [
        ("# Problem", ACCENT),
        ("  → DFT frequency calculations: slow (minutes to hours per molecule)", TEXT),
        ("  → Anharmonic VPT2: even more expensive — scales as O(N⁴)", TEXT),
        ("  → Bottleneck for any systematic benchmark campaign", TEXT),
        ("", DIM),
        ("# Solution", ACCENT),
        ("  → Hybrid ML-QM approach: ML computes dipole derivatives,", TEXT),
        ("    Gaussian 16 handles the VPT2 orchestration", TEXT),
        ("  → Real-time injection via ZMQ bridge — no Gaussian modifications", TEXT),
        ("  → 10-100× faster than pure DFT for dipole-heavy steps", TEXT),
        ("", DIM),
        ("# Impact", ACCENT),
        ("  → Rapid IR spectral predictions across chemical space", TEXT),
        ("  → Systematic benchmark: which ML model combination works best?", TEXT),
        ("  → Thesis question: energy surface vs dipole model quality", TEXT),
    ])


def slide_ir_theory(prs):
    two_col_slide(prs, "$ man ir_spectroscopy",
        left_lines=[
            ("# What is IR spectroscopy?", ACCENT),
            ("", DIM),
            ("  Molecules vibrate at specific", TEXT),
            ("  frequencies determined by:", TEXT),
            ("  → bond stiffness (k)", TEXT),
            ("  → atomic masses (m)", TEXT),
            ("", DIM),
            ("  ω = √(k/m)         (harmonic)", TEXT),
            ("", DIM),
            ("  IR light is absorbed when:", TEXT),
            ("  → photon frequency = vib. frequency", TEXT),
            ("  → dipole moment changes (Δμ ≠ 0)", TEXT),
            ("", DIM),
            ("# Energy units", ACCENT),
            ("  1 cm⁻¹  ≈ 0.124 meV", TEXT),
            ("  kT@300K ≈ 200 cm⁻¹", TEXT),
            ("  IR range: 400–4000 cm⁻¹", TEXT),
        ],
        right_lines=[
            ("# Why do we care?", ACCENT),
            ("", DIM),
            ("  IR spectrum = molecular fingerprint", TEXT),
            ("  → Identify functional groups", TEXT),
            ("  → Confirm synthesis products", TEXT),
            ("  → Study protein conformations", TEXT),
            ("", DIM),
            ("# The harmonic limit", ACCENT),
            ("", DIM),
            ("  V(x) ≈ V₀ + ½k·x²", TEXT),
            ("", DIM),
            ("  Works well for fundamentals.", TEXT),
            ("  Fails for overtones:", TEXT),
            ("  → predicts exactly 2×ω₀", TEXT),
            ("  → reality: slightly less", TEXT),
            ("  → asymmetric potential well", TEXT),
        ],
        split=5.0,
    )


def slide_vpt2(prs):
    two_col_slide(prs, "$ cat vpt2_theory.md",
        left_lines=[
            ("# Beyond the harmonic approximation", ACCENT),
            ("", DIM),
            ("  Real potential = anharmonic:", TEXT),
            ("  V(x) = ½kx² + αx³ + βx⁴ + …", TEXT),
            ("", DIM),
            ("  The cubic term (αx³) causes:", TEXT),
            ("  → overtone at ~2ω (not exactly)", TEXT),
            ("  → frequency shifts", TEXT),
            ("  → mode coupling", TEXT),
            ("", DIM),
            ("# VPT2 — Vibrational Perturbation Theory", ACCENT),
            ("  (2nd order)", TEXT),
            ("", DIM),
            ("  Treat anharmonicity as a perturbation.", TEXT),
            ("  Gaussian computes 3rd + 4th derivatives", TEXT),
            ("  of V and μ numerically (finite diff.).", TEXT),
            ("  → Overtones, combination bands", TEXT),
        ],
        right_lines=[
            ("# What Gaussian needs from ML", ACCENT),
            ("", DIM),
            ("  For each displaced geometry:", TEXT),
            ("  → energy  (from MACE-OMOL/MP)", TEXT),
            ("  → dipole  (from MACE4IR/Espaloma)", TEXT),
            ("", DIM),
            ("  Gaussian displaces the molecule", TEXT),
            ("  hundreds of times and asks the", TEXT),
            ("  external program each time.", TEXT),
            ("", DIM),
            ("# The dipole gap", ACCENT),
            ("", DIM),
            ("  MACE energy models: good forces,", TEXT),
            ("  bad dipoles (not trained for it).", TEXT),
            ("", DIM),
            ("  → Need a dedicated dipole model", TEXT),
            ("  → MACE4IR: trained on dipole data,", TEXT),
            ("    78 elements, E(3)-equivariant", TEXT),
        ],
        split=5.0,
    )


def slide_architecture(prs):
    content_slide(prs, "$ python3 main.py  # system architecture",  [
        ("  [molecule.xyz]", GREEN),
        ("       ↓", DIM),
        ("  ┌─────────────────────────────────────────────────────┐", DIM),
        ("  │  GEOMETRY OPTIMIZER  (MACE-OMOL-0 + ASE)           │", DIM),
        ("  └──────────────┬──────────────────────────────────────┘", DIM),
        ("                 ↓", DIM),
        ("  ┌──────────────┴──────────────┐", DIM),
        ("  │ DFT BASELINE                │   ┌─────────────────────────────┐", DIM),
        ("  │ B3LYP/6-31G(d,p)           │   │ ML FREQ CALC                │", DIM),
        ("  │ gaussian_freq()             │   │ load_mace_calculator()      │", DIM),
        ("  │ parse_gaussian()            │   │ get_dipole_calculator()     │", DIM),
        ("  └──────────────┬──────────────┘   │ setup_zmq_server()          │", DIM),
        ("                 │                  │ launch_gaussian()           │", DIM),
        ("                 └──────────────────┤ zmq_dipole_loop()           │", DIM),
        ("                                    └──────────────┬──────────────┘", DIM),
        ("                                                   ↓", DIM),
        ("  ┌─────────────────────────────────────────────────────────────┐", DIM),
        ("  │ ANALYSIS  eigenvector_match → kde_broaden → metrics → HTML │", DIM),
        ("  └─────────────────────────────────────────────────────────────┘", DIM),
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
    two_col_slide(prs, "$ cat zmq_bridge.md",
        left_lines=[
            ("# The bridge: how it works", ACCENT),
            ("", DIM),
            ("  Gaussian 16: 'External' keyword", TEXT),
            ("  → calls a helper script for each", TEXT),
            ("     displaced geometry", TEXT),
            ("", DIM),
            ("  ┌──────────────────────────────┐", DIM),
            ("  │  Gaussian (Fortran)          │", DIM),
            ("  │  ↓  writes geometry file     │", DIM),
            ("  │  ↓  calls gm_helper.py       │", DIM),
            ("  └────────────┬─────────────────┘", DIM),
            ("               │ ZMQ (IPC socket)", DIM),
            ("  ┌────────────┴─────────────────┐", DIM),
            ("  │  gm_main.py  (Python)        │", DIM),
            ("  │  → MACE energy + forces      │", DIM),
            ("  │  → dipole calculator         │", DIM),
            ("  │  ← returns fort.7 format     │", DIM),
            ("  └──────────────────────────────┘", DIM),
        ],
        right_lines=[
            ("# Key engineering challenges", ACCENT),
            ("", DIM),
            ("  LINGER=0:", TEXT),
            ("  → ZMQ socket must close cleanly", TEXT),
            ("  → Fixed: explicit SIGKILL timeout", TEXT),
            ("", DIM),
            ("  Absolute paths:", TEXT),
            ("  → Gaussian needs full path to", TEXT),
            ("    helper script in .gjf files", TEXT),
            ("  → Resolved at CLI startup via", TEXT),
            ("    MACE_HELPER_SCRIPT_PATH env var", TEXT),
            ("", DIM),
            ("  MACE module loading:", TEXT),
            ("  → Dipole model saved with different", TEXT),
            ("    class paths than runtime", TEXT),
            ("  → mace_loader.py handles remapping", TEXT),
            ("    via torch.load(pickle_module=…)", TEXT),
        ],
        split=5.1,
    )


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


def slide_results_water(prs):
    two_col_slide(prs, "$ cd results/water/ && cat report.md",
        left_lines=[
            ("# Molecule: H₂O (water)", ACCENT),
            ("  3 atoms · 3 normal modes", TEXT),
            ("  Ref: B3LYP/6-31G(d,p)", DIM),
            ("  Analysis: anharmonic (VPT2)", DIM),
            ("", DIM),
            ("# Calculator combinations (8 total)", ACCENT),
            ("", DIM),
            ("  mace_anicc + espaloma", TEXT),
            ("  → MAE  16.8 cm⁻¹  R² 0.9999", GREEN),
            ("", DIM),
            ("  mace_off + espaloma", TEXT),
            ("  → MAE  33.2 cm⁻¹  R² 0.9999", GREEN),
            ("", DIM),
            ("  mace_omol + espaloma", TEXT),
            ("  → MAE  39.0 cm⁻¹  R² 0.9997", TEXT),
            ("", DIM),
            ("  mace_mp + espaloma", TEXT),
            ("  → MAE 129.9 cm⁻¹  R² 0.9971", YELLOW),
        ],
        right_lines=[
            ("# Key observations", ACCENT),
            ("", DIM),
            ("  R² consistently > 0.997 across", GREEN),
            ("  all combinations — frequency", TEXT),
            ("  ordering is well-preserved.", TEXT),
            ("", DIM),
            ("  Best: MACE-ANI + Espaloma", ACCENT),
            ("  → MAE only 16.8 cm⁻¹", TEXT),
            ("  → Trained on organic H,C,N,O,S", TEXT),
            ("  → Small molecule specialist", TEXT),
            ("", DIM),
            ("  Worst: MACE-MP (large universal)", YELLOW),
            ("  → Jack of all trades,", TEXT),
            ("    master of none for organics", TEXT),
            ("", DIM),
            ("  Intensity (R²) low across all:", YELLOW),
            ("  → ML dipole magnitudes need work", TEXT),
            ("  → Active research direction", TEXT),
        ],
        split=5.1,
    )


def slide_results_aspirin(prs):
    two_col_slide(prs, "$ cd results/aspirin/ && cat report.md",
        left_lines=[
            ("# Molecule: Aspirin (C₉H₈O₄)", ACCENT),
            ("  acetylsalicylic acid", DIM),
            ("  19 atoms · 51 normal modes", TEXT),
            ("  Ref: B3LYP/6-31G(d,p)", DIM),
            ("  Analysis: harmonic", DIM),
            ("", DIM),
            ("# Results (mace_omol + espaloma)", ACCENT),
            ("", DIM),
            ("  MAE   41.3 cm⁻¹", TEXT),
            ("  RMSE  72.7 cm⁻¹", TEXT),
            ("  R²    0.9984", GREEN),
            ("  Speedup  ~18×", GREEN),
            ("", DIM),
            ("# 17× more modes than water", ACCENT),
            ("  → Method scales to real molecules", GREEN),
            ("  → R² stays high (0.9984)", TEXT),
            ("  → MAE slightly higher than water", TEXT),
            ("    — expected for complex system", TEXT),
        ],
        right_lines=[
            ("# Why aspirin?", ACCENT),
            ("", DIM),
            ("  Pharmacologically relevant molecule.", TEXT),
            ("  Covers C=O stretch (~1750 cm⁻¹),", TEXT),
            ("  aromatic ring modes, O-H stretch.", TEXT),
            ("", DIM),
            ("  A complete IR identification test.", TEXT),
            ("", DIM),
            ("# Scaling behaviour", ACCENT),
            ("", DIM),
            ("  H₂O:     3 modes  → ~1× speedup", DIM),
            ("  Aspirin: 51 modes → ~18× speedup", GREEN),
            ("", DIM),
            ("  Gaussian DFT scales badly with size.", TEXT),
            ("  ML cost grows much more slowly.", TEXT),
            ("  → Advantage grows with molecule.", TEXT),
            ("", DIM),
            ("  Intensity R²: ~0.05 (still poor)", YELLOW),
            ("  → Frequencies work well,", TEXT),
            ("    intensities are future work.", TEXT),
        ],
        split=5.1,
    )


def slide_status(prs):
    two_col_slide(prs, "$ git log --oneline --graph",
        left_lines=[
            ("# Current status (v1.1 in progress)", ACCENT),
            ("", DIM),
            ("  ✓ ZMQ bridge — stable", GREEN),
            ("  ✓ Anharmonic (VPT2) pipeline", GREEN),
            ("  ✓ CLI: mace-gaussian run/list", GREEN),
            ("  ✓ 4 calculator combinations", GREEN),
            ("  ✓ 131 tests, CI on every push", GREEN),
            ("  ✓ Aspirin, water, CH₄, acoh, BH₃·NH₃", GREEN),
            ("  ✓ HTML reports per molecule", GREEN),
            ("", DIM),
            ("  ⟳ mace_off + mace_anicc wiring", YELLOW),
            ("  ⟳ Batch runner (25-mol campaign)", YELLOW),
            ("  ⟳ PubChem fetcher", YELLOW),
        ],
        right_lines=[
            ("# Thesis question", ACCENT),
            ("", DIM),
            ("  Does energy surface quality or", TEXT),
            ("  dipole model quality dominate", TEXT),
            ("  IR accuracy?", TEXT),
            ("", DIM),
            ("  → Run systematic 25-molecule", TEXT),
            ("    benchmark across functional", TEXT),
            ("    groups and molecular sizes", TEXT),
            ("", DIM),
            ("# Long-term vision", ACCENT),
            ("", DIM),
            ("  Universal ML-QM spectroscopy tool", TEXT),
            ("  any molecule · any combination ·", TEXT),
            ("  one command:", DIM),
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

    slide_title(prs)           # 0 — neofetch header
    slide_motivation(prs)      # 1 — problem / solution / impact
    slide_ir_theory(prs)       # 2 — what is IR, harmonic limit (life scientists)
    slide_vpt2(prs)            # 3 — VPT2, anharmonicity, the dipole gap
    slide_architecture(prs)    # 4 — system architecture diagram
    slide_dipole_methods(prs)  # 5 — MACE4IR + Espaloma + energy models
    slide_zmq(prs)             # 6 — ZMQ bridge + engineering challenges
    slide_mode_matching(prs)   # 7 — eigenvector dot product
    slide_results_water(prs)   # 8 — water: 8 calculator combos
    slide_results_aspirin(prs) # 9 — aspirin: 51 modes, scaling
    slide_status(prs)          # 10 — current status + next steps
    slide_questions(prs)       # 11 — questions

    out = "/home/mot/mace_gaussian/presentation/presentation_v2.pptx"
    prs.save(out)

    print(f"✓ Saved: {out}")
    print(f"  12 slides  (~15 min)")
    print(f"  Slides: title · motivation · IR theory · VPT2 ·")
    print(f"          architecture · dipole methods · ZMQ bridge ·")
    print(f"          mode matching · water results · aspirin results ·")
    print(f"          status · questions")


if __name__ == "__main__":
    main()
