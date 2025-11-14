# MACE-Gaussian Interface Research Presentation

This directory contains scripts and assets for generating a research presentation about the MACE-Gaussian Interface project.

## Structure

```
presentation/
├── assets/
│   ├── code_snippets/    # Generated terminal-style code images
│   ├── flowcharts/       # Generated flowchart diagrams
│   └── plots/            # Exported analysis plots
├── examples/             # Test examples and style exploration
├── generate_assets.py    # Main script to generate all presentation assets
├── generate_pptx.py      # Script to create PowerPoint programmatically
└── requirements.txt      # Python dependencies
```

## Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Generate all assets
python generate_assets.py

# Create PowerPoint
python generate_pptx.py
```

## Asset Generation Methods

### 1. Terminal-Style Code Snippets
- Uses matplotlib with monospace fonts
- Dark terminal theme with syntax highlighting
- Multiple style options (VS Code dark, Monokai, etc.)

### 2. Flowcharts
- Mermaid diagrams (text-based, clean)
- Graphviz (programmatic layout)
- Custom matplotlib drawings

### 3. PowerPoint Generation
- Uses python-pptx library
- Programmatically creates slides
- Embeds generated images
- Professional styling

## Usage

```bash
# Test different code snippet styles
python examples/test_code_styles.py

# Test flowchart generation
python examples/test_flowcharts.py

# Generate final presentation
python generate_pptx.py --output presentation.pptx
```
