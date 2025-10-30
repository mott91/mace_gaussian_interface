#!/bin/bash
# Install the dual MACE packages required for this project
# These packages provide both standard MACE and dipole-enabled MACE calculations

set -e  # Exit on error

echo "Installing MACE packages..."

# Install standard MACE package
echo "Installing mace-torch from mace_ML_pkg..."
pip install -e ./mace_ML_pkg

# Install dipole-enabled MACE package
echo "Installing mace-dipole-core from mace_dipole_pkg..."
pip install -e ./mace_dipole_pkg

echo "MACE packages installed successfully!"
echo ""
echo "Verifying installation..."
python -c "import mace; print(f'mace-torch version: {mace.__version__}')"
python -c "import mace_dipole_core; print(f'mace-dipole-core version: {mace_dipole_core.__version__}')" || echo "Note: mace-dipole-core may not expose version"

echo ""
echo "Installation complete! You can now use the MACE calculators."
