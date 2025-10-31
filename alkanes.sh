#!/usr/bin/env zsh
# Hardcoded batch run for n-alkanes (methane \u2192 octane)

set -e  # stop if any command fails

alias py='python3'

echo "=== Starting alkane calculations ==="

py cli.py run C2H6_ase.xyz   # ethane
py cli.py run C3H8_ase.xyz   # propane
py cli.py run n_alkane_4C.xyz   # butane
py cli.py run n_alkane_5C.xyz   # pentane
py cli.py run n_alkane_6C.xyz   # hexane
py cli.py run n_alkane_7C.xyz   # heptane
py cli.py run n_alkane_8C.xyz   # octane

echo "=== All alkane calculations complete ==="
