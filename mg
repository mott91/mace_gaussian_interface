#!/usr/bin/env bash
# MACE-Gaussian (mg) wrapper script
# This script finds cli.py relative to its own location and calls it with Python

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Path to cli.py (should be in the same directory)
CLI_PATH="${SCRIPT_DIR}/cli.py"

# Check if cli.py exists
if [ ! -f "$CLI_PATH" ]; then
    echo "Error: cli.py not found at $CLI_PATH" >&2
    exit 1
fi

# Call Python with cli.py and pass all arguments
exec python3 "$CLI_PATH" "$@"
