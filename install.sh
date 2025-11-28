#!/bin/bash
# BOLDGenotyper Installation Script for Mac/Linux
# This script sets up the conda environment and installs the package

set -e  # Exit on error

echo "========================================="
echo "BOLDGenotyper Installation"
echo "========================================="
echo ""

# Check for conda
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Step 1: Creating conda environment..."
conda env create -f environment.yml

echo ""
echo "Step 2: Activating environment..."
eval "$(conda shell.bash hook)"
conda activate boldgenotyper

echo ""
echo "Step 3: Installing boldgenotyper package..."
pip install -e .

echo ""
echo "Step 4: Verifying installation..."
if boldgenotyper --version; then
    echo ""
    echo "========================================="
    echo "✓ Installation successful!"
    echo "========================================="
    echo ""
    echo "To use boldgenotyper, activate the environment:"
    echo "  conda activate boldgenotyper"
    echo ""
    echo "Then run:"
    echo "  boldgenotyper --help"
    echo ""
else
    echo ""
    echo "Installation completed, but verification failed."
    echo "Please activate the environment and try running boldgenotyper manually:"
    echo "  conda activate boldgenotyper"
    echo "  boldgenotyper --version"
    exit 1
fi
