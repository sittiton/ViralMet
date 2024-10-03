#!/bin/bash

set -e

echo "Starting ViralMet installation..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Conda first."
    exit 1
fi

# Remove existing environment if it exists
if conda env list | grep -q "viralmet"; then
    echo "Removing existing viralmet environment..."
    conda remove -n viralmet --all -y
fi

# Create new environment
echo "Creating new conda environment..."
conda env create -f environment.yml

# Activate environment and verify installation
echo "Verifying installation..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate viralmet

# Check critical components
commands=("fastqc" "trimmomatic" "bowtie2" "spades.py" "blastn" "diamond" "samtools" "mafft" "iqtree" "kraken2" "prokka")

for cmd in "${commands[@]}"; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Warning: $cmd is not properly installed"
    else
        echo "$cmd is properly installed"
    fi
done

echo "Installation complete. Please use 'conda activate viralmet' to start using the pipeline."
