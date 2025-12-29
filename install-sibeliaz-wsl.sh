#!/usr/bin/env bash
set -euo pipefail

echo "Installing progressiveMauve via conda in WSL..."

# Ensure conda exists
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: Conda not found inside WSL."
  echo "Install Miniconda or Anaconda first."
  exit 1
fi

# Initialise conda for non-interactive shells
eval "$(conda shell.bash hook)"

# Create env if it doesn't exist
if ! conda env list | grep -q '^mauve_env '; then
  conda create -y -n mauve_env
fi

conda activate mauve_env

# Channels (idempotent)
conda config --add channels conda-forge || true
conda config --add channels bioconda || true
conda config --add channels defaults || true

# Install
conda install -y progressivemauve

echo
echo "Done."
echo "Activate later with: conda activate mauve_env"
