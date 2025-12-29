# Dockerfile for SibeliaZ (AMD64 Linux)
# The image runs as amd64 and will be emulated on ARM64 Macs
# Build command: docker build -t sibeliaz-conda:latest .
# Or for cross-platform build: docker buildx build --platform linux/amd64,linux/arm64 -t sibeliaz-conda:latest --load .
# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:latest

# Prevent timezone and locale issues
ENV DEBIAN_FRONTEND=noninteractive

# Create and activate a conda environment with Python 3.11
# Install SibeliaZ and twopaco from bioconda channel
RUN conda create -y -n bioenv python=3.11 && \
    conda install -y -n bioenv -c bioconda sibeliaz twopaco && \
    conda clean --all -y

# Set up the conda environment to be default
SHELL ["conda", "run", "-n", "bioenv", "/bin/bash", "-c"]
ENV PATH="/opt/conda/envs/bioenv/bin:$PATH"

# Create workspace for mounted files
WORKDIR /workspace
VOLUME ["/workspace"]

# Default command - allows bash execution
CMD ["/bin/bash"]