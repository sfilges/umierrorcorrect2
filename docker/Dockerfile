################## BASE IMAGE ######################
FROM ubuntu:24.04 AS builder

################## METADATA ######################
LABEL base.image="ubuntu:24.04"
LABEL version="2"
LABEL software="umierrorcorrect"
LABEL software.version="1.0.0"
LABEL about.summary="UMI error correction pipeline for barcoded amplicon sequencing"
LABEL about.home="https://github.com/sfilges/umierrorcorrect"
LABEL about.tags="Genomics,UMI,Bioinformatics"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies and uv in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libcurl4-openssl-dev \
    bzip2 \
    xz-utils \
    gcc \
    g++ \
    python3-dev \
    pigz \
    bwa \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set up Python environment with uv
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1

# Install Python packages using uv (faster than pip)
RUN uv pip install --no-cache \
    umierrorcorrect \
    cutadapt

# Set working directory
WORKDIR /data

# Default command
CMD ["run_umierrorcorrect", "--help"]
