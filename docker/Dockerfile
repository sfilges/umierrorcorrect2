################## BASE IMAGE ######################
FROM ubuntu:24.04

################## METADATA ######################
ARG VERSION="0.32.2"
LABEL base.image="ubuntu:24.04" \
      version="${VERSION}" \
      software="umierrorcorrect2" \
      software.version="${VERSION}" \
      about.summary="Modernized UMI error correction pipeline for barcoded amplicon sequencing" \
      about.home="https://github.com/sfilges/umierrorcorrect2" \
      about.tags="Genomics,UMI,Bioinformatics"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install system dependencies
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
    python3-pip \
    pigz \
    bwa \
    mafft \
    samtools \
    fastp \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:0.5.24 /uv /usr/local/bin/uv

# Set up Python environment with uv
ENV UV_COMPILE_BYTECODE=1 \
    VIRTUAL_ENV=/opt/venv \
    PATH="/opt/venv/bin:$PATH"

WORKDIR /app

# Create a non-root user and set up directories
RUN groupadd -r umicorr && useradd -r -g umicorr -m -d /home/umicorr umicorr && \
    mkdir -p /data /opt/venv && chown -R umicorr:umicorr /data /opt/venv /app

# Switch to non-root user
USER umicorr

# Install dependencies first (better caching)
# We use uv sync --no-install-project to only install dependencies
COPY --chown=umicorr:umicorr pyproject.toml uv.lock ./
RUN uv venv $VIRTUAL_ENV && \
    uv sync --no-dev --no-install-project

# Copy the rest of the application and install it
COPY --chown=umicorr:umicorr . .
RUN uv pip install --no-cache .

WORKDIR /data

# Default command
ENTRYPOINT ["umierrorcorrect2"]