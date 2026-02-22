# Stage 1: Build stage
FROM ghcr.io/astral-sh/uv:0.5.24-python3.12-bookworm-slim AS builder

# Set environment variables for uv
ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    VIRTUAL_ENV=/opt/venv

WORKDIR /app

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies using uv
# We use cache mounts to keep the image small and speed up repeated builds
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=umierrorcorrect2,target=umierrorcorrect2 \
    --mount=type=bind,source=README.md,target=README.md \
    uv venv $VIRTUAL_ENV && \
    uv sync --no-dev --no-install-project

# Copy source and install the project
COPY . .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install .

# Stage 2: Final runtime stage
FROM python:3.12-slim-bookworm

ARG VERSION="0.33.2"
LABEL base.image="python:3.12-slim-bookworm" \
      version="${VERSION}" \
      software="umierrorcorrect2" \
      software.version="${VERSION}" \
      about.summary="Modernized UMI error correction pipeline for barcoded amplicon sequencing" \
      about.home="https://github.com/sfilges/umierrorcorrect2" \
      about.tags="Genomics,UMI,Bioinformatics"

# Prevent interactive prompts and Python bloat
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    PATH="/opt/venv/bin:$PATH" \
    VIRTUAL_ENV=/opt/venv

WORKDIR /app

# Install runtime dependencies (bioinformatics tools and shared libraries)
# Essential libraries for tools (bwa, samtools, fastp) and Python libs (pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
    pigz \
    bwa \
    mafft \
    samtools \
    fastp \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    libcurl4 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user and set up directories
RUN groupadd -r umicorr && useradd -r -g umicorr -m -d /home/umicorr umicorr && \
    mkdir -p /data && chown -R umicorr:umicorr /data /app

# Copy the pre-built virtual environment from the builder stage
COPY --from=builder --chown=umicorr:umicorr /opt/venv /opt/venv

# Switch to non-root user
USER umicorr

WORKDIR /data

# Default command
ENTRYPOINT ["umierrorcorrect2"]
