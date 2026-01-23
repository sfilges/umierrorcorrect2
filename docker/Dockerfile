################## BASE IMAGE ######################
# Build this image from the project root:
# docker build -f docker/Dockerfile -t umierrorcorrect2:0.30.2 .

FROM ubuntu:24.04

################## METADATA ######################
LABEL base.image="ubuntu:24.04"
LABEL version="0.30.2"
LABEL software="umierrorcorrect"
LABEL software.version="0.30.2"
LABEL about.summary="Modernized UMI error correction pipeline for barcoded amplicon sequencing"
LABEL about.home="https://github.com/sfilges/umierrorcorrect2"
LABEL about.tags="Genomics,UMI,Bioinformatics"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

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
    && rm -rf /var/lib/apt/lists/*

# Install fastp binary
RUN curl -L http://opengene.org/fastp/fastp > /usr/local/bin/fastp && \
    chmod a+x /usr/local/bin/fastp

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set working directory for application
WORKDIR /app

# Set up Python environment with uv
ENV UV_COMPILE_BYTECODE=1
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy project files and install
COPY . .
RUN uv venv $VIRTUAL_ENV && \
    uv pip install --no-cache .[fast] cutadapt

# Set working directory for data
WORKDIR /data

# Default command
ENTRYPOINT ["umierrorcorrect"]
CMD ["--help"]