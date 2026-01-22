# UMIErrorCorrect

[![PyPI version](https://badge.fury.io/py/umierrorcorrect.svg)](https://badge.fury.io/py/umierrorcorrect)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

A modern, high-performance pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI).

This package is a **complete modernization** of the original [UMIErrorCorrect](https://github.com/stahlberggroup/umierrorcorrect) published in *Clinical Chemistry* (2022).

## Key Features

- **High Performance**: Parallel processing of genomic regions and fastp-based preprocessing.
- **Modern Tooling**: Built with `typer`, `pydantic`, `loguru`, and `hatch`.
- **Easy Installation**: Fully PEP 621 compliant, installable via `pip` or `uv`.
- **Comprehensive**: From raw FASTQ to error-corrected VCFs and consensus statistics.
- **Robust**: Extensive test suite and type safety.

## Dependencies

- `fastp` for preprocessing
- `bwa` for alignment

## Installation

Use [uv](https://github.com/astral-sh/uv) for lightning-fast installation:

```bash
uv pip install umierrorcorrect
```

Or standard pip:

```bash
pip install umierrorcorrect
```

## Quick Start

Run the full pipeline on a single sample:

```bash
umierrorcorrect batch \
    -r1 sample_R1.fastq.gz \
    -r2 sample_R2.fastq.gz \
    -r hg38.fa \
    -o results/ \
    -ul 12 \
    -sl 16 \
    --fastp
```

For detailed instructions, see the **[User Guide](docs/USER_GUIDE.md)** or run:

```bash
umierrorcorrect --help
```

## Documentation

- [User Guide](docs/USER_GUIDE.md): Detailed usage instructions for all commands.
- [Implementation Details](docs/IMPLEMENTATION.md): Architecture and design overview.
- [Docker Guide](docs/docker.md): Running with containers.

## Citation

> Osterlund T., Filges S., Johansson G., Stahlberg A. *UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular Identifiers*, Clinical Chemistry, 2022. [doi:10.1093/clinchem/hvac136](https://doi.org/10.1093/clinchem/hvac136)
