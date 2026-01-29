# UMIErrorCorrect2

[![PyPI version](https://badge.fury.io/py/umierrorcorrect2.svg)](https://badge.fury.io/py/umierrorcorrect2)
[![CI](https://github.com/sfilges/umierrorcorrect2/actions/workflows/ci.yml/badge.svg)](https://github.com/sfilges/umierrorcorrect2/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/sfilges/umierrorcorrect2/branch/master/graph/badge.svg?token=)](https://codecov.io/gh/sfilges/umierrorcorrect2)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
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

### Mandatory

- [bwa](https://github.com/lh3/bwa) for alignment

### Optional

- [fastp](https://github.com/OpenGene/fastp) for preprocessing
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control
- [multiqc](https://seqera.io/multiqc/) for quality control / report aggregation

Fastp is **highly recommended**, but not mandatory, for preprocessing. If you do not have fastp installed or run with `--no-fastp`, the pipeline will use `cutadapt` for adapter trimming only.

The `--no-qc` flag disables quality control steps. If QC is enabled (default) but fastqc or multiqc are not installed, the pipeline will raise a warning but finish successfully.

## Installation

Use [uv](https://github.com/astral-sh/uv) for lightning-fast installation:

```bash
uv pip install umierrorcorrect2
```

Or standard pip:

```bash
pip install umierrorcorrect2
```

## Quick Start

The command-line tool is named `umierrorcorrect`. Run the full pipeline on a single sample:

```bash
umierrorcorrect batch \
    -r1 sample_R1.fastq.gz \
    -r2 sample_R2.fastq.gz \
    -r hg38.fa \
    -o results/
```

Run the pipeline on multiple samples in a folder (searches recursively for FASTQ files):

```bash
umierrorcorrect batch \
    -i folder_with_fastq_files/ \
    -r hg38.fa \
    -o results/
```

For detailed instructions, see the **[User Guide](docs/USER_GUIDE.md)** or run:

```bash
umierrorcorrect
```

## Documentation

- [User Guide](docs/USER_GUIDE.md): Detailed usage instructions for all commands.
- [Docker Guide](docs/DOCKER.md): Running with containers.
- [Implementation Details](docs/IMPLEMENTATION.md): Architecture and design overview.

## Citation

> Osterlund T., Filges S., Johansson G., Stahlberg A. *UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular Identifiers*, Clinical Chemistry, 2022. [doi:10.1093/clinchem/hvac136](https://doi.org/10.1093/clinchem/hvac136)
