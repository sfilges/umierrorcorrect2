# UMIErrorCorrect User Guide

Welcome to the **UMIErrorCorrect** user guide. This document provides detailed instructions on installing, configuring, and running the pipeline.

## Table of Contents

1. [Installation](#1-installation)
2. [Pipeline Overview](#2-pipeline-overview)
3. [Running the Pipeline (Batch Mode)](#3-running-the-pipeline-batch-mode)
4. [Input Modes](#4-input-modes)
5. [Configuration Options](#5-configuration-options)
6. [Output Files](#6-output-files)
7. [Running Individual Steps](#7-running-individual-steps)

---

## 1. Installation

### Prerequisites

- **Python 3.9+**
- **External Tools**:
  - `bwa` (Burrows-Wheeler Aligner) must be in your PATH.
  - `samtools` (optional but recommended)
  - `gzip` / `pigz`

### Methods

#### Option A: Using `uv` (Recommended)

`uv` is a modern Python package installer that is significantly faster than pip.

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create a virtual environment and install
uv venv
source .venv/bin/activate
uv pip install umierrorcorrect2
```

#### Option B: Standard pip

```bash
pip install umierrorcorrect2
```

#### Option C: Docker

```bash
docker pull ghcr.io/sfilges/umierrorcorrect:latest
```

---

## 2. Pipeline Overview

The pipeline consists of the following sequential steps:

1. **Preprocessing**:
    - Quality filtering (via `fastp`).
    - Adapter trimming (via `fastp` or `cutadapt`).
    - UMI extraction: Moves UMI sequences from read to header.
2. **Mapping**:
    - Aligns reads to the reference genome using `bwa mem`.
3. **UMI Clustering**:
    - Groups reads mapping to the same position.
    - Clusters UMIs within groups using adjacency-based clustering (Hamming/Levenshtein distance).
4. **Consensus Generation**:
    - Collapses reads in each UMI cluster into a single consensus sequence.
    - Calculates consensus quality scores.
5. **Variant Calling**:
    - Detects variants in the error-corrected consensus reads.

---

## 3. Running the Pipeline (Batch Mode)

The `run` command is the main entry point for running the full pipeline. It handles sample discovery, parallel execution, and logging.

### Basic Syntax

```bash
umierrorcorrect2 run [INPUT_MODE] -r reference.fa -o output_dir/ [OPTIONS]
```

### Essential Arguments

| Argument | Description | Required? |
| :--- | :--- | :--- |
| `-r`, `--reference` | Path to indexed reference genome FASTA | Yes |
| `-o`, `--output-dir` | Directory where results will be saved | Yes |
| `-ul`, `--umi-length` | Length of the UMI sequence (bp) | Yes (default: 19) |
| `-sl`, `--spacer-length` | Length of spacer between UMI and read (bp) | Yes (default: 16) |
| `-t`, `--threads` | Total number of CPU threads to use | No (default: 8) |

---

## 4. Input Modes

You must specify exactly one of the following input modes:

### A. Single Sample Mode

Process a specific pair of FASTQ files.

```bash
umierrorcorrect2 run \
    -r1 sample_R1.fastq.gz \
    -r2 sample_R2.fastq.gz \
    -r hg38.fa \
    -o results/
```

### B. Directory Mode

Automatically discover all FASTQ pairs in a directory.

```bash
umierrorcorrect2 run \
    -i /path/to/fastq_directory/ \
    -r hg38.fa \
    -o results/
```

### C. Sample Sheet Mode

Use a CSV/TSV file to define samples.

```bash
umierrorcorrect2 run \
    --sample-sheet samples.csv \
    -r hg38.fa \
    -o results/
```

**Sample Sheet Format**:
A CSV with columns: `sample_name`, `read1`, `read2` (optional).

```csv
sample_name,read1,read2
SampleA,data/A_R1.fq.gz,data/A_R2.fq.gz
SampleB,data/B_R1.fq.gz,data/B_R2.fq.gz
```

---

## 5. Configuration Options

### Preprocessing & QC

- `--fastp / --no-fastp`: Enable/disable `fastp` (default: Enabled).
- `--trim-adapters`: Trim adapters (default: Enabled).
- `--qc`: Generate FastQC and MultiQC reports (default: Enabled).
- `-q, --phred-score`: Min Q-score for `fastp` (default: 20).

### UMI Handling

- `-d, --edit-distance`: Max edit distance for UMI clustering (default: 1).
- `--dual-index`: Use if UMIs are present on both R1 and R2.
- `--reverse-index`: Use if UMI is on R2 instead of R1.

---

## 6. Output Files

The output directory will process the following structure:

```text
results/
├── batch_summary.tsv       # Summary of all processed samples
├── qc_reports/             # MultiQC and FastQC reports
│   └── multiqc_report.html
├── sample_name/            # Per-sample results
│   ├── sample_name.bam                  # Raw aligned reads
│   ├── sample_name_consensus_reads.bam  # Error-corrected consensus reads
│   ├── sample_name.cons                 # Consensus statistics (TSV)
│   ├── sample_name.vcf                  # Called variants
│   ├── sample_name_umi_stats.txt        # UMI Grouping statistics
│   └── logs/
└── ...
```

### Key Files Explained

- `_consensus_reads.bam`: The "clean" BAM file. Contains one read per UMI family (consensus). Use this for IGV visualization.
- `.cons`: A tabular file containing base counts (A, C, G, T, N, Indels) for every genomic position covered.
- `.vcf`: Variants called from the consensus data.

---

## 7. Running Individual Steps

For advanced users, pipeline steps can be run individually.

### Preprocess

```bash
umierrorcorrect2 preprocess -r1 input_R1.fq.gz -o output_dir/ -ul 12 -sl 16
```

### Mapping

```bash
umierrorcorrect2 mapping -r1 input_R1.fq.gz -r reference.fa -o output_dir/
```

### Consensus

```bash
umierrorcorrect2 consensus -b mapped.bam -o output_dir/
```

### Statistics

```bash
umierrorcorrect2 stats -o output_dir/ --cons-bam consensus.bam
```

### Variant Calling

```bash
umierrorcorrect2 variants -o output_dir/ -r reference.fa --cons sample.cons
```
